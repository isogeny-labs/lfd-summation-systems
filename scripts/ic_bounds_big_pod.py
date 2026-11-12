#!/usr/bin/env python3
"""Launcher and remote control for ONE rented Runpod CPU pod (IC bounds, k=14 runs).

Standard library only. Runpod calls and pod calls both go through curl (this
container's egress proxy answers urllib POSTs with 403). The Runpod API key
(env RUNPOD_API_KEY) reaches curl only as a config file on stdin. The pod holds
no account secret: only a random per-pod control token (env CTL_TOKEN on the pod),
kept locally in <state dir>/<name>.token (mode 0600).

  create --name N [--flavors cpu5m,cpu3m] [--vcpu 32] [--volume-gb 300] [--disk-gb 40]
         [--max-seconds 604800] [--dry-run | --yes-create]
  ls [--all] | status ID | stop ID --yes | terminate ID --yes | prices
  ping N | push N LOCAL [REMOTE] | run N JOB -- python3 ARGS... | jobs N
  tail N JOB [--bytes B] | signal N JOB TERM|KILL | get N REMOTE LOCAL [--max-bytes B] | sys N
  ctl-source [--out FILE]        write the control server source (for testing)

Global option --base-url URL overrides https://<pod-id>-8000.proxy.runpod.net (tests).
State dir: env POD_STATE_DIR (default: scratchpad/pod).
"""
import argparse, base64, hashlib, json, os, secrets, subprocess, sys, tempfile, time

API = "https://rest.runpod.io/v1"
GQL = "https://api.runpod.io/graphql"
IMAGE = "runpod/base:1.0.2-ubuntu2404"
STATE_DIR = os.environ.get(
    "POD_STATE_DIR",
    "/tmp/claude-0/-home-user/a6e79f2a-a5ec-5904-b377-343a14909fff/scratchpad/pod")

# --------------------------------------------------------------------------
# Control server (runs on the pod as /root/ctl.py). Standard library only.
# --------------------------------------------------------------------------
CTL_SOURCE = r'''#!/usr/bin/env python3
"""Pod control server: token-guarded file transfer and job runner confined to ROOT.

Env: CTL_TOKEN (required), CTL_ROOT (default /workspace), CTL_PORT (default 8000),
CTL_HOST (default 0.0.0.0). Every request needs "Authorization: Bearer <token>".
Mode "--run EXITFILE -- argv..." is the per-job wrapper (records the exit code).
"""
import hashlib, hmac, json, os, re, shutil, signal, subprocess, sys, tempfile, threading, time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

ROOT = os.path.realpath(os.environ.get("CTL_ROOT", "/workspace"))
PORT = int(os.environ.get("CTL_PORT", "8000"))
HOST = os.environ.get("CTL_HOST", "0.0.0.0")
TOKEN = os.environ.get("CTL_TOKEN", "")
JOBS_FILE = os.path.join(ROOT, "jobs.json")
CODE_DIR = os.path.join(ROOT, "code")
LOGS_DIR = os.path.join(ROOT, "logs")
EXIT_DIR = os.path.join(ROOT, "jobs")
UPLOAD_TOPS = ("code", "in")
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,99}$")
LOCK = threading.RLock()
JOBS = {}
PROCS = {}


class Bad(Exception):
    def __init__(self, code, msg):
        Exception.__init__(self, msg)
        self.code, self.msg = code, msg


def safe_path(rel):
    """Map a client path to a real path under ROOT, or raise Bad(400)."""
    if not isinstance(rel, str) or not rel or "\x00" in rel:
        raise Bad(400, "bad path")
    if rel.startswith("/"):
        if rel != ROOT and not rel.startswith(ROOT + "/"):
            raise Bad(400, "absolute path outside root")
        rel = os.path.relpath(rel, ROOT) if rel != ROOT else "."
    if any(p == ".." for p in rel.split("/")):
        raise Bad(400, "path traversal")
    full = os.path.realpath(os.path.join(ROOT, rel))
    if full != ROOT and not full.startswith(ROOT + os.sep):
        raise Bad(400, "path escapes root")
    return full


def job_name(n):
    if not isinstance(n, str) or not NAME_RE.match(n):
        raise Bad(400, "bad job name")
    return n


def proc_stat(pid):
    try:
        with open("/proc/%d/stat" % pid) as f:
            s = f.read()
        rest = s[s.rindex(")") + 2:].split()
        return rest[0], rest[19]
    except Exception:
        return None


def job_alive(j):
    st = proc_stat(j["pid"])
    return bool(st) and st[0] != "Z" and st[1] == j.get("pstart")


def save_jobs():
    tmp = JOBS_FILE + ".tmp"
    with open(tmp, "w") as f:
        json.dump(JOBS, f, indent=1)
    os.replace(tmp, JOBS_FILE)


def load_jobs():
    try:
        with open(JOBS_FILE) as f:
            JOBS.update(json.load(f))
    except Exception:
        pass


def exit_file(name):
    return os.path.join(EXIT_DIR, name + ".exit")


def reap_once():
    changed = False
    with LOCK:
        for name, j in JOBS.items():
            if j.get("state") != "running":
                continue
            p = PROCS.get(name)
            if p is not None:
                p.poll()
            code, have = None, False
            try:
                with open(exit_file(name)) as f:
                    code, have = int(f.read().strip()), True
            except Exception:
                pass
            if have or not job_alive(j):
                if not have:
                    try:  # the file may land just as the process ends
                        with open(exit_file(name)) as f:
                            code, have = int(f.read().strip()), True
                    except Exception:
                        pass
                if not have and j.get("signaled") == "KILL":
                    code, have = -9, True
                j["state"] = "finished"
                j["exit_code"] = code
                j["end"] = time.time()
                if not have:
                    j["note"] = "no exit code recorded"
                PROCS.pop(name, None)
                changed = True
        if changed:
            save_jobs()


def reaper():
    while True:
        try:
            reap_once()
        except Exception as e:
            sys.stderr.write("reaper: %r\n" % (e,))
        time.sleep(1.0)


def start_job(name, argv):
    job_name(name)
    if (not isinstance(argv, list) or not argv or not all(isinstance(a, str) and "\x00" not in a for a in argv)
            or len(argv) > 500):
        raise Bad(400, "argv must be a non-empty list of strings")
    if argv[0] != "python3":
        raise Bad(400, "only python3 jobs are accepted")
    with LOCK:
        j = JOBS.get(name)
        if j and j.get("state") == "running":
            reap_once()
            j = JOBS.get(name)
            if j and j.get("state") == "running":
                raise Bad(409, "job %s is running" % name)
        for d in (CODE_DIR, LOGS_DIR, EXIT_DIR):
            os.makedirs(d, exist_ok=True)
        try:
            os.unlink(exit_file(name))
        except FileNotFoundError:
            pass
        env = {k: v for k, v in os.environ.items() if k not in ("CTL_TOKEN", "CTL_B64")}
        log = open(os.path.join(LOGS_DIR, name + ".log"), "wb")
        try:
            p = subprocess.Popen(
                [sys.executable, os.path.abspath(__file__), "--run", exit_file(name), "--"] + argv,
                cwd=CODE_DIR, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                start_new_session=True, env=env, shell=False, close_fds=True)
        finally:
            log.close()
        st = proc_stat(p.pid)
        JOBS[name] = {"name": name, "pid": p.pid, "pstart": st[1] if st else None, "command": argv,
                      "start": time.time(), "state": "running", "exit_code": None}
        PROCS[name] = p
        save_jobs()
        return dict(JOBS[name])


def signal_job(name, sig):
    job_name(name)
    if sig not in ("TERM", "KILL"):
        raise Bad(400, "signal must be TERM or KILL")
    with LOCK:
        j = JOBS.get(name)
        if not j:
            raise Bad(404, "no such job")
        if j.get("state") != "running" or not job_alive(j):
            raise Bad(409, "job is not running")
        j["signaled"] = sig
        save_jobs()
        os.killpg(j["pid"], signal.SIGTERM if sig == "TERM" else signal.SIGKILL)
    return {"name": name, "signal": sig}


def jobs_view():
    reap_once()
    now = time.time()
    with LOCK:
        out = []
        for j in JOBS.values():
            v = dict(j)
            v["elapsed_s"] = round((j.get("end") or now) - j["start"], 1)
            out.append(v)
    return sorted(out, key=lambda v: v["start"])


def sys_info():
    def rd(p):
        with open(p) as f:
            return f.read()
    mem = {}
    for line in rd("/proc/meminfo").splitlines():
        k, _, v = line.partition(":")
        if k in ("MemTotal", "MemAvailable"):
            mem[k] = int(v.split()[0]) // 1024
    du = shutil.disk_usage(ROOT)
    info = {"uptime_s": float(rd("/proc/uptime").split()[0]), "nproc": os.cpu_count(),
            "affinity": len(os.sched_getaffinity(0)), "mem_total_mb": mem.get("MemTotal"),
            "mem_available_mb": mem.get("MemAvailable"),
            "workspace_disk_gb": {"total": round(du.total / 1e9, 1), "used": round(du.used / 1e9, 1),
                                  "free": round(du.free / 1e9, 1)},
            "loadavg": rd("/proc/loadavg").split()[:3], "python": sys.version.split()[0]}
    probe = ("import json,numpy\nr={'numpy':numpy.__version__}\n"
             "try:\n import threadpoolctl\n r['threadpoolctl']=threadpoolctl.__version__\n"
             " r['blas']=[[d.get('internal_api'),d.get('version'),d.get('num_threads')] for d in threadpoolctl.threadpool_info()]\n"
             "except Exception as e:\n r['blas_error']=repr(e)\n"
             "try:\n c=numpy.show_config(mode='dicts')['Build Dependencies']['blas']\n r['numpy_build_blas']=[c.get('name'),c.get('version')]\n"
             "except Exception as e:\n r['numpy_build_blas_error']=repr(e)\nprint(json.dumps(r))")
    try:
        o = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True, timeout=60)
        info["numpy_probe"] = json.loads(o.stdout) if o.returncode == 0 else {"error": o.stderr[-300:]}
    except Exception as e:
        info["numpy_probe"] = {"error": repr(e)}
    return info


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class H(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    server_version = "podctl/1"

    def log_message(self, fmt, *args):
        sys.stderr.write("%s %s\n" % (self.address_string(), fmt % args))

    def send_json(self, obj, code=200):
        b = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(b)))
        self.end_headers()
        self.wfile.write(b)

    def authed(self):
        h = self.headers.get("Authorization", "")
        want = "Bearer " + TOKEN
        return bool(TOKEN) and hmac.compare_digest(h.encode(), want.encode())

    def dispatch(self, method):
        if not self.authed():  # before routing and before reading any body
            self.close_connection = True
            return self.send_json({"error": "unauthorized"}, 401)
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            fn = ROUTES.get((method, u.path))
            if fn is None:
                raise Bad(404, "no such endpoint")
            fn(self, q)
        except Bad as e:
            self.close_connection = True
            self.send_json({"error": e.msg}, e.code)
        except (BrokenPipeError, ConnectionResetError):
            self.close_connection = True
        except Exception as e:
            self.close_connection = True
            try:
                self.send_json({"error": "internal: %r" % (e,)}, 500)
            except Exception:
                pass

    def do_GET(self): self.dispatch("GET")
    def do_PUT(self): self.dispatch("PUT")
    def do_POST(self): self.dispatch("POST")
    def do_DELETE(self): self.dispatch("DELETE")
    def do_HEAD(self): self.dispatch("HEAD")

    def read_json(self):
        n = int(self.headers.get("Content-Length") or 0)
        if n > 1 << 20:
            raise Bad(413, "body too large")
        try:
            return json.loads(self.rfile.read(n) or b"{}")
        except Exception:
            raise Bad(400, "bad json")


def r_ping(h, q):
    h.send_json({"ok": True, "time": time.time(), "root": ROOT})


def r_sys(h, q):
    h.send_json(sys_info())


def r_put_file(h, q):
    full = safe_path(q.get("path"))
    rel = os.path.relpath(full, ROOT)
    if rel.split(os.sep)[0] not in UPLOAD_TOPS or rel == ".":
        raise Bad(400, "uploads go under %s" % "/".join(UPLOAD_TOPS))
    if os.path.isdir(full):
        raise Bad(400, "is a directory")
    n = int(h.headers.get("Content-Length", "-1"))
    if n < 0:
        raise Bad(411, "Content-Length required")
    os.makedirs(os.path.dirname(full), exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".upload-", dir=os.path.dirname(full))
    ok = False
    try:
        sha, left = hashlib.sha256(), n
        with os.fdopen(fd, "wb") as f:
            while left:
                b = h.rfile.read(min(1 << 20, left))
                if not b:
                    raise Bad(400, "short body")
                f.write(b)
                sha.update(b)
                left -= len(b)
        os.chmod(tmp, 0o644)
        os.replace(tmp, full)  # atomic
        ok = True
    finally:
        if not ok:
            try:
                os.unlink(tmp)
            except FileNotFoundError:
                pass
    h.send_json({"path": rel, "size": n, "sha256": sha.hexdigest()})


def regular(full):
    if not os.path.isfile(full):
        raise Bad(404, "not a file")
    return os.path.getsize(full)


def r_stat(h, q):
    full = safe_path(q.get("path"))
    size = regular(full)
    h.send_json({"path": os.path.relpath(full, ROOT), "size": size, "sha256": sha256_file(full)})


def r_get_file(h, q):
    full = safe_path(q.get("path"))
    size = regular(full)
    mx = int(q.get("max", str(20 << 20)))
    if size > mx:
        raise Bad(413, "file is %d bytes, over the limit %d" % (size, mx))
    h.send_response(200)
    h.send_header("Content-Type", "application/octet-stream")
    h.send_header("Content-Length", str(size))
    h.end_headers()
    left = size
    with open(full, "rb") as f:
        while left:
            b = f.read(min(1 << 20, left))
            if not b:
                break
            h.wfile.write(b)
            left -= len(b)
    if left:
        h.close_connection = True


def r_run(h, q):
    d = h.read_json()
    h.send_json(start_job(d.get("name"), d.get("argv")))


def r_jobs(h, q):
    h.send_json({"jobs": jobs_view()})


def r_tail(h, q):
    name = job_name(q.get("name"))
    n = max(1, min(int(q.get("bytes", "4000")), 8 << 20))
    p = os.path.join(LOGS_DIR, name + ".log")
    if not os.path.isfile(p):
        raise Bad(404, "no log")
    size = os.path.getsize(p)
    with open(p, "rb") as f:
        f.seek(max(0, size - n))
        b = f.read()
    h.send_response(200)
    h.send_header("Content-Type", "text/plain; charset=utf-8")
    h.send_header("Content-Length", str(len(b)))
    h.end_headers()
    h.wfile.write(b)


def r_signal(h, q):
    d = h.read_json()
    h.send_json(signal_job(d.get("name"), d.get("signal")))


ROUTES = {("GET", "/ping"): r_ping, ("GET", "/sys"): r_sys, ("PUT", "/file"): r_put_file,
          ("GET", "/file"): r_get_file, ("GET", "/stat"): r_stat, ("POST", "/run"): r_run,
          ("GET", "/jobs"): r_jobs, ("GET", "/tail"): r_tail, ("POST", "/signal"): r_signal}


def run_wrapper(exitfile, argv):
    """Per-job wrapper: runs argv, records its exit code, survives a server restart."""
    noop = lambda *a: None  # handlers (unlike SIG_IGN) are reset in the exec'd child
    for s in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(s, noop)
    try:
        rc = subprocess.call(argv, shell=False)
    except OSError:
        rc = 127
    tmp = exitfile + ".tmp"
    with open(tmp, "w") as f:
        f.write(str(rc))
    os.replace(tmp, exitfile)
    sys.exit(0)


def main():
    if len(sys.argv) > 3 and sys.argv[1] == "--run" and sys.argv[3] == "--":
        run_wrapper(sys.argv[2], sys.argv[4:])
    if not TOKEN:
        sys.exit("CTL_TOKEN is not set")
    for d in (ROOT, CODE_DIR, LOGS_DIR, EXIT_DIR, os.path.join(ROOT, "in"), os.path.join(ROOT, "out")):
        os.makedirs(d, exist_ok=True)
    load_jobs()
    threading.Thread(target=reaper, daemon=True).start()
    srv = ThreadingHTTPServer((HOST, PORT), H)
    srv.daemon_threads = True
    signal.signal(signal.SIGTERM, lambda *a: os._exit(0))
    sys.stderr.write("ctl listening on %s:%d root=%s\n" % (HOST, srv.server_address[1], ROOT))
    srv.serve_forever()


if __name__ == "__main__":
    main()
'''

# --------------------------------------------------------------------------
# Pod start command (bash -c): lifetime cap, control server, pip, then image start.
# --------------------------------------------------------------------------
START = r'''( sleep "$POD_MAX_SECONDS"; echo "lifetime cap reached" >> /workspace/boot.log; kill -TERM 1 ) >/tmp/pod-lifetime.log 2>&1 &
export DEBIAN_FRONTEND=noninteractive
mkdir -p /workspace/code /workspace/in /workspace/out /workspace/logs
echo "$CTL_B64" | base64 -d > /root/ctl.py
( python3 -m pip install -q $PIP_PKGS \
  || python3 -m pip install -q --break-system-packages $PIP_PKGS \
  || { apt-get update -qq && apt-get install -y -qq python3-numpy python3-threadpoolctl; } ) >> /workspace/boot.log 2>&1
setsid python3 /root/ctl.py >> /workspace/ctl.log 2>&1 < /dev/null &
echo "control server started $(date -u +%FT%TZ)" >> /workspace/boot.log
[ -x /start.sh ] && exec /start.sh
exec sleep infinity'''


# --------------------------------------------------------------------------
# curl helpers
# --------------------------------------------------------------------------
def curl(method, url, token=None, *, data=None, upload=None, out=None, timeout=90, extra=()):
    """Run curl; the bearer token goes in a config on stdin. Returns (status, body bytes)."""
    own = out is None
    if own:
        fd, out = tempfile.mkstemp(prefix="podctl-")
        os.close(fd)
    cmd = ["curl", "-sS", "-m", str(timeout), "--connect-timeout", "30", "-K", "-",
           "--noproxy", "127.0.0.1,localhost", "-X", method, url, "-o", out, "-w", "%{http_code}"]
    datafile = None
    if data is not None:
        # the body can carry the pod's control token: keep it off the command line
        fd, datafile = tempfile.mkstemp(prefix="podctl-body-")
        with os.fdopen(fd, "w") as f:
            f.write(data)
        cmd += ["-H", "Content-Type: application/json", "--data-binary", "@" + datafile]
    if upload is not None:
        cmd += ["-T", upload, "-H", "Expect:"]
    cmd += list(extra)
    cfg = 'header = "Authorization: Bearer %s"\n' % token
    try:
        r = subprocess.run(cmd, input=cfg, capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit("curl failed (%d): %s" % (r.returncode, r.stderr.strip()[:300]))
        status = int(r.stdout.strip() or 0)
        body = b""
        if own:
            with open(out, "rb") as f:
                body = f.read()
        return status, body
    finally:
        for tmp in ([out] if own else []) + ([datafile] if datafile else []):
            try:
                os.unlink(tmp)
            except FileNotFoundError:
                pass


def rp(method, path, body=None, url=None):
    """Runpod API call (key from env, passed to curl on stdin only)."""
    key = os.environ.get("RUNPOD_API_KEY")
    if not key:
        raise SystemExit("RUNPOD_API_KEY is not set")
    st, b = curl(method, url or API + path, key, data=json.dumps(body) if body is not None else None)
    try:
        return json.loads(b or b"null")
    except ValueError:
        raise SystemExit("Runpod answered %d: %s" % (st, b[:300]))


# --------------------------------------------------------------------------
# state: token + record per pod name
# --------------------------------------------------------------------------
def check_name(n):
    import re
    if not re.match(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,62}$", n or ""):
        raise SystemExit("bad pod name %r (letters, digits, _ and - only)" % n)
    return n


def tok_path(n): return os.path.join(STATE_DIR, check_name(n) + ".token")
def rec_path(n): return os.path.join(STATE_DIR, check_name(n) + ".json")


def read_token(n):
    try:
        with open(tok_path(n)) as f:
            return f.read().strip()
    except FileNotFoundError:
        raise SystemExit("no token for %r under %s (not created by this tool?)" % (n, STATE_DIR))


def recorded():
    out = []
    if os.path.isdir(STATE_DIR):
        for fn in sorted(os.listdir(STATE_DIR)):
            if fn.endswith(".json"):
                try:
                    with open(os.path.join(STATE_DIR, fn)) as f:
                        rec = json.load(f)
                    if isinstance(rec, dict) and "id" in rec and "name" in rec:   # skip other JSON files
                        out.append(rec)
                except Exception:
                    pass
    return out


def require_recorded(pod_id):
    if not any(r.get("id") == pod_id for r in recorded()):
        raise SystemExit("pod %s is not recorded under %s by `create`; refusing" % (pod_id, STATE_DIR))


def base_url(a):
    if getattr(a, "base_url", None):
        return a.base_url.rstrip("/")
    try:
        with open(rec_path(a.name)) as f:
            return "https://%s-8000.proxy.runpod.net" % json.load(f)["id"]
    except FileNotFoundError:
        raise SystemExit("no record for %r under %s" % (a.name, STATE_DIR))


def call_pod(a, method, path, **kw):
    return curl(method, base_url(a) + path, read_token(a.name), **kw)


def need_ok(st, body, what):
    if st != 200:
        raise SystemExit("%s failed: HTTP %d %s" % (what, st, body.decode("utf-8", "replace")[:300]))
    return json.loads(body)


def sha256_local(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def q(s):
    from urllib.parse import quote
    return quote(s, safe="/")


# --------------------------------------------------------------------------
# Runpod commands
# --------------------------------------------------------------------------
def build_body(a, token):
    src = CTL_SOURCE.encode()
    body = {
        "name": a.name, "cloudType": "SECURE", "imageName": a.image,
        "containerDiskInGb": a.disk_gb, "volumeInGb": a.volume_gb, "volumeMountPath": "/workspace",
        "ports": ["8000/http"], "dockerEntrypoint": ["bash", "-c"], "dockerStartCmd": [START],
        "env": {"POD_MAX_SECONDS": str(a.max_seconds), "CTL_TOKEN": token, "PIP_PKGS": a.pip,
                "CTL_B64": base64.b64encode(src).decode()},
    }
    if getattr(a, "datacenters", None):   # restrict placement, e.g. to data centres where a run was fast
        body["dataCenterIds"] = a.datacenters.split(",")
    if a.gpu:   # GPU pod (same field names as the lab launcher)
        body.update({"computeType": "GPU", "gpuTypeIds": a.gpu.split(","), "gpuCount": a.gpu_count,
                     "minRAMPerGPU": a.min_ram_per_gpu, "minVCPUPerGPU": a.min_vcpu_per_gpu})
    else:
        body.update({"computeType": "CPU", "cpuFlavorIds": a.flavors.split(","), "vcpuCount": a.vcpu})
    return body


def cmd_create(a):
    check_name(a.name)
    if a.dry_run:
        body = build_body(a, "<redacted>")
        print(json.dumps(body, indent=1))
        return
    if not a.yes_create:
        raise SystemExit("refusing to create a pod without --yes-create (use --dry-run to inspect)")
    os.makedirs(STATE_DIR, mode=0o700, exist_ok=True)
    if os.path.exists(tok_path(a.name)):
        raise SystemExit("a pod named %r is already recorded; choose another name" % a.name)
    token = secrets.token_urlsafe(32)
    fd = os.open(tok_path(a.name), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, "w") as f:
        f.write(token + "\n")
    r = rp("POST", "/pods", build_body(a, token))
    if not isinstance(r, dict) or "id" not in r:
        os.unlink(tok_path(a.name))
        raise SystemExit("pod not created: %s" % json.dumps(r)[:600])
    rec = {"id": r["id"], "name": a.name, "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "max_seconds": a.max_seconds, "costPerHr": r.get("costPerHr")}
    with open(rec_path(a.name), "w") as f:
        json.dump(rec, f, indent=1)
    cost = r.get("costPerHr")
    print(json.dumps({"id": r["id"], "name": r.get("name"), "flavor": r.get("cpuFlavorId"),
                      "gpu": (r.get("gpu") or {}).get("displayName") if isinstance(r.get("gpu"), dict) else r.get("gpu"),
                      "machine": {k: (r.get("machine") or {}).get(k) for k in ("gpuTypeId", "location", "dataCenterId")} if isinstance(r.get("machine"), dict) else None,
                      "vcpu": r.get("vcpuCount"), "memoryInGb": r.get("memoryInGb"), "costPerHr": cost,
                      "desiredStatus": r.get("desiredStatus")}))
    if isinstance(cost, (int, float)):
        print("worst case at the lifetime cap (%ds): $%.2f" % (a.max_seconds, cost * a.max_seconds / 3600))
    print("proxy:", "https://%s-8000.proxy.runpod.net/" % r["id"])
    print("next: ping %s   (terminate when done: terminate %s --yes)" % (a.name, r["id"]))


def cmd_ls(a):
    if a.all:
        for p in rp("GET", "/pods"):
            print(p["id"], p.get("name"), p.get("desiredStatus"), p.get("vcpuCount"), p.get("memoryInGb"), p.get("costPerHr"))
        return
    for r in recorded():
        p = rp("GET", "/pods/" + r["id"])
        st = p.get("desiredStatus") if isinstance(p, dict) and "id" in p else "gone/unknown: %s" % json.dumps(p)[:80]
        print(r["id"], r["name"], st, r.get("created"), r.get("costPerHr"))


def cmd_status(a):
    p = rp("GET", "/pods/" + a.id)
    print(json.dumps({k: p.get(k) for k in ("id", "name", "desiredStatus", "vcpuCount", "memoryInGb",
                                           "costPerHr", "lastStatusChange")}) if isinstance(p, dict) else p)


def cmd_stop(a):
    require_recorded(a.id)
    if not a.yes:
        raise SystemExit("stop needs --yes")
    print(rp("POST", "/pods/%s/stop" % a.id))


def cmd_terminate(a):
    require_recorded(a.id)
    if not a.yes:
        raise SystemExit("terminate needs --yes")
    print(rp("DELETE", "/pods/" + a.id))


def cmd_prices(a):
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    for inst in a.instances.split(","):
        fl = inst.split("-")[0]
        qq = '{ cpuFlavors { id specifics(input:{instanceId:"%s"}) { stockStatus securePrice } } }' % inst
        r = rp("POST", "", {"query": qq}, url=GQL)
        for f in r["data"]["cpuFlavors"]:
            if f["id"] == fl:
                print(json.dumps({"queried": stamp, "instance": inst, **f["specifics"]}))


# --------------------------------------------------------------------------
# Remote control commands
# --------------------------------------------------------------------------
def cmd_ping(a):
    t = time.time()
    st, b = call_pod(a, "GET", "/ping", timeout=30)
    d = need_ok(st, b, "ping")
    print(json.dumps({"ok": d.get("ok"), "rtt_s": round(time.time() - t, 3), "root": d.get("root")}))


def cmd_push(a):
    if not os.path.isfile(a.local):
        raise SystemExit("no such file: " + a.local)
    name = a.remote or os.path.basename(a.local)
    if name.startswith("/") or ".." in name.split("/"):
        raise SystemExit("remote name must be relative, without ..")
    dest = a.dest or ("code" if name.endswith(".py") else "in")
    path = "%s/%s" % (dest, name)
    local_sha = sha256_local(a.local)
    st, b = call_pod(a, "PUT", "/file?path=" + q(path), upload=a.local, timeout=a.timeout)
    d = need_ok(st, b, "push")
    if d.get("sha256") != local_sha or d.get("size") != os.path.getsize(a.local):
        raise SystemExit("push verification FAILED: local %s, remote %s" % (local_sha, d.get("sha256")))
    print(json.dumps({"pushed": d["path"], "size": d["size"], "sha256": local_sha, "verified": True}))


def cmd_run(a):
    if not a.cmd or a.cmd[0] != "python3":
        raise SystemExit("usage: run NAME JOB -- python3 args... (only python3 jobs)")
    st, b = call_pod(a, "POST", "/run", data=json.dumps({"name": a.job, "argv": a.cmd}), timeout=60)
    d = need_ok(st, b, "run")
    print(json.dumps({k: d[k] for k in ("name", "pid", "state", "command")}))


def cmd_jobs(a):
    st, b = call_pod(a, "GET", "/jobs", timeout=60)
    for j in need_ok(st, b, "jobs")["jobs"]:
        print("%-24s %-9s exit=%-5s %8.0fs pid=%-7s %s" % (j["name"], j["state"], j.get("exit_code"),
              j.get("elapsed_s", 0), j["pid"], " ".join(j["command"])[:100]))


def cmd_tail(a):
    st, b = call_pod(a, "GET", "/tail?name=%s&bytes=%d" % (q(a.job), a.bytes), timeout=60)
    if st != 200:
        raise SystemExit("tail failed: HTTP %d %s" % (st, b.decode("utf-8", "replace")[:300]))
    sys.stdout.write(b.decode("utf-8", "replace"))


def cmd_signal(a):
    st, b = call_pod(a, "POST", "/signal", data=json.dumps({"name": a.job, "signal": a.sig}), timeout=60)
    print(json.dumps(need_ok(st, b, "signal")))


def cmd_get(a):
    st, b = call_pod(a, "GET", "/stat?path=" + q(a.remote), timeout=600)
    meta = need_ok(st, b, "stat")
    if meta["size"] > a.max_bytes:
        raise SystemExit("remote file is %d bytes, over --max-bytes %d" % (meta["size"], a.max_bytes))
    part = a.local + ".part"
    st, _ = call_pod(a, "GET", "/file?path=%s&max=%d" % (q(a.remote), a.max_bytes), out=part, timeout=a.timeout)
    if st != 200:
        msg = open(part, "rb").read()[:300] if os.path.exists(part) else b""
        if os.path.exists(part):
            os.unlink(part)
        raise SystemExit("get failed: HTTP %d %s" % (st, msg.decode("utf-8", "replace")))
    got = sha256_local(part)
    if got != meta["sha256"] or os.path.getsize(part) != meta["size"]:
        os.unlink(part)
        raise SystemExit("get verification FAILED: remote %s, local %s" % (meta["sha256"], got))
    os.replace(part, a.local)
    print(json.dumps({"saved": a.local, "size": meta["size"], "sha256": got, "verified": True}))


def cmd_sys(a):
    st, b = call_pod(a, "GET", "/sys", timeout=120)
    print(json.dumps(need_ok(st, b, "sys"), indent=1))


def cmd_ctl_source(a):
    if a.out:
        with open(a.out, "w") as f:
            f.write(CTL_SOURCE)
    else:
        sys.stdout.write(CTL_SOURCE)


# --------------------------------------------------------------------------
def main():
    argv = sys.argv[1:]
    tail_cmd = []
    if "--" in argv:
        i = argv.index("--")
        argv, tail_cmd = argv[:i], argv[i + 1:]
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base-url", default=None, help="override the pod URL (tests)")
    sub = p.add_subparsers(dest="cmd", required=True)

    def remote(name, **kw):
        sp = sub.add_parser(name, **kw)
        sp.add_argument("--base-url", default=argparse.SUPPRESS)
        sp.add_argument("name")
        return sp

    c = sub.add_parser("create")
    c.add_argument("--name", required=True)
    c.add_argument("--flavors", default="cpu5m,cpu3m")
    c.add_argument("--vcpu", type=int, default=32)
    c.add_argument("--volume-gb", type=int, default=300)
    c.add_argument("--disk-gb", type=int, default=40)
    c.add_argument("--max-seconds", type=int, default=604800)
    c.add_argument("--image", default=IMAGE)
    c.add_argument("--pip", default="numpy threadpoolctl", help="pip packages installed at boot (space separated)")
    c.add_argument("--gpu", help='GPU type id(s), comma separated, e.g. "NVIDIA B200"; omit for a CPU pod')
    c.add_argument("--gpu-count", type=int, default=1)
    c.add_argument("--datacenters", help="comma separated data centre ids (dataCenterIds), e.g. US-NC-1,EU-FR-1")
    c.add_argument("--min-ram-per-gpu", type=int, default=8)
    c.add_argument("--min-vcpu-per-gpu", type=int, default=2)
    c.add_argument("--dry-run", action="store_true")
    c.add_argument("--yes-create", action="store_true")
    l = sub.add_parser("ls"); l.add_argument("--all", action="store_true")
    sub.add_parser("status").add_argument("id")
    for n in ("stop", "terminate"):
        s = sub.add_parser(n); s.add_argument("id"); s.add_argument("--yes", action="store_true")
    pr = sub.add_parser("prices")
    pr.add_argument("--instances", default="cpu5m-32-256,cpu3m-32-256")
    remote("ping")
    s = remote("push"); s.add_argument("local"); s.add_argument("remote", nargs="?")
    s.add_argument("--dest", choices=["code", "in"]); s.add_argument("--timeout", type=int, default=3600)
    s = remote("run"); s.add_argument("job")
    remote("jobs")
    s = remote("tail"); s.add_argument("job"); s.add_argument("--bytes", type=int, default=4000)
    s = remote("signal"); s.add_argument("job"); s.add_argument("sig", choices=["TERM", "KILL"])
    s = remote("get"); s.add_argument("remote"); s.add_argument("local")
    s.add_argument("--max-bytes", type=int, default=20 << 20); s.add_argument("--timeout", type=int, default=3600)
    remote("sys")
    s = sub.add_parser("ctl-source"); s.add_argument("--out")

    a = p.parse_args(argv)
    handlers = {"create": cmd_create, "ls": cmd_ls, "status": cmd_status, "stop": cmd_stop,
                "terminate": cmd_terminate, "prices": cmd_prices, "ping": cmd_ping, "push": cmd_push,
                "jobs": cmd_jobs, "tail": cmd_tail, "signal": cmd_signal, "get": cmd_get, "sys": cmd_sys,
                "ctl-source": cmd_ctl_source}
    if a.cmd == "run":
        a.cmd = tail_cmd  # cmd_run reads a.cmd as the argv list
        cmd_run(a)
    else:
        handlers[a.cmd](a)


if __name__ == "__main__":
    main()
