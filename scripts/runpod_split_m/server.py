"""Runner for one RunPod CPU pod: receives a bundle, runs one job script, serves status and results.
Standard library only.  Auth: header X-Token must equal $RUNNER_TOKEN.  The job is killed after
$MAX_SECONDS.  Endpoints: PUT /upload/<name>, PUT /append/<name>, POST /start, GET /status, GET /fetch.
"""
import http.server, io, json, os, subprocess, tarfile, threading, time

TOKEN, WORK = os.environ['RUNNER_TOKEN'], '/work'
MAX_SECONDS = int(os.environ.get('MAX_SECONDS', '3600'))
os.makedirs(WORK + '/out', exist_ok=True)
ALLOWED = {n: WORK + '/' + n for n in ('bundle.tar.gz', 'job.sh', 'cell.env')}
state = dict(state='idle', started=None, ended=None, exit=None)
proc = None


def watchdog():
    global proc
    while proc.poll() is None:
        if time.time() - state['started'] > MAX_SECONDS:
            proc.kill()
            state['state'] = 'killed-time-limit'
        time.sleep(5)
    state['exit'], state['ended'] = proc.returncode, time.time()
    if state['state'] == 'running':
        state['state'] = 'done'


class H(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def reply(self, code, body, ctype='application/json'):
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def ok(self):
        if self.headers.get('X-Token') != TOKEN:
            self.reply(403, b'{}')
            return False
        return True

    def do_PUT(self):
        if not self.ok():
            return
        kind, _, name = self.path.strip('/').partition('/')
        path = ALLOWED.get(name)      # fixed destinations only; the request never forms a path
        if kind not in ('upload', 'append') or path is None:
            return self.reply(404, b'{}')
        data = self.rfile.read(int(self.headers['Content-Length']))
        with open(path, 'wb' if kind == 'upload' else 'ab') as f:
            f.write(data)
        self.reply(200, json.dumps(dict(name=os.path.basename(path), size=os.path.getsize(path))).encode())

    def do_POST(self):
        global proc
        if not self.ok():
            return
        if self.path != '/start' or state['state'] != 'idle':
            return self.reply(409, json.dumps(state).encode())
        log = open(WORK + '/out/job.log', 'wb')
        proc = subprocess.Popen(['bash', WORK + '/job.sh'], cwd=WORK, stdout=log, stderr=subprocess.STDOUT)
        state.update(state='running', started=time.time())
        threading.Thread(target=watchdog, daemon=True).start()
        self.reply(200, json.dumps(state).encode())

    def do_GET(self):
        if not self.ok():
            return
        if self.path == '/status':
            tail = b''
            if os.path.exists(WORK + '/out/job.log'):
                with open(WORK + '/out/job.log', 'rb') as f:
                    f.seek(max(0, os.path.getsize(WORK + '/out/job.log') - 1500))
                    tail = f.read()
            mem = open('/proc/meminfo').read().split('\n')[:3]
            body = dict(state, now=time.time(), max_seconds=MAX_SECONDS, tail=tail.decode(errors='replace'), mem=mem,
                        files=sorted(os.listdir(WORK + '/out')))
            return self.reply(200, json.dumps(body).encode())
        if self.path == '/fetch':
            buf = io.BytesIO()
            with tarfile.open(fileobj=buf, mode='w:gz') as t:
                t.add(WORK + '/out', arcname='out')
            return self.reply(200, buf.getvalue(), 'application/gzip')
        self.reply(404, b'{}')


http.server.ThreadingHTTPServer(('0.0.0.0', 8000), H).serve_forever()
