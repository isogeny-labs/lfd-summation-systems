"""End-to-end replay refusal controls, against a separately generated gate."""
import hashlib
from pathlib import Path
import pytest
from replay import replay

HERE=Path(__file__).resolve().parent
SRC=HERE/'sample'


def clone(tmp_path):
    for name in ['input.txt','derivations.txt','basis.txt']:
        (tmp_path/name).write_bytes((SRC/name).read_bytes())
    return [tmp_path/name for name in ['input.txt','derivations.txt','basis.txt']]


def test_closed_certificate_replays(tmp_path):
    inp,ledger,basis=clone(tmp_path)
    r=replay(inp,ledger,basis,hashlib.sha256(basis.read_bytes()).hexdigest(),closed=True)
    assert r['derived_rank']==5 and r['one_excluded'] and not r['contains_one']


def test_forged_output_fails_even_with_changed_digest(tmp_path):
    inp,ledger,basis=clone(tmp_path);lines=basis.read_text().splitlines();t=lines[0].split();t[-1]=str((int(t[-1])+1)%65537);lines[0]=' '.join(t);basis.write_text('\n'.join(lines)+'\n')
    with pytest.raises(ValueError,match='derived row mismatch'):
        replay(inp,ledger,basis,hashlib.sha256(basis.read_bytes()).hexdigest())


def test_future_reference_refused(tmp_path):
    inp,ledger,basis=clone(tmp_path);lines=ledger.read_text().splitlines();lines[1]='1 shift 999 0 1 0';ledger.write_text('\n'.join(lines)+'\n')
    with pytest.raises(ValueError,match='shift reference'):
        replay(inp,ledger,basis,hashlib.sha256(basis.read_bytes()).hexdigest())


def test_missing_row_cannot_pass_target_digest(tmp_path):
    inp,ledger,basis=clone(tmp_path);target=hashlib.sha256(basis.read_bytes()).hexdigest();basis.write_text('')
    with pytest.raises(ValueError,match='target digest'):replay(inp,ledger,basis,target)


def test_degree_D_source_is_not_shiftable(tmp_path):
    inp,ledger,basis=clone(tmp_path)
    # First generator has degree 2. Lower D to 2; its shifted source is then inadmissible.
    t=inp.read_text().split();t[4]='2';inp.write_text(' '.join(t)+'\n')
    with pytest.raises(ValueError,match='inadmissible source degree'):
        replay(inp,ledger,basis,hashlib.sha256(basis.read_bytes()).hexdigest())


@pytest.fixture(scope='module')
def rust_checker(tmp_path_factory):
    import subprocess
    out=tmp_path_factory.mktemp('rust-replay')/'checker'
    subprocess.run(['rustc','-O',str(HERE/'replay.rs'),'-o',str(out)],check=True)
    return out


def test_rust_checker_closed_sample(tmp_path,rust_checker):
    import subprocess,json
    inp,ledger,basis=clone(tmp_path)
    r=subprocess.run([str(rust_checker),str(inp),str(ledger),str(basis),'closed'],capture_output=True,text=True,check=True)
    j=json.loads(r.stdout);assert j['derived_rank']==5 and j['closed_check'] and not j['contains_one']


@pytest.mark.parametrize('corrupt',['basis','reference','degree'])
def test_rust_checker_refuses_corruption(tmp_path,rust_checker,corrupt):
    import subprocess
    inp,ledger,basis=clone(tmp_path)
    if corrupt=='basis':
        lines=basis.read_text().splitlines();t=lines[0].split();t[-1]=str((int(t[-1])+1)%65537);lines[0]=' '.join(t);basis.write_text('\n'.join(lines)+'\n')
    elif corrupt=='reference':
        lines=ledger.read_text().splitlines();lines[1]='1 shift 999 0 1 0';ledger.write_text('\n'.join(lines)+'\n')
    else:
        t=inp.read_text().split();t[4]='2';inp.write_text(' '.join(t)+'\n')
    r=subprocess.run([str(rust_checker),str(inp),str(ledger),str(basis)],capture_output=True,text=True)
    assert r.returncode!=0


def test_runner_refuses_existing_artifact_directory(tmp_path,monkeypatch):
    import run as runner
    monkeypatch.setattr(runner,'HERE',tmp_path)
    old=tmp_path/'runs'/'existing';old.mkdir(parents=True);(old/'result.json').write_text('preserve me')
    with pytest.raises(FileExistsError,match='refusing to overwrite'):
        runner.run_case('existing',17,2,2,0,3,[{0:1}])
    assert (old/'result.json').read_text()=='preserve me'
