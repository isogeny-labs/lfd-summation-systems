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


@pytest.fixture(scope='module',params=['baseline','indexed'])
def rust_checker(tmp_path_factory,request):
    import subprocess
    out=tmp_path_factory.mktemp('rust-replay')/'checker'
    subprocess.run(['rustc','-O',str(HERE/('replay.rs' if request.param=='baseline' else 'replay_indexed.rs')),'-o',str(out)],check=True)
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

@pytest.fixture(scope='module',params=['heap','indexed'])
def heap_producer(tmp_path_factory,request):
    import subprocess
    out=tmp_path_factory.mktemp('heap-producer')/'producer'
    subprocess.run(['rustc','-O',str(HERE/(request.param+'_producer.rs')),'-o',str(out)],check=True)
    return out


@pytest.mark.parametrize('compressed',[False,True])
def test_heap_certificate_matches_sample(tmp_path,heap_producer,compressed):
    import subprocess,gzip,json
    inp,ledger,basis=clone(tmp_path)
    expected=[ledger.read_bytes(),basis.read_bytes()]
    paths=[tmp_path/('new-'+name+('.gz' if compressed else '')) for name in ['derivations.txt','basis.txt']]
    r=subprocess.run([str(heap_producer),str(inp),*map(str,paths)],capture_output=True,text=True,check=True)
    assert json.loads(r.stdout)['status']=='closed'
    decoded=[gzip.decompress(p.read_bytes()) if compressed else p.read_bytes() for p in paths]
    assert decoded==expected
    for path,data in zip([ledger,basis],decoded):path.write_bytes(data)
    assert replay(inp,ledger,basis,hashlib.sha256(basis.read_bytes()).hexdigest(),closed=True)['one_excluded']


def test_rust_gzip_replay_and_truncation(tmp_path,rust_checker):
    import gzip,subprocess
    inp,ledger,basis=clone(tmp_path)
    paths=[]
    for p in [ledger,basis]:
        dest=p.with_suffix('.txt.gz');dest.write_bytes(gzip.compress(p.read_bytes()));paths.append(dest)
    cmd=[str(rust_checker),str(inp),*map(str,paths),'closed']
    assert subprocess.run(cmd,capture_output=True).returncode==0
    paths[0].write_bytes(paths[0].read_bytes()[:-8])
    assert subprocess.run(cmd,capture_output=True).returncode!=0


@pytest.mark.parametrize('name',['gate-k2-D4','gate-k3-D4','gate-m3-k2-D5','gate-m4-k2-D5','gate-m5-k3-D5'])
def test_heap_matches_baseline_on_chain_and_split_inputs(tmp_path,heap_producer,rust_checker,name):
    import subprocess,json
    baseline=tmp_path/'baseline-producer'
    subprocess.run(['rustc','-O',str(HERE/'sparse_producer.rs'),'-o',str(baseline)],check=True)
    inp=HERE/'inputs'/(name+'.txt')
    results=[]
    for tag,binary in [('baseline',baseline),('heap',heap_producer)]:
        ledger=tmp_path/(tag+'.ledger');basis=tmp_path/(tag+'.basis')
        result=json.loads(subprocess.run([str(binary),str(inp),str(ledger),str(basis)],capture_output=True,text=True,check=True).stdout)
        assert result['status'] in ['closed','unit-witness']
        results.append((ledger.read_bytes(),basis.read_bytes(),result['contains_one']))
    assert results[0]==results[1]
    cmd=[str(rust_checker),str(inp),str(tmp_path/'heap.ledger'),str(tmp_path/'heap.basis')]
    if not results[1][2]:cmd.append('closed')
    assert subprocess.run(cmd,capture_output=True).returncode==0


@pytest.mark.parametrize('corrupt',['input','compressed','decoded'])
def test_compressed_receipt_binding_refuses_mutation(tmp_path,corrupt):
    import gzip,json
    from verify_compressed import verify
    inp,ledger,basis=clone(tmp_path)
    artifacts={}
    for p in [ledger,basis]:
        dest=p.with_suffix('.txt.gz');dest.write_bytes(gzip.compress(p.read_bytes()))
        artifacts[dest.name]={'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'decoded_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    receipt={'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'artifacts':artifacts}
    if corrupt=='input':inp.write_text(inp.read_text()+' 99')
    elif corrupt=='compressed':(tmp_path/'basis.txt.gz').write_bytes(b'changed')
    else:receipt['artifacts']['basis.txt.gz']['decoded_sha256']='wrong'
    (tmp_path/'result.json').write_text(json.dumps(receipt))
    with pytest.raises(ValueError,match='digest mismatch'):verify(tmp_path)


def test_indexed_engine_refuses_oversized_domain_before_writing(tmp_path):
    import subprocess
    binary=tmp_path/'indexed'
    subprocess.run(['rustc','-O',str(HERE/'indexed_producer.rs'),'-o',str(binary)],check=True)
    inp=tmp_path/'input';inp.write_text('7340033 2 18 0 7 1 100 1000\n0\n')
    ledger=tmp_path/'ledger';basis=tmp_path/'basis'
    result=subprocess.run([str(binary),str(inp),str(ledger),str(basis)],capture_output=True,text=True)
    assert result.returncode!=0 and 'indexed-space cap exceeded' in result.stderr
    assert not ledger.exists() and not basis.exists()


@pytest.fixture(scope='module')
def restart_producer(tmp_path_factory):
    import subprocess
    binary=tmp_path_factory.mktemp('restart')/'producer'
    subprocess.run(['rustc','-O',str(HERE/'restart_producer.rs'),'-o',str(binary)],check=True)
    return binary


def test_capped_proof_resumes_to_checked_closure(tmp_path,restart_producer,rust_checker):
    import subprocess,json,gzip
    # Stop after the generator; restart must discover additional independent
    # shifts and prove closure rather than merely re-label the seed as exact.
    inp,_,_=clone(tmp_path);tokens=inp.read_text().split();original=tokens[:];tokens[6]='1';inp.write_text(' '.join(tokens)+'\n')
    ledger=tmp_path/'seed.ledger.gz';basis=tmp_path/'seed.basis.gz'
    first=json.loads(subprocess.run([str(restart_producer),str(inp),str(ledger),str(basis)],capture_output=True,text=True,check=True).stdout)
    assert first['status']=='capped' and first['derived_rank']==1
    assert subprocess.run([str(rust_checker),str(inp),str(ledger),str(basis)],capture_output=True).returncode==0
    inp.write_text(' '.join(original)+'\n');outledger=tmp_path/'resumed.ledger.gz';outbasis=tmp_path/'resumed.basis.gz'
    second=json.loads(subprocess.run([str(restart_producer),str(inp),str(outledger),str(outbasis),str(ledger),str(basis)],capture_output=True,text=True,check=True).stdout)
    assert second['status']=='closed' and second['derived_rank']==5
    assert gzip.decompress(outledger.read_bytes()).startswith(gzip.decompress(ledger.read_bytes()))
    assert gzip.decompress(outbasis.read_bytes()).startswith(gzip.decompress(basis.read_bytes()))
    assert subprocess.run([str(rust_checker),str(inp),str(outledger),str(outbasis),'closed'],capture_output=True).returncode==0


def test_forged_seed_is_refused_by_independent_final_replay(tmp_path,restart_producer,rust_checker):
    import subprocess,gzip
    inp,ledger,basis=clone(tmp_path)
    lines=basis.read_text().splitlines();row=lines[0].split();row[-1]=str((int(row[-1])+1)%65537);lines[0]=' '.join(row)
    seed=tmp_path/'forged.gz';seed.write_bytes(gzip.compress(('\n'.join(lines)+'\n').encode()))
    outledger=tmp_path/'output.ledger.gz';outbasis=tmp_path/'output.basis.gz'
    result=subprocess.run([str(restart_producer),str(inp),str(outledger),str(outbasis),str(ledger),str(seed)],capture_output=True)
    if result.returncode==0:
        assert subprocess.run([str(rust_checker),str(inp),str(outledger),str(outbasis),'closed'],capture_output=True).returncode!=0


def test_restart_wrapper_refuses_different_polynomial_system(tmp_path,monkeypatch):
    import sys
    import run_heap as runner
    import verify_compressed
    seed=tmp_path/'seed';seed.mkdir();original=(SRC/'input.txt').read_text();(seed/'input.txt').write_text(original)
    inp=tmp_path/'changed-input';tokens=original.split();tokens[-1]=str((int(tokens[-1])+1)%65537);inp.write_text(' '.join(tokens)+'\n')
    monkeypatch.setattr(verify_compressed,'verify',lambda path,**kwargs:None)
    monkeypatch.setattr(runner,'HERE',tmp_path)
    monkeypatch.setattr(sys,'argv',['run_heap.py',str(inp),'new-run','--engine','restart','--resume',str(seed)])
    with pytest.raises(ValueError,match='polynomial system or degree'):runner.main()
    assert not (tmp_path/'runs'/'new-run').exists()


def test_partial_source_cursor_resumes_exact_proof(tmp_path,restart_producer,rust_checker):
    import subprocess,json,gzip
    inp,_,_=clone(tmp_path);original=inp.read_text();tokens=original.split();tokens[6]='2';inp.write_text(' '.join(tokens)+'\n')
    ledger=tmp_path/'seed-ledger.gz';basis=tmp_path/'seed-basis.gz'
    seed=json.loads(subprocess.run([str(restart_producer),str(inp),str(ledger),str(basis)],capture_output=True,text=True,check=True).stdout)
    assert seed['status']=='capped' and seed['derived_rank']==2 and seed['products']%4!=0
    inp.write_text(original);outledger=tmp_path/'ledger.gz';outbasis=tmp_path/'basis.gz'
    result=json.loads(subprocess.run([str(restart_producer),str(inp),str(outledger),str(outbasis),str(ledger),str(basis),str(seed['products'])],capture_output=True,text=True,check=True).stdout)
    assert result['status']=='closed' and result['derived_rank']==5
    assert gzip.decompress(outledger.read_bytes())==(SRC/'derivations.txt').read_bytes()
    assert gzip.decompress(outbasis.read_bytes())==(SRC/'basis.txt').read_bytes()
    assert subprocess.run([str(rust_checker),str(inp),str(outledger),str(outbasis),'closed'],capture_output=True).returncode==0


def test_verified_prefix_cache_is_bound_to_current_artifacts(tmp_path,monkeypatch):
    import gzip,json
    import verify_compressed as checker
    inp,ledger,basis=clone(tmp_path);artifacts={}
    for p in [ledger,basis]:
        dest=p.with_suffix('.txt.gz');dest.write_bytes(gzip.compress(p.read_bytes()))
        artifacts[dest.name]={'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'decoded_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    receipt={'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'artifacts':artifacts,'derived_rank':5,'contains_one':False,'status':'closed'}
    (tmp_path/'result.json').write_text(json.dumps(receipt))
    output=tmp_path/'checker-home';output.mkdir();(output/'results').mkdir();(output/'replay.rs').write_bytes((HERE/'replay.rs').read_bytes())
    monkeypatch.setattr(checker,'HERE',output)
    assert checker.verify(tmp_path,closed=True)['closed_check']
    assert checker.verify(tmp_path,closed=True,reuse=True)['cache_reused']
    (tmp_path/'basis.txt.gz').write_bytes(gzip.compress(b'changed'))
    with pytest.raises(ValueError,match='compressed digest mismatch'):checker.verify(tmp_path,closed=True,reuse=True)


def test_false_cursor_cannot_pass_independent_closedness(tmp_path,restart_producer,rust_checker):
    import subprocess,json,gzip
    inp,_,_=clone(tmp_path);original=inp.read_text();tokens=original.split();tokens[6]='2';inp.write_text(' '.join(tokens)+'\n')
    ledger=tmp_path/'seed-ledger.gz';basis=tmp_path/'seed-basis.gz'
    subprocess.run([str(restart_producer),str(inp),str(ledger),str(basis)],capture_output=True,check=True)
    rows=[list(map(int,line.split())) for line in gzip.decompress(basis.read_bytes()).decode().splitlines()]
    eligible=sum(max(x.bit_count() for x in row[2::2])<3 for row in rows)
    inp.write_text(original);outledger=tmp_path/'ledger.gz';outbasis=tmp_path/'basis.gz'
    # Pretend all admissible source rows have been processed. Membership can
    # still be valid, but exactness must fail the independent closure gate.
    result=json.loads(subprocess.run([str(restart_producer),str(inp),str(outledger),str(outbasis),str(ledger),str(basis),str(4*eligible)],capture_output=True,text=True,check=True).stdout)
    assert result['status']=='closed' and result['derived_rank']==2
    assert subprocess.run([str(rust_checker),str(inp),str(outledger),str(outbasis),'closed'],capture_output=True).returncode!=0
