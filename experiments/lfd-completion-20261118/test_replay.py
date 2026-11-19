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


@pytest.fixture(scope='module')
def priority_producer(tmp_path_factory):
    import subprocess
    binary=tmp_path_factory.mktemp('priority')/'producer'
    subprocess.run(['rustc','-O',str(HERE/'priority_producer.rs'),'-o',str(binary)],check=True)
    return binary


@pytest.mark.parametrize('name',['gate-k2-D4','gate-k3-D4','gate-m3-k2-D5','gate-m4-k2-D5','gate-m5-k3-D5'])
def test_degree_priority_schedule_keeps_checked_reading(tmp_path,priority_producer,rust_checker,name):
    import subprocess,json
    inp=HERE/'inputs'/(name+'.txt');ledger=tmp_path/'ledger.gz';basis=tmp_path/'basis.gz'
    result=json.loads(subprocess.run([str(priority_producer),str(inp),str(ledger),str(basis)],capture_output=True,text=True,check=True).stdout)
    known=json.loads((HERE/'results'/(name+'.json')).read_text())
    assert result['status'] in ['closed','unit-witness'] and result['contains_one']==known['contains_one']
    cmd=[str(rust_checker),str(inp),str(ledger),str(basis)]
    if not result['contains_one']:
        assert result['derived_rank']==known['derived_rank'];cmd.append('closed')
    assert subprocess.run(cmd,capture_output=True).returncode==0


def test_priority_schedule_rejects_fifo_cursor(tmp_path,monkeypatch):
    import run_heap as runner
    import sys
    monkeypatch.setattr(sys,'argv',['run_heap.py','unused','tag','--engine','priority','--resume','unused','--resume-cursor'])
    with pytest.raises(ValueError,match='FIFO restart engine'):runner.main()


@pytest.mark.parametrize('compiled',[False,True])
def test_compact_falsification_controls(tmp_path,compiled):
    import subprocess,sys,json
    out=tmp_path/'falsification.json'
    cmd=[sys.executable,str(HERE/'compact_falsify.py'),'--trials','64','--out',str(out)]
    if compiled:cmd.append('--compiled')
    subprocess.run(cmd,capture_output=True,text=True,check=True)
    assert json.loads(out.read_text())['all_agree']


@pytest.mark.parametrize('backend',['python','rust'])
@pytest.mark.parametrize('corrupt',[None,'row','future','degree','digest'])
def test_compact_membership_replay_refuses_invalid_proofs(tmp_path,corrupt,backend):
    import json
    from compact_candidate import Ring,Compact
    from compact_replay import verify
    if backend=='rust':
        from compact_replay_rust import verify
    # A high chain bit squares to 1; the auxiliary variable exercises the
    # independent bound on reducer multipliers.
    inp=tmp_path/'input';inp.write_text('17 1 2 1 2 10 100 10000\n1\n1 2 1\n')
    C=Compact(Ring(17,1,2,1,2));assert C.run([{2:1}])['contains_one']
    proof={'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'trace':C.trace}
    if corrupt=='row':proof['trace'][-1]['row'][0]=(proof['trace'][-1]['row'][0][0],2)
    elif corrupt=='future':proof['trace'][-1]['source']=['carry',999,1]
    elif corrupt=='degree':proof['trace'][-1]['reductions'].append([0,3<<2,1])
    elif corrupt=='digest':proof['input_sha256']='wrong'
    path=tmp_path/'proof';path.write_text(json.dumps(proof))
    if corrupt:
        with pytest.raises(ValueError):verify(inp,path)
    else:
        result=verify(inp,path);assert result['contains_one'] and not result['closed_check']


@pytest.mark.parametrize('implementation',['prototype','independent'])
@pytest.mark.parametrize('kind',['closed-control','missing-carry','missing-pair'])
def test_compact_condition_falsifier(tmp_path,monkeypatch,kind,implementation):
    import json,shutil
    import importlib
    checker=importlib.import_module('compact_conditions' if implementation=='prototype' else 'compact_conditions_independent')
    from compact_falsify import text_input
    from compact_candidate import Ring
    home=tmp_path/'checker';home.mkdir();(home/'results').mkdir();shutil.copyfile(HERE/'compact_candidate.py',home/'compact_candidate.py');monkeypatch.setattr(checker,'HERE',home)
    d=tmp_path/'case';d.mkdir()
    if kind=='missing-pair':
        R=Ring(17,1,1,3,3);x,y,z=[1<<(1+6*j) for j in range(3)];gens=[{x+y:1,0:1},{x+z:1,0:1}]
    else:
        R=Ring(17,1,2,0,2);gens=[{2:1,0:16}] if kind=='closed-control' else [{2:1}]
    inp=d/'input.txt';inp.write_text(text_input(R,gens));(d/'candidate.json').write_text(json.dumps({'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'polynomials':[sorted(g.items()) for g in gens]}))
    result=checker.check(d)
    if kind=='closed-control':assert result['status']=='finite-conditions-passed'
    else:
        assert result['status']=='conditions-failed'
        assert result['missing']['kind']==('carries' if kind=='missing-carry' else 'pairs')


@pytest.fixture(scope='module')
def compiled_compact(tmp_path_factory):
    import subprocess
    binary=tmp_path_factory.mktemp('compiled-compact')/'producer'
    subprocess.run(['rustc','-O',str(HERE/'compact_producer.rs'),'-o',str(binary)],check=True)
    return binary


def test_compiled_compact_retains_valid_capped_prefix(tmp_path,compiled_compact):
    import subprocess,json
    from compact_replay import verify
    inp=tmp_path/'input';tokens=(SRC/'input.txt').read_text().split();tokens[6]='2';inp.write_text(' '.join(tokens)+'\n');proof=tmp_path/'proof'
    result=json.loads(subprocess.run([str(compiled_compact),str(inp),str(proof)],capture_output=True,text=True,check=True).stdout)
    assert result['status']=='capped' and result['compact_rows']==2
    data=json.loads(proof.read_text());data['input_sha256']=hashlib.sha256(inp.read_bytes()).hexdigest();proof.write_text(json.dumps(data))
    checked=verify(inp,proof);assert checked['verified_rows']==2 and not checked['closed_check']


def test_compiled_compact_high_chain_wrap_witness(tmp_path,compiled_compact):
    import subprocess,json
    from compact_replay import verify
    inp=tmp_path/'input';inp.write_text(f'7340033 2 18 0 2 10 100 10000\n1\n1 {1<<35} 1\n');proof=tmp_path/'proof'
    result=json.loads(subprocess.run([str(compiled_compact),str(inp),str(proof)],capture_output=True,text=True,check=True).stdout)
    assert result['status']=='unit-witness' and result['contains_one']
    data=json.loads(proof.read_text());data['input_sha256']=hashlib.sha256(inp.read_bytes()).hexdigest();proof.write_text(json.dumps(data))
    assert verify(inp,proof)['contains_one']
    from compact_replay_rust import verify as rust_verify
    assert rust_verify(inp,proof)['contains_one']


def test_compact_runner_refuses_existing_tag(tmp_path,monkeypatch):
    import run_compact as runner
    import sys
    from types import SimpleNamespace
    existing=tmp_path/'runs'/'existing';existing.mkdir(parents=True);(existing/'result.json').write_text('keep this')
    monkeypatch.setattr(runner,'HERE',tmp_path);monkeypatch.setattr(runner.shutil,'disk_usage',lambda path:SimpleNamespace(free=4*2**30))
    monkeypatch.setattr(sys,'argv',['run_compact.py','unused','existing','--engine','rust'])
    with pytest.raises(FileExistsError):runner.main()
    assert (existing/'result.json').read_text()=='keep this'


def test_compiled_compact_three_auxiliary_variables(tmp_path,compiled_compact):
    import subprocess,json
    from compact_candidate import Ring,Compact
    from compact_falsify import text_input,direct
    from compact_replay import verify
    R=Ring(17,1,2,3,3);u,v,w=[1<<(2+6*j) for j in range(3)]
    gens=[{2*w:1,u:16},{w+u:1,v:16},{2*v:1,0:16}]
    inp=tmp_path/'input';inp.write_text(text_input(R,gens));proof=tmp_path/'proof'
    out=json.loads(subprocess.run([str(compiled_compact),str(inp),str(proof)],capture_output=True,text=True,check=True).stdout)
    C=Compact(R);py=C.run(gens);assert not C.capped and out['status']=='provisional-complete'
    data=json.loads(proof.read_text());assert data['trace']==json.loads(json.dumps(C.trace))
    data['input_sha256']=hashlib.sha256(inp.read_bytes()).hexdigest();proof.write_text(json.dumps(data));assert not verify(inp,proof)['contains_one']
    binary=tmp_path/'reference';subprocess.run(['rustc','-O',str(HERE/'sparse_producer.rs'),'-o',str(binary)],check=True)
    reference=direct(binary,tmp_path,R,gens);assert reference['status']=='closed' and reference['derived_rank']==C.candidate_rank()


@pytest.mark.parametrize('bad',[True,'0 999'])
def test_compact_text_conversion_rejects_token_injection(tmp_path,bad):
    import json
    from compact_replay_rust import convert
    inp=tmp_path/'input';inp.write_bytes((SRC/'input.txt').read_bytes());proof=tmp_path/'proof'
    proof.write_text(json.dumps({'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'trace':[{'source':['input',bad],'inverse':1,'reductions':[],'row':[[0,1]]}]}))
    with pytest.raises(ValueError,match='unsigned integer'):convert(inp,proof,tmp_path/'converted')


def test_compact_rust_replay_binds_converted_bytes(tmp_path,monkeypatch):
    import json
    import compact_replay_rust as checker
    from compact_candidate import Ring,Compact
    inp=tmp_path/'input';inp.write_text('17 1 2 0 2 10 100 10000\n1\n1 2 1\n');C=Compact(Ring(17,1,2,0,2));C.run([{2:1}]);proof=tmp_path/'proof'
    proof.write_text(json.dumps({'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'trace':C.trace}))
    def mutate(input_path,text,binary):
        text.write_text('changed');return {'status':'passed','verified_rows':2,'contains_one':True}
    monkeypatch.setattr(checker,'run_backend',mutate)
    with pytest.raises(ValueError,match='converted text changed'):checker.verify(inp,proof)


@pytest.mark.parametrize('frontier',[False,True])
def test_coprime_pruning_readings_against_direct_closure(tmp_path,frontier):
    import subprocess,sys,json
    out=tmp_path/'pruned.json'
    subprocess.run([sys.executable,str(HERE/'compact_falsify.py'),'--trials','64','--compiled','--prune-coprime','--out',str(out)]+(['--frontier'] if frontier else []),capture_output=True,text=True,check=True)
    assert json.loads(out.read_text())['all_agree']


def test_coprime_product_pair_is_skipped(tmp_path,compiled_compact):
    import subprocess,json
    from compact_replay import verify
    inp=tmp_path/'input';inp.write_text('17 1 1 2 2 10 100 10000\n2\n2 0 16 2 1\n2 0 16 128 1\n');proof=tmp_path/'proof'
    result=json.loads(subprocess.run([str(compiled_compact),str(inp),str(proof),'2000000','coprime'],capture_output=True,text=True,check=True).stdout)
    assert result['status']=='provisional-complete' and result['coprime_pairs_skipped']==1 and result['pairs_processed']==0
    data=json.loads(proof.read_text());data['input_sha256']=hashlib.sha256(inp.read_bytes()).hexdigest();proof.write_text(json.dumps(data))
    assert verify(inp,proof)['verified_rows']==2
