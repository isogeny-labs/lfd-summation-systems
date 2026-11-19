"""Bind compact JSON proofs and replay a deterministic text conversion in Rust."""
import argparse,hashlib,json,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
CONVERTER_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
BACKEND_SOURCE_BYTES=(HERE/'compact_replay.rs').read_bytes()

def convert(inp,proof,dest):
    inp=Path(inp);proof=Path(proof);body=proof.read_bytes();data=json.loads(body);input_sha=hashlib.sha256(inp.read_bytes()).hexdigest()
    if data.get('input_sha256')!=input_sha:raise ValueError('input digest')
    trace=data['trace']
    if not isinstance(trace,list):raise ValueError('trace shape')
    def number(value):
        if type(value) is not int or not 0<=value<2**64:raise ValueError('unsigned integer required')
        return str(value)
    with Path(dest).open('w') as stream:
        stream.write(str(len(trace))+'\n')
        for ident,node in enumerate(trace):
            source=node['source'];name=source[0]
            if name=='input' and len(source)==2:kind=0;i=source[1];j=0
            elif name=='carry' and len(source)==3:kind=1;i,j=source[1:]
            elif name=='pair' and len(source)==3:kind=2;i,j=source[1:]
            else:raise ValueError('source kind/shape')
            reductions=node['reductions'];rows=node['row']
            values=[ident,kind,i,j,node['inverse'],len(reductions)]
            for item in reductions:
                if len(item)!=3:raise ValueError('reduction shape')
                values.extend(item)
            values.append(len(rows))
            for pair in rows:
                if len(pair)!=2:raise ValueError('row shape')
                values.extend(pair)
            stream.write(' '.join(map(number,values))+'\n')
    converted=Path(dest).read_bytes()
    return {'input_sha256':input_sha,'proof_sha256':hashlib.sha256(body).hexdigest(),'text_sha256':hashlib.sha256(converted).hexdigest(),'text_bytes':len(converted),'nodes':len(trace)}

def run_backend(inp,text,binary):
    result=subprocess.run([str(binary),str(inp),str(text)],capture_output=True,text=True,timeout=1800)
    if result.returncode:raise ValueError('compact Rust replay failed: '+result.stderr[-2000:])
    return json.loads(result.stdout)

def verify(inp,proof):
    inp=Path(inp);proof=Path(proof);source=HERE/'compact_replay.rs'
    with tempfile.TemporaryDirectory(prefix='compact-check-',dir=proof.parent) as tmp:
        tmp=Path(tmp);text=tmp/'proof.txt';binding=convert(inp,proof,text);binary=tmp/'checker'
        source_bytes=BACKEND_SOURCE_BYTES;source_sha=hashlib.sha256(source_bytes).hexdigest();snapshot=tmp/'compact_replay.rs';snapshot.write_bytes(source_bytes);subprocess.run(['rustc','-O',str(snapshot),'-o',str(binary)],check=True)
        binary_sha=hashlib.sha256(binary.read_bytes()).hexdigest();result=run_backend(inp,text,binary)
        if hashlib.sha256(text.read_bytes()).hexdigest()!=binding['text_sha256']:raise ValueError('converted text changed during replay')
    if hashlib.sha256(proof.read_bytes()).hexdigest()!=binding['proof_sha256']:raise ValueError('proof changed during replay')
    if hashlib.sha256(inp.read_bytes()).hexdigest()!=binding['input_sha256']:raise ValueError('input changed during replay')
    if result['verified_rows']!=binding['nodes']:raise ValueError('row count mismatch')
    result.update(binding,verifier_source_sha256=source_sha,verifier_binary_sha256=binary_sha,converter_source_sha256=CONVERTER_SHA256)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('proof',type=Path);ap.add_argument('--out',type=Path);a=ap.parse_args();result=verify(a.input,a.proof)
    text=json.dumps(result,indent=2)+'\n'
    if a.out:a.out.write_text(text)
    print(text)
