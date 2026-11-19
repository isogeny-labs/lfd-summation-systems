#!/usr/bin/env python3
"""Losslessly archive this task's terminal raw files after verifying their hashes."""
import gzip,hashlib,json,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent


def digest(stream):
    h=hashlib.sha256()
    for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def main():
    report=[]
    for receipt in (HERE/'runs').glob('*/result.json'):
        record=json.loads(receipt.read_text());archives={}
        # A terminal receipt is required; live workers have no result.json.
        for name in ['derivations.txt','basis.txt']:
            src=receipt.parent/name;dest=receipt.parent/(name+'.gz')
            expected=record['artifact_sha256'][name]
            if src.exists():
                with src.open('rb') as f:assert digest(f)==expected
                if not dest.exists():
                    temporary=dest.with_suffix('.gz.tmp')
                    with src.open('rb') as inp,temporary.open('wb') as raw:
                        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,compresslevel=1,mtime=0) as out:shutil.copyfileobj(inp,out,1<<20)
                    temporary.replace(dest)
            with gzip.open(dest,'rb') as decoded:assert digest(decoded)==expected,'archive differs from raw receipt'
            with dest.open('rb') as encoded:compressed=digest(encoded)
            archives[name]=dict(file=dest.name,decoded_sha256=expected,archive_sha256=compressed,bytes=dest.stat().st_size)
            if src.exists():src.unlink()  # Reversible: the exact bytes have just been verified in the archive.
        record['archives']=archives
        record['archive_tool_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        receipt.write_text(json.dumps(record,indent=2)+'\n')
        (HERE/'results'/f'{receipt.parent.name}.json').write_bytes(receipt.read_bytes())
        (HERE/'inputs'/f'{receipt.parent.name}.txt').write_bytes((receipt.parent/'input.txt').read_bytes())
        report.append(dict(run=receipt.parent.name,status=record['status'],archive_bytes=sum(a['bytes'] for a in archives.values())))
    print(json.dumps(report))
if __name__=='__main__':main()
