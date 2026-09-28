#!/usr/bin/env python3
"""Download explicitly selected, versioned references and verify upstream checksums."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import urllib.request


def digest_file(path,algorithm):
    h=hashlib.new(algorithm)
    with Path(path).open('rb') as fh:
        for block in iter(lambda:fh.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def download(item,out,opener=urllib.request.urlopen):
    name=item['filename']
    if Path(name).name!=name or name in ('.','..'):raise ValueError('Invalid reference filename')
    if not item['url'].startswith('https://'):raise ValueError('Reference URL must use HTTPS')
    algorithm=item['checksum_algorithm'];expected=item['checksum'].lower()
    if algorithm not in ('md5','sha256'):raise ValueError('Unsupported checksum')
    if len(expected)!={'md5':32,'sha256':64}[algorithm] or any(c not in '0123456789abcdef' for c in expected):raise ValueError('Invalid checksum')
    out=Path(out);out.mkdir(parents=True,exist_ok=True);dest=out/name
    if dest.exists():
        if digest_file(dest,algorithm)!=expected:raise ValueError(f'Existing file checksum mismatch: {dest}; move it aside before retrying')
        return dest
    fd,tmp=tempfile.mkstemp(prefix=name+'.',suffix='.part',dir=out)
    try:
        with os.fdopen(fd,'wb') as target,opener(item['url'],timeout=120) as response:
            for block in iter(lambda:response.read(1024*1024),b''):target.write(block)
        if digest_file(tmp,algorithm)!=expected:raise ValueError(f'Download checksum mismatch: {name}')
        # Hard-link creation is atomic and refuses an existing target, including a concurrent download.
        os.link(tmp,dest)
    finally:
        Path(tmp).unlink(missing_ok=True)
    return dest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,default=Path(__file__).with_name('gencode_v48.json'))
    p.add_argument('--out',type=Path,default=Path(__file__).with_name('downloads'))
    choices=p.add_mutually_exclusive_group(required=True)
    choices.add_argument('--list',action='store_true');choices.add_argument('--id',action='append');choices.add_argument('--all',action='store_true')
    args=p.parse_args();manifest=json.loads(args.manifest.read_text());items=manifest['files']
    if args.list:
        for item in items:print(f"{item['id']}\t{item['filename']}\t{item['regions']}")
        return
    wanted=set(args.id or [i['id'] for i in items]);known={i['id'] for i in items}
    if wanted-known:p.error('Unknown reference IDs: '+', '.join(sorted(wanted-known)))
    for item in items:
        if item['id'] in wanted:
            print('Downloading/verifying '+item['filename'],flush=True)
            print(download(item,args.out),flush=True)


if __name__=='__main__':main()
