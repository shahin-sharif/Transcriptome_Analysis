import hashlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location('refs',Path(__file__).resolve().parents[1]/'references'/'DownloadReferences.py')
refs=importlib.util.module_from_spec(spec);spec.loader.exec_module(refs)

class DownloadTests(unittest.TestCase):
    def test_verified_atomic_download_and_existing_file(self):
        data=b'>synthetic\nACGT\n'
        item=dict(filename='synthetic.fa',url='https://example.invalid/synthetic.fa',checksum_algorithm='sha256',checksum=hashlib.sha256(data).hexdigest())
        with tempfile.TemporaryDirectory() as tmp:
            p=refs.download(item,tmp,opener=lambda *a,**k:io.BytesIO(data))
            self.assertEqual(p.read_bytes(),data)
            refs.download(item,tmp,opener=lambda *a,**k:self.fail('Existing verified file must not be downloaded'))
            p.write_bytes(b'wrong')
            with self.assertRaises(ValueError):refs.download(item,tmp)
            self.assertEqual(p.read_bytes(),b'wrong')
    def test_failed_checksum_leaves_no_partial_file(self):
        item=dict(filename='test.fa',url='https://example.invalid/test',checksum_algorithm='md5',checksum='0'*32)
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):refs.download(item,tmp,opener=lambda *a,**k:io.BytesIO(b'bad'))
            self.assertEqual(list(Path(tmp).iterdir()),[])

if __name__=='__main__':unittest.main()
