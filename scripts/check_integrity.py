"""Verify exact evidence bytes, including the complete manifest file set."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def verify(proof_dir=None):
    proof = proof_dir or ROOT / 'proof'
    manifest = json.loads((proof / 'SHA256.json').read_text(encoding='utf-8'))
    actual = {p.relative_to(proof).as_posix() for p in proof.rglob('*')
              if p.is_file() and p != proof / 'SHA256.json'}
    if actual != set(manifest):
        raise RuntimeError(f'Manifest file-set mismatch: missing={set(manifest)-actual}, extra={actual-set(manifest)}')
    for name, expected in manifest.items():
        if hashlib.sha256((proof / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f'SHA256 mismatch: {name}')
    return {'status': 'PASS', 'proof_files_verified': len(manifest)}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
