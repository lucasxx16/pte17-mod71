from pathlib import Path
import json,hashlib
here=Path(__file__).resolve().parent
source=here.parent.parent/'D'
hashes={}
for name in ['check_multiplicity4.py','IntegerLatticeEnumerator.cs','lattice_certificate_input.json','lattice_certificate_input.txt','verify_lattice_certificate.py','run_integer_lattice.ps1']:
    a=(source/name).read_bytes();b=(here/name).read_bytes()
    assert a==b
    hashes[name]=hashlib.sha256(a).hexdigest()
original=json.loads((source/'multiplicity4_result.json').read_text())
replay=json.loads((here/'multiplicity4_result.json').read_text())
original_time=original.pop('elapsed_seconds');replay_time=replay.pop('elapsed_seconds')
assert original==replay
def read_lattice(path):
    return dict(line.split('=',1) for line in path.read_text().splitlines() if '=' in line)
original_lattice=read_lattice(source/'lattice_enumeration_result.txt')
replay_lattice=read_lattice(here/'lattice_enumeration_result.txt')
original_lattice_time=original_lattice.pop('elapsed_seconds');replay_lattice_time=replay_lattice.pop('elapsed_seconds')
assert original_lattice==replay_lattice
assert not (here/'lattice_candidate_solutions.txt').read_text()
out={'status':'PASS','all_source_copies_identical':True,'m4_all_semantic_fields_identical':True,
     'lattice_all_semantic_fields_identical':True,'candidate_file_empty':True,
     'original_m4_seconds':original_time,'replay_m4_seconds':replay_time,
     'original_lattice_seconds':original_lattice_time,'replay_lattice_seconds':replay_lattice_time,
     'hashes':hashes}
(here/'completed_runs_comparison.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
