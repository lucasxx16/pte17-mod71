"""Reproduce in an isolated copy, leaving the committed evidence unchanged."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

from check_integrity import ROOT, verify


def read_result(path):
    result = json.loads(path.read_text(encoding='utf-8'))
    result.pop('elapsed_seconds', None)
    return result


def read_lattice(path):
    result = dict(line.split('=', 1) for line in path.read_text().splitlines() if '=' in line)
    result.pop('elapsed_seconds', None)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true', help='Certificates and controls only; no exhaustive searches')
    parser.add_argument('--powershell', help='PowerShell 7+ executable; defaults to pwsh')
    args = parser.parse_args()
    integrity = verify()
    shell = args.powershell or shutil.which('pwsh')
    if not shell:
        parser.error('PowerShell 7+ with C# Add-Type is required; use --powershell PATH')
    version_check = subprocess.run(
        [str(shell), '-NoProfile', '-Command', 'if ($PSVersionTable.PSVersion.Major -lt 7) { exit 1 }'],
        capture_output=True,
    )
    if version_check.returncode:
        parser.error('PowerShell 7+ is required for the checked C# audit; Windows PowerShell 5.1 is insufficient')
    run_name = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8]
    work = ROOT / '.runs' / run_name
    proof = work / 'proof'
    shutil.copytree(ROOT / 'proof', proof)
    logs = work / 'logs'
    logs.mkdir()
    summary = {'status': 'RUNNING', 'full_enumerations': not args.quick,
               'integrity': integrity, 'steps': [], 'python': sys.version}
    summary_path = work / 'run_summary.json'

    def save():
        summary_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')

    def run(label, folder, filename, powershell=False):
        if powershell:
            command = [str(shell), '-NoProfile', '-File', filename]
        else:
            command = [sys.executable, '-E', filename]
        print(f'Running {label}', flush=True)
        started = time.monotonic()
        log = logs / (label + '.log')
        with log.open('w', encoding='utf-8') as stream:
            completed = subprocess.run(command, cwd=folder, stdout=stream, stderr=subprocess.STDOUT)
        record = {'step': label, 'exit_code': completed.returncode,
                  'elapsed_seconds': time.monotonic() - started,
                  'log': log.relative_to(work).as_posix()}
        summary['steps'].append(record)
        save()
        if completed.returncode:
            raise RuntimeError(f'{label} failed; inspect {log}')
        print(f'PASS {label} ({record["elapsed_seconds"]:.1f}s)', flush=True)

    d = proof / 'round1' / 'D'
    c = proof / 'round1' / 'C' / 'audit_D'
    e = proof / 'round1' / 'E' / 'audit_D'
    print(f'Run directory: {work}', flush=True)
    save()
    try:
        run('reference_controls', proof / 'reference', 'verify_local.py')
        run('D_exact_certificate', d, 'verify_lattice_certificate.py')
        run('C_exact_certificate', c, 'check_certificate_independent.py')
        run('E_exact_certificate', e, 'independent_certificate_audit.py')
        run('C_positive_path', c, 'positive_path_control.py')
        run('E_all_isqrt_inputs', e, 'run_isqrt_audit.ps1', powershell=True)
        if not args.quick:
            run('D_full_multiplicity4', d, 'check_multiplicity4.py')
            run('E_different_full_multiplicity4', e, 'independent_m4_audit.py')
            run('C_screen_samples', c, 'check_m4_screen_samples.py')
            run('C_all_survivors_Hasse', c, 'check_m4_survivors_independent.py')
            run('D_full_integer_lattice', d, 'run_integer_lattice.ps1', powershell=True)
            original = ROOT / 'proof' / 'round1' / 'D'
            if read_result(d / 'multiplicity4_result.json') != read_result(original / 'multiplicity4_result.json'):
                raise RuntimeError('Full multiplicity-4 output differs from the saved evidence')
            if read_lattice(d / 'lattice_enumeration_result.txt') != read_lattice(original / 'lattice_enumeration_result.txt'):
                raise RuntimeError('Full lattice output differs from the saved evidence')
            if (d / 'lattice_candidate_solutions.txt').read_bytes():
                raise RuntimeError('Unexpected lattice candidate')
            run('E_compare_saved_replay_with_fresh_D', e, 'compare_completed_runs.py')
            summary['full_semantic_outputs_match_saved_evidence'] = True
        summary['committed_evidence_unchanged'] = verify()
        summary['status'] = 'PASS'
    except Exception as exc:
        summary['status'] = 'FAIL'
        summary['error'] = str(exc)
        raise
    finally:
        save()
    print(json.dumps({'status': summary['status'], 'full_enumerations': not args.quick,
                      'steps_passed': len(summary['steps']), 'summary': str(summary_path)}, indent=2), flush=True)


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Do not disable assertions with -O')
    main()
