"""Export complete actual-record ledgers and reconcile their readback to MP06.

These are derived audit tables, never repaired/native CRISPResso annotations.
Raw sequence, raw quality and every declared local coding are retained. Empty
coding means unresolved, not a negative label. No read is omitted for quality,
alignment failure or disagreement. Standard library; no engine execution.
"""
from pathlib import Path
from collections import Counter, defaultdict
import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import itertools
import json

ROOT = Path(__file__).resolve().parents[2]
FRONT = ROOT / 'extensions/correspondence_frontier'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


mech = module('integrity_mechanism', FRONT / 'analyze_mechanism.py')
ext = module('integrity_external', FRONT / 'analyze_external.py')
FAMILIES = {'F3': ('T', 'M_fw', 'P_fw'),
            'F4': ('T', 'M_fw', 'P_fw', 'PE_default'),
            'F9': ('T', 'M_fw', 'P_fw', 'M_rc', 'P_rc', 'MP_fw',
                   'MP_rc', 'PE_default', 'PE_zero')}
CFTR_OPS = FAMILIES['F9']
HEK3_OPS = ('T', 'M_fw', 'P_fw', 'M_rc', 'P_rc')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def write_rows(path, ops, rows):
    fields = ['run', 'read_id', 'raw_sequence', 'raw_quality', 'Q30'] + list(ops)
    # Zero timestamp and empty original filename make the gzip reproducible.
    with path.open('wb') as b, gzip.GzipFile(filename='', mode='wb', fileobj=b,
                                           mtime=0) as z:
        with io.TextIOWrapper(z, encoding='utf-8', newline='') as text:
            writer = csv.DictWriter(text, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def cftr_rows():
    refs = [mech.raw.load(ROOT / 'raw_replay/metadata' / n)
            for n in ('REFERENCE.json', 'REFERENCE_PROGRAMMED.json')]
    for rid, records in itertools.groupby(mech.cor.records(),
                                         key=lambda x: x[0]['run_accession']):
        maps, weights, seen = {}, {}, defaultdict(Counter)
        for mode in CFTR_OPS[3:]:
            folder = ((FRONT / 'inputs/PE_default') if mode == 'PE_default'
                      else FRONT / 'inputs/cftr/runs' / mode)
            maps[mode], _, weights[mode] = mech.read_mode(
                folder / ('CRISPResso_on_' + rid), not mode.endswith('_fw'), refs)
        original = iter(mech.raw.fastq(ROOT / 'raw_replay/inputs/fastq' /
                                      (rid + '.fastq.gz')))
        for _, name, sequence, quality, old in records:
            raw_name, raw_sequence, raw_quality = next(original)
            assert (name, sequence) == (raw_name, raw_sequence)
            key = min(sequence, mech.rc(sequence))
            codes = dict(zip(CFTR_OPS[:3], (old['target'],
                         old['mutant_projection'], old['programmed_projection'])))
            for mode, mapping in maps.items():
                value = mapping.get(key) if quality else None
                codes[mode] = value[0] if value is not None else None
                if value is not None:
                    seen[mode][key] += 1
            yield dict(run=rid, read_id=name, raw_sequence=sequence,
                       raw_quality=raw_quality, Q30=int(quality),
                       **{k: v if v is not None else '' for k, v in codes.items()})
        assert next(original, None) is None
        for mode in maps:
            assert seen[mode] == weights[mode], (rid, mode, 'allele abundance')


def hek3_rows():
    base = FRONT / 'inputs/hek3'
    spec = json.loads((base / 'RULE_FREEZE.json').read_text())
    for inp in json.loads((base / 'FILTERED_INPUTS.json').read_text()):
        rid = inp['run']
        maps, weights, seen = {}, {}, defaultdict(Counter)
        for mode in HEK3_OPS[1:]:
            ref = spec['reference'] if mode.startswith('M') else spec['programmed_reference']
            maps[mode], weights[mode], _ = ext.tablemap(
                base / 'runs' / mode / ('CRISPResso_on_' + rid) /
                'Alleles_frequency_table.zip', ref, spec['local_lo'],
                spec['local_hi'], mode.endswith('_rc'))
        n = qpass = 0
        for name, sequence, qualities in mech.raw.fastq(base / 'inputs' /
                                                        (rid + '.fastq.gz')):
            quality = sum(ord(c) - 33 for c in qualities) >= 30 * len(qualities)
            key = min(sequence, mech.rc(sequence))
            codes = {'T': ext.exact_target(sequence, spec['left_anchor'],
                     spec['right_anchor']) if quality else None}
            for mode, mapping in maps.items():
                codes[mode] = mapping.get(key) if quality else None
                if quality and key in mapping:
                    seen[mode][key] += 1
            n += 1
            qpass += quality
            yield dict(run=rid, read_id=name, raw_sequence=sequence,
                       raw_quality=qualities, Q30=int(quality),
                       **{k: v if v is not None else '' for k, v in codes.items()})
        assert (n, qpass) == (inp['N'], inp['Q30'])
        assert all(seen[m] == weights[m] for m in maps), rid


def readback(path, cohort):
    """Read actual saved bytes, rejoin every raw record and rebuild all counts.

    Classification below does not call the aggregate-analysis classification
    routines; codings come from the saved ledger, not the expected count tables.
    """
    counts, nraw, quality = defaultdict(Counter), Counter(), Counter()
    cftr = cohort == 'cftr'
    ops = CFTR_OPS if cftr else HEK3_OPS
    base = FRONT / 'inputs/hek3'
    spec = json.loads((base / 'RULE_FREEZE.json').read_text())
    exact = mech.raw.load(ROOT / 'raw_replay/metadata/REFERENCE.json')['exact_programmed_core'] + 'AT'
    previous, originals, ids = None, None, set()
    with gzip.open(path, 'rt', newline='') as f:
        for row in csv.DictReader(f):
            rid = row['run']
            if rid != previous:
                if originals is not None:
                    assert next(originals, None) is None, previous
                folder = ROOT / 'raw_replay/inputs/fastq' if cftr else base / 'inputs'
                originals = iter(mech.raw.fastq(folder / (rid + '.fastq.gz')))
                previous, ids = rid, set()
            assert row['read_id'] not in ids
            ids.add(row['read_id'])
            assert next(originals) == (row['read_id'], row['raw_sequence'], row['raw_quality'])
            q = sum(ord(c) - 33 for c in row['raw_quality']) >= 30 * len(row['raw_quality'])
            assert int(row['Q30']) == q
            assert q or not any(row[k] for k in ops)
            nraw[rid] += 1
            quality[rid] += q
            codes = {k: row[k] or None for k in ops}

            def sig(s):
                if cftr:
                    p = ''.join(mech.raw.CODONS[s[i:i+3]] for i in range(0, len(s), 3)) if len(s) % 3 == 0 else None
                    state = (int(p == 'GTIKENIIFGVSYD'), int(s == exact))
                    assert not state[1] or state[0]
                    return state
                return int(s == spec['programmed_local']), int(s == spec['unedited_local'])

            universe = {(0, 0), (1, 0), (1, 1)} if cftr else {(0, 0), (1, 0), (0, 1)}
            families = FAMILIES if cftr else {'F5': HEK3_OPS}
            events = ('native', 'exact_programmed', 'native_nonexact') if cftr else ('programmed_local', 'unedited_local')
            for family, required in families.items():
                strings = [codes[k] for k in required]
                complete = all(s is not None for s in strings)
                states = {sig(s) for s in strings} if complete else universe
                for cert in ('target_relative', 'string_witness'):
                    domain = states if cert == 'target_relative' or complete and len(set(strings)) == 1 else universe
                    for j, event in enumerate(events):
                        labels = {int(s == (1, 0)) if j == 2 else s[j] for s in domain}
                        label = 'U' if len(labels) > 1 else 'E' if labels == {1} else 'NE'
                        counts[rid, family, cert, event][label] += 1
            if not cftr:
                for op, coding in codes.items():
                    for j, event in enumerate(events):
                        label = 'U' if coding is None else 'E' if sig(coding)[j] else 'NE'
                        counts[rid, op, 'single_projection', event][label] += 1
    assert next(originals, None) is None
    expected = FRONT / ('expected/mechanism/FAMILY_CERTIFICATE_COUNTS.csv' if cftr
                        else 'expected/external/HEK3_EVENT_BOUNDS.csv')
    checked = 0
    with expected.open() as f:
        for row in csv.DictReader(f):
            key = row['run'], row['family' if cftr else 'operator'], row['certificate'], row['event']
            actual = counts.pop(key)
            assert nraw[row['run']] == int(row['N'])
            assert all(actual[k] == int(row[k]) for k in ('E', 'NE', 'U')), key
            checked += 1
    assert not counts
    assert sum(nraw.values()) == (177006 if cftr else 218930)
    return {'status': 'PASS', 'raw_records_rejoined_in_order': sum(nraw.values()),
            'Q30': sum(quality.values()), 'unique_record_ids_per_run': True,
            'all_frozen_count_rows_reconstructed': checked,
            'per_run_N': dict(nraw), 'bytes': path.stat().st_size,
            'sha256': digest(path), 'file': path.name}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--cohort', choices=('all', 'cftr', 'hek3'), default='all')
    args = p.parse_args()
    out = args.output.resolve()
    assert not out.exists() and not out.is_relative_to(ROOT)
    out.mkdir(parents=True)
    for name, record in json.loads((FRONT / 'INPUT_MANIFEST.json').read_text())['files'].items():
        path = FRONT / name
        assert path.stat().st_size == record['bytes'] and digest(path) == record['sha256'], name
    mech.raw.verify_inputs()
    result = {'status': 'RUNNING', 'scope': 'Complete derived record audit from authentic raw reads and complete native allele tables; no repaired native annotations, new engine run, molecular accuracy or human-independent reproduction.'}
    dump(out / 'LEDGER_VERIFICATION.json', result)
    for cohort in ('cftr', 'hek3') if args.cohort == 'all' else (args.cohort,):
        path = out / (cohort.upper() + '_COMPLETE_RECORD_LEDGER.csv.gz')
        write_rows(path, CFTR_OPS if cohort == 'cftr' else HEK3_OPS,
                   cftr_rows() if cohort == 'cftr' else hek3_rows())
        result[cohort] = readback(path, cohort)
        dump(out / 'LEDGER_VERIFICATION.json', result)
        print(cohort, json.dumps(result[cohort]), flush=True)
    result['status'] = 'PASS'
    dump(out / 'LEDGER_VERIFICATION.json', result)


if __name__ == '__main__':
    main()
