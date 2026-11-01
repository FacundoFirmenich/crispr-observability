"""Reconstruct actual coding catalogs from original FASTQ and frozen allele tables.

No alignment engine, model fitting, generated sequence, or synthetic biological
record is introduced. Frozen M8 scientific outputs are read and never modified.
"""
from pathlib import Path
from collections import Counter, defaultdict
from fractions import Fraction
import csv, gzip, hashlib, io, json, re, sys, zipfile

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R/'src'))
from observed_contracts.core import peptide, NATIVE, PROPERTIES, TARGETS, factorization, grouped_withdrawal

OPS = ['target', 'mutant_projection', 'programmed_projection']
LEFT, RIGHT = 'TGGATTATGCCT', 'ATGAATATAGATACAGAAGC'
COMP = str.maketrans('ACGTN', 'TGCAN')

def get(name): return json.loads((R/'inputs'/name).read_text())
def put(name, obj): (R/'results'/name).write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rc(s): return s.translate(COMP)[::-1]

def extract_target(seq):
    hits = []
    for s in [seq, rc(seq)]:
        for a in re.finditer('(?='+LEFT+')', s):
            for b in re.finditer('(?='+RIGHT+')', s):
                lo = a.start()+len(LEFT)
                if lo <= b.start() <= lo+350:
                    hits.append(s[lo:b.start()]+RIGHT[:2])
    return hits[0] if len(hits) == 1 and not set(hits[0])-set('ACGT') else None

def project(reference_alignment, observed_alignment, ref):
    assert len(reference_alignment) == len(observed_alignment)
    assert reference_alignment.replace('-', '') == ref['sequence']
    lo, hi = ref['core_start_0'], ref['coding_end_exclusive_0']
    pos, pieces, flanks = -1, [], {}
    for a, b in zip(reference_alignment, observed_alignment):
        if a != '-':
            pos += 1
            if pos in [lo-1, hi]: flanks[pos] = b
            keep = lo <= pos < hi
        else:
            keep = lo <= pos+1 < hi
        if keep and b != '-': pieces.append(b)
    s = ''.join(pieces)
    if not s or set(s)-set('ACGT') or any(flanks.get(k) not in ('A','C','G','T') for k in [lo-1,hi]):
        return None
    return s

def allele_map(path, ref):
    values, abundances = {}, Counter()
    with zipfile.ZipFile(path) as z:
        with z.open(z.namelist()[0]) as f:
            for row in csv.DictReader(io.TextIOWrapper(f), delimiter='\t'):
                s = row['Aligned_Sequence'].replace('-', '')
                value = project(row['Reference_Sequence'], row['Aligned_Sequence'], ref)
                assert s not in values or values[s] == value
                values[s] = value
                abundances[s] += int(row['#Reads'])
    return values, abundances

def add_catalog(counters, examples, op, run, coding, read_id):
    if coding is None: return
    for key in [run, 'ALL']:
        counters[op][key][coding] += 1
        examples[op][key].setdefault(coding, {'run': run, 'read_id': read_id, 'operator': op, 'coding': coding})

def main():
    for f in json.loads((R/'INPUT_MANIFEST.json').read_text()):
        p = R/f['path']
        assert p.stat().st_size == f['bytes'] and sha(p) == f['sha256'], f['path']
    design = get('RUN_DESIGN.json')
    refs = [get('REFERENCE_MUTANT.json'), get('REFERENCE_PROGRAMMED.json')]
    exact = refs[0]['exact_programmed_core']+'AT'
    counters = {o:defaultdict(Counter) for o in OPS+['witness']}
    examples = {o:defaultdict(dict) for o in OPS+['witness']}
    witness_counts = {}
    correspondence, details, raw_receipts = [], [], []
    all_n = 0
    export = R/'results/OBSERVED_RECORD_CODINGS.jsonl.gz'
    with export.open('wb') as raw_out, gzip.GzipFile(filename='',mode='wb',fileobj=raw_out,mtime=0) as gz, io.TextIOWrapper(gz,encoding='utf-8',newline='\n') as out:
        for run in design['runs']:
            rid = run['run_accession']
            maps = [allele_map(R/'inputs'/kind/f'{rid}.zip', ref) for kind, ref in zip(['mutant_alleles','programmed_alleles'], refs)]
            observed_abundances = [Counter(),Counter()]
            cnt, wc = Counter(), {'native_event':Counter(), 'exact_programmed_event':Counter()}
            digest, n = hashlib.sha256(), 0
            with gzip.open(R/'inputs/raw'/f'{rid}.fastq.gz','rt') as f:
                while header := f.readline():
                    s, plus, q = f.readline().strip(), f.readline(), f.readline().strip()
                    assert plus.startswith('+') and len(s) == len(q)
                    name = header.split()[0][1:]
                    n += 1
                    digest.update((s+'\n'+q+'\n').encode())
                    quality = sum(ord(c)-33 for c in q) >= 30*len(q)
                    codes = [None,None,None]
                    if quality:
                        codes[0] = extract_target(s)
                        for i,(mp,_) in enumerate(maps):
                            keys = {v for v in [s, rc(s)] if v in mp}
                            assert len(keys) <= 1
                            if keys:
                                key = next(iter(keys))
                                codes[i+1] = mp[key]
                                observed_abundances[i][key] += 1
                    cnt['raw'] += 1
                    if not quality: cnt['quality_unknown'] += 1
                    resolved = all(c is not None for c in codes)
                    common = resolved and len(set(codes)) == 1
                    if resolved:
                        cnt['all_three_resolved'] += 1
                        ps = [peptide(c) for c in codes]
                        ns = [p == NATIVE for p in ps]
                        es = [c == exact for c in codes]
                        same_peptide = all(p is not None for p in ps) and len(set(ps)) == 1
                        cnt['same_coding'] += int(common)
                        cnt['same_valid_peptide'] += int(same_peptide)
                        cnt['same_native_label'] += int(len(set(ns)) == 1)
                        cnt['same_exact_label'] += int(len(set(es)) == 1)
                        if not common:
                            cnt['different_coding'] += 1
                            cnt['different_coding_same_valid_peptide'] += int(same_peptide)
                            cnt['different_coding_same_native_label'] += int(len(set(ns)) == 1)
                            details.append({'run':rid,'read_id':name,'codings':dict(zip(OPS,codes)),'peptides':ps,'native_labels':ns,'exact_labels':es,'same_valid_peptide':same_peptide})
                    if common:
                        cnt['witness'] += 1
                        add_catalog(counters, examples, 'witness', rid, codes[0], name)
                        for event, truth in [('native_event',peptide(codes[0]) == NATIVE),('exact_programmed_event',codes[0] == exact)]: wc[event]['E' if truth else 'NE'] += 1
                    else:
                        for event in wc: wc[event]['U'] += 1
                        cnt['missing_correspondence' if quality and not resolved else 'coding_disagreement' if quality else 'quality_failure'] += 1
                    for o,c in zip(OPS,codes): add_catalog(counters, examples, o, rid, c, name)
                    out.write(json.dumps({'run':rid,'read_id':name,'raw_sequence_sha256':hashlib.sha256(s.encode()).hexdigest(),'common_Q30':quality,'codings':dict(zip(OPS,codes)),'witness':common},sort_keys=True,separators=(',',':'))+'\n')
            assert n == run['raw_reads'] and digest.hexdigest() == run['sequence_quality_stream_sha256']
            assert all(a == expected for a,(_,expected) in zip(observed_abundances,maps)), rid
            for event,c in wc.items(): witness_counts[rid,event] = {'N':n,'E':c['E'],'NE':c['NE'],'U':c['U']}
            raw_receipts.append({'run':rid,'records':n,'sequence_quality_stream_sha256':digest.hexdigest(),'matches_frozen_stream':True,'both_allele_abundances_preserved':True})
            correspondence.append({'run':rid,'preparation':run['source_replica'],'condition':run['condition'],**dict(cnt)})
            all_n += n
    assert all_n == 177006
    frozen = get('M8_WITNESS_COUNTS.json')
    for r in frozen:
        event = 'native_event' if r['event'].startswith('native') else 'exact_programmed_event'
        assert witness_counts[r['run_accession'],event] == {k:r[k] for k in ['N','E','NE','U']}
    catalogs = {op:{run:[{'coding':s,'count':n,'example':examples[op][run][s]} for s,n in sorted(counts.items())] for run,counts in byrun.items()} for op,byrun in counters.items()}
    put('CODING_CATALOGS.json', {'exact_programmed_coding':exact,'translation':'Fixed standard table and declared M8 local frame; no contextual Sec reinterpretation','catalogs':catalogs})
    atlas = []
    for op in OPS:
        for run,catalog in catalogs[op].items():
            for rep in PROPERTIES:
                for target in TARGETS:
                    atlas.append({'operator':op,'run':run,**factorization(catalog,rep,target,exact)})
    blocks = []
    for contrast in get('M8_WITNESS_CONTRASTS.json'):
        tr, cr = contrast['treated'], contrast['comparator']
        event = 'native_event' if contrast['event'].startswith('native') else 'exact_programmed_event'
        c = grouped_withdrawal(catalogs['witness'][tr],catalogs['witness'][cr],witness_counts[tr,event],witness_counts[cr,event],event,exact)
        assert Fraction(c['original_lower_exact']) == Fraction(contrast['L_exact'])
        blocks.append({'treated':tr,'control':cr,'preparation':contrast['preparation'],**c})
    put('CORRESPONDENCE_BY_RUN.json', correspondence)
    put('CODING_DISAGREEMENT_RECORDS.json', details)
    put('FACTORIZATION_ATLAS.json', atlas)
    put('SEQUENCE_FAMILY_WITHDRAWAL.json', blocks)
    summary = {'status':'PASS_RECONSTRUCTION_AND_FINITE_AUDIT','records':all_n,'raw_streams':len(raw_receipts),'counts_match_M8':len(frozen),'per_reference_abundances_checked':2*len(raw_receipts),'correspondence_totals':dict(sum((Counter({k:v for k,v in r.items() if k not in ['run','preparation','condition'] and type(v)==int}) for r in correspondence),Counter())),'factorization_cases':len(atlas),'cases_with_counterexamples':sum(not a['factors_on_observed_catalog'] for a in atlas),'whole_sequence_family_certificates':len(blocks),'minimum_families_distribution':dict(Counter(b['minimum_sequence_families'] for b in blocks)),'original_raw_streams':raw_receipts,'no_new_biological_observations':True,'external_alignment_engine_rerun':False,'external_genealogy':False}
    put('AUDIT_SUMMARY.json', summary)
    print(json.dumps(summary,ensure_ascii=False))

if __name__ == '__main__': main()
