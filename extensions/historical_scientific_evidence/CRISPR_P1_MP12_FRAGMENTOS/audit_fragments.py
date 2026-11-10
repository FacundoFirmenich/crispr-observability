"""Audit actual FASTQ fragment linkage in the original studies and known deposits."""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from acquire import HEK3, CFTR, CFTR_EXTRA


HERE = Path(__file__).resolve().parent
SOURCES = HERE / 'sources'
MP07 = HERE.parent / 'CRISPR_P1_MP07_RECUPERACION'
I2 = MP07 / 'i2'
DUPE = ('SRR29536570', 'SRR29536571', 'SRR29536572')
CLUSTER = re.compile(rb'([^ :]+:\d+:[^ :]+:\d+:\d+:\d+:\d+)/(1|2)')


def digest(path, algorithm='sha256'):
    with path.open('rb') as f:
        return hashlib.file_digest(f, algorithm).hexdigest()


def raw_path(run):
    if run in HEK3:
        folder = I2 / 'extensions/correspondence_frontier/inputs/hek3/inputs'
    elif run in CFTR:
        folder = I2 / 'raw_replay/inputs/fastq'
    elif run in CFTR_EXTRA:
        folder = I2 / 'extensions/reproduction/geometry/additional_fastq'
    else:
        folder = SOURCES
    return folder / (run + '.fastq.gz')


def fastq(path):
    with gzip.open(path, 'rb') as stream:
        for lines in zip(stream, stream, stream, stream, strict=True):
            header, seq, plus, quality = (line.rstrip(b'\r\n') for line in lines)
            assert header.startswith(b'@') and plus.startswith(b'+') and seq and len(seq) == len(quality)
            read_id, *tokens = header[1:].split()
            assert re.fullmatch(re.split('[_.]',path.name,maxsplit=1)[0].encode()+rb'\.\d+', read_id)
            key = next((m.groups() for token in tokens if (m := CLUSTER.fullmatch(token))), None)
            yield read_id.decode(), key, seq, quality


def read_report(project):
    with (SOURCES / (project+'.tsv')).open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    lookup = {row['run_accession']: row for row in rows}
    assert len(rows) == len(lookup)
    return lookup


def source_evidence(project, rows, targets):
    matched = {r: rows[r] for r in targets}
    samples = {row['sample_accession'] for row in matched.values()}
    experiments = {row['experiment_accession'] for row in matched.values()}
    peers = [r for r,row in rows.items() if r not in matched and
             (row['sample_accession'] in samples or row['experiment_accession'] in experiments)]
    alias = re.compile(r'^Fig2b_HEK3_10_rep[123](?:[_.]|$)' if project == 'PRJNA565979'
                       else r'^Fig4a[_-]', re.I)
    return {'runs':len(rows), 'layout_platform':[dict(layout=a, instrument=b, count=n)
            for (a,b),n in sorted(Counter((r['library_layout'],r['instrument_model'])
                                          for r in rows.values()).items())],
            'other_same_sample_or_experiment':sorted(peers),
            'matching_aliases':{r:row['run_alias'] for r,row in sorted(rows.items()) if alias.search(row['run_alias'])},
            'paired_aliases':{r:row['run_alias'] for r,row in sorted(rows.items()) if row['library_layout']=='PAIRED'}}


def case_sets():
    sets = defaultdict(set)
    name = MP07/'complete_record_ledgers/CFTR_COMPLETE_RECORD_LEDGER.csv.gz'
    with gzip.open(name, 'rt', newline='') as f:
        for row in csv.DictReader(f):
            strings = (row['T'],row['M_fw'],row['P_fw'])
            if all(strings) and len(set(strings)) > 1:
                sets['CFTR_three_extractor_disagreement'].add(row['read_id'])
    assert len(sets['CFTR_three_extractor_disagreement']) == 83
    for locus in ('hek3','cftr'):
        with (MP07/'panel_MP07'/f'{locus}_witnesses.jsonl').open(encoding='utf-8') as f:
            sets[locus.upper()+'_path_panel'] = {row['read_id'] for row in map(json.loads,f)}
    return sets


def main():
    projects = {'PRJNA565979':read_report('PRJNA565979'), 'PRJNA1055086':read_report('PRJNA1055086')}
    evidence = {p:source_evidence(p, rows, HEK3 if p=='PRJNA565979' else CFTR)
                for p, rows in projects.items()}
    stats, unavailable, all_ids = {}, [], set()
    by_cluster = defaultdict(list)
    duplicate_streams = defaultdict(list)
    sequences = Counter()
    for run in (*HEK3,*CFTR,*CFTR_EXTRA,*DUPE):
        row = projects['PRJNA565979' if run in HEK3 else 'PRJNA1055086'][run]
        xml = ET.parse(SOURCES / (row['experiment_accession']+'.xml'))
        assert [el.tag for el in xml.findall('.//LIBRARY_LAYOUT/*')] == [row['library_layout']]
        assert row['library_layout']=='SINGLE' and len(row['fastq_ftp'].split(';'))==1 and not row['submitted_ftp']
        path = raw_path(run)
        if not path.is_file():
            assert run in DUPE
            unavailable.append(run)
            continue
        counts, flowcells, sequence_quality = Counter(), Counter(), hashlib.sha256()
        n = bases = 0
        for rid, key, seq, quality in fastq(path):
            assert rid not in all_ids
            all_ids.add(rid)
            n += 1
            bases += len(seq)
            payload = seq+b'\n'+quality+b'\n'
            sequence_quality.update(payload)
            if run == 'SRR10285861':
                sequences[seq] += 1
            if key is None:
                counts['unkeyed'] += 1
            else:
                cluster, end = key
                counts['R'+end.decode()] += 1
                flowcells[b':'.join(cluster.split(b':')[:4]).decode()] += 1
                by_cluster[cluster].append((run,rid,end,hashlib.sha256(payload).digest()))
        assert (n,bases) == (int(row['read_count']), int(row['base_count'])), run
        md5 = digest(path, 'md5')
        signature = sequence_quality.hexdigest()
        duplicate_streams[signature].append(run)
        stats[run] = {'reads':n,'bases':bases,'sample':row['sample_accession'],'experiment':row['experiment_accession'],
                      'headers':dict(counts),'flowcell_lanes':dict(flowcells),'sha256':digest(path),
                      'local_md5':md5,'archive_md5':row['fastq_md5'],'sequence_quality_sha256':signature}
        if run in DUPE:
            assert md5 == row['fastq_md5'], run

    comparisons, cross_end_ids, repeated_ids = Counter(), set(), set()
    ledger = HERE/'CLUSTER_LINKS.csv.gz'
    with ledger.open('wb') as base, gzip.GzipFile(filename='',fileobj=base,mode='wb',mtime=0) as packed:
        with io.TextIOWrapper(packed,encoding='utf-8',newline='') as stream:
            writer=csv.writer(stream)
            writer.writerow(('cluster','run_a','read_a','end_a','run_b','read_b','end_b','comparison'))
            for key, hits in sorted(by_cluster.items()):
                for (r1,id1,e1,d1),(r2,id2,e2,d2) in itertools.combinations(sorted(hits),2):
                    kind = ('opposite_ends' if e1!=e2 else 'identical_same_end' if d1==d2
                            else 'different_basecalls_or_qualities_same_end')
                    comparisons[r1,r2,kind] += 1
                    if kind=='opposite_ends':
                        cross_end_ids.update((id1,id2))
                    elif kind=='identical_same_end':
                        repeated_ids.update((id1,id2))
                    writer.writerow((key.decode(),r1,id1,e1.decode(),r2,id2,e2.decode(),kind))
    cases=case_sets()
    assert all(ids <= all_ids for ids in cases.values())
    out={'protocol_sha256':digest(HERE/'PROTOCOL.md'),'projects':evidence,'raw':stats,'unavailable':unavailable,
         'raw_total':len(all_ids),'distinct_clusters':len(by_cluster),
         'cluster_comparisons':[dict(run_a=a,run_b=b,type=k,n=n) for (a,b,k),n in sorted(comparisons.items())],
         'cluster_ledger_sha256':digest(ledger),
         'exact_duplicate_streams':[sorted(runs) for runs in duplicate_streams.values() if len(runs)>1],
         'opposite_end_record_ids':sorted(cross_end_ids),
         'targeted_case_audit':{name:{'n':len(ids),'same_end_repeat':len(ids & repeated_ids),
                                      'opposite_end':len(ids & cross_end_ids)} for name,ids in sorted(cases.items())},
         'manual_hek3_sequence_key':{'total':sum(sequences.values()),'distinct':len(sequences),
                                     'singleton_records':sum(n==1 for n in sequences.values()),
                                     'largest_sequence_family':max(sequences.values())}}
    (HERE/'FRAGMENT_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:out[k] for k in ('raw_total','distinct_clusters','unavailable','cluster_comparisons',
                                        'exact_duplicate_streams','targeted_case_audit','manual_hek3_sequence_key')},indent=2))


if __name__=='__main__':
    main()
