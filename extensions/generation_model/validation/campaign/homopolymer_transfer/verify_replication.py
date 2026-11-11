"""Verify a completed panel and compare it with the previous two donors."""

import argparse
import collections
import gzip
import hashlib
import json
import zipfile
from pathlib import Path


def sha(data):
    return hashlib.sha256(data).hexdigest()


def classify(reads):
    scores=[r['scores'] for r in reads]
    if len(scores)<10 or 5*sum(min(s)<=59 for s in scores)<4*len(scores):
        return 'U'
    supports=[sum(s[i]<=s[1-i] for s in scores) for i in (0,1)]
    winners=[str(i) for i,n in enumerate(supports) if 2*n>len(scores)]
    return winners[0] if len(winners)==1 else 'U'


def tally(rows):
    before,after,removed=(collections.Counter() for _ in range(3))
    for row in rows:
        def category(value):
            return 'U' if value=='U' else 'C' if value==row['expected'] else 'D'
        initial,final=category(row['primary']),category(row['exact_majority'])
        before[initial]+=1
        after[final]+=1
        if initial!='U' and final=='U':
            removed[initial]+=1
    return {'initial':{k:before[k] for k in 'CDU'},'final':{k:after[k] for k in 'CDU'},
            'withdrawn':{k:removed[k] for k in 'CD'}}


def verify(donors, panel):
    assert set(donors)==set(panel['donors'])
    seen=set()
    out={}
    for donor,data in donors.items():
        specs=panel['donors'][donor]['strata']
        assert set(data['strata'])=={key for key,rows in specs.items() if rows}
        rows_by_length=collections.defaultdict(list)
        names=set()
        pairs=0
        for key,rows in data['strata'].items():
            expected={r['pos']:r for r in specs[key]}
            assert len(rows)==len(expected) and {r['pos'] for r in rows}==set(expected)
            for row in rows:
                pos=row['pos']
                assert pos not in seen
                seen.add(pos)
                assert all(row[k]==v for k,v in expected[pos].items())
                reads=row['reads']
                ids=[r['id_sha256'] for r in reads]
                assert len(ids)==len(set(ids))==min(20,row['eligible'])
                names.update(ids)
                pairs+=len(ids)
                primary=classify(reads)
                assert row['primary']==primary
                gate=primary if primary!='U' and 2*sum(r['flank']==primary for r in reads)>len(reads) else 'U'
                assert row['exact_majority']==gate
            rows_by_length[key.split('_')[1]].extend(rows)
        out[donor]={'pairs':pairs,'distinct_reads':len(names),'length':{key:tally(rows) for key,rows in rows_by_length.items()},
                    'strata':{key:tally(rows) for key,rows in data['strata'].items()}}
    return out,seen


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--reference-archive',type=Path,required=True)
    args=p.parse_args()
    root=args.repo/'results'
    new=json.loads((root/'HOMOPOLYMER_HG003_HG004_RESULT.json').read_bytes())
    old=json.loads((root/'HOMOPOLYMER_HG006_HG007_RESULT.json').read_bytes())
    previous=(root/'HOMOPOLYMER_FROZEN_PANEL.json').read_bytes()
    panel=(root/'HOMOPOLYMER_HG003_HG004_FROZEN_PANEL.json').read_bytes()
    aliases=(root/'HOMOPOLYMER_ALIAS_READBACK.json').read_bytes()
    inventory=(root/'HOMOPOLYMER_PUBLIC_SOURCE_INVENTORY.json').read_bytes()
    runner=(args.repo/'campaign/homopolymer_transfer/run_panel.py').read_bytes()
    with zipfile.ZipFile(args.reference_archive) as z:
        reference=gzip.decompress(z.read('inputs/chr20_sequence.txt.gz'))
    assert new['binding']=={'inventory':sha(inventory),'panel':sha(panel),'runner':sha(runner),
                            'aliases':sha(aliases),'reference_sha256':sha(reference)}
    assert old['binding']['panel']==sha(previous)
    assert json.loads(panel)['excluded_panel_sha256']==sha(previous)
    a,positions_a=verify(old['donors'],json.loads(previous))
    b,positions_b=verify(new['donors'],json.loads(panel))
    assert len(positions_a)==128 and len(positions_b)==113 and not positions_a&positions_b
    result={'total_distinct_loci':len(positions_a|positions_b),'donors':a|b}
    (root/'HOMOPOLYMER_REPLICATION_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    for donor,data in result['donors'].items():
        print(donor,'pairs',data['pairs'],'distinct reads',data['distinct_reads'],data['length'])


if __name__=='__main__':
    main()
