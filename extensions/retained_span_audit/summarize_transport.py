"""Reproduce MP14 aggregate counts and explicit zero-cell reporting."""
import pathlib,json,gzip,csv,collections,hashlib
R=pathlib.Path(__file__).resolve().parent/'TRANSPORT_MP14'
summary=json.loads((R/'SUMMARY.json').read_text());protocol=json.loads((R/'PROTOCOL.json').read_text())
published=json.loads((R/'PUBLISHED_RECONCILIATION.json').read_text())
assert set(published['clones'])==set(protocol['groups'])|{'10'}
protocol_sha=hashlib.sha256((R/'PROTOCOL.json').read_bytes()).hexdigest()
c=collections.Counter();admit=collections.Counter();cross=collections.Counter();totals=collections.Counter();rows=[]
for name in summary['result_files']:
    j=json.loads((R/name).read_text());assert j['protocol_sha256']==protocol_sha
    ledger=R/f"clone{j['clone']}_pairs.jsonl.gz"
    assert hashlib.sha256(ledger.read_bytes()).hexdigest()==j['ledger_sha256']
    with gzip.open(ledger,'rt') as f:assert sum(1 for _ in f)==sum(j['totals'].values())
    totals.update(j['totals'])
    local={(x['library'],x['padding'],x['target'],x['delta'],x['call']):x['pairs'] for x in j['counts']};c.update(local)
    for x in j['admission']:admit[x['library'],x['admitted_mates']]+=x['pairs']
    for x in j['r1_vs_r2_literal']:cross[x['library'],x['padding'],x['target'],x['delta'],x['r1'],x['r2_literal']]+=x['pairs']
    for t in ['DNA','cDNA']:
        for p in [0,9]:
            for a in ['long40','short18']:
                for d in protocol['deltas']:
                    assert sum(local.get((t,p,a,d,k),0) for k in ['E','NE','U'])==j['totals'][t]
                    for call in ['E','NE','U','mate_conflict']:
                        rows.append({'clone':j['clone'],'library':t,'padding':p,'target':a,'delta':d,'call':call,'pairs':local.get((t,p,a,d,call),0)})
assert len(rows)==3808 and sum(totals.values())==416656
aggregate={'clones':17,'totals':dict(totals),
 'primary_counts':[{'library':t,'target':a,**{k:c[t,9,a,36,k] for k in ['E','NE','U','mate_conflict']}} for t in ['DNA','cDNA'] for a in ['long40','short18']],
 'admission':[{'library':t,'admitted_mates':m,'pairs':v} for (t,m),v in sorted(admit.items())],
 'primary_opposed_r1_r2_event_comparisons':sum(cross[t,9,a,36,x,y] for t in totals for a in ['long40','short18'] for x,y in [('E','NE'),('NE','E')]),
 'primary_joint_determinate_event_comparisons':sum(cross[t,9,a,36,x,y] for t in totals for a in ['long40','short18'] for x in ['E','NE'] for y in ['E','NE']),
 'full_strata':952,'published_clone_catalog_reconciled':True,
 'limitations':['same-study computational holdout','no author-declared main target alleles outside clone10','mate agreement is not independently adjudicated biological accuracy','raw-denominator admission loss retained']}
(R/'AGGREGATE.json').write_text(json.dumps(aggregate,indent=2)+'\n')
with (R/'COMPLETE_COUNTS.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('VALIDATED',aggregate['clones'],sum(totals.values()),aggregate['full_strata'],len(rows))
