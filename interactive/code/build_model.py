"""Build the public, prose-free numerical companion from frozen source tables."""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import csv, json, hashlib, argparse

HERE = Path(__file__).resolve().parents[1]
PREPS = {'1':'CFB-26-16','2':'CFB-26-02','3':'CFB-06-04'}
EVENTS = {
 'CTT_plus_full_SE2':'exact_programmed',
 'CTT_full_SE2_on_witnessed_observed_coding':'exact_programmed',
 'native_local_peptide_all_observed_sequences':'native_peptide',
 'native_local_peptide_on_witnessed_observed_coding':'native_peptide',
}
OPS = [
 ('broad_edge','Edge anchor','none'),('broad_current','Inner anchor','none'),
 ('target','Target motifs','none'),('compatible_union','Compatible union','none'),
 ('crispresso_projected','Mutant projection','mutant'),
 ('programmed_reference','Programmed projection','programmed'),
 ('reference_agreement','Reference-label agreement','both'),
 ('triple_observed_coding_witness','Common-string witness','both'),
]
OP_REF = {x[0]:x[2] for x in OPS}

def partition(N,E,NE,U):
    if any(type(x) is not int or x < 0 for x in [N,E,NE,U]) or N == 0 or E+NE+U != N:
        raise ValueError('A positive integer N and nonnegative E, NE, U with E+NE+U=N are required.')
    return Fraction(E,N), Fraction(E+U,N)

def endpoints(N,E,NE,U):
    lo,hi=partition(N,E,NE,U)
    return {'N':N,'E':E,'NE':NE,'U':U,'L':float(lo),'H':float(hi),
            'L_pp':float(100*lo),'H_pp':float(100*hi),'L_exact':str(lo),'H_exact':str(hi)}

def contrast(t,c):
    tl,th=partition(t['N'],t['E'],t['NE'],t['U'])
    cl,ch=partition(c['N'],c['E'],c['NE'],c['U'])
    lo,hi=tl-ch,th-cl
    return {'L':float(lo),'H':float(hi),'L_pp':float(100*lo),'H_pp':float(100*hi),
            'L_exact':str(lo),'H_exact':str(hi)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sources',type=Path,default=HERE/'sources')
    ap.add_argument('--output',type=Path,default=HERE/'evidence_model.json');args=ap.parse_args()
    src=args.sources; source=[]; tables={}
    for path in sorted((src/'tables').glob('*.csv')):
        rows=list(csv.DictReader(path.open(newline='')));tables[path.name]=rows
        source.append({'id':path.stem,'path':'tables/'+path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':len(rows),'kind':'frozen_CSV'})
    for fn,pub in [('WITNESS_BOUNDS.json','raw_replay/expected/WITNESS_BOUNDS.json'),('JOINT_STABILITY.json','engine_replay/original/JOINT_STABILITY.json'),('fig1_coordinate_geometry.csv','figures/inputs/fig1_coordinate_geometry.csv')]:
        p=src/'numerical'/fn;source.append({'id':p.stem,'path':pub,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'kind':'frozen_numerical_artifact'})
    def ref(fn,i): return {'id':Path(fn).stem,'path':'tables/'+fn,'row':i+2}
    streams=[]
    for i,r in enumerate(tables['Table_S1_Primary_Run_Manifest.csv']):
        streams.append({'run':r['run_accession'],'preparation':r['preparation'],'condition':r['condition'],'N':int(r['raw_reads']),'replica':int(r['source_replica']),'sha256':r['sequence_quality_stream_sha256'],'source':ref('Table_S1_Primary_Run_Manifest.csv',i)})
    runmap={r['run']:r for r in streams};obs=[]
    for fn in ['Table_S5_Projection_Counts.csv','Table_S7_Coding_Witness_Counts.csv']:
        for i,r in enumerate(tables[fn]):
            run=r['run_accession'];op=r['observer'];event=EVENTS[r['event']]
            v={'id':'obs:'+op+':'+event+':'+run,'run':run,'preparation':runmap[run]['preparation'],'condition':runmap[run]['condition'],'event':event,'sourceEvent':r['event'],'operator':op,'reference':OP_REF[op],'source':ref(fn,i)}
            v.update(endpoints(*(int(r[k]) for k in ['N','E','NE','U'])))
            assert Fraction(v['L_exact'])==Fraction(r['L_exact']) and Fraction(v['H_exact'])==Fraction(r['H_exact'])
            obs.append(v)
    bykey={(v['preparation'],v['condition'],v['event'],v['operator']):v for v in obs}
    contrasts=[]
    for fn in ['Table_S6_Projection_Contrasts.csv','Table_S8_Coding_Witness_Contrasts.csv']:
        for i,r in enumerate(tables[fn]):
            prep=PREPS[r['preparation']];event=EVENTS[r['event']];op=r['observer'];control=r['control']
            t=bykey[(prep,'PE treated',event,op)];c=bykey[(prep,control,event,op)]
            v={'id':'contrast:'+op+':'+event+':'+prep+':'+control,'preparation':prep,'event':event,'operator':op,'reference':OP_REF[op],'control':control,'treated':t['run'],'comparator':c['run'],'treatedCounts':{k:t[k] for k in ['N','E','NE','U']},'controlCounts':{k:c[k] for k in ['N','E','NE','U']},'source':ref(fn,i)}
            v.update(contrast(t,c));assert Fraction(v['L_exact'])==Fraction(r['L_exact']) and Fraction(v['H_exact'])==Fraction(r['H_exact']);contrasts.append(v)
    geom=[]
    for i,r in enumerate(tables['Table_1_Native_Contrast_Intervals.csv']):
        lo,hi=Fraction(r['L_exact']),Fraction(r['H_exact']);geom.append({'id':'geometry:'+r['preparation']+':'+r['operator'],'preparation':r['preparation'],'operator':r['operator'],'event':'native_peptide','L':float(lo),'H':float(hi),'L_pp':float(100*lo),'H_pp':float(100*hi),'L_exact':str(lo),'H_exact':str(hi),'source':ref('Table_1_Native_Contrast_Intervals.csv',i)})
    software=[]
    for i,r in enumerate(tables['Table_S11_Direct_Software_and_Event_Comparison.csv']):
        den=int(r['aligned_all_references']);mod=int(r['modified']);unmod=int(r['unmodified']);assert mod+unmod+int(r['discarded'])==den
        software.append({'run':r['run_accession'],'preparation':r['preparation'],'condition':r['condition'],'reference':'mutant' if r['configuration']=='observed_mutant_reference' else 'programmed','N':int(r['raw_N']),'qualityPass':int(r['engine_input_Q30']),'aligned':den,'modified':mod,'unmodified':unmod,'softwareModifiedPercent':float(r['modified_percent_reported']),'softwareStatistic':'MODIFIED / aligned; different event and denominator from native-peptide raw hull','source':ref('Table_S11_Direct_Software_and_Event_Comparison.csv',i)})
    sensitivity=[]
    for i,r in enumerate(tables['Table_S12_Alignment_Parameter_Sensitivity.csv']):
        run=r['run_accession'];v={'id':'grid:'+r['job_id']+':'+EVENTS[r['event']],'job_id':r['job_id'],'run':run,'preparation':runmap[run]['preparation'],'condition':runmap[run]['condition'],'event':EVENTS[r['event']],'sourceEvent':r['event'],'operator':'alignment_projection','reference':r['reference'],'setting':r['setting'],'rangeScope':r['range_scope'],'qualityPass':int(r['Q30_pass']),'aligned':int(r['aligned']),'softwareModifiedPercent':float(r['software_MODIFIED_percent_aligned_reported']),'source':ref('Table_S12_Alignment_Parameter_Sensitivity.csv',i)}
        v.update(endpoints(*(int(r[k]) for k in ['raw_N','E','NE','U'])));assert Fraction(v['L_exact'])==Fraction(r['L_exact']) and Fraction(v['H_exact'])==Fraction(r['H_exact']);sensitivity.append(v)
    gindex={(r['preparation'],r['condition'],r['event'],r['reference'],r['setting']):r for r in sensitivity};scontrasts=[]
    for r in sensitivity:
        if r['condition']!='PE treated': continue
        for control in ['mock','untreated']:
            c=gindex[(r['preparation'],control,r['event'],r['reference'],r['setting'])]
            v={k:r[k] for k in ['preparation','event','reference','setting','rangeScope']};v.update({'id':'gridcontrast:'+r['id']+':'+control,'control':control,'treated':r['run'],'comparator':c['run'],'operator':'alignment_projection','source':{'id':'Table_S12_Alignment_Parameter_Sensitivity','path':'tables/Table_S12_Alignment_Parameter_Sensitivity.csv','treatedRow':r['source']['row'],'controlRow':c['source']['row']}});v.update(contrast(r,c));scontrasts.append(v)
    env=[]
    for i,r in enumerate(tables['Table_S12_PE_Mock_Endpoint_Envelopes.csv']):
        sub=[v for v in scontrasts if v['preparation']==PREPS[r['preparation']] and v['event']==EVENTS[r['event']] and v['control']=='mock' and v['rangeScope']=='finite_OAT_grid']
        assert len(sub)==14 and all(v['L']>0 for v in sub)
        assert abs(min(v['L_pp'] for v in sub)-float(r['L_min_pp']))<1e-10
        assert abs(max(v['H_pp'] for v in sub)-float(r['H_max_pp']))<1e-10
        env.append({'preparation':PREPS[r['preparation']],'event':EVENTS[r['event']],'profiles':14,'L_min_pp':float(r['L_min_pp']),'H_max_pp':float(r['H_max_pp']),'L_min_reference':r['L_min_reference'],'L_min_setting':r['L_min_setting'],'H_max_reference':r['H_max_reference'],'H_max_setting':r['H_max_setting'],'allPositive':True,'source':ref('Table_S12_PE_Mock_Endpoint_Envelopes.csv',i)})
    w=json.loads((src/'numerical/WITNESS_BOUNDS.json').read_text());domains=[]
    for r in w['domains']:
        c=r['counts'];v={'run':r['run_accession'],'preparation':runmap[r['run_accession']]['preparation'],'condition':r['condition'],'N':r['N_raw'],'witness':c.get('identical_observed_coding_witness',0),'qualityFailed':c.get('common_quality_fail',0),'missingCorrespondence':c.get('missing_correspondence_witness',0),'codingDisagreement':c.get('different_coding',0),'equalLabelsDifferentCoding':c.get('equal_native_labels_different_coding',0),'source':{'id':'WITNESS_BOUNDS','path':'raw_replay/expected/WITNESS_BOUNDS.json','section':'domains'}}
        assert sum(v[k] for k in ['witness','qualityFailed','missingCorrespondence','codingDisagreement'])==v['N'];domains.append(v)
    stability=[]
    for r in json.loads((src/'numerical/JOINT_STABILITY.json').read_text()):
        stability.append({'run':r['run'],'preparation':runmap[r['run']]['preparation'],'event':EVENTS[r['event']],'scope':r['scope'],'profiles':r['profiles'],'conflicts':r['counts'].get('E_NE_conflict',0),'counts':r['counts'],'source':{'id':'JOINT_STABILITY','path':'engine_replay/original/JOINT_STABILITY.json'}})
    conflicts=Counter(r['run'] for r in tables['Table_S9_Reference_Changes.csv']);assert sum(conflicts.values())==137
    external=[]
    for i,r in enumerate(tables['Table_3_External_Catalog_Accounting.csv']):
        external.append({'case':r['case'],'inputUnit':r['input_unit'],'inputUnits':int(r['input_units']),'tableRows':int(r['table_rows']),'assignmentWeight':int(r['assignment_weight']),'modified':int(r['MODIFIED']),'unmodified':int(r['UNMODIFIED']),'unresolvedTemplateWeight':int(r['unresolved_template_weight']),'jointResolvedSignatures':int(r['joint_resolved_signatures']),'outsideAligned':int(r['raw_outside_aligned_table']),'source':ref('Table_3_External_Catalog_Accounting.csv',i)})
    trunc=[]
    for fn in ['Table_S2_M1_All_Counts.csv','Table_S4_Original_Ambiguity_Veto.csv']:
        for i,r in enumerate(tables[fn]):
            if r['run_accession'] not in runmap or r['policy']!='common_raw_meanQ30':continue
            v={k:r[k] for k in ['observer','trim','policy']};v.update({'run':r['run_accession'],'preparation':runmap[r['run_accession']]['preparation'],'condition':runmap[r['run_accession']]['condition'],'event':EVENTS[r['event']],'variant':'original_ambiguity_veto' if fn.startswith('Table_S4') else 'nominal','source':ref(fn,i)});v.update(endpoints(*(int(r[k]) for k in ['N','E','NE','U'])));trunc.append(v)
    summary={'rawRecords':sum(r['N'] for r in streams),'witnessRetained':sum(r['witness'] for r in domains),'unknown':sum(r['N']-r['witness'] for r in domains),'qualityFailed':sum(r['qualityFailed'] for r in domains),'missingCorrespondence':sum(r['missingCorrespondence'] for r in domains),'codingDisagreement':sum(r['codingDisagreement'] for r in domains),'equalLabelsDifferentCoding':sum(r['equalLabelsDifferentCoding'] for r in domains),'referenceNativeLabelChanges':sum(conflicts.values()),'alignmentJobs':len({r['job_id'] for r in sensitivity}),'observationRows':len(obs),'contrastRows':len(contrasts),'sensitivityEventRows':len(sensitivity),'sensitivityContrastRows':len(scontrasts),'externalTableRows':sum(r['tableRows'] for r in external)}
    assert summary['rawRecords']==177006 and summary['witnessRetained']==170725 and summary['unknown']==6281
    assert summary['codingDisagreement']==83 and summary['equalLabelsDifferentCoding']==18 and summary['alignmentJobs']==144
    model={'schemaVersion':'1.0','release':'MP03D','title':'CRISPR observability — evidence explorer','generatedDate':'2026-10-28','scope':{'description':'Offline numerical and epistemic companion built from frozen source data; manuscript and supplement text are withheld from this public artifact.','interpretation':'Conditional finite-record identification hulls under correct accepted labels; no sampling confidence coverage is claimed.','notClaims':['Molecular classification accuracy','Independent human replication','New biological cohort','Haplotype-specific function','Universal observer ranking']},'catalogs':{'preparations':[{'id':v,'replica':int(k)} for k,v in PREPS.items()],'events':[{'id':'native_peptide','label':'Native local peptide','definition':'Observed local DNA translates in silico to GTIKENIIFGVSYD; no protein-function measurement.'},{'id':'exact_programmed','label':'Exact programmed coding','definition':'Exact CTT plus full SE2 programmed coding string.'}],'operators':[{'id':o,'label':l,'reference':r} for o,l,r in OPS],'references':[{'id':'none','label':'Motif operator'},{'id':'mutant','label':'Mutant reference'},{'id':'programmed','label':'Programmed reference'},{'id':'both','label':'Both references required'}]},'summary':summary,'sources':source,'streams':streams,'observations':obs,'contrasts':contrasts,'geometryContrasts':geom,'softwareMetrics':software,'sensitivityObservations':sensitivity,'sensitivityContrasts':scontrasts,'sensitivityEnvelopes':env,'witnessDomains':domains,'jointStability':stability,'referenceConflicts':[{'run':k,'preparation':runmap[k]['preparation'],'changes':v,'source':{'id':'Table_S9_Reference_Changes','path':'tables/Table_S9_Reference_Changes.csv'}} for k,v in conflicts.items()],'externalCases':external,'truncationObservations':trunc,'geometry':[{'label':r['label'],'start':int(r['start']),'end':int(r['end']),'read_end':int(r['read_end'])} for r in csv.DictReader((src/'numerical/fig1_coordinate_geometry.csv').open())],'graph':graph(),'calculationContract':{'partition':'E + NE + U = N; all are nonnegative integers and N > 0.','hull':'[E/N, (E+U)/N]','contrast':'[E_t/N_t - (E_c+U_c)/N_c, (E_t+U_t)/N_t - E_c/N_c]','completion':'(E+k)/N for integer k from0 throughU is an attainable hypothetical completion, not an observed estimate. If a slider uses q in[0,1], set k=round(q*U) and displayk.','counterfactualWarning':'Changing q or counts does not add reads or biological evidence. Published records remain immutable.','assumptions':['Accepted labels are correct for the declared event.','Every unknown may independently be event or non-event.','Two-arm worst-case completion is independent between declared raw universes.']}}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(model,ensure_ascii=False,separators=(',',':'))+'\n')
    report={'status':'PASS','model':args.output.name,'sha256':hashlib.sha256(args.output.read_bytes()).hexdigest(),'sources':len(source),'summary':summary,'checks':['All144 observation partitions and rational extrema','All96 source contrast extrema reconstructed from E/NE/U','All288 alignment event partitions and rational extrema','All192 alignment two-arm hulls reconstructed','All6 envelope rows match the14-profile finite OAT family','Nine disjoint raw-domain partitions and retained totals','Every source file fingerprint retained','No full manuscript or supplement prose used']}
    (args.output.parent/'MODEL_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

def graph():
    nodes=[]
    def n(id,type,label,status,mechanism,ids,scope,falsification):
        nodes.append({'id':id,'type':type,'label':label,'status':status,'mechanism':mechanism,'evidenceIds':ids,'scope':scope,'falsification':falsification})
    n('E_geometry','evidence','Finite-read geometry','observed_configuration','The shifted distal anchor ends at225 while the documented read ends at224.',['fig1_coordinate_geometry'],'Illustrated coordinate configuration, not every product.','Check the source coordinates and original read length.')
    n('E_anchor','evidence','Anchor contrast data','replayed','Three motif contracts inspect the same raw universes.',['Table_1_Native_Contrast_Intervals'],'Common original-read Q30 mask; PE versus mock.','Reconstruct E/NE/U under the frozen operators.')
    n('E_strings','evidence','Coding correspondence counts','replayed','Raw-domain partitions separate equal strings, missing correspondence and83 different coding strings.',['WITNESS_BOUNDS'],'Nine primary streams; original target plus two baseline projections.','Audit orientation, boundaries and observed string extraction.')
    n('E_counts','evidence','Frozen event partitions','replayed','Accepted events, non-events and unknowns sum to each raw universe.',['Table_S5_Projection_Counts','Table_S7_Coding_Witness_Counts'],'DNA sequence events; two event definitions retained.','Any nonconserved count or raw identity requires correction.')
    n('E_contrasts','evidence','Two-arm conditional hulls','recomputed','Compatible extreme completions define each treated-minus-control hull.',['Table_S6_Projection_Contrasts','Table_S8_Coding_Witness_Contrasts'],'Declared preparation, event, operator and raw denominators.','Exact rational replay must recover both endpoints.')
    n('E_grid','evidence','144 alignment jobs','executed_retrospective','Finite gap-penalty and score profiles alter alignments under two references.',['Table_S12_Alignment_Parameter_Sensitivity'],'126 baseline/OAT jobs plus18 separate stock-profile jobs.','Repeat the pinned engine and compare record identities and counts.')
    n('E_labels','evidence','Distinct-record conflicts','audited','Some read identifiers receive both E and NE across profiles.',['JOINT_STABILITY'],'Scopes overlap; counts must not be added across families.','Reconstruct per-identifier signatures within each declared family.')
    n('E_external','evidence','External assignment catalogs','audited','572 assignment-weight units lack a resolved template in P23H/WT.',['Table_3_External_Catalog_Accounting'],'Catalog rows and source weights, not fabricated read identities.','An authenticated unique template mapping would change the resolved domain.')
    n('A_quality','assumption','Common input and quality policy','declared','Operator contrasts require a fixed raw universe and shared original-read mask.',['Table_S1_Primary_Run_Manifest'],'Counts from incompatible universes cannot be pooled.','A changed stream or mask invalidates that comparison contract.')
    n('A_labels','assumption','Accepted labels are correct','conditional','Hull coverage depends on the accepted event/non-event assignments.',['Table_S7_Coding_Witness_Counts'],'This premise is not established by agreement alone.','Independent molecular truth can reveal accepted-label errors.')
    n('A_completion','assumption','Worst-case unknown completion','mathematical','Unknown records may independently belong to either event class.',[],'Finite-record extrema; attainable fractions are discrete.','Restrictions linking unknown labels would change the identified set.')
    n('A_domain','assumption','Target requires a resolved template','declared','Complete-template equality requires a unique source reference assignment.',['Table_3_External_Catalog_Accounting'],'Ambiguous assignments remain outside the resolved target domain.','A unique authenticated reference-to-template mapping refines the domain.')
    n('A_validation_scope','assumption','Validation scope is computational','declared','These sources record replay and retrospective sensitivity on the same observations; a separate molecular truth or prospective cohort is not supplied.',[],'Statement of this released evidence scope, not a universal absence claim.','An added authenticated truth channel or prospective evaluation changes the evidence scope.')
    n('C_anchor','claim','Observation can depend on product geometry','supported_local','Moving a required anchor beyond the read end can change unresolved mass.',['Table_1_Native_Contrast_Intervals'],'Local finite-read contracts, not a universal ranking.','Revise if the claimed missing complete anchor is present under the frozen orientation.')
    n('C_strings','claim','Same label need not mean same sequence','supported_local','Eighteen equal native labels occur among unequal observed coding strings.',['WITNESS_BOUNDS'],'Sequence correspondence on a record, not unique physical molecule identity.','Revise if independent boundary/orientation checking removes the unequal strings.')
    n('C_hull','claim','Unknown mass constrains the event fraction','conditional_exact','E/N and(E+U)/N are the sharp finite-record extreme completions.',['Table_S7_Coding_Witness_Counts'],'Accepted-label premise; not a confidence interval.','Partition errors or wrong accepted labels require a different calculation.')
    n('C_witness','claim','Common-string witness supports positive contrasts','supported_conditional','A common coding string is required before classifying the witnessed event.',['Table_S8_Coding_Witness_Contrasts'],'Original three-operator witness; three preparations, raw denominators.','Frozen raw replay must recover the published counts and endpoints.')
    n('C_grid','claim','Stable sign can coexist with changing labels','supported_finite_grid','Positive computed conditional PE-minus-mock lower endpoints coexist with163 native-label conflicts across14 OAT profiles.',['Table_S12_PE_Mock_Endpoint_Envelopes','JOINT_STABILITY'],'Retrospective finite grid; original witness is not redefined.','A pinned rerun that changes these extrema or distinct identifiers requires correction.')
    n('C_domain','claim','Reference ambiguity restricts the target domain','supported_catalog','An unresolved reference assignment cannot certify complete-template equality.',['Table_3_External_Catalog_Accounting'],'P23H/WT assignment catalog; source weights retained.','A verified unique template assignment would refine this local domain.')
    n('B_accuracy','boundary','Molecular accuracy remains uncalibrated','open','Computational agreement does not supply an independent truth channel.',[],'No comparative accuracy ranking or full-allele function.','Requires authenticated orthogonal molecular evidence.')
    n('B_transport','boundary','Future transport remains unvalidated','open','A retrospective configuration family is not prospective validation.',[],'New accessions and rules must be specified before outcome inspection.','Requires independent, prospectively frozen evaluation.')
    edges=[('H_anchor',['E_geometry','E_anchor','A_quality'],'C_anchor','Geometry + shared policy + counts'),('H_strings',['E_strings','A_quality'],'C_strings','String comparison + authenticated raw domain'),('H_hull',['E_counts','A_labels','A_completion'],'C_hull','Counts + accepted-label premise + completion model'),('H_witness',['E_strings','E_counts','E_contrasts','A_quality','A_labels','A_completion'],'C_witness','Correspondence + partitions + shared input + two-arm assumptions'),('H_grid',['E_grid','E_labels','A_quality','A_labels','A_completion'],'C_grid','Conditional parameter extrema + per-record conflicts + hull assumptions'),('H_domain',['E_external','A_domain'],'C_domain','Catalog identity + declared target domain'),('H_accuracy',['C_hull','A_validation_scope'],'B_accuracy','Identified bounds do not supply truth calibration'),('H_transport',['C_grid','A_validation_scope'],'B_transport','Local support does not certify future transport')]
    return {'semantics':'Navigation of evidence dependencies. AND means every listed premise is required for this support path. Toggling a premise changes displayed support, not historical evidence or scientific truth.','nodes':nodes,'hyperedges':[{'id':i,'type':'AND','premises':p,'conclusion':c,'label':l,'scope':'Declared finite domain only.'} for i,p,c,l in edges]}

if __name__=='__main__': main()
