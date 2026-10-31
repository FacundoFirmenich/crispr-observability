"""Exact reference-context checks. These do not identify the samples' isoforms."""
from pathlib import Path
import json,hashlib,itertools,csv,argparse
ROOT=Path(__file__).resolve().parent
S=ROOT/'sources'
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'results');O=ap.parse_args().output;O.mkdir(parents=True,exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
codons=dict(zip((''.join(c)for c in itertools.product('TCAG',repeat=3)),
 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
def translate(s):
 assert len(s)%3==0 and set(s)<=set('ACGT')
 return ''.join(codons[s[i:i+3]]for i in range(0,len(s),3))
u=json.loads((S/'uniprot_cftr_full.raw.json').read_text());c=json.loads((S/'ensembl_cftr_stable_cds.raw.json').read_text());t=json.loads((S/'ensembl_cftr_stable_transcript.raw.json').read_text());af=json.loads((S/'alphafold_cftr_metadata.raw.json').read_text());g=json.loads((S/'ensembl_cftr_gene.raw.json').read_text())
protein=u['sequence']['value'];cds=c['seq'];peptide='GTIKENIIFGVSYD'
assert translate(cds)==protein+'*'
hits=[i for i in range(len(protein))if protein.startswith(peptide,i)];assert hits==[499]
start=hits[0]*3;native=cds[start:start+3*len(peptide)]
reference=json.loads((ROOT/'frozen_sources/raw_replay/metadata/REFERENCE.json').read_text())
programmed=reference['exact_programmed_core']+'AT';observed_modal_local=reference['mutant_core']+'AT'
assert translate(native)==translate(programmed)==peptide
assert native!=programmed
assert cds[1520:1523]=='CTT'
local_deletion=native[:1520-start]+native[1523-start:]
assert local_deletion==observed_modal_local
rows=[]
for i in range(len(peptide)):
 a=native[i*3:i*3+3];b=programmed[i*3:i*3+3];position=500+i
 rows.append({'canonical_residue_1':position,'amino_acid':peptide[i],'canonical_cds_first_1':1498+3*i,'canonical_codon':a,'programmed_codon':b,'synonymous_difference':a!=b,'deleted_bases_overlap_codon':any(1521<=j<=1523 for j in range(1498+3*i,1501+3*i))})
reference_set=[]
for model in af:
 seq=model['sequence'];loc=[i+1 for i in range(len(seq))if seq.startswith(peptide,i)]
 if loc:reference_set.append({'uniprot_accession':model['uniprotAccession'],'alphafold_entry':model['entryId'],'sequence_length':len(seq),'native_peptide_start_1':loc,'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'model_version':model['latestVersion'],'model_sequence_matches_canonical_uniprot':seq==protein})
assert len(reference_set)==3 and len({r['sequence_length']for r in reference_set})==3
record={'status':'PASS','analysis_type':'PUBLIC_REFERENCE_AUDIT_NOT_SAMPLE_OR_FUNCTION_VALIDATION','canonical_uniprot':u['primaryAccession'],'canonical_sequence_version':u['entryAudit']['sequenceVersion'],'canonical_protein_length':len(protein),'ensembl_gene':g['id'],'ensembl_gene_version':g['version'],'canonical_transcript':t['id']+'.'+str(t['version']),'assembly':t['assembly_name'],'cds_length_including_stop':len(cds),'cds_translation_exactly_canonical_protein_plus_terminal_stop':True,'local_peptide':peptide,'canonical_local_peptide_coordinates_1_inclusive':[500,513],'canonical_local_cds_coordinates_1_inclusive':[1498,1539],'canonical_local_DNA':native,'programmed_local_DNA':programmed,'modal_untreated_reference_local_DNA':observed_modal_local,'canonical_and_programmed_DNA_differ_at_bases':sum(a!=b for a,b in zip(native,programmed)),'synonymous_codons_changed':[r['canonical_residue_1']for r in rows if r['synonymous_difference']],'coding_deletion_validated':{'expression':'c.1521_1523delCTT','deleted_reference_bases':cds[1520:1523],'protein_event':'p.Phe508del','overlaps_canonical_codons':[r['canonical_residue_1']for r in rows if r['deleted_bases_overlap_codon']],'exactly_matches_preserved_untreated_modal_local_core':True},'reference_set_local_peptide_counterexample':reference_set,'interpretation':'The same declared local peptide is present in three published sequence-reference records of unequal complete lengths. This is a finite reference-set non-identification example, not evidence that any of those isoforms occurred or were expressed in the study samples. AlphaFold copies reference sequences and is not an independent molecular truth channel.','limitations':['No sample isoform, haplotype, transcript, protein expression or function assigned.','No newly generated biomolecule or new measured biological observation.','Public reference lookups are a post hoc annotation audit; original counts, policies and conclusions are frozen.'],'source_files_sha256':{p.name:sha(p)for p in [S/'uniprot_cftr_full.raw.json',S/'ensembl_cftr_stable_cds.raw.json',S/'ensembl_cftr_stable_transcript.raw.json',S/'alphafold_cftr_metadata.raw.json',ROOT/'frozen_sources/raw_replay/metadata/REFERENCE.json']}}
(O/'MOLECULAR_REFERENCE_AUDIT.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
with (O/'CANONICAL_AND_PROGRAMMED_CODONS.csv').open('w',newline='')as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with (O/'LOCAL_PEPTIDE_REFERENCE_SCOPE.csv').open('w',newline='')as f:
 w=csv.DictWriter(f,fieldnames=list(reference_set[0]));w.writeheader();w.writerows(reference_set)
for name,sequence,kind in [('canonical_local_DNA',native,'dna'),('programmed_local_DNA',programmed,'dna'),('declared_native_local_peptide',peptide,'protein')]:
 p=O/(name+'.fasta');p.write_text('>'+name+' | reference_context_only\n'+sequence+'\n')
print(json.dumps(record,ensure_ascii=False))
