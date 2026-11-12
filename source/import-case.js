import {Vector3,Box3} from 'three';

const schema='observability-case/1';
const maxAtoms=30000, maxBases=360;
const elements={C:'#82baff',N:'#acccff',O:'#f6aab8',P:'#d8c5ff',S:'#f7c1ad',H:'#e7eaff'};

function finiteCoordinates(values){
  if(!Array.isArray(values)||values.length!==3||!values.every(n=>Number.isFinite(n)&&Math.abs(n)<=1e6))
    throw new RangeError('Invalid or unsupported Cartesian coordinate.');
  return values;
}

function validatePortable(input){
  if(input?.schema!==schema||typeof input.name!=='string'||!input.name.trim()||input.name.length>240)
    throw new Error('Use a valid observability-case/1 file with a case name.');
  if(input.kind==='sequence'){
    if(typeof input.sequence!=='string'||!(/^[ACGTN]+$/).test(input.sequence)||input.sequence.length<2||
      input.sequence.length>maxBases)
      throw new RangeError(`Use 2–${maxBases} DNA bases in A, C, G, T or N.`);
    return {schema,kind:'sequence',name:input.name,format:input.format||'CASE JSON',sequence:input.sequence};
  }
  if(input.kind!=='structure'||input.units!=='angstrom'||!Array.isArray(input.atoms)||
    !input.atoms.length||input.atoms.length>maxAtoms||!Array.isArray(input.bonds)||input.bonds.length>6*maxAtoms)
    throw new Error('Use an ångström structure with a valid atom array and explicit bond array.');
  const atoms=input.atoms.map(atom=>{
    if(!atom||typeof atom.element!=='string'||!/^[A-Z][A-Z]?$/.test(atom.element)||
      typeof atom.label!=='string'||atom.label.length>100||typeof atom.chain!=='string'||
      !atom.chain||atom.chain.length>32||atom.chain==='*')
      throw new Error('Invalid atom metadata.');
    return {xyz:finiteCoordinates(atom.xyz),element:atom.element,label:atom.label,chain:atom.chain};
  });
  const keySet=new Set();
  const bonds=[];
  for(const pair of input.bonds){
    if(!Array.isArray(pair)||pair.length!==2||
      !pair.every(n=>Number.isInteger(n)&&n>=0&&n<atoms.length)||pair[0]===pair[1])
      throw new Error('An explicit bond references an invalid atom index.');
    const [a,b]=[Math.min(...pair),Math.max(...pair)],key=`${a}:${b}`;
    if(!keySet.has(key)){keySet.add(key);bonds.push([a,b]);}
  }
  return {schema,kind:'structure',name:input.name,format:input.format||'CASE JSON',units:'angstrom',atoms,bonds};
}

function readPdb(text){
  const atoms=[],rawBonds=[],indices=new Map();
  let firstModelClosed=false;
  for(const line of text.split(/\r?\n/)){
    if(line.startsWith('ENDMDL'))firstModelClosed=true;
    if(!firstModelClosed&&/^(ATOM  |HETATM)/.test(line)&&[' ','A'].includes(line[16])){
      const fields=[line.slice(30,38),line.slice(38,46),line.slice(46,54)];
      const serial=Number(line.slice(6,11));
      if(fields.some(field=>!field.trim())||!Number.isInteger(serial)||serial<=0)
        throw new Error('Invalid PDB atom record.');
      const atom=line.slice(12,16).trim(),res=line.slice(17,20).trim();
      const chain=line.slice(21,22).trim()||'Unassigned';
      const element=line.slice(76,78).trim().toUpperCase()||atom.replace(/^\d+/,'')[0];
      indices.set(serial,atoms.length);
      atoms.push({xyz:fields.map(Number),element,chain,label:`${res} ${line.slice(22,27).trim()} · ${atom} · ${chain}`});
    }
    if(line.startsWith('CONECT')){
      const fields=line.slice(6,31).match(/.{1,5}/g)?.map(v=>v.trim()).filter(Boolean).map(Number)||[];
      for(const target of fields.slice(1))rawBonds.push([fields[0],target]);
    }
  }
  const bonds=rawBonds.filter(([a,b])=>a!==b&&indices.has(a)&&indices.has(b)).map(pair=>pair.map(n=>indices.get(n)));
  return {atoms,bonds};
}

function readXyz(text){
  const lines=text.trimStart().split(/\r?\n/);
  const count=Number(lines[0]);
  if(!Number.isInteger(count)||count<1||count>maxAtoms)throw new RangeError('Invalid XYZ atom count.');
  const rows=lines.slice(2).filter(line=>line.trim());
  if(rows.length!==count)throw new Error('Use a single-frame XYZ file with exactly the declared atom count.');
  const atoms=rows.map((row,i)=>{
    const [symbol,x,y,z]=row.trim().split(/\s+/);
    const element=symbol.toUpperCase();
    return {element,xyz:[x,y,z].map(Number),chain:'1',label:`${element} · atom ${i+1}`};
  });
  return {atoms,bonds:[]};
}

export function readCase(filename,content){
  const ext=filename.split('.').at(-1).toLowerCase();
  const name=filename.slice(0,200);
  if(ext==='json')return validatePortable(JSON.parse(content));
  if(['fa','fasta','fna'].includes(ext)){
    const [head,...lines]=content.trim().split(/\r?\n/);
    if(!head?.startsWith('>')||lines.some(line=>line.startsWith('>')))
      throw new Error('Use a single-record DNA FASTA file.');
    const sequence=lines.join('').replace(/\s+/g,'').toUpperCase();
    return validatePortable({schema,kind:'sequence',name:head.slice(1,200).trim()||name,sequence,format:'FASTA'});
  }
  if(ext!=='pdb'&&ext!=='xyz')throw new Error('Supported inputs are PDB, XYZ, DNA FASTA and observability-case/1 JSON.');
  return validatePortable({schema,kind:'structure',name,units:'angstrom',
    format:ext.toUpperCase(),...(ext==='pdb'?readPdb(content):readXyz(content))});
}

export function importedGeometry(model,chain){
  const points=[],links=[],selected=new Map();
  const box=new Box3();
  model.atoms.forEach(({xyz,label,element,chain:group},id)=>{
    if(chain!=='*'&&chain!==group)return;
    const p=new Vector3(...xyz);
    box.expandByPoint(p);
    const point={p,xyz,label,r:element==='H'?.27:.5,color:elements[element]||'#bfa9f2'};
    selected.set(id,point);
    points.push(point);
  });
  for(const [a,b] of model.bonds){
    const p=selected.get(a),q=selected.get(b);
    if(p&&q&&p.p.distanceToSquared(q.p)>1e-8)links.push({a:p.p,b:q.p,r:.16,color:'#9cadd7'});
  }
  const center=box.getCenter(new Vector3());
  points.forEach(point=>point.p.sub(center));
  const frame=Math.max(7,box.getSize(new Vector3()).length()*.78);
  return {points,links,frame};
}
