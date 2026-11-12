// Read the atom_site loop of a textual PDBx/mmCIF coordinate file.
function* tokens(source){
  let i=0;
  while(i<source.length){
    const char=source[i];
    if(/\s/.test(char)){i++;continue;}
    if(char==='#'){while(i<source.length&&source[i]!=='\n')i++;continue;}
    if(char===';'&&(i===0||source[i-1]==='\n')){
      const stop=source.indexOf('\n;',i+1);
      if(stop<0)throw new Error('Unterminated mmCIF text field.');
      yield {value:source.slice(i+1,stop),plain:false};
      i=stop+2;continue;
    }
    if(char==='"'||char==="'"){
      const begin=++i;
      while(i<source.length&&(source[i]!==char||(i+1<source.length&&!/\s/.test(source[i+1]))))i++;
      if(i>=source.length)throw new Error('Unterminated mmCIF quoted value.');
      yield {value:source.slice(begin,i++),plain:false};
      continue;
    }
    const begin=i++;
    while(i<source.length&&!/\s/.test(source[i]))i++;
    yield {value:source.slice(begin,i),plain:true};
  }
}
const absent=value=>value===undefined||value==='.'||value==='?';
const tag=t=>t?.plain&&t.value.startsWith('_');
const control=t=>t?.plain&&/^(?:loop_|stop_|save_|data_|global_)/i.test(t.value);

export function readMmcif(source){
  const lex=tokens(source);
  let token=lex.next().value;
  while(token){
    if(!token.plain||token.value.toLowerCase()!=='loop_'){token=lex.next().value;continue;}
    const columns=[];
    while(tag(token=lex.next().value))columns.push(token.value.toLowerCase());
    if(!columns.some(key=>key.startsWith('_atom_site.')))continue;
    const ix=name=>columns.indexOf('_atom_site.'+name);
    const axes=['cartn_x','cartn_y','cartn_z'].map(ix);
    const element=ix('type_symbol');
    if(axes.includes(-1)||element<0)throw new Error('mmCIF atom_site requires Cartn_x/y/z and type_symbol.');
    const atoms=[];
    let model;
    while(token&&!tag(token)&&!control(token)){
      const row=[token.value];
      for(let k=1;k<columns.length;k++){
        token=lex.next().value;
        if(!token||tag(token)||control(token))throw new Error('Incomplete mmCIF atom_site row.');
        row.push(token.value);
      }
      const value=(...names)=>names.map(n=>row[ix(n)]).find(v=>!absent(v));
      const presentModel=value('pdbx_pdb_model_num')||'1';
      model ??= presentModel;
      const alt=value('label_alt_id','auth_alt_id');
      if(presentModel===model&&(!alt||alt==='A')){
        const xyz=axes.map(k=>absent(row[k])?NaN:Number(row[k]));
        const chain=value('auth_asym_id','label_asym_id')||'Unassigned';
        const atom=value('auth_atom_id','label_atom_id')||row[element];
        const res=value('auth_comp_id','label_comp_id')||'';
        const pos=value('auth_seq_id','label_seq_id')||'';
        atoms.push({xyz,element:row[element].toUpperCase(),chain,label:`${res} ${pos} · ${atom} · ${chain}`.trim()});
        if(atoms.length>30000)throw new Error('This build supports at most 30,000 atoms.');
      }
      token=lex.next().value;
    }
    return {atoms,bonds:[]};
  }
  throw new Error('No PDBx/mmCIF atom_site coordinate loop found.');
}
