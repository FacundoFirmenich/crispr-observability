const lessons=[
  {title:'Find the base outside the read',mode:'edge',preset:'edge',
   instructions:['Load the archived witness. The anchor [210, 222) shifts to [213, 225).',
     'Click a pink base beyond the boundary, then press Read +1 nt.',
     'Check the count and click one of the nearby read ends in the sweep.'],
   expect:'At read end 224, 11 of 12 anchor bases are observed. At 225, all 12 fit. The last included index is one less than the read end.',
   pass:r=>r?.kind==='coordinate_containment'&&r.inputs.start===210&&r.inputs.end===222&&r.inputs.shift===3&&r.inputs.L===225},
  {title:'Compare spans with the same label',mode:'compare',preset:'same-label',
   instructions:['Load the second archived record. Compare the programmed and mutant reference intervals.',
     'Read the two extracted strings. Under Interval B, choose Mutant reference · optimal 1.',
     'Toggle Display reverse complement to check that the overlap stays constant.'],
   expect:'The default outputs share the non-event label despite lengths 72 and 39 nt. The stored alternative optimal span restores string and coordinate agreement.',
   pass:r=>r?.record==='SRR29623271.868'&&r.inputs.selectors[1]==='1:0'&&r.result.identicalString},
  {title:'Inspect a modeled duplex',mode:'atom',preset:'atom',
   instructions:['Load the model built from the archived 42-base span. Press 12-pair detail.',
     'Orbit, zoom, and click an atom to inspect its coordinates in ångströms.',
     'Isolate a strand or change to Space-filling. Use Export PNG to save the current camera and crop.'],
   expect:'The crop [129, 141) contains 492 heavy atoms when both strands are shown. A crop or visibility control changes the display selection.',
   pass:r=>r?.kind==='ideal_duplex_inspection'&&r.inputs.start===129&&r.inputs.end===141},
  {title:'Load and replay your own case',mode:'import',
   instructions:['Use Open a case or replay. Load a structure for atom distances or a FASTA for span comparison.',
     'For a structure, choose two atoms and vary the distance threshold. For FASTA, change the interval endpoints.',
     'Download the experiment JSON, then reopen it with the same file control.'],
   expect:'The replay restores the case and settings and recomputes the numerical result. Imported case data stays on the device.',
   pass:r=>r?.kind==='imported_structure_inspection'&&Boolean(r.result.measurement)}
];

const controlNotes={
  edge:[
    ['Anchor start/end','These coordinates specify the unshifted required interval. Start is included and end is excluded.'],
    ['Product displacement','Adds the same integer to both endpoints. It translates the interval and preserves its width.'],
    ['Read end and sweep','The observed range is [0, read end). Containment changes when read end reaches the shifted anchor end. Pin one result before varying a parameter.'],
    ['Read boundary in 3D','The plane marks the exclusive read end. Turning it off changes the drawing only.']],
  compare:[
    ['Record and rule selectors','An archived rule selects a saved span. An optimal alternative selects a different saved placement for that rule.'],
    ['Interval endpoints','Editing an endpoint switches that interval to Custom. Values are always entered in forward record coordinates.'],
    ['Reverse complement','Changes the displayed orientation and maps [a, b) to [N − b, N − a). The number of shared positions is invariant.']],
  atom:[
    ['Start and end','Restrict both strands to the selected positions within [114, 156).'],
    ['Visible strands and representation','The complement is computed. Ball and stick shows covalent links. Space-filling changes glyph radii and hides the links.'],
    ['Fit view and export','Fit view resets the camera. Export PNG samples the current camera, including the selected zoom and orientation.']],
  import:[
    ['Chain and atom numbers','The numbers index currently visible atoms, starting at 1. Click two atoms or enter their visible indices.'],
    ['Distance threshold','Compares Euclidean distance in Å with the selected threshold. An equal distance is included.'],
    ['Case and replay','A case JSON stores normalized input. Replay JSON also carries the controls and camera and recomputes the result.']]
};

export function setupGuide($,load){
  let index=0,current,currentMode='edge';
  function render(){
    const lesson=lessons[index];
    $('guide-step').textContent=`${index+1} of ${lessons.length} · ${lesson.title}`;
    $('guide-instructions').replaceChildren(...lesson.instructions.map(line=>{
      const li=document.createElement('li');li.textContent=line;return li;
    }));
    $('guide-expect').textContent=lesson.expect;
    $('guide-previous').disabled=index===0;
    $('guide-next').disabled=index===lessons.length-1;
    $('guide-preset').hidden=!lesson.preset;
    $('guide-result').textContent=lesson.pass(current)?'Checkpoint reached in the current experiment.':
      currentMode===lesson.mode?'Try the controls and compare the current result with the expectation.':'Open this experiment to follow the checkpoint.';
  }
  $('guide-toggle').addEventListener('click',()=>{
    const open=$('guide-body').hidden;
    $('guide-body').hidden=!open;
    $('guide-toggle').setAttribute('aria-expanded',String(open));
    $('guide-toggle').textContent=open?'Close walkthrough':'Open walkthrough';
    if(open)render();
  });
  $('guide-previous').addEventListener('click',()=>{index--;render();});
  $('guide-next').addEventListener('click',()=>{index++;render();});
  $('guide-preset').addEventListener('click',()=>load(lessons[index].preset));
  function update(result,mode){
    current=result;
    if(currentMode!==mode||!$('control-help').firstChild){
      $('control-help').replaceChildren(...controlNotes[mode].map(([label,explanation])=>{
        const p=document.createElement('p'),title=document.createElement('strong');
        title.textContent=label;
        p.append(title,document.createTextNode(' — '+explanation));
        return p;
      }));
    }
    currentMode=mode;
    if(!$('guide-body').hidden)render();
  }
  return {update};
}
