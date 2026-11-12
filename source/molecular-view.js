import * as T from 'three';
import {TrackballControls} from 'three/addons/controls/TrackballControls.js';

const ink={blue:'#6cc9ff',violet:'#b1a2ff',rose:'#f5a8bc',pale:'#d6dcff',muted:'#48577e'};
const yAxis=new T.Vector3(0,1,0);
const v=(x,y,z)=>new T.Vector3(x,y,z);
const spheres=new T.SphereGeometry(1,12,8);
const smallSpheres=new T.SphereGeometry(1,8,5);
const rods=new T.CylinderGeometry(1,1,1,8);
const surface=new T.MeshStandardMaterial({roughness:.38,metalness:.08});
const selectionMaterial=new T.MeshBasicMaterial({color:ink.pale,transparent:true,opacity:.64,wireframe:true});

export class MolecularView {
  constructor(canvas, onPick) {
    this.canvas=canvas;
    this.renderer=new T.WebGLRenderer({canvas,antialias:true,powerPreference:'low-power'});
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio,1.75));
    this.renderer.setClearColor('#131c34');
    this.renderer.toneMapping=T.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure=1.45;
    this.scene=new T.Scene();
    this.camera=new T.PerspectiveCamera(39,1,.1,4000);
    this.scene.add(new T.HemisphereLight('#cbe3ff','#514064',2.15));
    const key=new T.DirectionalLight('#e4efff',2.8);
    key.position.set(-20,30,25);
    const fill=new T.DirectionalLight('#c5abff',1.8);
    fill.position.set(10,-20,-12);
    this.scene.add(key,fill);
    this.meshes=new T.Group();
    this.scene.add(this.meshes);
    this.controls=new TrackballControls(this.camera,canvas);
    this.controls.staticMoving=true;
    this.controls.rotateSpeed=2.6;
    this.controls.panSpeed=.55;
    this.controls.noPan=false;
    this.controls.minDistance=3;
    this.controls.maxDistance=1700;
    this.controls.addEventListener('change',()=>this.draw());
    this.scene.add(this.camera);
    this.marker=new T.Mesh(spheres,selectionMaterial);
    this.marker.visible=false;
    this.scene.add(this.marker);
    this.measurementMaterial=new T.LineDashedMaterial({color:'#f5c5dd',dashSize:.8,gapSize:.45,depthTest:false});
    this.measurement=new T.Line(new T.BufferGeometry(),this.measurementMaterial);
    this.measurement.visible=false;
    this.measurement.renderOrder=4;
    this.scene.add(this.measurement);
    this.frame=1;
    this.ray=new T.Raycaster();
    let start;
    canvas.addEventListener('pointerdown',e=>{start=[e.clientX,e.clientY];});
    canvas.addEventListener('pointerup',e=>{
      if(!start||e.button!==0||Math.hypot(e.clientX-start[0],e.clientY-start[1])>6||!this.points)return;
      const rect=canvas.getBoundingClientRect();
      this.ray.setFromCamera(new T.Vector2((e.clientX-rect.left)*2/rect.width-1,1-(e.clientY-rect.top)*2/rect.height),this.camera);
      const hit=this.ray.intersectObject(this.points,false)[0];
      if(hit)onPick(this.points.userData.rows[hit.instanceId]);
    });
    canvas.addEventListener('keydown',e=>{
      const directions={ArrowLeft:[0,1,.16],ArrowRight:[0,1,-.16],ArrowUp:[1,0,.16],ArrowDown:[1,0,-.16]};
      const dir=directions[e.key];
      if(!dir)return;
      e.preventDefault();
      const axis=v(dir[0],dir[1],0).applyQuaternion(this.camera.quaternion).normalize();
      const q=new T.Quaternion().setFromAxisAngle(axis,dir[2]);
      this.camera.position.sub(this.controls.target).applyQuaternion(q).add(this.controls.target);
      this.camera.up.applyQuaternion(q);
      this.camera.lookAt(this.controls.target);
      this.draw();
    });
    this.resize=new ResizeObserver(()=>{
      const rect=canvas.parentElement.getBoundingClientRect();
      const width=Math.round(rect.width),height=Math.round(rect.height);
      if(!width||!height||(width===this.width&&height===this.height))return;
      const first=this.width===undefined;
      this.width=width; this.height=height;
      this.camera.aspect=width/height;
      this.camera.updateProjectionMatrix();
      this.renderer.setSize(width,height,false);
      this.controls.handleResize();
      if(first)this.fit();
      this.draw();
    });
    this.resize.observe(canvas.parentElement);
    this.renderer.setAnimationLoop(()=>this.controls.update());
  }

  draw() {this.renderer.render(this.scene,this.camera);}

  cameraState() {
    return {position:this.camera.position.toArray(),up:this.camera.up.toArray(),target:this.controls.target.toArray()};
  }

  restoreCamera(state) {
    if(!state||!['position','up','target'].every(k=>Array.isArray(state[k])&&state[k].length===3&&
        state[k].every(n=>Number.isFinite(n)&&Math.abs(n)<=1e7)))throw new Error('Invalid camera state.');
    const position=new T.Vector3(...state.position),target=new T.Vector3(...state.target),up=new T.Vector3(...state.up);
    if(position.distanceTo(target)<.1||up.length()<.1)throw new Error('Invalid camera basis.');
    this.camera.position.copy(position);this.camera.up.copy(up.normalize());
    this.controls.target.copy(target);
    this.camera.lookAt(target);
    this.controls.update();
    this.draw();
  }

  renderPng(scale) {
    if(![1,2,4].includes(scale)||!this.width||!this.height)throw new Error('Choose a render scale of 1×, 2× or 4×.');
    const factor=Math.min(scale,4096/Math.max(this.width,this.height));
    const w=Math.max(1,Math.round(this.width*factor)),h=Math.max(1,Math.round(this.height*factor));
    const ratio=this.renderer.getPixelRatio(),aspect=this.camera.aspect;
    let dataUrl;
    try{
      this.renderer.setPixelRatio(1);
      this.renderer.setSize(w,h,false);
      this.camera.aspect=w/h;
      this.camera.updateProjectionMatrix();
      this.draw();
      dataUrl=this.canvas.toDataURL('image/png');
      if(!dataUrl.startsWith('data:image/png;base64,'))throw new Error('The browser could not encode the PNG.');
    }finally{
      this.renderer.setPixelRatio(ratio);
      this.renderer.setSize(this.width,this.height,false);
      this.camera.aspect=aspect;
      this.camera.updateProjectionMatrix();
      this.draw();
    }
    return {dataUrl,width:w,height:h};
  }

  fit() {
    const xScale=Math.max(1,1.25/this.camera.aspect);
    this.camera.up.set(0,1,0);
    this.controls.target.set(0,0,0);
    this.camera.position.set(this.frame*.11,this.frame*.18,this.frame*1.75*xScale);
    this.camera.lookAt(this.controls.target);
    this.controls.update();
    this.draw();
  }

  select(p,radius=1.9) {
    this.marker.visible=Boolean(p);
    if(p){this.marker.position.copy(p);this.marker.scale.setScalar(radius);}
    this.draw();
  }

  setMeasurement(a,b){
    this.measurement.visible=Boolean(a&&b);
    if(a&&b){
      this.measurement.geometry.dispose();
      this.measurement.geometry=new T.BufferGeometry().setFromPoints([a,b]);
      this.measurement.computeLineDistances();
    }
    this.draw();
  }

  setModel(points,links,frame,cut,reset) {
    for(const item of [...this.meshes.children]) {
      this.meshes.remove(item);
      if(item.isInstancedMesh)item.dispose();
      else {item.material.dispose();item.geometry.dispose();}
    }
    const dummy=new T.Object3D(), shade=new T.Color();
    this.points=new T.InstancedMesh(points.length>3500?smallSpheres:spheres,surface,points.length);
    points.forEach((item,i)=>{
      dummy.position.copy(item.p);dummy.scale.setScalar(item.r);dummy.quaternion.identity();
      dummy.updateMatrix();
      this.points.setMatrixAt(i,dummy.matrix);
      this.points.setColorAt(i,shade.set(item.color));
    });
    this.points.userData.rows=points;
    const sticks=new T.InstancedMesh(rods,surface,links.length);
    links.forEach(({a,b,r,color},i)=>{
      const delta=b.clone().sub(a);
      dummy.position.copy(a).lerp(b,.5);
      dummy.scale.set(r,delta.length(),r);
      dummy.quaternion.setFromUnitVectors(yAxis,delta.normalize());
      dummy.updateMatrix();
      sticks.setMatrixAt(i,dummy.matrix);
      sticks.setColorAt(i,shade.set(color));
    });
    this.meshes.add(sticks,this.points);
    if(cut!==null) {
      const sheet=new T.Mesh(new T.PlaneGeometry(22,22),
        new T.MeshBasicMaterial({color:ink.rose,opacity:.18,transparent:true,side:T.DoubleSide,depthWrite:false}));
      sheet.rotation.y=Math.PI/2;sheet.position.x=cut;
      const edges=new T.LineSegments(new T.EdgesGeometry(sheet.geometry),new T.LineBasicMaterial({color:ink.rose}));
      edges.rotation.y=sheet.rotation.y;edges.position.x=cut;
      this.meshes.add(sheet,edges);
    }
    this.frame=frame;
    this.marker.visible=false;
    this.measurement.visible=false;
    if(reset)this.fit();else this.draw();
  }
}

export function scaffold(lo,hi,style) {
  const points=[],links=[];
  const center=(lo+hi)/2;
  const x=i=>2.45*(i-center);
  for(let i=lo;i<hi;i++){
    const theta=(i-lo)*Math.PI/5;
    const color=style(i);
    const p=v(x(i+.5),6.3*Math.sin(theta),6.3*Math.cos(theta));
    const q=v(x(i+.5),6.3*Math.sin(theta+2.65),6.3*Math.cos(theta+2.65));
    const group=color===ink.muted?[color,color]:[color,ink.violet];
    points.push({p,r:.86,color:group[0],index:i,strand:0},{p:q,r:.86,color:group[1],index:i,strand:1});
    links.push({a:p,b:q,r:.21,color});
    if(i>lo)for(const [offset,point] of [[4,p],[3,q]]){
      links.push({a:points.at(-offset).p,b:point,r:.34,color:group[offset===4?0:1]});
    }
  }
  return {points,links,x,frame:(hi-lo)*1.31};
}

export function atomic(model,first,last,strand,representation) {
  const colors=[
    {C:'#76b4ff',N:'#9be1ff',O:'#f0b9b1',P:'#e0eeff'},
    {C:'#a998ff',N:'#d6bcff',O:'#f2b6cc',P:'#eedcff'}
  ];
  const byId=new Map(), points=[],links=[];
  const rangeCenter=(first+last-1)*3.38/2;
  model.points.forEach(([x,y,z,s,raw,element,atom],id)=>{
    const index=s?model.seq.length-1-raw:raw;
    if(index<first||index>=last||(strand!=='both'&&Number(strand)!==s))return;
    const r=representation==='space' ? {C:1.7,N:1.55,O:1.52,P:1.8}[element] : element==='P'?.62:.47;
    const point={p:v(z-rangeCenter,x,y),r,color:colors[s][element],index,strand:s,atom,element,xyz:[x,y,z]};
    byId.set(id,point);points.push(point);
  });
  if(representation!=='space'){
    for(const [a,b] of model.bonds){
      const p=byId.get(a),q=byId.get(b);
      if(p&&q)links.push({a:p.p,b:q.p,r:.19,color:p.strand?ink.violet:ink.blue});
    }
  }
  return {points,links,frame:Math.max(12,(last-first)*1.77)};
}

export {ink};
