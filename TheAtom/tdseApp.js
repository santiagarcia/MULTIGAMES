// tdseApp.js
// App wiring: UI, rendering, perf panel, periodic table, and interaction with TDSECore.
'use strict';

import { TDSECore } from './tdseCore.js';

const HARTREE = 27.2114; // eV

// DOM elements
const canvas = document.getElementById('canvas');
const hud = document.getElementById('hud');
const ctx = canvas.getContext('2d');

// Controls
const btnPause = document.getElementById('btnPause');
const btnReset = document.getElementById('btnReset');
const rngFade = document.getElementById('rngFade');
const rngNuc = document.getElementById('rngNuc');
const btnPickEl = document.getElementById('btnPickEl');
const curElement = document.getElementById('curElement');
const curZ = document.getElementById('curZ');
const rngTN = document.getElementById('rngTN');
const rngTL = document.getElementById('rngTL');
const rngTS = document.getElementById('rngTS');
const btnTDSEReset = document.getElementById('btnTDSEReset');
const btnBuildGS = document.getElementById('btnBuildGS');
const btnLight = document.getElementById('btnLight');
const rngPE = document.getElementById('rngPE');
const rngE0 = document.getElementById('rngE0');
const selGauge = document.getElementById('selGauge');
const selMap = document.getElementById('selMap');
// Overlays
const ptableOverlay = document.getElementById('ptableOverlay');
const ptableDiv = document.getElementById('ptable');
const btnClosePT = document.getElementById('btnClosePT');
// Accuracy panel
const accNorm = document.getElementById('accNorm');
const accEnergy = document.getElementById('accEnergy');
const accEmitted = document.getElementById('accEmitted');
const accDx = document.getElementById('accDx');
const accDy = document.getElementById('accDy');
// Perf panel
const perfFrame = document.getElementById('perfFrame');
const perfSolve = document.getElementById('perfSolve');
const perfGrid = document.getElementById('perfGrid');
const perfDt = document.getElementById('perfDt');
const perfGC = document.getElementById('perfGC');
// Value labels
const valFade = document.getElementById('valFade');
const valNuc = document.getElementById('valNuc');
const valTN = document.getElementById('valTN');
const valTL = document.getElementById('valTL');
const valTS = document.getElementById('valTS');
const valPE = document.getElementById('valPE');
const valE0 = document.getElementById('valE0');
const rngEta = document.getElementById('rngEta');
const rngM = document.getElementById('rngM');
const rngRcap = document.getElementById('rngRcap');
const rngEps = document.getElementById('rngEps');
const rngDt = document.getElementById('rngDt');
const valEta = document.getElementById('valEta');
const valM = document.getElementById('valM');
const valRcap = document.getElementById('valRcap');
const valEps = document.getElementById('valEps');
const valDt = document.getElementById('valDt');

// App state
let running = true;
let fpsEMA=0, lastT=performance.now();
let gcCount=0, lastFrameDelta=0;
let orbitYaw = 0, orbitPitch=0, zoom=1.0;
let W=0,H=0,CX=0,CY=0;
let off, ofx; // offscreen for psi^2 image
let imgData; // reused ImageData
const spark = document.getElementById('spark');
const spx = spark?.getContext('2d');
const sparkN = 280;
const ringDip = new Float32Array(sparkN);
const ringAcc = new Float32Array(sparkN);
let ringIdx=0;

// Protons visualization
let protonPts = [];

// Core solver
let core = new TDSECore({ N:+rngTN.value|0, L:+rngTL.value, Z:1, eps:+(rngEps?.value ?? 0.3), dt: 0.001 });

function resize(){
  const dpr = Math.min(devicePixelRatio||1, 2);
  W = canvas.clientWidth; H = canvas.clientHeight;
  canvas.width=W*dpr; canvas.height=H*dpr; ctx.setTransform(dpr,0,0,dpr,0,0);
  CX=W/2; CY=H/2;
}

function uiSync(){
  valFade.textContent = (+rngFade.value).toFixed(2);
  valNuc.textContent = String(+rngNuc.value|0);
  valTN.textContent = String(+rngTN.value|0);
  valTL.textContent = String(+rngTL.value|0);
  valTS.textContent = (+rngTS.value).toFixed(2)+'×';
  valPE.textContent = (+rngPE.value).toFixed(1)+' eV';
  valE0.textContent = (+rngE0.value).toFixed(3);
  if(valEta) valEta.textContent = (+rngEta.value).toFixed(2);
  if(valM) valM.textContent = String(+rngM.value|0);
  if(valRcap) valRcap.textContent = String(+rngRcap.value|0);
  if(valEps) valEps.textContent = (+rngEps.value).toFixed(2);
  if(valDt) valDt.textContent = (+rngDt.value).toFixed(4);
  buildProtons();
}

function nucleusGradient(){
  const r = +rngNuc.value;
  const g = ctx.createRadialGradient(CX, CY, 0, CX, CY, r);
  g.addColorStop(0, 'rgba(124,192,255,0.95)');
  g.addColorStop(1, 'rgba(124,192,255,0.05)');
  return g;
}

function buildProtons(){
  const Z = core.Z|0; if(Z<=0){ protonPts=[]; return; }
  const rN = +rngNuc.value * 0.6;
  const pts=[]; let placed=0,k=0; const maxRings=8;
  while(placed<Z && k<=maxRings){
    if(k===0){ pts.push([0,0]); placed++; k++; continue; }
    const ringCount = Math.min(6*k, Z-placed);
    for(let i=0;i<ringCount;i++){
      const ang = (i/ringCount)*Math.PI*2;
      const rr = (k/maxRings)*rN;
      pts.push([rr*Math.cos(ang), rr*Math.sin(ang)]);
    }
    placed += ringCount; k++;
  }
  protonPts = pts;
}

function mapColor(t){
  const mode = selMap.value;
  if(mode==='gray') { const v=Math.round(255*Math.max(0,Math.min(1,t))); return [v,v,v]; }
  if(mode==='viridis') return viridisColor(t);
  return turboColor(t);
}
function turboColor(t){
  t = t<0?0:(t>1?1:t);
  const r = Math.round(255*(0.135 - 0.616*t + 2.300*t*t - 2.349*t*t*t + 0.779*t*t*t*t));
  const g = Math.round(255*(0.091 - 0.279*t + 0.899*t*t - 0.818*t*t*t + 0.254*t*t*t*t + 0.003*t*t*t*t*t));
  const b = Math.round(255*(0.106 + 1.116*t - 1.975*t*t + 1.437*t*t*t - 0.330*t*t*t*t));
  return [Math.max(0,Math.min(255,r)), Math.max(0,Math.min(255,g)), Math.max(0,Math.min(255,b))];
}
function viridisColor(t){
  const stops = [
    [0.267,0.004,0.329], [0.283,0.141,0.458], [0.254,0.265,0.530],
    [0.207,0.372,0.553], [0.164,0.471,0.558], [0.133,0.658,0.517], [0.477,0.821,0.318]
  ];
  t = t<0?0:(t>1?1:t);
  const p=t*(stops.length-1); const i=Math.floor(p); const f=p-i; const a=stops[i]; const b=stops[Math.min(stops.length-1,i+1)];
  return [Math.round(255*(a[0]+(b[0]-a[0])*f)), Math.round(255*(a[1]+(b[1]-a[1])*f)), Math.round(255*(a[2]+(b[2]-a[2])*f))];
}

function draw(clear){
  const fade = +rngFade.value;
  if(clear){ ctx.clearRect(0,0,W,H); }
  else { ctx.fillStyle = `rgba(0,0,0,${fade})`; ctx.fillRect(0,0,W,H); }

  // density image
  const N=core.N; if(!off){ off=document.createElement('canvas'); ofx=off.getContext('2d'); }
  off.width=N; off.height=N;
  if(!imgData || imgData.width!==N || imgData.height!==N){ imgData = ofx.createImageData(N,N); }
  const img = imgData;
  let maxv=0; const re=core.re, im=core.im;
  for(let k=0;k<N*N;k++){ const v=re[k]*re[k]+im[k]*im[k]; if(v>maxv) maxv=v; }
  const invMax = maxv>0? 1/maxv:1;
  const data=img.data;
  for(let j=0;j<N;j++){
    for(let i=0;i<N;i++){
      const idx=j*N+i; const v=(re[idx]*re[idx]+im[idx]*im[idx])*invMax; const [R,G,B]=mapColor(1-Math.exp(-2*v));
      const p=idx*4; data[p]=R; data[p+1]=G; data[p+2]=B; data[p+3]=255;
    }
  }
  ofx.putImageData(img,0,0);

  ctx.save();
  ctx.translate(CX,CY);
  ctx.rotate(orbitYaw);
  const yScale=Math.max(0.2, Math.cos(orbitPitch));
  ctx.scale(zoom, zoom*yScale);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(off, -W/2, -H/2, W, H);
  // protons
  ctx.globalCompositeOperation='lighter';
  for(const p of protonPts){
    const [x,y]=p; const r=3;
    const grad=ctx.createRadialGradient(x,y,0,x,y,r*3); grad.addColorStop(0,'rgba(255,120,120,0.9)'); grad.addColorStop(1,'rgba(255,120,120,0)');
    ctx.fillStyle=grad; ctx.beginPath(); ctx.arc(x,y,r*2,0,Math.PI*2); ctx.fill();
    ctx.fillStyle='rgba(255,180,180,0.95)'; ctx.beginPath(); ctx.arc(x,y,r,0,Math.PI*2); ctx.fill();
  }
  ctx.restore();

  // nucleus glow
  ctx.save(); ctx.globalCompositeOperation='lighter'; ctx.fillStyle=nucleusGradient(); ctx.beginPath(); ctx.arc(CX,CY,+rngNuc.value,0,Math.PI*2); ctx.fill(); ctx.restore();

  // sparkline
  if(spx){
    spx.clearRect(0,0,spark.width,spark.height);
    // normalize recent window
    let minD=1e9,maxD=-1e9,minA=1e9,maxA=-1e9;
    for(let i=0;i<sparkN;i++){ const d=ringDip[i], a=ringAcc[i]; if(d<minD)minD=d; if(d>maxD)maxD=d; if(a<minA)minA=a; if(a>maxA)maxA=a; }
    const h=spark.height, w=spark.width;
    const yD=(v)=>{ return h - (v-minD)/Math.max(1e-6,(maxD-minD)) * h; };
    const yA=(v)=>{ return h - (v-minA)/Math.max(1e-6,(maxA-minA)) * h; };
    spx.lineWidth=1;
    // dipole in cyan
    spx.strokeStyle='#7cc0ff'; spx.beginPath();
    for(let i=0;i<sparkN;i++){ const x=i+0.5; const idx=(ringIdx+i)%sparkN; const y=yD(ringDip[idx]); if(i===0) spx.moveTo(x,y); else spx.lineTo(x,y); }
    spx.stroke();
    // accel in orange
    spx.strokeStyle='#ffb36b'; spx.beginPath();
    for(let i=0;i<sparkN;i++){ const x=i+0.5; const idx=(ringIdx+i)%sparkN; const y=yA(ringAcc[idx]); if(i===0) spx.moveTo(x,y); else spx.lineTo(x,y); }
    spx.stroke();
  }
}

function step(dt){
  // throttle substeps to keep solver time under ~12ms
  const t0=performance.now();
  const speed = +rngTS.value;
  let tRemain = dt*speed;
  const base = core.dt;
  let steps=0; let solveMs=0;
  while(tRemain>0){
    const d = Math.min(base, tRemain);
    // field
    const Eph = +rngPE.value; const omega = Eph/HARTREE; const E0=+rngE0.value;
    const Et = core.light.on ? (E0*Math.cos(omega*core.t)) : 0;
    const gauge = selGauge.value; // length uses E(t), velocity uses A(t) where A=-E/omega*sin(omega t)
    const field = (gauge==='length' || omega<1e-6) ? Et : (-E0/omega)*Math.sin(omega*core.t);
    const s0=performance.now();
    core.step(d, field, gauge);
    const s1=performance.now();
    solveMs += (s1-s0);
    tRemain -= d; steps++;
    if((performance.now()-t0)>14) break; // backpressure
  }
  perfSolve.textContent = `${solveMs.toFixed(2)} ms`;
  perfGrid.textContent = String(core.N);
  perfDt.textContent = core.dt.toFixed(4);
}

function tick(t){
  const dt = Math.min(0.05, (t-lastT)/1000); lastT=t;
  const f0=performance.now();
  if(running) step(dt);
  draw(false);
  const f1=performance.now();
  const frameMs=(f1-f0);
  // rough GC detection: spikes larger than 50ms vs previous delta
  if(frameMs - lastFrameDelta > 50) gcCount++;
  lastFrameDelta = frameMs;
  perfFrame.textContent = `${frameMs.toFixed(2)} ms`;
  fpsEMA = fpsEMA*0.9 + (1000/frameMs)*0.1;
  const I = 13.6 * core.Z * core.Z;
  const Eph = +rngPE.value;
  hud.textContent = `fps ${fpsEMA.toFixed(0)} | TDSE N ${core.N} | Z ${core.Z} | Eph ${Eph.toFixed(1)}eV vs I ${I.toFixed(1)}eV | P_out ${(core.emitted*100).toFixed(2)}%`;
  accNorm.textContent = core.norm.toFixed(6);
  accEnergy.textContent = core.energy.toFixed(3);
  accEmitted.textContent = (core.emitted*100).toFixed(2)+'%';
  if(accDx) accDx.textContent = (core.dipole?.x ?? 0).toFixed(3);
  if(accDy) accDy.textContent = (core.dipole?.y ?? 0).toFixed(3);
  if(perfGC) perfGC.textContent = String(gcCount);
  // push dipole and accel traces
  ringDip[ringIdx] = core.dipole.x;
  ringAcc[ringIdx] = core.acc?.x ?? 0;
  ringIdx = (ringIdx+1)%sparkN;
  requestAnimationFrame(tick);
}

// Events
btnPause.addEventListener('click', ()=>{ running=!running; btnPause.textContent = running? 'Pause':'Play'; btnPause.classList.toggle('primary', running); });
btnReset.addEventListener('click', ()=>{ core._initGuess(); core.normalize(); core.t=0; core.emitted=0; });
rngFade.addEventListener('input', uiSync);
rngNuc.addEventListener('input', uiSync);
selMap.addEventListener('change', uiSync);
rngTN.addEventListener('input', ()=>{ uiSync(); core.regrid({N:+rngTN.value|0}); buildProtons(); });
rngTL.addEventListener('input', ()=>{ uiSync(); core.regrid({L:+rngTL.value}); buildProtons(); });
rngTS.addEventListener('input', ()=>{ uiSync(); });
btnTDSEReset.addEventListener('click', ()=>{ core._initGuess(); core.normalize(); });
btnBuildGS.addEventListener('click', async ()=>{
  btnBuildGS.disabled=true; btnBuildGS.textContent='Building…';
  await core.buildGroundState(600, 0.002, (p)=>{ accEnergy.textContent = p.energy.toFixed(3); accNorm.textContent = p.norm.toFixed(6); });
  btnBuildGS.disabled=false; btnBuildGS.textContent='Build ground state';
});
btnLight.addEventListener('click', ()=>{ core.light.on = !core.light.on; btnLight.textContent = core.light.on? 'Light: On' : 'Light: Off'; });
btnPickEl.addEventListener('click', ()=>{ buildPTable(); ptableOverlay.style.display='grid'; });
btnClosePT.addEventListener('click', ()=>{ ptableOverlay.style.display='none'; });

if(rngEta) rngEta.addEventListener('input', ()=>{ valEta.textContent=(+rngEta.value).toFixed(2); core.cap.eta=+rngEta.value; });
if(rngM) rngM.addEventListener('input', ()=>{ valM.textContent=String(+rngM.value|0); core.cap.m=+rngM.value|0; });
if(rngRcap) rngRcap.addEventListener('input', ()=>{ valRcap.textContent=String(+rngRcap.value|0); core.cap.R=+rngRcap.value; core._buildCAP(); });
if(rngEps) rngEps.addEventListener('input', ()=>{ valEps.textContent=(+rngEps.value).toFixed(2); core.eps=+rngEps.value; core._buildPotential(); });
if(rngDt) rngDt.addEventListener('input', ()=>{ valDt.textContent=(+rngDt.value).toFixed(4); core.dt = +rngDt.value; });

// Orbit controls (right-click) and wheel zoom
canvas.addEventListener('contextmenu', (e)=> e.preventDefault());
let dragging=false,lastX=0,lastY=0;
canvas.addEventListener('pointerdown', (e)=>{ if(e.button===2){ dragging=true; lastX=e.clientX; lastY=e.clientY; canvas.setPointerCapture(e.pointerId);} });
canvas.addEventListener('pointermove', (e)=>{ if(dragging){ const dx=e.clientX-lastX, dy=e.clientY-lastY; orbitYaw += dx*0.005; orbitPitch = Math.max(-1.2, Math.min(1.2, orbitPitch + dy*0.005)); lastX=e.clientX; lastY=e.clientY; }});
canvas.addEventListener('pointerup', (e)=>{ if(e.button===2){ dragging=false; try{ canvas.releasePointerCapture(e.pointerId);}catch{} } });
canvas.addEventListener('wheel', (e)=>{ const z=zoom*Math.exp(-e.deltaY*0.0015); zoom=Math.max(0.2, Math.min(5,z)); e.preventDefault(); }, {passive:false});

// Periodic table
const PTABLE_1_36 = [
  {Z:1,s:'H',name:'Hydrogen',g:1,p:1},{Z:2,s:'He',name:'Helium',g:18,p:1},
  {Z:3,s:'Li',name:'Lithium',g:1,p:2},{Z:4,s:'Be',name:'Beryllium',g:2,p:2},{Z:5,s:'B',name:'Boron',g:13,p:2},{Z:6,s:'C',name:'Carbon',g:14,p:2},{Z:7,s:'N',name:'Nitrogen',g:15,p:2},{Z:8,s:'O',name:'Oxygen',g:16,p:2},{Z:9,s:'F',name:'Fluorine',g:17,p:2},{Z:10,s:'Ne',name:'Neon',g:18,p:2},
  {Z:11,s:'Na',name:'Sodium',g:1,p:3},{Z:12,s:'Mg',name:'Magnesium',g:2,p:3},{Z:13,s:'Al',name:'Aluminium',g:13,p:3},{Z:14,s:'Si',name:'Silicon',g:14,p:3},{Z:15,s:'P',name:'Phosphorus',g:15,p:3},{Z:16,s:'S',name:'Sulfur',g:16,p:3},{Z:17,s:'Cl',name:'Chlorine',g:17,p:3},{Z:18,s:'Ar',name:'Argon',g:18,p:3},
  {Z:19,s:'K',name:'Potassium',g:1,p:4},{Z:20,s:'Ca',name:'Calcium',g:2,p:4},{Z:21,s:'Sc',name:'Scandium',g:3,p:4},{Z:22,s:'Ti',name:'Titanium',g:4,p:4},{Z:23,s:'V',name:'Vanadium',g:5,p:4},{Z:24,s:'Cr',name:'Chromium',g:6,p:4},{Z:25,s:'Mn',name:'Manganese',g:7,p:4},{Z:26,s:'Fe',name:'Iron',g:8,p:4},{Z:27,s:'Co',name:'Cobalt',g:9,p:4},{Z:28,s:'Ni',name:'Nickel',g:10,p:4},{Z:29,s:'Cu',name:'Copper',g:11,p:4},{Z:30,s:'Zn',name:'Zinc',g:12,p:4},{Z:31,s:'Ga',name:'Gallium',g:13,p:4},{Z:32,s:'Ge',name:'Germanium',g:14,p:4},{Z:33,s:'As',name:'Arsenic',g:15,p:4},{Z:34,s:'Se',name:'Selenium',g:16,p:4},{Z:35,s:'Br',name:'Bromine',g:17,p:4},{Z:36,s:'Kr',name:'Krypton',g:18,p:4}
];
function buildPTable(){
  ptableDiv.innerHTML='';
  ptableDiv.style.gridTemplateColumns='repeat(18, 34px)';
  PTABLE_1_36.forEach(el=>{
    const d=document.createElement('div');
    d.className='el'; d.title=`${el.Z} ${el.name}`; d.textContent=el.s; d.style.gridColumn=String(el.g); d.style.gridRow=String(el.p);
    d.addEventListener('click', ()=>{
      core.Z=el.Z; curElement.textContent=el.s; curZ.textContent=String(el.Z); ptableOverlay.style.display='none';
      core._buildPotential(); core._initGuess(); core.normalize(); core.emitted=0; buildProtons();
    });
    ptableDiv.appendChild(d);
  });
}

// Bootstrap
function bootstrap(){
  resize(); uiSync(); buildPTable(); draw(true);
  requestAnimationFrame(tick);
}

window.addEventListener('resize', ()=>{ resize(); }, {passive:true});

bootstrap();
