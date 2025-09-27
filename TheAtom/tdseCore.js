// tdseCore.js
// Numerics core for 2D soft-core hydrogenic TDSE with CAP, ITP, and CN-ADI time-stepping.
// Vanilla JS modules, strict mode, typed arrays, no per-frame allocations in hot paths.

'use strict';

/** Utility: ensure value is finite during dev. No-op in prod if desired. */
export function assertFinite(x, label) {
  if (!Number.isFinite(x)) throw new Error(`Non-finite ${label ?? ''}`);
}

/** Colormap helpers kept in app for rendering; core stays headless. */

/** TDSE core state and operations */
export class TDSECore {
  /**
   * Construct the solver core.
   * @param {object} opts - { N, L, Z, eps, cap: {eta,m,R}, gauge }
   */
  constructor(opts={}) {
    this.N = opts.N ?? 128;
    this.L = opts.L ?? 24; // domain length in a0
    this.dx = this.L / this.N;
    this.Z = opts.Z ?? 1;
    this.eps = opts.eps ?? 0.3; // soft-core epsilon (a0)
    this.gauge = opts.gauge ?? 'length'; // 'length' | 'velocity'

    this.cap = {
      eta: opts.cap?.eta ?? 0.5,
      m: opts.cap?.m ?? 2,
      R: opts.cap?.R ?? 8,
      startFrac: opts.cap?.startFrac ?? 0.7 // start CAP at startFrac*half-width in radius
    };

    // time stepping
    this.dt = opts.dt ?? 0.001; // a.u.
    this.light = { E0: 0.02, omega: 20/27.2114, on: false, A0: 0.0 }; // omega in a.u., A0 derived

    // allocate typed arrays
    const NN = this.N*this.N;
    this.re = new Float32Array(NN);
    this.im = new Float32Array(NN);
    this.V = new Float32Array(NN);
    this.capMask = new Float32Array(NN); // exp(-eta*(r/R)^m*dt) style factor per step equivalent, but we pre-store profile and apply per-step factor

    // CN-ADI buffers: scratch lines for tridiagonal solves
    this.lineR = new Float32Array(this.N);
    this.lineI = new Float32Array(this.N);
    this.a = new Float32Array(this.N);
    this.b = new Float32Array(this.N);
    this.c = new Float32Array(this.N);
    this.tmp = new Float32Array(this.N);
  // Complex coefficients for CN-ADI tri-diagonal solves
  this.aR = new Float32Array(this.N); this.aI = new Float32Array(this.N);
  this.bR = new Float32Array(this.N); this.bI = new Float32Array(this.N);
  this.cR = new Float32Array(this.N); this.cI = new Float32Array(this.N);
  this.tmpCR = new Float32Array(this.N); this.tmpCI = new Float32Array(this.N);

    // observables
    this.norm = 0;
    this.energy = 0;
    this.emitted = 0;
    this.t = 0;

    // precompute grids
    this._buildPotential();
    this._buildCAP();

    // init with simple 1s guess; app can call buildGroundState() for ITP
    this._initGuess();
    this._updateObservables();
  }

  /** Compute soft-core Coulomb potential and cache grid coordinates */
  _buildPotential(){
    const N=this.N, dx=this.dx, half=this.L/2, eps=this.eps, Z=this.Z;
    this.x = new Float32Array(N);
    for(let i=0;i<N;i++) this.x[i] = -half + (i+0.5)*dx;
    // Precompute potential and its gradients for acceleration observable
    this.dVdx = new Float32Array(N*N);
    this.dVdy = new Float32Array(N*N);
    for(let j=0;j<N;j++){
      const y = -half + (j+0.5)*dx;
      for(let i=0;i<N;i++){
        const x = this.x[i];
        const r2 = x*x + y*y;
        this.V[j*N+i] = -Z / Math.sqrt(r2 + eps*eps);
      }
    }
    // central differences for dV; handle edges with one-sided diffs
    for(let j=0;j<N;j++){
      const jm = j>0? j-1:j, jp = j<N-1? j+1:j;
      for(let i=0;i<N;i++){
        const im = i>0? i-1:i, ip = i<N-1? i+1:i;
        const idx=j*N+i;
        this.dVdx[idx] = (this.V[j*N+ip] - this.V[j*N+im])/(2*dx);
        this.dVdy[idx] = (this.V[jp*N+i] - this.V[jm*N+i])/(2*dx);
      }
    }
  }

  /** Build CAP radial profile; multiplied into wave each step as exp(-eta*(r/R)^m*dt) */
  _buildCAP(){
    const N=this.N, dx=this.dx, half=this.L/2;
    const {eta,m,R,startFrac} = this.cap;
    // We store the profile f(r)=(max(0,(r/Rcap))^m). Real factor per step: exp(-eta * f(r) * dt)
    this.capProfile = new Float32Array(N*N);
    const Rcap = R; // user in a0 units
    for(let j=0;j<N;j++){
      const y = -half + (j+0.5)*dx;
      for(let i=0;i<N;i++){
        const x = -half + (i+0.5)*dx;
        const r = Math.hypot(x,y);
        const rr = Math.max(0, r - startFrac*half);
        const val = Math.pow(Math.max(0, rr / Math.max(1e-6,Rcap)), m);
        this.capProfile[j*N+i] = val;
      }
    }
  }

  /** Simple 1s-like initial guess exp(-Z r) */
  _initGuess(){
    const N=this.N, dx=this.dx, half=this.L/2;
    for(let j=0;j<N;j++){
      const y = -half + (j+0.5)*dx;
      for(let i=0;i<N;i++){
        const x = -half + (i+0.5)*dx;
        const r = Math.hypot(x,y);
        this.re[j*N+i] = Math.exp(-this.Z*r);
        this.im[j*N+i] = 0;
      }
    }
    this.normalize();
  }

  /** Reset arrays when N/L/Z/eps change (keeps light/gauge/CAP params) */
  regrid({N=this.N,L=this.L,Z=this.Z,eps=this.eps,dt=this.dt}={}){
    this.N = N; this.L = L; this.dx = L/N; this.Z=Z; this.eps=eps; this.dt=dt;
    const NN = N*N;
    this.re = new Float32Array(NN);
    this.im = new Float32Array(NN);
    this.V = new Float32Array(NN);
    this.capMask = new Float32Array(NN);
    this.lineR = new Float32Array(N);
    this.lineI = new Float32Array(N);
    this.a = new Float32Array(N);
    this.b = new Float32Array(N);
    this.c = new Float32Array(N);
    this.tmp = new Float32Array(N);
  this.aR = new Float32Array(N); this.aI = new Float32Array(N);
  this.bR = new Float32Array(N); this.bI = new Float32Array(N);
  this.cR = new Float32Array(N); this.cI = new Float32Array(N);
  this.tmpCR = new Float32Array(N); this.tmpCI = new Float32Array(N);
    this._buildPotential();
    this._buildCAP();
    this._initGuess();
    this.t=0; this.emitted=0;
    this._updateObservables();
  }

  /** L2 normalize wavefunction to 1.0 */
  normalize(){
    const N=this.N, dx=this.dx, area=dx*dx; let s=0; const re=this.re, im=this.im;
    for(let k=0;k<N*N;k++){ const r=re[k], ii=im[k]; s += r*r + ii*ii; }
    const inv = 1/Math.sqrt(Math.max(1e-12, s*area));
    for(let k=0;k<N*N;k++){ re[k]*=inv; im[k]*=inv; }
  }

  /** Compute observables: norm, energy (discrete), and dipole components */
  _updateObservables(){
    const N=this.N, dx=this.dx, area=dx*dx, re=this.re, im=this.im, V=this.V;
    let norm=0, Ex=0, Ey=0, Epot=0, Ekin=0, ax=0, ay=0;
    // Kinetic via 5-point Laplacian inner product: <psi|(-1/2 nabla^2)|psi>
    // We avoid allocations by on-the-fly laplacian.
    for(let j=0;j<N;j++){
      const jm = j>0? j-1 : j;
      const jp = j<N-1? j+1 : j;
      for(let i=0;i<N;i++){
        const imn = i>0? i-1 : i;
        const ip = i<N-1? i+1 : i;
        const idx=j*N+i;
        const rr = re[idx], ii = im[idx];
        norm += rr*rr + ii*ii;
        const lapRe = (re[j*N+ip] + re[j*N+imn] + re[jp*N+i] + re[jm*N+i] - 4*rr);
        const lapIm = (im[j*N+ip] + im[j*N+imn] + im[jp*N+i] + im[jm*N+i] - 4*ii);
        Ekin += -0.5 * (rr*lapRe + ii*lapIm) / (dx*dx);
        Epot += V[idx]*(rr*rr + ii*ii);
        const x = -this.L/2 + (i+0.5)*dx;
        const y = -this.L/2 + (j+0.5)*dx;
        Ex += x*(rr*rr + ii*ii);
        Ey += y*(rr*rr + ii*ii);
        // acceleration from force: a = -<∇V> plus field term in length gauge
        ax += this.dVdx ? (this.dVdx[idx]*(rr*rr+ii*ii)) : 0;
        ay += this.dVdy ? (this.dVdy[idx]*(rr*rr+ii*ii)) : 0;
      }
    }
    this.norm = norm*area;
    this.energy = Ekin*area + Epot*area;
    this.dipole = {x:Ex*area, y:Ey*area};
    // store acceleration; field term added in step() based on lastField/gauge
    this._ax_base = -ax*area; this._ay_base = -ay*area;
  }

  /** Imaginary-time propagation to ground state using CN-ADI with dt-> -i*dt (i.e., diffusion). */
  async buildGroundState(steps=500, dtImag=0.002, onProgress){
    // For simplicity here: perform repeated diffusion + renorm using a stabilized scheme.
    // We can reuse our ADI line solves: replace i with -1.
    for(let s=0;s<steps;s++){
      this._stepCN_ADI( dtImag, /*imag*/true, /*field*/0, /*gauge*/'length');
      this.normalize();
      this._updateObservables();
      if(onProgress && (s%10===0)) onProgress({step:s, energy:this.energy, norm:this.norm});
    }
  }

  /** Advance real time one step with CN-ADI. Gauge can be 'length' or 'velocity'. */
  step(dt, field, gauge){
    this._stepCN_ADI(dt, /*imag*/false, field, gauge ?? this.gauge);
    // CAP: apply after unitary step; update emitted = norm loss
    const before = this._normNoArea();
    const f = this.capProfile; const re=this.re, im=this.im; const eta=this.cap.eta;
    for(let k=0;k<re.length;k++){
      const fac = Math.exp(-eta * f[k] * dt);
      re[k]*=fac; im[k]*=fac;
    }
    const after = this._normNoArea();
    this.emitted += Math.max(0, (before-after)*this.dx*this.dx);
    this.t += dt;
    this._updateObservables();
    this.lastField = field||0; this.lastGauge = gauge;
    const fieldAccel = (gauge==='length') ? (-this.lastField) : 0; // a_x includes -E(t) in length gauge
    this.acc = { x: this._ax_base + fieldAccel, y: this._ay_base };
  }

  _normNoArea(){ let s=0; const re=this.re, im=this.im; for(let k=0;k<re.length;k++){ const r=re[k], i=im[k]; s+=r*r+i*i; } return s; }

  /** Core CN-ADI step: split along x then y with tridiagonal solves. */
  _stepCN_ADI(dt, imag, field, gauge){
    const N=this.N, re=this.re, im=this.im;
    // First half-step implicit in x, explicit in y (V and T_y on RHS)
    for(let j=0;j<N;j++){
      // Assemble RHS and complex tri-diagonal for row j
      this._buildRowX(j, dt, imag, field, gauge);
      this._thomasSolveComplex(this.aR, this.aI, this.bR, this.bI, this.cR, this.cI, this.lineR, this.lineI, this.tmpCR, this.tmpCI);
      // write back half-step result to arrays (overwriting current row)
      for(let i=0;i<N;i++){ const idx=j*N+i; re[idx]=this.lineR[i]; im[idx]=this.lineI[i]; }
    }
    // Second half-step implicit in y, explicit in x (V and T_x on RHS)
    for(let i=0;i<N;i++){
      this._buildColY(i, dt, imag, field, gauge);
      this._thomasSolveComplex(this.aR, this.aI, this.bR, this.bI, this.cR, this.cI, this.lineR, this.lineI, this.tmpCR, this.tmpCI);
      for(let j=0;j<N;j++){ const idx=j*N+i; re[idx]=this.lineR[j]; im[idx]=this.lineI[j]; }
    }
  }

  /** Assemble tridiagonal system for one row along x */
  _buildRowX(j, dt, imag, field, gauge){
    const N=this.N, dx=this.dx, re=this.re, im=this.im, V=this.V;
    const r=this.lineR, ii=this.lineI;
    // RHS seed: current row
    for(let i=0;i<N;i++){ const idx=j*N+i; r[i]=re[idx]; ii[i]=im[idx]; }
    // Explicit Y kinetic + potential + length-gauge potential term
    const jm = j>0? j-1:j, jp=j<N-1? j+1:j;
    const kfac = -0.5/(dx*dx);
    const E = (gauge==='length') ? (field||0) : 0;
    for(let i=0;i<N;i++){
      const idx=j*N+i;
      const lapYRe = (re[jp*N+i] - 2*re[idx] + re[jm*N+i]);
      const lapYIm = (im[jp*N+i] - 2*im[idx] + im[jm*N+i]);
      const x = -this.L/2 + (i+0.5)*dx;
      const Hr = kfac*lapYRe + V[idx]*re[idx] + (E ? (-E*x*re[idx]) : 0);
      const Hi = kfac*lapYIm + V[idx]*im[idx] + (E ? (-E*x*im[idx]) : 0);
      const ar = r[i], ai=ii[i];
      // (1 - i dt/2 H_y) psi: apply as real/imag mix
      if(!imag){
        r[i]  = ar + (dt*0.5)*Hi;
        ii[i] = ai - (dt*0.5)*Hr;
      }else{
        // imaginary-time: (1 - dt/2 H_y)
        r[i]  = ar - (dt*0.5)*Hr;
        ii[i] = ai - (dt*0.5)*Hi;
      }
    }
    // Velocity gauge term i A p_x on RHS (A along x)
    if(gauge==='velocity' && !imag){
      const A = field||0;
      for(let i=1;i<N-1;i++){
        const dRe = (re[j*N+(i+1)] - re[j*N+(i-1)])/(2*dx);
        const dIm = (im[j*N+(i+1)] - im[j*N+(i-1)])/(2*dx);
        const ar = r[i], ai=ii[i];
        // -i A p_x psi => add A d/dx terms with coupling
        r[i]  = ar - (0.5*dt)*(-A)*dIm;
        ii[i] = ai + (0.5*dt)*(-A)*dRe;
      }
    }
    // Build complex tri-diagonal for implicit x: (I + s T_x), s = i dt/2 (unitary) or s = dt/2 (imag)
    const gamma = dt/(4*dx*dx);
    if(!imag){
      // unitary: coefficients purely imaginary
      for(let i=0;i<N;i++){ this.aR[i]=0; this.bR[i]=1; this.cR[i]=0; this.aI[i]=-gamma; this.bI[i]=2*gamma; this.cI[i]=-gamma; }
      // Neumann boundaries: double magnitude on the single neighbor
      this.aI[0]=0; this.cI[0]=-2*gamma; this.bI[0]=2*gamma;
      this.cI[N-1]=0; this.aI[N-1]=-2*gamma; this.bI[N-1]=2*gamma;
    }else{
      // imaginary-time: coefficients purely real
      for(let i=0;i<N;i++){ this.aR[i]=-gamma; this.bR[i]=1+2*gamma; this.cR[i]=-gamma; this.aI[i]=0; this.bI[i]=0; this.cI[i]=0; }
      this.aR[0]=0; this.cR[0]=-2*gamma; this.bR[0]=1+2*gamma;
      this.cR[N-1]=0; this.aR[N-1]=-2*gamma; this.bR[N-1]=1+2*gamma;
    }
  }

  /** Assemble tridiagonal system for one column along y */
  _buildColY(i, dt, imag, field, gauge){
    const N=this.N, dx=this.dx, re=this.re, im=this.im, V=this.V;
    const r=this.lineR, ii=this.lineI;
    // RHS seed: current column
    for(let j=0;j<N;j++){ const idx=j*N+i; r[j]=re[idx]; ii[j]=im[idx]; }
    // Explicit X kinetic + potential + length-gauge term
    const kfac = -0.5/(dx*dx);
    const E = (gauge==='length') ? (field||0) : 0;
    for(let j=0;j<N;j++){
      const jm = j>0? j-1:j, jp=j<N-1? j+1:j;
      const idx=j*N+i;
      const ip = Math.min(N-1,i+1), imn = Math.max(0,i-1);
      const lapXRe = (re[j*N+ip] - 2*re[idx] + re[j*N+imn]);
      const lapXIm = (im[j*N+ip] - 2*im[idx] + im[j*N+imn]);
      const x = -this.L/2 + (i+0.5)*dx;
      const Hr = kfac*lapXRe + V[idx]*re[idx] + (E ? (-E*x*re[idx]) : 0);
      const Hi = kfac*lapXIm + V[idx]*im[idx] + (E ? (-E*x*im[idx]) : 0);
      const ar=r[j], ai=ii[j];
      if(!imag){
        r[j]  = ar + (0.5*dt)*Hi;
        ii[j] = ai - (0.5*dt)*Hr;
      }else{
        r[j]  = ar - (0.5*dt)*Hr;
        ii[j] = ai - (0.5*dt)*Hi;
      }
    }
    // No velocity-gauge term along y for A along x
    // Build complex tri-diagonal for implicit y
    const gamma = dt/(4*dx*dx);
    if(!imag){
      for(let j=0;j<N;j++){ this.aR[j]=0; this.bR[j]=1; this.cR[j]=0; this.aI[j]=-gamma; this.bI[j]=2*gamma; this.cI[j]=-gamma; }
      this.aI[0]=0; this.cI[0]=-2*gamma; this.bI[0]=2*gamma;
      this.cI[N-1]=0; this.aI[N-1]=-2*gamma; this.bI[N-1]=2*gamma;
    }else{
      for(let j=0;j<N;j++){ this.aR[j]=-gamma; this.bR[j]=1+2*gamma; this.cR[j]=-gamma; this.aI[j]=0; this.bI[j]=0; this.cI[j]=0; }
      this.aR[0]=0; this.cR[0]=-2*gamma; this.bR[0]=1+2*gamma;
      this.cR[N-1]=0; this.aR[N-1]=-2*gamma; this.bR[N-1]=1+2*gamma;
    }
  }

  /** Thomas algorithm solver (in-place on d) */
  _thomasSolve(a,b,c,d,tmp){
    const n=d.length;
    tmp[0] = c[0]/b[0];
    d[0] = d[0]/b[0];
    for(let i=1;i<n;i++){
      const denom = b[i] - a[i]*tmp[i-1];
      tmp[i] = (i<n-1) ? c[i]/denom : 0;
      d[i] = (d[i] - a[i]*d[i-1]) / denom;
    }
    for(let i=n-2;i>=0;i--){ d[i] = d[i] - tmp[i]*d[i+1]; }
  }

  /** Complex Thomas algorithm: solves (a x_{i-1} + b x_i + c x_{i+1} = d) for complex arrays */
  _thomasSolveComplex(aR,aI,bR,bI,cR,cI,dR,dI,tCR,tCI){
    const n=dR.length;
    // helper complex ops
    const inv = (xr,xi)=>{ const den=xr*xr+xi*xi || 1e-30; return [xr/den, -xi/den]; };
    const mul = (ar,ai,br,bi)=>[ar*br - ai*bi, ar*bi + ai*br];
    const sub = (ar,ai,br,bi)=>[ar-br, ai-bi];
    // Forward elimination: modify c and d in-place into tC' and d'
    // i=0: normalize
    let [br,bi] = [bR[0], bI[0]];
    let [ivr,ivi] = inv(br,bi);
    // c0' = c0 / b0
    let m = mul(cR[0],cI[0], ivr,ivi); tCR[0]=m[0]; tCI[0]=m[1];
    // d0' = d0 / b0
    m = mul(dR[0],dI[0], ivr,ivi); dR[0]=m[0]; dI[0]=m[1];
    for(let i=1;i<n;i++){
      // denom = b[i] - a[i]*c'[i-1]
      const ac = mul(aR[i],aI[i], tCR[i-1], tCI[i-1]);
      const denom = sub(bR[i],bI[i], ac[0], ac[1]);
      ;[ivr,ivi] = inv(denom[0], denom[1]);
      // c'[i] = c[i] / denom
      m = mul(cR[i],cI[i], ivr,ivi); tCR[i] = (i<n-1)? m[0]:0; tCI[i] = (i<n-1)? m[1]:0;
      // d'[i] = (d[i] - a[i]*d'[i-1]) / denom
      const ad = mul(aR[i],aI[i], dR[i-1], dI[i-1]);
      const num = sub(dR[i],dI[i], ad[0], ad[1]);
      m = mul(num[0],num[1], ivr,ivi); dR[i]=m[0]; dI[i]=m[1];
    }
    // Back substitution
    for(let i=n-2;i>=0;i--){
      // d[i] = d[i] - c'[i]*d[i+1]
      const cd = mul(tCR[i],tCI[i], dR[i+1], dI[i+1]);
      dR[i] = dR[i] - cd[0];
      dI[i] = dI[i] - cd[1];
    }
  }
}
