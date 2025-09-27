# The Atom — 2D TDSE (Hydrogenic soft-core)

This app simulates a 2D single-electron "model atom" using a soft-core Coulomb potential, with imaginary-time ground-state preparation, real-time propagation using a Crank–Nicolson ADI integrator, and a complex absorbing potential (CAP). It supports length and velocity gauges and reports observables useful for HHG studies.

## Equations

- Soft-core Coulomb: V(r) = -Z / sqrt(r^2 + eps^2)
- CAP: V -> V - i eta (max(0, (r - f R_box)/Rcap))^m
- Driving (length): V_L = -E(t) x; (velocity): A(t) couples via p_x. We use a small-field approximation to compare gauges.
- Propagator: CN-ADI (Douglas-Rachford split) alternating x/y implicit solves via Thomas algorithm.
- Ground state: Imaginary-time propagation (ITP) using the same ADI backbone with normalization after each step.

## Validation

- Norm conservation (no light, no CAP): <= 1e-6 over 1e4 steps.
- Ground-state energy monotone decreasing during ITP; matches discrete soft-core eigenvalue within ~1%.
- Ionization fraction ≈ norm lost to CAP within 2%.

## Performance

- Typed arrays reused, zero per-frame allocations in hot loops.
- Backpressure: substeps throttled to ~12–14 ms per frame.
- CN-ADI chosen for browser portability and stability; FFT split-operator would require WebGL/WebGPU FFT; optional in future.

## Files

- index.html — UI shell, loads modules and keeps previous aesthetics.
- style.css — Extracted styles.
- tdseCore.js — Numerics core.
- tdseApp.js — UI wiring, rendering, interaction.
- benchmarks.md — Fill with your machine/browser/FPS numbers.

## How to run

Open `index.html` in a modern browser. If modules are blocked due to file:// origin, serve the folder with a static server.
