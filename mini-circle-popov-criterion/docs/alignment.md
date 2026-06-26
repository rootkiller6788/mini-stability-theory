# Course Alignment

## Primary Courses

### MIT 6.241J — Dynamic Systems and Control
- Frequency-domain methods (Nyquist, Bode)
- Circle and Popov criteria for nonlinear stability
- KYP lemma and passivity theory
- Lyapunov stability theory

### Stanford AA203 — Optimal and Learning-based Control
- Nonlinear systems analysis
- Absolute stability and sector conditions
- LMI-based stability verification
- Multiplier methods

### Berkeley EE222 — Nonlinear Systems
- Lur'e systems and absolute stability
- Circle criterion (time-varying and time-invariant)
- Popov criterion and multiplier theory
- Passivity and small-gain theorems

### Caltech CDS110 — Introductory Control Theory
- Classical Nyquist stability criterion
- Gain and phase margins
- Introduction to nonlinear stability
- Sector-bounded nonlinearities

### ETH 227-0216 — Nonlinear Control Systems
- Absolute stability theory
- Circle and Popov criteria
- KYP lemma and positive real systems
- Lyapunov-based methods

## Implementation Coverage

| Topic | Implementation |
|-------|---------------|
| Transfer function evaluation |  — polynomial at s=jw |
| Nyquist plot generation |  — freq sweep |
| Critical disk computation |  |
| Winding number / encirclement |  |
| Circle criterion test |  |
| Half-plane case |  |
| Linear (point) case |  |
| Popov plot computation |  |
| Popov line test |  |
| Optimal multiplier search |  |
| Loop transformation |  |
| Criterion comparison |  |
| State-space to TF |  (Leverrier) |
| KYP LMI solver |  |
| Lyapunov equation solver |  |
| Passivity verification |  |
| Lur'e system construction |  + deep copy |
| Unified stability test |  |
| Gain/phase margins |  |
| Robust margin search |  |
| H-infinity norm |  |
| Sector verification |  |
| Sector transformations |  |
| Standard nonlinearities | saturation, deadzone, relay, quantizer, cubic |
| Passivity equivalence |  |
| Max sector search |  |
