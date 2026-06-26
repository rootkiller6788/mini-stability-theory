# Course Alignment -- LaSalle Invariance Principle
## MIT 6.241J
- **LaSalle (1960)**: Trajectories approach largest invariant set in {Vdot=0}
- **Barbashin-Krasovskii (1952)**: GAS requires V radially unbounded
- Connection: lasalle_invariant_set computes maximal invariant set
## Stanford AA203
- **Region of attraction**: Largest sublevel set in domain
- **Barrier certificates**: SOS for ROA estimation
## Key Formulas
LaSalle: Vdot<=0 => x(t)->M (largest invariant set in {Vdot=0})
BK: V unbounded + Vdot<0 => GAS
RA = {x0: lim phi(t,x0)=0 as t->inf}
## Exercises
1. Find maximal invariant set for 2D system
2. Estimate region of attraction via sublevel sets