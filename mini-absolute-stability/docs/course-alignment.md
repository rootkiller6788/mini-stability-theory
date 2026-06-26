# Course Alignment -- Absolute Stability (Lure Systems)
## MIT 6.241J: Dynamic Systems and Control
- **Lure problem**: Nonlinearity in sector [alpha,beta] with linear forward path G(s)
- **Circle criterion**: Re[(1+beta*G(jw))*(1+alpha*G(jw))^{-1}] > 0
- Connection: abs_stab_circle_test verifies circle condition
## Stanford AA203
- **Passivity**: Absolute stability via passivity of nonlinearity
- **KYP lemma**: LMI equivalent of frequency condition
## Key Formulas
Circle: Nyquist avoids critical circle centered at -(1/alpha+1/beta)/2
Popov: Re[(1+jw*eta)*G(jw)] + 1/k > 0 for all w>=0
## Exercises
1. Verify circle criterion for a Lure system
2. Compare circle vs. Popov criteria on benchmarks