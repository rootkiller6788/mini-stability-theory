# Course Alignment -- Circle & Popov Criteria
## MIT 6.241J
- **Circle criterion**: Nyquist must not intersect critical circle D(alpha,beta)
- **Popov criterion**: Nyquist lies to right of Popov line through -1/k+j*0
- Connection: popov_line_test checks Popov condition
## Stanford EE363
- **LMI formulation**: KYP lemma converts frequency conditions to LMIs
- **S-procedure**: Lossless for single nonlinearity
## Key Formulas
Popov: Re[(1+jw*eta)*G(jw)] + 1/k > 0
Critical circle: center=-(1/alpha+1/beta)/2, radius=|1/alpha-1/beta|/2
## Exercises
1. Plot Popov line for a nonlinear system
2. Compute critical circle from sector bounds