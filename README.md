# FCIQMC_utils
tiny tools for FCIQMC analysis

## requirements
pyblock

## draw.py
\>\>\> draw.draw(fname, estimator, rolling, cutoff)
\# draw step-energy curve for NECI, MNECI or DNECI outputs
  Args:
    fname: None or str. FCIMCStats or fciqmc_stats file.
    estimator: str. 'trial' for trial wf estimator or 'projE' for HF estimator
    rolling: int. rolling window to smooth step-energy curve
    cutoff: int. start step for drawing

\>\>\> draw.diagnostic(fname, cutoff)
\# draw numerator and denominator for trial and projection (HF) estimator
    
## get_etot.py
\>\>\> get_etot.etot(fname, cutoff)
\# perform block analysis and print energy/std err for NECI, MNECI or DNECI
