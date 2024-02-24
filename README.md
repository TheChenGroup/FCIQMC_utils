# FCIQMC_utils
tiny tools for FCIQMC analysis

## requirements
pyblock

## draw.py
\>\>\> draw.draw(fname, estimator, rolling, cutoff) <br>
\# draw step-energy curve for NECI, MNECI or DNECI outputs <br>
  Args: <br>
    fname: None or str. FCIMCStats or fciqmc_stats file. <br>
    estimator: str. 'trial' for trial wf estimator or 'projE' for HF estimator <br>
    rolling: int. rolling window to smooth step-energy curve <br>
    cutoff: int. start step for drawing <br>

\>\>\> draw.diagnostic(fname, cutoff) <br>
\# draw numerator and denominator for trial and projection (HF) estimator <br>
    
## get_etot.py
\>\>\> get_etot.etot(fname, cutoff) <br>
\# perform block analysis and print energy/std err for NECI, MNECI or DNECI <br>
