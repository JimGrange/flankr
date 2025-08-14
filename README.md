
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flankr v1.2.0

An R package implementing computational models of Eriksen flanker task
performance. Full details of how to use the package is available in the
following paper:

Grange, J.A. (2016). flankr: An R package for implementing computational
models of attentional selectivity. *Behavior Research Methods, 48,*
528–541.

- PDF Link:
  <https://link.springer.com/article/10.3758/s13428-015-0615-y>

## Updates for version 1.2.0

- 50% further efficiency in DSTP simulation speed
- 24% further efficiency in SSP simulation speed
- Note that the way random seeds are handled in both
  $\texttt{simulateDSTP}$ and $\texttt{simulateSSP}$ is slightly
  different to that in version 1.0.0 (initial release). Therefore, there
  may be very slight differences between simulation data (and therefore
  potentially very slight differences in best-fitting parameter values)
  between versions.
