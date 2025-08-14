
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flankr v1.2.0 <a><img src='images/logo.png' align="right" height="250"/></a>

$\texttt{flankr}$ is an R package implementing computational models of
Eriksen flanker task performance. The package allows simulation of the
models as well as fitting the models to participant data. Additional
utility functions allow plotting of the best-fitting model parameters
against observed data, as well as providing Bayesian Information
Criterion values for model competition.

Current models implemented in $\texttt{flankr}$ are: - The Shrinking
Spotlight Model (SSP) of White et al. (2011) - The Dual-Stage Two-Phase
Model (DSTP) of Hübner et al. (2010)

## Installation

The development version can be installed from
[GitHub](https://github.com/) with:

``` r
require(devtools)
#> Loading required package: devtools
#> Loading required package: usethis
devtools::install_github("JimGrange/flankr")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo JimGrange/flankr@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/t3/59rg8brd3rzc4ph50b6vbcwh0000gp/T/Rtmp9BISlS/remotes95c952cc77b4/JimGrange-flankr-0fe15fa/DESCRIPTION’ ... OK
#> * preparing ‘flankr’:
#> * checking DESCRIPTION meta-information ... OK
#> * cleaning src
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘flankr_1.1.0.tar.gz’
#> Installing package into '/private/var/folders/t3/59rg8brd3rzc4ph50b6vbcwh0000gp/T/RtmpwsUSUw/temp_libpath1330a4bf794d'
#> (as 'lib' is unspecified)
```

## User guide

Full details of how to use the package is available in the following
paper:

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

## References

Hübner, R., Steinhauser, M., & Lehle, C. (2010). A dual-stage two-phase
model of selective attention. *Psychological Review, 117(3)*, 759–784.
<https://doi.org/10.1037/a0019471>

White, C. N., Ratcliff, R., & Starns, J. S. (2011). Diffusion models of
the flanker task: Discrete versus gradual attentional selection.
*Cognitive Psychology, 63(4)*, 210–238.
<https://doi.org/10.1016/j.cogpsych.2011.08.001>
