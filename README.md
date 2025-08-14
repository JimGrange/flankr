
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flankr v1.2.0 <a><img src='images/logo.png' align="right" height="250"/></a>

$\texttt{flankr}$ is an R package implementing computational models of
Eriksen flanker task performance.

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
#> * checking for file ‘/private/var/folders/t3/59rg8brd3rzc4ph50b6vbcwh0000gp/T/RtmpyDMdof/remotes9297396fa29f/JimGrange-flankr-0fe15fa/DESCRIPTION’ ... OK
#> * preparing ‘flankr’:
#> * checking DESCRIPTION meta-information ... OK
#> * cleaning src
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘flankr_1.1.0.tar.gz’
#> Installing package into '/private/var/folders/t3/59rg8brd3rzc4ph50b6vbcwh0000gp/T/RtmpwsUSUw/temp_libpath1330a205df2db'
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
