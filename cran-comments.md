## Resubmission
This is a resubmission. In this version I have addressed all comments from CRAN
review as follows:

* Removed "This package", "Functions for", package name, title or similar" from
DESCRIPTION
* Added doi describing the methods in DESCRIPTION
* Added missing \value to .Rd files regarding exported methods 
* Replaced all dontrun{} with donttest{}. In addition, I've provided toy 
examples that run under 5 seconds, and kept donttest{} for the rest
* Changed all messages printed to console from print() to message()
* When par changed in function I've added an immediate call of on.exit()
* Removed setting of seed to specific number within functions


## Test environments
* local OS Sonoma v.14.3.1, R 4.5.1
* local Windows 10 Home, R 4.3.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs. 


## Win-builder test (devel and release)
There were no ERRORs or WARNINGS. 

There was 1 NOTE:

### Checking CRAN incoming feasibility

* New submission
  * This is a new release to CRAN
* Possibly misspelled words in DESCRIPTION: Attentional (2:45), DSTP (11:6), 
SSP (11:37), attentional (11:52)
  * These have all been checked and are correct.
