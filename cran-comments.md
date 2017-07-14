## Release Summary

This is a resubmission to deal with the following issues:

* The Description field in the DESCRIPTION file has been updated to better
distinguish it from existing packages.

* In a vignette, a link to a URL has been modified to fix a message about it 
being possibly invalid.

## Test environments

* ubuntu 12.04 (on travis-ci), R 3.4.0
* local Antegros Linux, R 3.4.1
* local Windows 10 build, R 3.4.1
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'James Melville <jlmelville@gmail.com>'

This is the first version of the package and my first submission to CRAN.

There was a message about possibly mis-spelled words in DESCRIPTION:

  BFGS (8:64, 9:20, 9:28)
  Broyden (8:30)
  Goldfarb (8:47)
  Shanno (8:56)

Those words are spelled correctly.

## Downstream dependencies

None.
