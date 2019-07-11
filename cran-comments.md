## Release Summary

This is a patch update to fix a user-reported bug.

## Test environments

* ubuntu 14.04 (on travis-ci), R 3.4.4, R 3.6.0, R-devel
* ubuntu 16.04 (on rhub), R 3.6.1
* fedora 30 (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.5.3, R 3.6.1
* local Windows 10 build, R 3.5.3
* Windows Server 2008 (on rhub) R-devel
* Windows Server 2012 (on appveyor) R 3.6.1
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'James Melville <jlmelville@gmail.com>'

There was a message about possibly mis-spelled words in DESCRIPTION:

  BFGS (8:64, 9:20, 9:28)
  Broyden (8:30)
  Goldfarb (8:47)
  Shanno (8:56)

Those words are spelled correctly.

## Downstream dependencies

None.
