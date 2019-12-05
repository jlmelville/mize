## Release Summary

This is a patch update to maintain compatibility with R-devel.

## Test environments

* ubuntu 16.04 (on travis-ci), R 3.5.3, R 3.6.1, R-devel
* ubuntu 16.04 (on rhub), R 3.6.1
* fedora 30 (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.5.3, R 3.6.1
* local Windows 10 build, R 3.6.1
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

## CRAN checks

There are two errors, for r-devel-linux-x86_64-debian-clang and 
r-devel-linux-x86_64-debian-gcc. This release is designed to expressly fix
these errors.

## Downstream dependencies

There are 2 downstream dependencies. 

* All completed R CMD CHECK without issues.
