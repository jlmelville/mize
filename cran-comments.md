## Release Summary

This is a patch update for compatibility with R-devel, fixing errors which have
appeared in the CRAN checks.

## Test environments

* local ubuntu 25.10 R 4.5.2, devel
* ubuntu 24.04 with clang and gcc12 (on rhub), devel
* fedora 42 with gcc15 (on rhub), devel
* ubuntu 24.04 (on github actions), R 4.4.3, R 4.5.2, devel
* mac OS X Sequoia (on github actions) R 4.5.2
* local Windows 11 build, R 4.5.2, devel
* Windows Server 2022 (on github actions), R 4.4.3, R 4.5.2
* win-builder (devel)

## R CMD check results

There were 4 ERRORs and 3 NOTEs.

The ERRORs (for r-devel-linux-x86_64-debian-clang, 
R-devel-linux-x86_64-debian-gcc, r-devel-linux-x86_64-fedora-clang and	
r-devel-windows-x86_64) are:

> Running ‘testthat.R’ [10s/12s]
   Running the tests in ‘tests/testthat.R’ failed.

The package has been updated to prevent these failures.

The NOTEs (for r-oldrel-macos-arm64, r-oldrel-macos-x86_64 and 
r-oldrel-windows-x86_64) are: 

> checking LazyData ... NOTE
  'LazyData' is specified without a 'data' directory

The DESCRIPTION file has been updated to fix this.

## Downstream dependencies

There are 3 downstream dependencies. 

* All completed R CMD CHECK without issues. There were no new NOTEs, WARNINGs or
ERRORs.
