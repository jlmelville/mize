## Release Summary

This is a patch update to fix a minor bug.

## Test environments

* ubuntu 16.04 (on travis-ci), R 3.6.3, R 4.0.0, R-devel
* ubuntu 16.04 (on rhub), R 3.6.1
* fedora 32 (on rhub), R-devel
* mac OS X High Sierra (on travis-ci), R 3.6.3, R 4.0.2
* local Windows 10 build, R 4.0.2
* Windows Server 2008 (on rhub) R-devel
* Windows Server 2012 (on appveyor) R 4.0.2
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

> * checking CRAN incoming feasibility ... NOTE
> Maintainer: 'James Melville <jlmelville@gmail.com>'
>
> There was a message about possibly mis-spelled words in DESCRIPTION:
>
>  BFGS (8:64, 9:20, 9:28)
>  Broyden (8:30)
>  Goldfarb (8:47)
>  Shanno (8:56)

Those words are spelled correctly.

For winbuilder checks only, there was a NOTE about a URL in a vignette
(<https://doi.org/10.1137/030601880>)

> Found the following (possibly) invalid URLs:
>  URL: https://doi.org/10.1137/030601880
>    From: inst/doc/convergence.html
>          inst/doc/mize.html
>    Status: Error
>    Message: libcurl error code 56:
>      	Recv failure: Connection was reset
      	
I am unable to reproduce this and can access the URL without problems via
curl and a web browser.

## CRAN checks

There are no NOTEs, ERRORs or WARNINGs.

## Downstream dependencies

There are 2 downstream dependencies. 

* All completed R CMD CHECK without issues.
