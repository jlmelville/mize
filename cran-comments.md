## Release Summary

This is the initial CRAN submission of mize.

## Test environments

* ubuntu 12.04 (on travis-ci), R 3.4.0
* local Windows 10 build, R 3.4.1
* win-builder (devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'James Melville <jlmelville@gmail.com>'

This is the first submission of the package and my first submission to CRAN.

There was a message about possibly mis-spelled words in DESCRIPTION:

  BFGS (8:64, 9:20, 9:28)
  Broyden (8:30)
  Goldfarb (8:47)
  Shanno (8:56)

Those words are spelled correctly.

There was a message about an invalid URL:

Found the following (possibly) invalid URLs:
  URL: http://www.cs.umd.edu/users/oleary/software/
    From: inst/doc/mize.html
    Status: Error
    Message: libcurl error code 35:
    	Unknown SSL protocol error in connection to www.cs.umd.edu:443

The URL resolves correctly in a web browser and via 
curl::curl(url = "http://www.cs.umd.edu/users/oleary/software/") in R 3.4.1
on my Windows 10 test environment.

## Downstream dependencies

None.
