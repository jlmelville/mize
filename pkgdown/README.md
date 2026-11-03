# Maintaining the website

Build the site locally from the package root:

```r
pkgdown::build_site()
```

The explicit `development: mode: devel` setting builds the current GitHub
version in `pkgdown-site/dev`, even though the package version has three
components. To serve the existing output locally, install the optional
`nanonext` package and run `pkgdown::preview_site()`. Stop the server with
`pkgdown::stop_preview()`.

The website workflow builds the CRAN-matched released source at the site root
and the current GitHub source under `/dev`. It publishes the `pkgdown-site`
directory to the `gh-pages` branch. A GitHub release can precede CRAN
acceptance; once CRAN accepts the release, rerun the workflow manually to
update the stable site.
