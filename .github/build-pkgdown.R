build_sites <- function() {
  root <- normalizePath(".")
  config <- yaml::read_yaml(file.path(root, "_pkgdown.yml"))
  destination <- file.path(root, config$destination)
  scratch <- tempfile("mize-release-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

  packages <- available.packages(
    repos = c(CRAN = "https://cloud.r-project.org"),
    type = "source"
  )
  version <- packages["mize", "Version"]
  tag <- paste0("v", version)
  archive <- file.path(scratch, "release.tar.gz")
  download.file(
    paste0(
      "https://github.com/jlmelville/mize/archive/refs/tags/",
      tag,
      ".tar.gz"
    ),
    archive,
    mode = "wb"
  )
  untar(archive, exdir = scratch)
  release <- file.path(scratch, paste0("mize-", version))
  stopifnot(
    read.dcf(file.path(release, "DESCRIPTION"))[1, "Version"] == version
  )

  released_config <- yaml::read_yaml(file.path(release, "_pkgdown.yml"))
  for (key in c("url", "template", "llm-docs")) {
    released_config[[key]] <- config[[key]]
  }
  released_config$navbar$structure$right <- config$navbar$structure$right
  released_config$navbar$components$versions <- config$navbar$components$versions
  released_config$development$mode <- "release"
  released_config$repo$branch <- tag
  yaml::write_yaml(released_config, file.path(release, "_pkgdown.yml"))
  file.copy(file.path(root, "pkgdown"), release, recursive = TRUE)

  pkgdown::build_site_github_pages(
    release,
    dest_dir = destination,
    install = TRUE,
    new_process = TRUE
  )
  pkgdown::build_site_github_pages(
    root,
    dest_dir = destination,
    install = FALSE,
    new_process = TRUE
  )
}

build_sites()
