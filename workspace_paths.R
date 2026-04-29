find_workspace_root <- function() {
  src_path <- NULL

  for (frame in sys.frames()) {
    if (!is.null(frame$ofile)) {
      src_path <- frame$ofile
    }
  }

  if (!is.null(src_path)) {
    return(normalizePath(dirname(src_path), winslash = "/", mustWork = TRUE))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

workspace_root <- find_workspace_root()
read_dir <- file.path(workspace_root, "read")
raw_data_dir <- file.path(workspace_root, "projectfile")
output_dir <- file.path(workspace_root, "output")
code_dir <- file.path(workspace_root, "projectmd")

workspace_dirs <- function() {
  list(
    root = workspace_root,
    read = read_dir,
    projectfile = raw_data_dir,
    output = output_dir,
    projectmd = code_dir
  )
}

read_path <- function(...) {
  file.path(read_dir, ...)
}

raw_data_path <- function(...) {
  file.path(raw_data_dir, ...)
}

output_path <- function(...) {
  path <- file.path(output_dir, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

code_path <- function(...) {
  file.path(code_dir, ...)
}