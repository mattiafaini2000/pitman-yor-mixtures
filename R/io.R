# Resolve inputs explicitly; generated files must stay inside the checkout.
project_root <- function(root = getwd()) {
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

project_input_path <- function(path, root = project_root()) {
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(root, path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

checked_project_path <- function(path, root = project_root()) {
  root <- project_root(root)
  path <- gsub("\\\\", "/", path)
  if (grepl("(^|/)\\.\\.(/|$)", path)) stop("Output paths cannot contain '..'.")
  if (!grepl("^(/|[A-Za-z]:/)", path)) path <- file.path(root, path)
  ancestor <- path
  while (!file.exists(ancestor) && !dir.exists(ancestor)) {
    next_ancestor <- dirname(ancestor)
    if (identical(ancestor, next_ancestor)) stop("Cannot resolve output directory.")
    ancestor <- next_ancestor
  }
  resolved <- normalizePath(ancestor, winslash = "/", mustWork = TRUE)
  inside <- function(candidate) {
    if (.Platform$OS.type == "windows") {
      candidate <- tolower(candidate)
      root <- tolower(root)
    }
    identical(candidate, root) || startsWith(candidate, paste0(root, "/"))
  }
  if (!inside(resolved) || !inside(path)) stop("Outputs must stay inside the project.")
  path
}

project_output_path <- function(path, root = project_root()) {
  path <- checked_project_path(path, root)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

ensure_project_directory <- function(path, root = project_root()) {
  path <- checked_project_path(path, root)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_project_csv <- function(object, path, root = project_root()) {
  utils::write.csv(object, project_output_path(path, root), row.names = FALSE)
}
