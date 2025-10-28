# Cleanup script for Quarto
files_to_remove <- c(
  list.files(pattern = "\\.(bst|spl|cls|ipynb|sty)$")
  # list.files("site_libs", full.names = TRUE, recursive = TRUE)
)

# Iterate over files and delete them if they exist
for (file in files_to_remove) {
  if (file.exists(file)) {
    file.remove(file)
    cat("Deleted:", file, "\n")
  } else {
    cat("File not found:", file, "\n")
  }
}

# Code to delete entire folder site_libs
fol_to_del <- c("site_libs", list.files(pattern = "_files"))
for (i in fol_to_del) {
  if (dir.exists(i)) {
    unlink(i, recursive = TRUE)
    cat("Deleted:", i, "\n")
  } else {
    cat("Folder not found:", i, "\n")
  }
}
