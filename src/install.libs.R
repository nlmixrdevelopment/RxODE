## Copy binaries/libraries to the appropriate directory.
libs <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH));
dir.create(libs,recursive=TRUE);
files <- file.path(".", sprintf("RxODE%s",.Platform$dynlib.ext));
message("Installing ", files, " to ", libs);
file.copy(files, libs, overwrite = TRUE);

if (.Platform$dynlib.ext != ".dll"){
    libs <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH),sprintf("libRxODE%s",.Platform$dynlib.ext));
    message("Installing ", files, " to ", libs);
    file.copy(files, libs, overwrite = TRUE);
}


