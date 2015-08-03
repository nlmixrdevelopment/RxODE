# Copy libs, headers, and binaries to to
# R_PACKAGE_DIR/{libs,include,bin} as per architecture

bin  <- file.path(R_PACKAGE_DIR, paste0("bin", R_ARCH))
libs <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
incl <- file.path(R_PACKAGE_DIR, "include")

if(!file.exists(bin))
   dir.create(bin, recursive = TRUE, showWarnings = FALSE)

if(!file.exists(libs))
   dir.create(libs, recursive = TRUE, showWarnings = FALSE)

if(!file.exists(incl))
   dir.create(incl, recursive = TRUE, showWarnings = FALSE)

# binary files
files <- file.path(".", "tran.exe")
message("Installing ", files, " to ", bin)
file.copy(files, bin, overwrite = TRUE)

# libs
files <- file.path("ode", "libodeaux.a")
message("Installing ", files, " to ", libs)
file.copy(files, libs, overwrite = TRUE)

# header files
files <- file.path("ode", "dop853.h")
message("Installing ", files, " to ", incl)
file.copy(files, incl, overwrite = TRUE)

