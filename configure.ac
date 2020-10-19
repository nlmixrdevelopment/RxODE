## This is based on https://stat.ethz.ch/pipermail/r-devel/2010-September/058488.html
AC_INIT(RxODE, 0.9 , matthew.fidler@gmail.com)
AC_CONFIG_SRCDIR([src/init.c])

# find R home and set CC/CFLAGS
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
RBIN="${R_HOME}/bin/R"
CC=`"${RBIN}" CMD config CC`;
CFLAGS=`"${RBIN}" CMD config CFLAGS`
LIBS="${PKG_LIBS}"

# Checks for programs.
AC_PROG_CC
# Check for OpenMP

CXX14STD=""

ISYSTEM="I"
if test "${GCC}" = yes; then
  ISYSTEM="isystem"
fi



AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(OPENMP_CXXFLAGS)
AC_SUBST(CXX14STD)
AC_SUBST(ISYSTEM)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

echo "unlink(\"R/RxODE_md5.R\"); writeLines(sprintf(\"RxODE.md5 <- \\\"%s\\\"\\n\", digest::digest(sapply(c(paste0(\"src/\", list.files(\"src\", pattern = \"[.](c|h|cpp|f)$\")), paste0(\"inst/include/\", list.files(\"inst/include\")), paste0(\"R/\", list.files(\"R/\", pattern = \"[.]R\"))), digest::digest, file = TRUE))), \"R/RxODE_md5.R\")" > build.R

${R_HOME}/bin/R CMD BATCH build.R
rm build.R
rm build.Rout

echo --------[[begin src/Makevars]]--------
cat src/Makevars
echo --------[[end src/Makevars]]--------

echo --------[[begin R/RxODE_md5.R]]--------
cat R/RxODE_md5.R
echo --------[[end R/RxODE_md5.R]]--------

