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

sh inst/tools/stripper/stripper

echo --------[[begin src/Makevars]]--------
cat src/Makevars
echo --------[[end src/Makevars]]--------

