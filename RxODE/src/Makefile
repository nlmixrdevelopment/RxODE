CC = `"$(R_HOME)/bin/R" CMD config CC`
FC = `"$(R_HOME)/bin/R" CMD config FC`
CFLAGS = `"$(R_HOME)/bin/R" CMD config CFLAGS`
FFLAGS = `"$(R_HOME)/bin/R" CMD config FFLAGS`

# Run make in subdirs
all: libodeaux libdparse libtran

libodeaux:
	echo "make libodeaux.a in ode/ ..."
	(cd ode; $(MAKE) CC="$(CC)" FC="$(FC)" CFLAGS="$(CFLAGS) -fPIC" FFLAGS="$(FFLAGS) -fPIC" libodeaux.a)

libdparse:
	echo "make libdparse.a in dparser/ ..."
	(cd dparser; $(MAKE) -f Makefile.org CC="$(CC)" CFLAGS="$(CFLAGS)" libdparse.a)

libtran:
	echo "make tran ..."
	($(CC) tran*.c -Idparser -Ldparser -ldparse -o tran.exe)

clean:
	echo "cleaning ..."
	(cd dparser;rm -f *.o *.a)
	(cd ode;rm -f *.o *.a)
	rm -f *.exe
