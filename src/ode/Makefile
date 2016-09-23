#CC = gcc
#FC = gfortran

.PHONY: all

LIBODEAUX = libodeaux.a

LIB_SRCS = opkda2.f dlsoda.f opkda1_abbr.f  vode.f 
LIB_OBJS = $(LIB_SRCS:%.f=%.o) dop853.o

all: $(LIBODEAUX)

$(LIBODEAUX): $(LIB_OBJS)
	ar cruv $@ $^
	ranlib $@

dop853.o: dop853.c


 
