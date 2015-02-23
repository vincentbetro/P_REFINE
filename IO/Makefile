
include ../Makefile.options

SRCS=\
	P_CGNS.cpp \

OBJECTS = $(SRCS:.cpp=.o)

io.a:	$(OBJECTS)
		$(LINKSTATIC) lib$@ $(OBJECTS)
		$(RANLIB) lib$@

io.so:	$(OBJECTS)
		$(LINKSHARED) $@ $(OBJECTS)

clean:
	/bin/rm -f *.o
	/bin/rm -f libio.a
	/bin/rm -f libio.so

include ../suffix.rules
