
include ../Makefile.options

SRCS=\
	adapt_mesh.cpp \
	refine_mesh.cpp \
	refine_elements.cpp \

OBJECTS = $(SRCS:.cpp=.o)

adapt.a:	$(OBJECTS)
		$(LINKSTATIC) lib$@ $(OBJECTS)
		$(RANLIB) lib$@

adapt.so:	$(OBJECTS)
		$(LINKSHARED) $@ $(OBJECTS)

clean:
	/bin/rm -f *.o
	/bin/rm -f libadapt.a
	/bin/rm -f libadapt.so

include ../suffix.rules
