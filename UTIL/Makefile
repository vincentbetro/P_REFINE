
include ../Makefile.options

SRCS=\
	List.cpp \
	sort.cpp \

OBJECTS = $(SRCS:.cpp=.o)

utility.a:		$(OBJECTS)
		$(LINKSTATIC) lib$@ $(OBJECTS)
		$(RANLIB) lib$@

utility.so:	$(OBJECTS)
		$(LINKSHARED) $@ $(OBJECTS)

clean:
	/bin/rm -f *.o
	/bin/rm -f libutility.a
	/bin/rm -f libutility.so

include ../suffix.rules

