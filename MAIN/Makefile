
include ../Makefile.options

LIBRARY_PATH += -L../IO -L../UTIL -L../ADAPT

LIBS += -ladapt -lio -lutility -lSpaceField -loctree -lhdf5 -lcgns
#LIBS += -lstencil -lgeometry -lsmooth -lio -lutility -lcgns /usr/local/hdf5-1.6.4-icc/lib/libhdf5.a 

LDFLAGS = $(LIBRARY_PATH) $(LIBS)

C_SRCS = P_REFINE.cpp \

SRCS = $(C_SRCS)
OBJECTS = $(C_SRCS:.cpp=.o)

all: P_REFINE

P_REFINE: io adapt $(OBJECTS)
	$(LINK) -o $@ $(OBJECTS) $(LDFLAGS)
	cp $@ ../$@.$(EXE_SUFFIX)
	cp $@ ../../$@.$(EXE_SUFFIX)

io:
	cd ../IO; $(MAKE)

utility:
	cd ../UTIL; $(MAKE)

adapt:
	cd ../ADAPT; $(MAKE)

clean:
	/bin/rm -f *.o
	/bin/rm -f P_REFINE

include ../suffix.rules
