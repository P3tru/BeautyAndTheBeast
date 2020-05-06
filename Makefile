# Basic Makefile

### Compilers
CC  = gcc
CXX = g++

DEBUG_LEVEL    = -g
EXTRA_CCFLAGS  = -W -Wall
CPPFLAGS       = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS        = $(CPPFLAGS)

RM = rm -f
MV = mv

SRCDIR := src
INCDIR := include

### ROOT
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

### RAT
RATLIBS  := -L$(RATROOT)/lib -lRATEvent

### BOOST
BOOSTCFLAGS := -I/data/snoplus/home/zsoldos/.local/boost-1.71.0
BOOSTLIBS   := -L/data/snoplus/home/zsoldos/.local/boost-1.71.0/lib -lboost_system -lboost_filesystem

### Python
PYTHONCFLAGS := $(shell python3-config --cflags)
PYTHONLIBS   := $(shell python3-config --ldflags)

CPPFLAGS  += -I$(INCDIR) $(ROOTCFLAGS) -I$(RATROOT)/include
CPPFLAGS  +=  $(BOOSTCFLAGS)
CPPFLAGS  +=  $(PYTHONCFLAGS)
EXTRALIBS  = $(ROOTLIBS)
EXTRALIBS += $(RATLIBS)
EXTRALIBS += -lz
EXTRALIBS += $(BOOSTLIBS)
EXTRALIBS += $(PYTHONLIBS)

SRCS = $(wildcard $(SRCDIR)/*.cc)
OBJS = $(subst .cc,.o,$(SRCS))

.PHONY: all clean 
.DEFAULT_GOAL = FlattenHits

help:
	@grep -h -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

all: FlattenHits

libcnpy.so: libcnpy.o
	@echo "Compiling library $@"
	$(CXX) $(CPPFLAGS) -shared $^ -o lib/$@ $(EXTRALIBS)
	$(RM) $^
libcnpy.o:
	$(CXX) $(CPPFLAGS) -fPIC -c src/cnpy.cpp -o $@ $(EXTRALIBS)


FlattenHits: libcnpy.so FlattenHits.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o FlattenHits FlattenHits.cc $(OBJS) $(EXTRALIBS) -L$(PWD)/lib -lcnpy
	$(RM) FlattenHits.o $(OBJS)

clean:
	$(RM) $(OBJS) libcnpy.o FlattenHits
