GPUARCH ?= 35
NVCC ?= nvcc
BOOSTINC ?= /usr/local/sage/local/include

METHOD ?= SCF
LMAX ?= 6
NMAX ?= 10

ifeq ($(BOOSTINC), NOBOOST)
	BOOST = -DNOBOOST
else
	BOOST = -I$(BOOSTINC)
endif

CFLAGS := -D$(METHOD) -DLMAX=$(LMAX) -DNMAX=$(NMAX) $(BOOST) -O3
LIBS = -lm -lmpi
CXXFLAGS += $(CFLAGS)
LDFLAGS  += -lm $(MUSE_LD_FLAGS)
CUDAFLAGS += -arch=sm_$(GPUARCH) $(CFLAGS)

CODELIB = libetics.a

# object files needed for the library and the standalone
CODEOBJS1 = integrate.o mathaux.o scf.o ic.o

# object files only needed for the standalone
CODEOBJS2 = io.o main.o

AR = ar ruv
RANLIB = ranlib

all: $(CODELIB) standalone

clean:
	$(RM) *.o *.a etics

$(CODELIB): $(CODEOBJS1)
	$(NVCC) $(CUDAFLAGS) -dlink $(CODEOBJS1) -o dlink.o
	$(AR) $@ $(CODEOBJS1) dlink.o
	$(RANLIB) $@

library: $(CODELIB)

standalone: $(CODEOBJS1) $(CODEOBJS2)
	$(NVCC) $(LIBS) $(LINKFLAGS) $(CUDAFLAGS) -o etics $(CODEOBJS1) $(CODEOBJS2)

.SUFFIXES: .cu .cpp .o

.cu.o: $<
	$(NVCC) $(CUDAFLAGS) -dc -o $@ $<

.cpp.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
