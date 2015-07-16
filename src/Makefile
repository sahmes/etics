GPUARCH ?= 35
NVCC ?= nvcc

METHOD ?= SCF
LMAX ?= 6
NMAX ?= 10

CFLAGS := -D$(METHOD) -DLMAX=$(LMAX) -DNMAX=$(NMAX) -O3 -I/usr/local/sage/local/include
LIBS = -lm -lmpi
CXXFLAGS += $(CFLAGS)
LDFLAGS  += -lm $(MUSE_LD_FLAGS)
CUDAFLAGS += -arch=sm_$(GPUARCH) $(CFLAGS)

CODELIB = libetics.a

# object files needed for the library and the standalone
CODEOBJS1 = integrate.o mathaux.o scf.o

# object files only needed for the standalone
CODEOBJS2 = io.o ic.o main.o

AR = ar ruv
RANLIB = ranlib

all: $(CODELIB)

clean:
	$(RM) *.o *.a etics

$(CODELIB): $(CODEOBJS1)
	$(NVCC) $(CUDAFLAGS) -dlink $(CODEOBJS1) -o dlink.o
	$(AR) $@ $(CODEOBJS1) dlink.o
	$(RANLIB) $@

standalone: $(CODEOBJS1) $(CODEOBJS2)
	$(NVCC) $(LIBS) $(LINKFLAGS) $(CUDAFLAGS) -o etics $(CODEOBJS1) $(CODEOBJS2)

.SUFFIXES: .cu .cpp .o

.cu.o: $<
	$(NVCC) $(CUDAFLAGS) -dc -o $@ $<

.cpp.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
