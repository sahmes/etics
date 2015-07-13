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

CODEOBJS = integrate.o io.o main.o mathaux.o scf.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB)

clean:
	$(RM) -f *.o *.a etics

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

standalone: $(CODEOBJS)
	 $(NVCC) $(LIBS) $(LINKFLAGS) $(CUDAFLAGS) -o etics $(CODEOBJS)

.SUFFIXES: .cu .cpp .o

.cu.o: $<
	$(NVCC) $(CUDAFLAGS) -dc -o $@ $<

.cpp.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
