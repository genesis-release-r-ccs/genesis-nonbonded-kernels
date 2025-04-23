# HOWTO compile
# make FC=frtpx/frt/gfortran/flang/armflang OMP=omp/noomp SIMD=dir/nodir KIND=Generic/Oct PERF=y/n
SEARCH := frt frtpx
is_FJ := $(words $(filter $(FC), $(SEARCH)))
#  select your path
TOPDIR_FJ=
TOPDIR_GCC=
TOPDIR_ACFL=
TOPDIR_LLVM=
# Fugaku
#MARCH=armv8.2-a+sve 
# Graviton 3
#MARCH=armv8.4-a+sve
# Graviton 4
#MARCH=armv9-a+sve

ifneq ($(is_FJ), 0)
	PERFHELPDIR=$(TOPDIR_FJ)
	OMPPARA=-Kopenmp
	ifeq ($(SIMD), dir)
		KOCL=-Kocl
	endif
	ifeq ($(Kfast), y) 
  	  FCFLAGS = -m $(KOCL) -Kfast,parallel,simd=2,nounroll,swp_strong,noalias=s $(OMPPARA)
        else
	  FCFLAGS = -m $(KOCL) -Knounroll,noalias=s -O3 $(OMPPARA) 
	endif
	LDFLAGS = $(OMPPARA)  -Kparallel -Nfjomplib 
	PPFLAGS = -P -Cpp -Cfpp
	ifeq ($(FC),frt)
		FPP=frt
	else
		FPP=frtpx
	endif
else ifeq ($(FC),gfortran)
        Kfast=n
	PERFHELPDIR=$(TOPDIR_GCC)
	OMPPARA=-fopenmp
	FCFLAGS = -O3 -march=$(MARCH) -Wall $(OMPPARA) -ffree-line-length-none -fallow-argument-mismatch
	LDFLAGS = $(OMPPARA)
	PPFLAGS=-traditional-cpp -traditional
	FPP=cpp
else ifeq ($(FC),flang)
        Kfast=n
	PERFHELPDIR=$(TOPDIR_LLVM)
	OMPPARA=-fopenmp
	FCFLAGS = -O3 -march=$(MARCH) -Wall $(OMPPARA) 
	LDFLAGS = $(OMPPARA) 
	PPFLAGS=-traditional-cpp -traditional
	FPP=cpp
else ifeq ($(FC),armflang)
        Kfast=n
	PERFHELPDIR=$(TOPDIR_ACFL)
	OMPPARA=-fopenmp
	FCFLAGS = -O3 -march=$(MARCH) -Wall $(OMPPARA) 
	LDFLAGS = $(OMPPARA) -L$(PERFHELPDIR) -lperf_helper
	PPFLAGS=-traditional-cpp -traditional
	FPP=cpp
endif
ifeq ($(PERF), y) 
	DIRHEADER=perf
	DEFS=-DHAVE_PERF
	PERFLINK = 
	PERF_HELPER_MOD=$(PERFHELPDIR)/perf_helper_mod.mod
	INCLUDES=-I$(PERFHELPDIR)
else
	DIRHEADER=normal
	DEFS=
	PERFLINK =
	PERF_HELPER_MOD=
	INCLUDES=
endif
ifeq ($(SIMD),dir)
	DEFS+=-DDIR
endif
ifeq ($(OMP),omp)
	DEFS+=-DOMP_PARALLEL
endif

KERNEL=kernel
ifeq ($(Kfast), y) 
  DIR=$(DIRHEADER)_$(FC)_org_$(OMP)_$(SIMD)
else
  DIR=$(DIRHEADER)_$(FC)_$(OMP)_$(SIMD)
endif

OBJS= sub_data_io.o initialize_data.o kernel_$(KIND).o main.o

.SUFFIXES : .o .fpp .f90

.fpp.f90:
ifneq ($(is_FJ), 0)
	cp $*.fpp $*.f90
	$(FPP) $(PPFLAGS) $(DEFS) $*.f90
	mv $*.cpp.f90 $*.f90
else
	$(FPP) $(PPFLAGS) $(DEFS) $*.fpp $*.f90
	$(FC)  $(INCLUDES) $(FCFLAGS) -c $*.f90
endif
.f90.o:
	$(FC) $(INCLUDES) $(FCFLAGS) -c $*.f90

$(KERNEL) : $(OBJS) $(PERF_HELPER_MOD)
	mkdir -p $(DIR)
	$(FC) -o $(DIR)/$(KERNEL) $(OBJS) $(LDFLAGS) $(PERFLINK)


.PHONY: clean
clean:
	rm -f *.o *.mod *.bak *.s *.f90 *.lst
$(PERF_HELPER_MOD): 
	cp $(PERF_HELPER_MOD) .
assemble: kernel_$(KIND).f90 $(OBJS)
	mkdir -p $(DIR)
	$(FC) $(INCLUDES) $(FCFLAGS) -S kernel_$(KIND).f90
	mv kernel_$(KIND).s $(DIR)
