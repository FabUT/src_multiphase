
#F90=mpiifort -O2 -fp-model precise -assume byterecl -g -traceback -fp-stack-check -ftrapuv -check pointers -check bounds

#####F90=mpiifort -O2 -g -traceback -check bounds -ftrapuv

#F90=mpiifort -O0 -fp-model precise -assume byterecl -g -traceback -fp-stack-check -ftrapuv -check pointers -check bounds

F90=gfortran -g -fbacktrace -fbounds-check -ffree-line-length-none -ffpe-trap=invalid


#
#opt= -mcmodel=large -fbacktrace -fbounds-check 

### FOR PETSC
PDIR=/d/ftholin/Program/TARANIS/LIBS/petsc-3.9.1
opt= #-I$(PDIR)/include -I$(PDIR)/arch-linux2-c-debug/include

LIBS=\
-Wl,-rpath,/d/ftholin/Program/TARANIS/LIBS/petsc-3.9.1/arch-linux2-c-debug/lib -Wl,-rpath,/d/ftholin/Program/TARANIS/LIBS/petsc-3.9.1/arch-linux2-c-debug/lib -L/d/ftholin/Program/TARANIS/LIBS/petsc-3.9.1/arch-linux2-c-debug/lib -Wl,-rpath,/d/ftholin/Program/TARANIS/LIBS/petsc-3.9.1/arch-linux2-c-debug/lib -Wl,-rpath,/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib/debug_mt -L/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib/debug_mt -Wl,-rpath,/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib -L/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib -Wl,-rpath,/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/ipp/lib/intel64 -L/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/ipp/lib/intel64 -Wl,-rpath,/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64_lin -L/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64_lin -Wl,-rpath,/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64_lin -L/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64_lin -Wl,-rpath,/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64/gcc4.4 -L/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64/gcc4.4 -Wl,-rpath,/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/daal/lib/intel64_lin -L/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/daal/lib/intel64_lin -Wl,-rpath,/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64_lin/gcc4.4 -L/opt/tools/intel/studio/2017/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64_lin/gcc4.4 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -Wl,-rpath,/opt/intel/mpi-rt/2017.0.0/intel64/lib/debug_mt -Wl,-rpath,/opt/intel/mpi-rt/2017.0.0/intel64/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lHYPRE -llapack -lblas -lX11 -lstdc++ -ldl -lmpifort -lmpi -lmpigi -lrt -lpthread -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lgcc_s -lirc_s -lstdc++ -ldl


###

list=mod_data.o mod_input.o mod_mat_tab.o mod_Fermi_Dirac.o air_ETL.o mod_euler.o mod_elec.o mod_output.o 

#multi:  main.f90 $(list)
#	$(F90) $(opt) $(list) $(LIBS) $< -o $@

multi:  main.f90 $(list)
	$(F90) $(opt) $(list) $< -o $@

mod_input.o mod_input.mod: mod_input.f90
	$(F90) $(opt) -c $< 

mod_euler.o mod_euler.mod: mod_euler.f90
	$(F90) $(opt) -c $< 

air_ETL.o air_ETL.mod: air_ETL.f90
	$(F90) $(opt) -c $< 

mod_mat_tab.o mod_mat_tab.mod: mod_mat_tab.f90
	$(F90) $(opt) -c $< 

mod_elec.o mod_elec.mod: mod_elec.f90
	$(F90) $(opt) -c $< 

mod_data.o mod_data.mod: mod_data.f90
	$(F90) $(opt) -c $< 

mod_snes.o mod_snes.mod: mod_snes.F90
	$(F90) $(opt) $(PETSC_LIB) -c mod_snes.F90

mod_Fermi_Dirac.o mod_Fermi_Dirac.mod: mod_Fermi_Dirac.f90
	$(F90) $(opt) -c $<

mod_output.o mod_output.mod: mod_output.f90
	$(F90) $(opt) -c $< 

clean:
	rm *.o *.mod
