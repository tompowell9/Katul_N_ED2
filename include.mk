#------------------------------------------------------------------------------------------#
#    This is the file you need to adjust depending on your system and needs.               #
#------------------------------------------------------------------------------------------#



#------ Define make (gnu makes works best) ------------------------------------------------#
MAKE = /usr/bin/make

# RAPP root directory
RAPP_ROOT=/home/dmedvigy/code/ed2/RAPP

#------------------------------------------------------------------------------------------#
#    HDF5 libraries. You don't have to include them if you don't want to use anything in   #
# HDF5, in which case you can set up USE_HDF5 to 0, and leave HDF5_INCS and HDF5_LIBS in   #
# blank.                                                                                   #
#------------------------------------------------------------------------------------------#
USE_HDF5=1
HDF5_DIR=/home/dmedvigy/hdf5-1.8.10-patch1/hdf5
HDF5_INCS=-I$(HDF5_DIR)/include -I$(HDF5_DIR)/lib
HDF5_LIBS=-L$(HDF5_DIR)/lib -lz -lm -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 
#------------------------------------------------------------------------------------------#
#    NetCDF libraries. You don't have to include them if you don't want to use anything in #
# netCDF, in which case you can set up USE_NCDF to 0, and leave NCDF_INCS and NCDF_LIBS in #
# blank.                                                                                   #
#------------------------------------------------------------------------------------------#
USE_NCDF=1
NCDF_DIR=/home/dmedvigy/netcdf-4.2.1.1/netcdf
NCDF_INCS=-I$(NCDF_DIR)/include
NCDF_LIBS=-L$(NCDF_DIR)/lib -lnetcdff -lnetcdf

#------ Defining the compiler and library paths in case they are not in LD_LIBRARY_PATH ---#
CMACH=PC_LINUX1
F_COMP=ifort
C_COMP=icc
LOADER=ifort
LIBS=

##################################### COMPILER OPTIONS #####################################

#------------------------------------------------------------------------------------------#
# A. Pickiest - Use this whenever you change arguments on functions and subroutines.       #
#               This will perform the same tests as B but it will also check whether all   #
#               arguments match between subroutine declaration and subroutine calls.       #
#               WARNING: In order to really check all interfaces you must compile with     #
#                        this option twice:                                                #
#               1. Compile (./compile.sh)                                                  #
#               2. Prepare second compilation(./2ndcomp.sh)                                #
#               3. Compile one more time (./compile.sh)                                    #
#               If the compilation fails either at step 1 or 3, then your code has inter-  #
#                  face problems. If it successfully compiles, then you can switch to B.   #
#------------------------------------------------------------------------------------------#
#USE_INTERF=0
#F_OPTS= -FR -O0 -recursive -Vaxlib  -check all -g -fpe0 -ftz -gen-interfaces            \
         -warn interfaces -debug extended -debug inline_debug_info -debug-parameters all \
         -traceback -ftrapuv
#C_OPTS= -O0 -DLITTLE  -g -traceback -debug extended
#LOADER_OPTS= -FR -O0 -Vaxlib  -check all -g -fpe0 -ftz -gen-interfaces \
             -warn interfaces -debug extended -debug inline_debug_info  \
             -debug-parameters all -traceback -ftrapuv
#C_LOADER_OPTS=-v -g -traceback
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# B. Pickiest with no interface - This will compile fast but the run will be slow due to   #
#    the -O0 option. However, by setting -O0 you will take full advantage of the intel     #
#    debugger.                                                                             #
#    Ideally, you should compile your code with this option whenever you make any changes. #
#    Note, however, that if you change arguments you should first try A.                   #
#------------------------------------------------------------------------------------------#
#USE_INTERF=1
#F_OPTS= -FR -O0 -recursive -Vaxlib  -check all -g -fpe0 -ftz  -debug extended \
         -debug inline_debug_info -debug-parameters all -traceback -ftrapuv
#C_OPTS= -O0 -DLITTLE  -g -traceback -debug extended
#LOADER_OPTS= -FR -O0 -Vaxlib  -check all -g -fpe0 -ftz -debug extended        \
              -debug inline_debug_info -debug-parameters all -traceback -ftrapuv
#C_LOADER_OPTS=-v -g -traceback
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# C. Fast debugging - This will check pretty much the same as B, but with higher optimiza- #
#    tion. However, by setting -O2 you won't be able to see all variables on the debugger. #
#    You should use this if your run seems to be working fine at the beginning, but it is  #
#    failing or giving instabilities, or funny results after a long time.                  #
#------------------------------------------------------------------------------------------#
#USE_INTERF=1
#F_OPTS= -FR -O2 -recursive -Vaxlib  -check all -g -fpe0 -ftz  -debug extended \
#        -debug inline_debug_info -debug-parameters all -traceback -ftrapuv
#C_OPTS= -O2 -DLITTLE  -g -traceback -debug extended
#LOADER_OPTS= -FR -O2 -Vaxlib  -check all -g -fpe0 -ftz -debug extended        \
#             -debug inline_debug_info -debug-parameters all -traceback -ftrapuv
#C_LOADER_OPTS=-v -g -traceback
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# D. Fast check - This will check pretty much the same as C, but it will not set up        #
#    anything for the debugger. Use this only if you really don't want 249544to deal with idb    #
#    or if you have a good idea of which problem you are dealing with.                     #
#------------------------------------------------------------------------------------------#
#USE_INTERF=1
#F_OPTS= -FR -O2 -recursive -Vaxlib -check all -fpe0 -ftz  -traceback -ftrapuv
#C_OPTS= -O2 -DLITTLE -traceback
#LOADER_OPTS= -FR -O2 -Vaxlib  -check all -fpe0 -ftz -traceback -ftrapuv
#C_LOADER_OPTS=-v -traceback
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
# E. Fast - This is all about performance, use only when you are sure that the model has   #
#           no code problem, and you want results asap. This will not check for any        #
#           problems, which means that this is an option suitable for end users, not de-   #
#           velopers.                                                                      #
#------------------------------------------------------------------------------------------#
USE_INTERF=1
F_OPTS=-O3 -static -xHost -traceback
C_OPTS=-O3 -static -DUNDERSCORE -DLITTLE
LOADER_OPTS=
C_LOADER_OPTS=-v -traceback
#------------------------------------------------------------------------------------------#

############################################################################################

#------ Archive command -------------------------------------------------------------------#
ARCHIVE=ar rs
