
#=================================================
#		Define compiler names
#=================================================

FC      = gfortran
FCFLAGS = -x f95-cpp-input
FPIC    = -fPIC
LIBTOOL = ar src

#=================================================
#		Define variables and options
#=================================================

# Case chosen to wrapping
CASE = 3
ifeq ($(CASE), 1) 
	DIR_NAME     = simple_case
	MODULE_NAME  = module

	SOURCE       = build/wrappers/test 

	WRAPPER_FILE = build/wrappers/test

	KIND_MAP     =
	move         = 
	string       =
else ifeq ($(CASE), 2)
	DIR_NAME     = example_cylinder
	MODULE_NAME  = calcul_volume

	SOURCE       = build/wrappers/DNAD \
	       	       build/wrappers/cyldnad

	WRAPPER_FILE = build/wrappers/DNAD \
	       	       build/wrappers/cyldnad

	KIND_MAP     = -k kind_map
	move         =
	string       =
else ifeq ($(CASE), 3)
	DIR_NAME     = dassflow1d_test
	MODULE_NAME  = dassflow1d

	SOURCE       = build/wrappers/m_common \
		       build/wrappers/m_linear_algebra 
#\
#		       build/wrappers/m_mesh \
#		       build/wrappers/m_obs \
#		       build/wrappers/m_sw_mono

	WRAPPER_FILE = build/wrappers/m_linear_algebra 
#\
#		       build/wrappers/m_mesh \
#		       build/wrappers/m_obs \
#		       build/wrappers/m_sw_mono

	KIND_MAP     = -k kind_map
	move         = -M
	string       = -S 128
endif

OBJ_TMP        = $(addsuffix .o, ${SOURCE})
F90_TMP        = $(addsuffix .f90, ${SOURCE})
F90_WRAPPER    = $(addsuffix .f90, ${WRAPPER_FILE})
PACKAGE	       = -P # Create a package instead of a module
#PACKAGE       =

#=================================================
#		Python wrappers
#=================================================

build/wrappers/_${MODULE_NAME}.so: build/wrappers/libsrc.a $(F90_TMP)
	cd build/wrappers ; f90wrap -m ${MODULE_NAME} $(patsubst build/wrappers/%.f90, %.f90, ${F90_WRAPPER}) ${KIND_MAP} ${string} ${PACKAGE} ${move}
	cd build/wrappers ; f2py-f90wrap --fcompiler=${FC} --opt="-O3" --build-dir . -c -m _${MODULE_NAME} -L. -lsrc f90wrap*.f90
#	cp ./src/cases_test/dassflow1d_test/wrappers/finish_to_gen_wrappers.pl ./build/wrappers/
#	cd ./build/wrappers ; perl finish_to_gen_wrappers.pl

#=================================================
#	Copy and compile Fortran files
#=================================================

files: copy compile

copy:
	echo "# ==================================================" ; \
	echo " Copy source files " ; \
	echo "# ==================================================" ; \
	cp ./src/cases_test/${DIR_NAME}/*.f90 ./build/wrappers
	cp ./src/cases_test/${DIR_NAME}/add_files/* ./build/wrappers | true

compile:
	echo "# ==================================================" ; \
	echo " Compile fortran files " ; \
	echo "# ==================================================" ; \
	for file_name in ${SOURCE}; do \
	  echo "# ==================================================" ; \
	  $(FC) $(FPIC) $(FCFLAGS) -B ./build/wrappers -c $(addsuffix .f90, $$file_name) -o $(addsuffix .o, $$file_name); \
	done
	mv *.mod ./build/wrappers
	echo "# ==================================================" ; \
	echo " Create library source " ; \
	echo "# ==================================================" ; \
	${LIBTOOL} build/wrappers/libsrc.a ${OBJ_TMP}
	echo "# ==================================================" ; \
	echo " Done!							*" ; \

#==================================================
#		Run Python file
#==================================================

run:
	echo "# =================================================" ;\
	echo " Run Python file " ; \
	echo "# =================================================" ; \
	python ./bin/${DIR_NAME}.py 

#==================================================
#		Save wrapped files
#==================================================

save:
	echo "# =================================================" ;\
	echo " Save wrapped files " ; \
	echo "# =================================================" ; \
	mkdir -p ./build/wrapped_files/${DIR_NAME}
	rm -rf ./build/wrapped_files/${DIR_NAME}
	cp -R ./build/wrappers ./build/wrapped_files/${DIR_NAME} 

#==================================================
#			Clean
#==================================================

cleanstep:
	@echo "================================================================================"
	@echo "  Clean step wrapper files                                                      *"
	@echo "================================================================================"
	@rm -rf build/wrappers/f90wrap*
	@rm -rf build/wrappers/*.py	
	@rm -rf build/wrappers/${MODULE_NAME}/
	@rm -rf build/wrappers/src.linux*/
	@rm -rf build/wrappers/*.so

cleancomp:
	@echo "================================================================================"
	@echo "  Clean compilation files                                                      *"
	@echo "================================================================================"
	@rm -rf build/wrappers/*.o
	@rm -rf build/wrappers/*.mod
	@rm -rf build/wrappers/f90wrap*
	@rm -rf build/wrappers/*.py	
	@rm -rf build/wrappers/${MODULE_NAME}/
	@rm -rf build/wrappers/src.linux*/
	@rm -rf build/wrappers/*.so
	@rm -rf build/wrappers/*.a

clean:
	@echo "================================================================================"
	@echo "  Clean all wrapper files                                                     *"
	@echo "================================================================================"
	@rm -rf build/wrappers/*

cleanall: clean cleansave
	@echo "================================================================================"
	@echo "  Clean all files                                                      *"
	@echo "================================================================================"

cleansave:
	@echo "================================================================================"
	@echo "  Clean saving files                                                      *"
	@echo "================================================================================"
	@rm -rf build/wrapped_files/*


