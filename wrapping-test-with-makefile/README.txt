
- Copy source files and compile Fortran files:
	$ make files

- Wrapper code:
	$ make

- Save wrapped files in ./build/wrapped_files/
	$ make save

- Run Python file:
	$ source ./scripts/env.sh
	$ make run
(!!!It is first necessary to modify the path of PYTHONPATH in ./scripts/env.sh file).

**Clean files before recompiling the /src files.

- Clean step wrapper files:
	$ make cleanstep

- Clean compilation files:
	$ make cleancomp

- Clean all wrapper files:
	$ make clean

- Clean saving files:
	$ make cleansave

- Clean all wrapper and saving files:
	$ make cleanall

**Modify the value of variable ${CASE} in Makefile for each case of source files. 

