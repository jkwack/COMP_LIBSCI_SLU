
# Path
PATH_JK  = .

# Source file
SOURCE_LAPACK = $(PATH_JK)/SA_LApack.f90

# Compilers
CC = cc 
FF = ftn

slu:
	rm -rf use_SLU
	rm -rf c_fortran_dgssv_transposed.o
	$(CC) -c c_fortran_dgssv_transposed.c
	$(FF) -cpp -fopenmp JK_gen_solve_w_SuperLU.f90 c_fortran_dgssv_transposed.o -o use_SLU

cleanall::
	rcsclean -q
	rm -f *.o
	rm -f $(PATH_JK)/*.o
	rm -f $(PATH_JK)/use_SLU


