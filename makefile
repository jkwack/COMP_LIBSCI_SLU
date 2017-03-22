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
	rm -f use_SLU


