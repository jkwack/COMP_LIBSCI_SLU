Description: Create a tri-diagonal matrix and then solve it with SuperLU and LApack on BlueWaters
	1) Read M via command line arguments
	2) Generate a CSR matrix
	3) Solve it with SuperLU and LApack
	4) Compute condition numbers of the matrix
	5) Generate a CSR matrix for the reduced system
	6) LU factorization with SuperLU and LApack
	7) Compute condition numbers of the matrix for the reduced system

 written by JaeHyuk Kwack (jkwack2@illinois.edu)
            Scientific and Engineering Applications Support
            National Center for Supercomputing Applications
            University of Illinois at Urbana-Champaign

How to build:
	source source_me	:: Loading required modules
	make slu		:: generating an executable
	pat_build -w use_SLU	:: generating an CrayPat-compatible executable

How to run:	
	qsub Job_SLU  
