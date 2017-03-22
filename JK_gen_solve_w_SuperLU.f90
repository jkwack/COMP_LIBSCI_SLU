!
!  Description: Create a tri-diagonal matrix and then solve it with SuperLU and LApack
!
!  Filename:   JK_gen_solve_w_SuperLU.f90
!  Update description :
!         1) Read M via command line arguments
!         2) Generate a CSR matrix
!         3) Solve it with SuperLU and LApack
!         4) Compute condition numbers of the matrix
!         5) Generate a CSR matrix for the reduced system
!         6) LU factorization with SuperLU and LApack
!         7) Compute condition numbers of the matrix for the reduced system
!
!         written by JaeHyuk Kwack (jkwack2@illinois.edu)
!                    Scientific and Engineering Applications Support
!                    National Center for Supercomputing Applications
!                    University of Illinois at Urbana-Champaign
! -----------------------------------------------------------------------

program main
use omp_lib
implicit none
#if defined(CRAYPAT)
include "pat_apif.h"
#endif
character(len=32)                   :: arg
integer                             :: M = 128
integer                             :: nR, nCoI
integer,dimension(:),allocatable    :: CoI,RoP
real*8,dimension(:),allocatable     :: Amat
real*8,dimension(:,:),allocatable   :: Bvec,X_SLU,X_LA
real*8                              :: t0,t1
integer                             :: i,j,k,id,ierr
integer                             :: idF, idT
real*8,parameter                    :: pi=3.141592653589793238
real*8,parameter                    :: h = 2.0d0*pi/8192.0d0

! For SuperLU
integer                             :: iopt,info,nrhs=1
integer*8                           :: f_factors

! For LApack
integer,dimension(:),allocatable    :: ipiv
real*8,dimension(:,:),allocatable   :: S
character                           :: TRANS = 'N'

! For condition numbers
character,dimension(2)              :: NORM =(/'1','I'/)
real*8,dimension(2)                 :: ANORM, RCOND
integer,dimension(:),allocatable    :: IWORK
real*8,dimension(:),allocatable     :: WORK
real*8                              :: dlange

interface
!
   subroutine Gen_CSR_1(M,Amat,Bvec,CoI,RoP,nR,nCoI,h)
   implicit none
   integer,intent(in)                              :: M
   real*8,dimension(:),allocatable,intent(inout)   :: Amat
   real*8,dimension(:,:),allocatable,intent(inout) :: Bvec
   integer,dimension(:),allocatable,intent(inout)  :: CoI,RoP
   integer,intent(inout)                           :: nR,nCoI
   real*8,intent(in)                               :: h
   end subroutine Gen_CSR_1
!
   subroutine Gen_CSR_2(M,Amat,X,CoI,RoP,nR,nCoI)
   implicit none
   integer,intent(in)                              :: M
   real*8,dimension(:),allocatable,intent(inout)   :: Amat
   real*8,dimension(:,:),intent(inout)             :: X
   integer,dimension(:),allocatable,intent(inout)  :: CoI,RoP
   integer,intent(inout)                           :: nR,nCoI
   end subroutine Gen_CSR_2
!
end interface

do i=1,iargc()
   call getarg(i,arg)
   if (i==1) read(arg,*) M
enddo

!
! Creating a CSR matrix and residual vectors
!
call Gen_CSR_1(M,Amat,Bvec,CoI,RoP,nR,nCoI,h)

!
! Using SuperLU
!
allocate(X_SLU(nR,4));     X_SLU(:,:) = Bvec(:,:)
! LU decomposition with SuperLU
iopt=1
#if defined(CRAYPAT)
call PAT_record(PAT_STATE_ON,ierr)
call PAT_region_begin(1,'SLU_LU_decom',ierr)
#endif
t0 = omp_get_wtime()
call c_fortran_dgssv_transposed(iopt,nR,nCoI,nrhs,Amat,CoI,RoP,X_SLU,nR,f_factors,info)
t1 = omp_get_wtime()
#if defined(CRAYPAT)
call PAT_region_end(1,ierr)
call PAT_record(PAT_STATE_OFF,ierr)
#endif
write(*,'(a)') "================================================================================="
write(*,'(a60,f12.6,a)') " Wall time for LU decomposition by SuperLU = ",t1-t0," seconds"
write(*,'(a)') "================================================================================="
! Triangular solve with SuperLU
iopt=2
#if defined(CRAYPAT)
call PAT_record(PAT_STATE_ON,ierr)
call PAT_region_begin(2,'SLU_LU_solve',ierr)
#endif
t0 = omp_get_wtime()
do i=1,4
   call c_fortran_dgssv_transposed(iopt,nR,nCoI,nrhs,Amat,CoI,RoP,X_SLU(:,i),nR,f_factors,info)
enddo
t1 = omp_get_wtime()
#if defined(CRAYPAT)
call PAT_region_end(2,ierr)
call PAT_record(PAT_STATE_OFF,ierr)
#endif
write(*,'(a)') "================================================================================="
write(*,'(a60,f12.6,a)') " Wall time for solving four vectors by SuperLU = ",t1-t0," seconds"
write(*,'(a)') "================================================================================="
! Free all the storage for SuperLU
iopt=3
call c_fortran_dgssv_transposed(iopt,nR,nCoI,nrhs,Amat,CoI,RoP,X_SLU,nR,f_factors,info)
! Print out results by SuperLU
write(*,*) 'x1 ='
write(*,'(4(f15.6,7x))') X_SLU(1,1:4)
write(*,'(4(f15.6,7x))') X_SLU(2,1:4)
write(*,*) 'xM ='
write(*,'(4(f15.6,7x))') X_SLU(nR-1,1:4)
write(*,'(4(f15.6,7x))') X_SLU(nR,1:4)

!
! Using LA pack
!
! Converting a CSR matrix to a dense matrix
allocate(S(nR,nR))
call JK_CSR_to_Dense(S,nR,Amat,CoI,RoP,nCoI)
! LU decomposition with LApack
allocate(ipiv(nR),X_LA(nR,4));         X_LA(:,:) = Bvec(:,:)
#if defined(CRAYPAT)
call PAT_record(PAT_STATE_ON,ierr)
call PAT_region_begin(3,'LA_LU_decom',ierr)
#endif
t0 = omp_get_wtime()
call DGETRF( nR, nR, S, nR, ipiv, info )
t1 = omp_get_wtime()
#if defined(CRAYPAT)
call PAT_region_end(3,ierr)
call PAT_record(PAT_STATE_OFF,ierr)
#endif
write(*,'(a)') "================================================================================="
write(*,'(a60,f12.6,a)') " Wall time for LU decomposition by LApack = ",t1-t0," seconds"
write(*,'(a)') "================================================================================="
! Triangular solve with LApack
#if defined(CRAYPAT)
call PAT_record(PAT_STATE_ON,ierr)
call PAT_region_begin(4,'LA_LU_solve',ierr)
#endif
t0 = omp_get_wtime()
call DGETRS( TRANS, nR, 4, S, nR, ipiv, X_LA, nR, info )
t1 = omp_get_wtime()
#if defined(CRAYPAT)
call PAT_region_end(4,ierr)
call PAT_record(PAT_STATE_OFF,ierr)
#endif
write(*,'(a)') "================================================================================="
write(*,'(a60,f12.6,a)') " Wall time for solving by LApack = ",t1-t0," seconds"
write(*,'(a)') "================================================================================="
! Print out results by LApack
write(*,*) 'x1 ='
write(*,'(4(f15.6,7x))') X_LA(1,1:4)
write(*,'(4(f15.6,7x))') X_LA(2,1:4)
write(*,*) 'xM ='
write(*,'(4(f15.6,7x))') X_LA(nR-1,1:4)
write(*,'(4(f15.6,7x))') X_LA(nR,1:4)

!
! Computing Condition Numbers with LApack
!
allocate(WORK(4*nR),IWORK(nR))
t0 = omp_get_wtime()
do i=1,2
   ANORM(i) = dlange( NORM(i),nR,NR,S,nR,WORK)
   call DGECON( NORM(i), nR, S, nR, ANORM(i), RCOND(i), WORK, IWORK, info )
enddo
t1 = omp_get_wtime()
write(*,'(a)') "================================================================================="
write(*,'(a60,f12.6,a)') " Wall time for DGECON by LApack = ",t1-t0," seconds"
write(*,'(a)') "================================================================================="
! Print out condition numbers
write(*,'(a30,e15.6,2(a,e15.6))') 'Based on 1-norm, ANORM = ',ANORM(1),', RCOND = ',RCOND(1),', COND = ',1.0d0/RCOND(1)
write(*,'(a30,e15.6,2(a,e15.6))') 'Based on inf-norm, ANORM = ',ANORM(2),', RCOND = ',RCOND(2),', COND = ',1.0d0/RCOND(2)

!
! Generating the reduced system
!
write(*,'(//,a)') "================================================================================="
write(*,'(a)') "================================================================================="
write(*,'(a)') "                       Generating the reduced system "
write(*,'(a)') "================================================================================="
write(*,'(a)') "================================================================================="
deallocate(Amat,Bvec,CoI,RoP)
call Gen_CSR_2(M,Amat,X_LA,CoI,RoP,nR,nCoI)

!
! Using SuperLU
!
! LU decomposition with SuperLU
iopt=1
#if defined(CRAYPAT)
call PAT_record(PAT_STATE_ON,ierr)
call PAT_region_begin(5,'SLU_LU_decom',ierr)
#endif
t0 = omp_get_wtime()
call c_fortran_dgssv_transposed(iopt,nR,nCoI,nrhs,Amat,CoI,RoP,X_SLU,nR,f_factors,info)
t1 = omp_get_wtime()
#if defined(CRAYPAT)
call PAT_region_end(5,ierr)
call PAT_record(PAT_STATE_OFF,ierr)
#endif
write(*,'(a)') "================================================================================="
write(*,'(a60,f12.6,a)') " Wall time for LU decomposition by SuperLU = ",t1-t0," seconds"
write(*,'(a)') "================================================================================="

!
! Using LA pack
!
! Converting a CSR matrix to a dense matrix
call JK_CSR_to_Dense(S,nR,Amat,CoI,RoP,nCoI)
! LU decomposition with LApack
#if defined(CRAYPAT)
call PAT_record(PAT_STATE_ON,ierr)
call PAT_region_begin(6,'LA_LU_decom',ierr)
#endif
t0 = omp_get_wtime()
call DGETRF( nR, nR, S, nR, ipiv, info )
t1 = omp_get_wtime()
#if defined(CRAYPAT)
call PAT_region_end(6,ierr)
call PAT_record(PAT_STATE_OFF,ierr)
#endif
write(*,'(a)') "================================================================================="
write(*,'(a60,f12.6,a)') " Wall time for LU decomposition by LApack = ",t1-t0," seconds"
write(*,'(a)') "================================================================================="

!
! Computing Condition Numbers with LApack
!
t0 = omp_get_wtime()
do i=1,2
   ANORM(i) = dlange( NORM(i),nR,NR,S,nR,WORK)
   call DGECON( NORM(i), nR, S, nR, ANORM(i), RCOND(i), WORK, IWORK, info )
enddo
t1 = omp_get_wtime()
write(*,'(a)') "================================================================================="
write(*,'(a60,f12.6,a)') " Wall time for DGECON by LApack = ",t1-t0," seconds"
write(*,'(a)') "================================================================================="
! Print out condition numbers
write(*,'(a30,e15.6,2(a,e15.6))') 'Based on 1-norm, ANORM = ',ANORM(1),', RCOND = ',RCOND(1),', COND = ',1.0d0/RCOND(1)
write(*,'(a30,e15.6,2(a,e15.6))') 'Based on inf-norm, ANORM = ',ANORM(2),', RCOND = ',RCOND(2),', COND = ',1.0d0/RCOND(2)

end program main
!
!
!
!
!
subroutine Gen_CSR_1(M,Amat,Bvec,CoI,RoP,nR,nCoI,h)
implicit none
integer,intent(in)                              :: M
real*8,dimension(:),allocatable,intent(inout)   :: Amat
real*8,dimension(:,:),allocatable,intent(inout) :: Bvec
integer,dimension(:),allocatable,intent(inout)  :: CoI,RoP
integer,intent(inout)                           :: nR,nCoI
real*8,intent(in)                               :: h
real*8,dimension(2,6)                           :: Blocked_Row
integer                                         :: i,j,k,id
integer                                         :: idF, idT

!
! Creating a CSR matrix and residual vectors
!
nR = 2*M
if (M==8192) then
   nCoI = 6*nR
else
   nCoI = 6*nR - 8
endif
allocate(Amat(nCoI),CoI(nCoI),RoP(nR+1),Bvec(nR,4))
!                            A                B                C
Blocked_Row(1,:) = (/  51.0d0,  9.0d0*h,108.0d0,   0.0d0*h, 51.0d0, -9.0d0*h/)
Blocked_Row(2,:) = (/-138.0d0,-18.0d0*h,  0.0d0,-108.0d0*h,138.0d0,-18.0d0*h/)
print "(2x,a,i6,a,i6,a,i8)",' M = ',M,", nR=",nR,", nCoI=",nCoI
id = 1;        RoP(1) = 1
do i=1,M
   do j=1,2
      if (M == 8192) then
         RoP((i-1)*2+j+1) = RoP((i-1)*2+j)+6
      else
         if (i == 1 .or. i==M) then
            RoP((i-1)*2+j+1) = RoP((i-1)*2+j)+4
         else
            RoP((i-1)*2+j+1) = RoP((i-1)*2+j)+6
         endif
      endif
      if (M==8192) then
         if (i==1) then
            Amat(id:id+5) = Blocked_Row(j,[3,4,5,6,1,2])
            idF = 2*(i-1)+1
            CoI(id:id+5) = (/idF,idF+1,idF+2,idF+3,2*M-1,2*M/)
         elseif (i==M) then
            Amat(id:id+5) = Blocked_Row(j,[5,6,1,2,3,4])
            idF = 2*(i-1)-1
            CoI(id:id+5) = (/1,2,idF,idF+1,idF+2,idF+3/)
         else
            Amat(id:id+5) = Blocked_Row(j,:)
            idF = 2*(i-1)-1
            CoI(id:id+5) = (/idF,idF+1,idF+2,idF+3,idF+4,idF+5/)
         endif
         id = id + 6
      else
         if (i==1) then
            Amat(id:id+3) = Blocked_Row(j,[3,4,5,6])
            idF = 2*(i-1)+1
            CoI(id:id+3) = (/idF,idF+1,idF+2,idF+3/)
            id = id + 4
         elseif (i==M) then
            Amat(id:id+3) = Blocked_Row(j,[1,2,3,4])
            idF = 2*(i-1)-1
            CoI(id:id+3) = (/idF,idF+1,idF+2,idF+3/)
            id = id + 4
         else
            Amat(id:id+5) = Blocked_Row(j,:)
            idF = 2*(i-1)-1
            CoI(id:id+5) = (/idF,idF+1,idF+2,idF+3,idF+4,idF+5/)
            id = id + 6
         endif
      endif

   enddo
enddo

Bvec(1:2,1:2) = Blocked_Row(1:2,1:2)
Bvec(nR-1:nR,3:4) = Blocked_Row(1:2,5:6)
!
end subroutine Gen_CSR_1
!
!
!
!
!
subroutine Gen_CSR_2(M,Amat,X,CoI,RoP,nR,nCoI)
implicit none
integer,intent(in)                              :: M
real*8,dimension(:),allocatable,intent(inout)   :: Amat
real*8,dimension(:,:),intent(inout)             :: X
integer,dimension(:),allocatable,intent(inout)  :: CoI,RoP
integer,intent(inout)                           :: nR,nCoI
real*8,dimension(4,6)                           :: Blocked_Row
integer                                         :: i,j,k,id
integer                                         :: idF, idT

!
! Creating a CSR matrix for the reduced system
!
nR = 2*M
nCoI = 6*nR
allocate(Amat(nCoI),CoI(nCoI),RoP(nR+1))
!
Blocked_Row(1,1:2) = (/  1.0d0,  0.0d0  /)
Blocked_Row(2,1:2) = (/  0.0d0,  1.0d0  /)
Blocked_Row(1:2,3:6) = X(1:2,1:4)
Blocked_Row(3,5:6) = (/  1.0d0,  0.0d0  /)
Blocked_Row(4,5:6) = (/  0.0d0,  1.0d0  /)
Blocked_Row(3:4,1:4) = X(nR-1:nR,1:4)
!do i=1,4
!   write(*,'(6(f15.6,x))') Blocked_Row(i,:)
!enddo

print "(2x,a,i6,a,i6,a,i8)",' M = ',M,", nR=",nR,", nCoI=",nCoI
id = 1;        RoP(1) = 1;          k=0
do i=1,M
   do j=1,2
      RoP((i-1)*2+j+1) = RoP((i-1)*2+j)+6
      k = k + 1
      if (k>4) k = k - 4
      if (i==1) then
         Amat(id:id+5) = Blocked_Row(k,[3,4,5,6,1,2])
         idF = 2*(i-1)+1
         CoI(id:id+5) = (/idF,idF+1,idF+2,idF+3,2*M-1,2*M/)
      elseif (i==M) then
         Amat(id:id+5) = Blocked_Row(k,[5,6,1,2,3,4])
         idF = 2*(i-1)-1
         CoI(id:id+5) = (/1,2,idF,idF+1,idF+2,idF+3/)
      else
         Amat(id:id+5) = Blocked_Row(k,:)
         idF = 2*(i-1)-1
         CoI(id:id+5) = (/idF,idF+1,idF+2,idF+3,idF+4,idF+5/)
      endif
      id = id + 6
   enddo
enddo
!write(*,'(6(f15.6,x))') Amat(:)
!write(*,'(6(i7,x))') CoI(:)
!write(*,'(20(i5,x))') RoP(:)

end subroutine Gen_CSR_2
!
!
!
!
!
subroutine JK_CSR_to_Dense(SDn,n,S,CoIWp,RoPWp,nnz)
implicit none
integer,intent(in)                  :: n,nnz
real*8,dimension(n,n),intent(inout) :: SDn
real*8,dimension(nnz),intent(in)    :: S
integer,dimension(nnz),intent(in)   :: CoIWp
integer,dimension(n+1),intent(in)   :: RoPWp
integer                             :: i,j,istart,iend,irow,icol

SDn(:,:) = 0.0d0
do i=1,n
   istart = RoPWp(i)
   iend = RoPWp(i+1)-1
   irow = i
   do j=istart,iend
      icol = CoIWp(j)
      SDn(irow,icol) = S(j)
   enddo
enddo

end subroutine JK_CSR_to_Dense
