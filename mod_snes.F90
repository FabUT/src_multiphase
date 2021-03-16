module mod_snes

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscts.h>
!!!!!#include <petsc/finclude/petscviewer.h>

use petsc
use petscsys
use petscis
use petscvec
use petscmat
use petscpc
use petscksp
use petscsnes
!!!!!use petscviewer

!use mod_data, only : iout, rang, nb_procs, imatrix, PR, Pi, nod, elem, edge, mesh, Input, Line, Ncolmax,&
!                      output_folder, folder_loc, folder_list


use mod_data, only : PR, iout
!
implicit none
!
!!!  Variables PETSC :

PetscViewer    ::  lab
SNES  ::         snes  ! nonlinear solver context 
KSP   ::         ksp   ! linear solver context 
PC    ::         pc    ! preconditioner context 
Vec   ::         x,r,b ! solution, residual vectors 
Mat   ::         J     ! Jacobian matrix  


procedure(F_test), pointer :: Function_snes => NULL() 
procedure(J_test), pointer :: Jacobian_snes => NULL() 

integer :: N_snes
integer, allocatable :: ix(:)
real(PR), allocatable :: Jac(:,:)
real(PR), allocatable :: F(:)
real(PR), allocatable :: x0(:)
real(PR), allocatable :: b0(:)
logical :: restart_snes=.false. !!! si taille du problème change

contains


subroutine init_SNES(N,x_in,b_in,F_function,J_function)

   implicit none
   integer, intent(in) :: N !!! taille du problème
   real(PR), intent(in) :: b_in(1:N) !!! terme source
   real(PR), intent(inout) :: x_in(1:N) !!! solution
   procedure(F_test) :: F_function !!! fonction
   procedure(J_test) :: J_function !!! jacobienne
   integer :: i, ierr

   if(restart_snes)then
    deallocate(ix)
    deallocate(x0)
    deallocate(b0)
    deallocate(Jac)
    deallocate(F)
    call Vecdestroy(x,ierr)
    call Vecdestroy(r,ierr)
    call Vecdestroy(b,ierr)
    call Matdestroy(J,ierr)
    call snesdestroy(snes,ierr)
   endif

   write(iout,*) '  Initialization de SNES'

   N_snes=N
   allocate(ix(1:N)) ; forall(i=1:N) ix(i)=i-1
   allocate(x0(1:N)) ; x0=0.0_PR
   allocate(b0(1:N)) ; b0=0.0_PR
   allocate(Jac(1:N,1:N)) ; Jac=0.0_PR
   allocate(F(1:N)) ; F=0.0_PR

   call PetscInitialize('option.dat',ierr)

   !!!--- Create nonlinear solver context
   call SNESCreate(PETSC_COMM_SELF,snes,ierr)

   if(ierr.ne.0) call crash('snescreate')

   call SNESSetFromOptions(snes,ierr)

   if(ierr.ne.0) call crash('snesiptions')

   !!!--- Create vectors for solution and nonlinear function
   !call VecCreate(PETSC_COMM_WORLD,x,ierr)
   !call VecSetSizes(x,PETSC_DECIDE,2,ierr)
   !call VecSetFromOptions(x,ierr)
   !call VecDuplicate(x,r,ierr)

   call VecCreateSeq(PETSC_COMM_SELF,N,x,ierr)
   call VecSet(x,0.0_PR,ierr)
   call VecDuplicate(x,r,ierr)
   call VecDuplicate(x,b,ierr)
   if(ierr.ne.0) call crash('vectors')


   !!!--- Create Jacobian matrix data structure
   call MatCreate(PETSC_COMM_SELF,J,ierr)
   if(ierr.ne.0) call crash('matcreate')

   call MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,N,N,ierr)
   if(ierr.ne.0) call crash('matsetsize')

   ! !!!call MatSetFromOptions(J,ierr)
   call MatSetUp(J,ierr)
   if(ierr.ne.0) call crash('matsetup')

   !do i=1,N
   ! call MatSetValues(J,1,ix(i),N,A(i)%I(1:N),A(i)%C(1:N),INSERT_VALUES,ierr)    
   !enddo

   restart_snes=.false.

   !!!---definition de la fonction et de la jacobienne :

   Function_snes => F_function
   Jacobian_snes => J_function

   call SNESSetFunction(snes,r,FormFunction,0,ierr)
   if(ierr.ne.0) call crash('function')

   call SNESSetJacobian(snes,J,J,FormJacobian,0,ierr)
   if(ierr.ne.0) call crash('Jacobian')

   write(iout,*)  '   fini :) !'

end subroutine init_SNES

subroutine solve_SNES(N,x_in,b_in)

   implicit none
   integer, intent(in) :: N
   real(PR), intent(in) :: b_in(1:N)
   real(PR), intent(inout) :: x_in(1:N)
   integer :: i, ierr

   !write(iout,*) '  Start SNES resolution !'

   !!!---set initial solution
   call VecSetValues(x,N,ix(1:N),x_in(1:N),INSERT_VALUES,ierr)
   !if(ierr.ne.0) call crash('SetValue x')

   !!!!---set source term
   !call VecSetValues(b,N,ix(1:N),b_in(1:N),INSERT_VALUES,ierr)
   !if(ierr.ne.0) call crash('SetValue b')

   call SNESSolve(snes,b,x,ierr)
   !if(ierr.ne.0) call crash('snessolve')

   call VecGetvalues(x,N,ix(1:N),x_in(1:N),ierr)

   !call VecView(x,PETSC_VIEWER_STDOUT_WORLD,ierr)
   !call VecView(r,PETSC_VIEWER_STDOUT_WORLD,ierr)

   !write(iout,*) '   Snes fini !'

end subroutine solve_SNES

subroutine FormFunction(snes,xx,ff,dummy)
implicit none
     SNES    :: snes !!! context snes
     Vec     :: xx,ff !!! input/ function vectors
     integer :: dummy !!! user defined function context
     integer :: ierr, i, N

     N=N_snes

     call VecGetvalues(xx,N,ix(1:N),x0(1:N),ierr)
     if(ierr.ne.0) call crash('GetValues FormFunction')

     call Function_snes(N,x0(1:N),F(1:N))

     call VecSetValues(ff,N,ix(1:N),f(1:N),INSERT_VALUES,ierr)
     if(ierr.ne.0) call crash('SetValues FormFunction')

     !call VecView(ff,PETSC_VIEWER_STDOUT_WORLD,ierr)

end subroutine FormFunction

subroutine FormJacobian(snes,xx,JJ,B,dummy)
implicit none
     SNES    :: snes  !!! context snes
     Vec     :: xx    !!! input vector
     Mat     :: JJ,B  !!! Jacobian matrix
     integer :: dummy !!! user defined function context
     integer :: ierr, i, N

     N=N_snes

     call VecGetvalues(xx,N,ix(1:N),x0(1:N),ierr)

     call Jacobian_snes(N,x0(1:N),Jac(1:N,1:N))

     do i=1,N
        call MatSetValues(JJ,1,ix(i),N,ix(1:N),Jac(i,1:N),INSERT_VALUES,ierr)    
        if(ierr.ne.0) call crash('MatsetValues FormJacobian')
     enddo
     do i=1,N
        call MatSetValues(B,1,ix(i),N,ix(1:N),Jac(i,1:N),INSERT_VALUES,ierr)    
        if(ierr.ne.0) call crash('MatsetValues FormJacobian2')
     enddo

     !call MatSetValues(JJ,N,ix(1:N),N,ix(1:N),J(1:N,1:N),INSERT_VALUES,ierr)    
     !call MatSetValues(B, N,ix(1:N),N,ix(1:N),J(1:N,1:N),INSERT_VALUES,ierr)

     call MatAssemblyBegin(JJ,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyEnd(JJ,MAT_FINAL_ASSEMBLY,ierr)

     call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)

     !call PetscViewerASCIIOpen(PETSC_COMM_SELF,'matrix.dat',lab,ierr)
     !call MatView(JJ,lab,ierr) ; if(ierr.ne.0) call crash('MatView')

end subroutine FormJacobian

subroutine F_test(N,x,f)

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)

     !!---exemple 1 : x(1)=1 ; x(2)=2
     !f(1)=x(1)**2+x(1)*x(2)-3.0_PR
     !f(2)=x(1)*x(2)+x(2)**2-6.0_PR
     !!---exemple 2: result: x=1.07716 ; y=0.0897632
     f(1)=sin(3.0_PR*x(1))+x(2)
     f(2)=x(1)*x(2)-12.0_PR*x(2)**2
    
end subroutine F_test

subroutine J_test(N,x,J)

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)

     !!!---exemple 1
     !J(1,1)=2.0_PR*x(1)+x(2) ; J(1,2)=x(1)
     !J(2,1)=x(2)             ; J(2,2)=x(1)+2.0_PR*x(2)

     !!!---exemple 2:
     J(1,1)=3.0_PR*cos(3.0_PR*x(1)) ; J(1,2)=1.0_PR
     J(2,1)=x(2)                    ; J(2,2)=x(1)-24.0_PR*x(2)

end subroutine J_test


subroutine crash(message)

   implicit none
   character(len=*), intent(in) :: message
   integer :: errcode, ierr
 
   write(*,*) '-----  PROBLEM mod_para :c --------'
   write(*,*) message
   write(*,*) '-----------------------------------'

   !call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)

   !if(ierr.ne.0) write(iout,*) 'X( PROBLEM during MPI_Abort !'

end subroutine crash


end module mod_snes
