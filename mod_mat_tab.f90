module mod_mat_tab

implicit none

integer, parameter, public :: PR=selected_real_kind(8)
real(PR), parameter, public :: Pi=2.0_PR*Asin(1.0_PR)

integer :: outtabE=0
integer :: outtabr=0
integer :: outtabT=0
integer :: outtabP=0

!!!---pour interpolations:
integer :: ir1,ir2,iT1,iT2,iTsat1, iTsat2
real(PR) :: zNE,zNW,zSE,zSW

type mat_tab
   character(len=2) :: nom
   integer :: Z
   real(PR) :: A
   real(PR) :: rhos, sig0, lamb0
   integer :: mode=1
   integer :: N_nu, N_T, N_sat, N_rho
   real(PR) :: rho0, rhof, T0, Tf, DT, DTm1, nu0, nuf, Dlogr, Dlogrm1, logr0, logrf
   real(PR) :: Emin, Emax, Pmin, Pmax, Tsatmin, Tsatmax, Esatmin, Esatmax
   real(PR), allocatable :: nu(:), T(:), rho(:)   
   real(PR), allocatable :: P(:,:), E(:,:), Cv(:,:), Cp(:,:), g(:,:)
   real(PR), allocatable :: c(:,:), Zb(:,:), sig(:,:), lamb(:,:), Seeb(:,:), pinf(:,:)
   real(PR), allocatable :: Tm(:)
   real(PR), allocatable :: T_sat(:), P_sat(:)
   real(PR), allocatable :: nu_sat(:,:), rho_sat(:,:), E_sat(:,:), Cv_sat(:,:), Cp_sat(:,:), g_sat(:,:)
   real(PR), allocatable :: c_sat(:,:), Zb_sat(:,:), sig_sat(:,:), lamb_sat(:,:), Seeb_sat(:,:)
   real(PR), allocatable :: pinf_sat(:,:)
end type mat_tab


type(mat_tab), allocatable :: mat(:)


contains

subroutine read_table(im,file_in)

   implicit none
   integer, intent(in) :: im
   character(len=*) :: file_in
   real(PR) :: P1, T1, E1, Cv1, Cp1, pinf, nu
   integer :: i,j,id, ierr 
   character(len=50) :: dump


   print*, ' START READ TABLE ... '
   id=12
   open(unit=id ,file=file_in,status='old',action='read',iostat=ierr)
   if(ierr.ne.0)then
     print*, " PROBLEM OPENING", file_in, " :( "
     stop
   else
     print*, " FILE", file_in, " OPENED :) "
   endif

   read(id,'(A)') dump

   read(id,'(A)') dump
   i=index(dump,'for') ; read(dump(i+4:),*) mat(im)%nom
   print*, 'nom=', mat(im)%nom

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%mode
   print*, 'mode=', mat(im)%mode

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%Z
   print*, 'Z=', mat(im)%Z

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%A
   print*, 'A=', mat(im)%A

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%rhos
   print*, 'rhos=', mat(im)%rhos

   if(mat(im)%mode.eq.1)then
      read(id,'(A)') dump
      i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%sig0
      print*, 'sig0=', mat(im)%sig0

      read(id,'(A)') dump
      i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%lamb0
      print*, 'lamb0=', mat(im)%lamb0
   endif


   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%N_nu
   print*, 'N_nu=', mat(im)%N_nu 
   mat(im)%N_rho=mat(im)%N_nu

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%rho0
   print*, 'rho0=', mat(im)%rho0

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%rhof
   print*, 'rhof=', mat(im)%rhof

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%N_T
   print*, 'N_T=', mat(im)%N_T

   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%T0
   print*, 'T0=', mat(im)%T0
 
   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%Tf
   print*, 'Tf=', mat(im)%Tf
 
   if(mat(im)%mode.eq.1)then
     allocate(mat(im)%T(1:mat(im)%N_T))
     read(id,*) dump
     do i=1,mat(im)%N_T
       read(id,'(A)') dump
       read(dump(3:),*) mat(im)%T(i)
       !print*, dump 
     enddo
   endif

   read(id,'(A)') dump
   !print*, '1', dump

   read(id,'(A)') dump
   !print*, '2', dump

   read(id,'(A)') dump
   !print*, '3', dump

  allocate(mat(im)%P(1:mat(im)%N_nu,1:mat(im)%N_T))    ; mat(im)%P(:,:)=0.0_PR
  allocate(mat(im)%E(1:mat(im)%N_nu,1:mat(im)%N_T))    ; mat(im)%E(:,:)=0.0_PR
  allocate(mat(im)%Cv(1:mat(im)%N_nu,1:mat(im)%N_T))   ; mat(im)%Cv(:,:)=0.0_PR
  allocate(mat(im)%g(1:mat(im)%N_nu,1:mat(im)%N_T))    ; mat(im)%g(:,:)=0.0_PR
  allocate(mat(im)%c(1:mat(im)%N_nu,1:mat(im)%N_T))    ; mat(im)%c(:,:)=0.0_PR
  allocate(mat(im)%Zb(1:mat(im)%N_nu,1:mat(im)%N_T))   ; mat(im)%Zb(:,:)=0.0_PR

  if(mat(im)%mode.eq.1)then
     allocate(mat(im)%nu(1:mat(im)%N_nu))             ; mat(im)%nu(:)=0.0_PR
     allocate(mat(im)%rho(1:mat(im)%N_rho))           ; mat(im)%rho(:)=0.0_PR
     allocate(mat(im)%sig(1:mat(im)%N_nu,1:mat(im)%N_T))  ; mat(im)%sig(:,:)=0.0_PR
     allocate(mat(im)%lamb(1:mat(im)%N_nu,1:mat(im)%N_T)) ; mat(im)%lamb(:,:)=0.0_PR
     allocate(mat(im)%Seeb(1:mat(im)%N_nu,1:mat(im)%N_T)) ; mat(im)%Seeb(:,:)=0.0_PR
     allocate(mat(im)%Cp(1:mat(im)%N_nu,1:mat(im)%N_T))   ; mat(im)%Cp(:,:)=0.0_PR
  else
     allocate(mat(im)%pinf(1:mat(im)%N_nu,1:mat(im)%N_T))  ; mat(im)%pinf(:,:)=0.0_PR
  endif



   do j=1,mat(im)%N_T
     read(id,'(A)') dump
     read(id,'(A)') dump
     read(id,'(A)') dump

     if(mat(im)%mode.eq.1)then
        do i=1,mat(im)%N_nu
           read(id,*) mat(im)%nu(i), P1, mat(im)%P(i,j), E1, mat(im)%E(i,j), mat(im)%Cv(i,j), mat(im)%Cp(i,j),&
                      mat(im)%c(i,j), mat(im)%g(i,j), mat(im)%Zb(i,j),&
                      mat(im)%sig(i,j), mat(im)%lamb(i,j), mat(im)%Seeb(i,j), pinf
        enddo
     else
        do i=1,mat(im)%N_nu
           read(id,*) mat(im)%P(i,j), mat(im)%E(i,j), mat(im)%Cv(i,j),&
                      mat(im)%c(i,j), mat(im)%g(i,j), mat(im)%Zb(i,j),&
                      mat(im)%pinf(i,j) 
        enddo
     endif
     read(id,'(A)') dump
     !print*, '7', dump 
     read(id,'(A)') dump
     !print*, '8', dump 
   enddo

   !!!---lecture de la courbe de saturation
   read(id,'(A)') dump
   read(id,'(A)') dump
   i=index(dump,'=') ; read(dump(i+1:),*) mat(im)%N_sat
 
   allocate(mat(im)%nu_sat(1:2,1:mat(im)%N_sat))   ; mat(im)%nu_sat(:,:)=0.0_PR
   allocate(mat(im)%rho_sat(1:2,1:mat(im)%N_sat))   ; mat(im)%rho_sat(:,:)=0.0_PR
   allocate(mat(im)%T_sat(1:mat(im)%N_sat))        ; mat(im)%T_sat(:)=0.0_PR
   allocate(mat(im)%P_sat(1:mat(im)%N_sat))        ; mat(im)%P_sat(:)=0.0_PR
   allocate(mat(im)%E_sat(1:2,1:mat(im)%N_sat))    ; mat(im)%E_sat(:,:)=0.0_PR
   allocate(mat(im)%Cv_sat(1:2,1:mat(im)%N_sat))   ; mat(im)%Cv_sat(:,:)=0.0_PR
   allocate(mat(im)%c_sat(1:2,1:mat(im)%N_sat))    ; mat(im)%c_sat(:,:)=0.0_PR
   allocate(mat(im)%g_sat(1:2,1:mat(im)%N_sat))    ; mat(im)%g_sat(:,:)=0.0_PR
   allocate(mat(im)%Zb_sat(1:2,1:mat(im)%N_sat))   ; mat(im)%Zb_sat(:,:)=0.0_PR

   if(mat(im)%mode.eq.1)then
      allocate(mat(im)%Cp_sat(1:2,1:mat(im)%N_sat))   ; mat(im)%Cp_sat(:,:)=0.0_PR
      allocate(mat(im)%sig_sat(1:2,1:mat(im)%N_sat))  ; mat(im)%sig_sat(:,:)=0.0_PR
      allocate(mat(im)%lamb_sat(1:2,1:mat(im)%N_sat)) ; mat(im)%lamb_sat(:,:)=0.0_PR
      allocate(mat(im)%Seeb_sat(1:2,1:mat(im)%N_sat)) ; mat(im)%Seeb_sat(:,:)=0.0_PR
   else
      allocate(mat(im)%pinf_sat(1:2,1:mat(im)%N_sat)) ; mat(im)%pinf_sat(:,:)=0.0_PR
   endif

   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump

   !!!---lecture de la courbe de bulle

   if(mat(im)%mode.eq.1)then
      do i=1, mat(im)%N_sat
        read(id,*) mat(im)%T_sat(i), mat(im)%P_sat(i), mat(im)%nu_sat(1,i),&   !, mat(im)%nu_sat(2,i),&
                                               mat(im)%E_sat(1,i),&    !, mat(im)%E_sat(2,i),&
                                               mat(im)%Cv_sat(1,i),&   !, mat(im)%Cv_sat(2,i),&
                                               mat(im)%Cp_sat(1,i),&   !, mat(im)%Cp_sat(2,i),&
                                               mat(im)%c_sat(1,i),&    !, mat(im)%c_sat(2,i),&
                                               mat(im)%g_sat(1,i),&    !, mat(im)%g_sat(2,i),&
                                               mat(im)%Zb_sat(1,i),&   !, mat(im)%Zb_sat(2,i),&
                                               mat(im)%sig_sat(1,i),&  !, mat(im)%sig_sat(2,i),&
                                               mat(im)%lamb_sat(1,i),& !, mat(im)%lamb_sat(2,i),&
                                               mat(im)%Seeb_sat(1,i),& !, mat(im)%Seeb_sat(2,i),&
                                               pinf            
      enddo
   else
      do i=1, mat(im)%N_sat
        read(id,*) mat(im)%T_sat(i), mat(im)%P_sat(i), mat(im)%nu_sat(1,i),&   !, mat(im)%nu_sat(2,i),&
                                               mat(im)%E_sat(1,i),&    !, mat(im)%E_sat(2,i),&
                                               mat(im)%Cv_sat(1,i),&   !, mat(im)%Cv_sat(2,i),&
                                               mat(im)%c_sat(1,i),&    !, mat(im)%c_sat(2,i),&
                                               mat(im)%g_sat(1,i),&    !, mat(im)%g_sat(2,i),&
                                               mat(im)%Zb_sat(1,i),&   !, mat(im)%Zb_sat(2,i),&
                                               mat(im)%pinf_sat(1,i)            
      enddo
   endif


   read(id,'(A)') dump !!! il ya N_sat+1 lignes
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump

   !!!---lecture de la courbe de goute

   if(mat(im)%mode.eq.1)then
      do i=1, mat(im)%N_sat
        read(id,*) mat(im)%T_sat(i), mat(im)%P_sat(i), mat(im)%nu_sat(2,i),&   !, mat(im)%nu_sat(2,i),&
                                               mat(im)%E_sat(2,i),&    !, mat(im)%E_sat(2,i),&
                                               mat(im)%Cv_sat(2,i),&   !, mat(im)%Cv_sat(2,i),&
                                               mat(im)%Cp_sat(2,i),&   !, mat(im)%Cp_sat(2,i),&
                                               mat(im)%c_sat(2,i),&    !, mat(im)%c_sat(2,i),&
                                               mat(im)%g_sat(2,i),&    !, mat(im)%g_sat(2,i),&
                                               mat(im)%Zb_sat(2,i),&   !, mat(im)%Zb_sat(2,i),&
                                               mat(im)%sig_sat(2,i),&  !, mat(im)%sig_sat(2,i),&
                                               mat(im)%lamb_sat(2,i),& !, mat(im)%lamb_sat(2,i),&
                                               mat(im)%Seeb_sat(2,i),& !, mat(im)%Seeb_sat(2,i),&
                                               pinf            
      enddo
   else
      do i=1, mat(im)%N_sat
        read(id,*) mat(im)%T_sat(i), mat(im)%P_sat(i), mat(im)%nu_sat(2,i),&   !, mat(im)%nu_sat(2,i),&
                                               mat(im)%E_sat(2,i),&    !, mat(im)%E_sat(2,i),&
                                               mat(im)%Cv_sat(2,i),&   !, mat(im)%Cv_sat(2,i),&
                                               mat(im)%c_sat(2,i),&    !, mat(im)%c_sat(2,i),&
                                               mat(im)%g_sat(2,i),&    !, mat(im)%g_sat(2,i),&
                                               mat(im)%Zb_sat(2,i),&   !, mat(im)%Zb_sat(2,i),&
                                               mat(im)%pinf_sat(2,i)            
      enddo
   endif

   read(id,'(A)') dump !!! il ya N_sat+1 lignes
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump
   read(id,'(A)') dump

   !!!---lecture de la courbe de fusion

   allocate(mat(im)%Tm(1:mat(im)%N_nu)) ; mat(im)%Tm(:)=0.0_PR
   do i=1, mat(im)%N_nu
     read(id,*) nu, mat(im)%Tm(i)
   enddo


   close(id)

   print*, ' READ FINI ! :) '

   mat(im)%rho_sat(:,:)=1.0_PR/mat(im)%nu_sat(:,:)

   if(mat(im)%mode.eq.1)then
      mat(im)%rho(:)=1.0_PR/mat(im)%nu(:)
      open(unit=1,file='density_tab.dat',status='replace')
      do i=1,mat(im)%N_rho
         write(1,*) mat(im)%rho(i)
      enddo
      close(1)
   endif



    mat(im)%nu0=1.0_PR/mat(im)%rhof ; mat(im)%nuf=1.0_PR/mat(im)%rho0
    !mat(im)%T0=mat(im)%T(1) ; mat(im)%Tf=mat(im)%T(mat(im)%N_T)
    mat(im)%DT=(mat(im)%Tf-mat(im)%T0)/real(mat(im)%N_T-1,PR) !!! N_T value=N_T-1 intervalles
    mat(im)%DTm1=1.0_PR/mat(im)%DT
    mat(im)%logr0=log(mat(im)%rho0) 
    mat(im)%logrf=log(mat(im)%rhof) 
    !mat(im)%logr0=log(mat(im)%rho(1))
    !mat(im)%logrf=log(mat(im)%rho(mat(im)%N_rho))
    mat(im)%Dlogr=(mat(im)%logrf-mat(im)%logr0)/real(mat(im)%N_rho-1,PR)
    mat(im)%Dlogrm1=1.0_PR/mat(im)%Dlogr
    mat(im)%Emin=minval(mat(im)%E(:,:))
    mat(im)%Emax=maxval(mat(im)%E(:,:))
    mat(im)%Pmin=minval(mat(im)%P(:,:))
    mat(im)%Pmax=maxval(mat(im)%P(:,:))
    mat(im)%Esatmax=maxval(mat(im)%E_sat(1,:))
    mat(im)%Esatmin=minval(mat(im)%E_sat(1,:))
    mat(im)%Tsatmax=maxval(mat(im)%T_sat(:))
    mat(im)%Tsatmin=minval(mat(im)%T_sat(:))


    if(mat(im)%mode.eq.2)then
     allocate(mat(im)%T(1:mat(im)%N_T))
     do i=1,mat(im)%N_T
       mat(im)%T(i)=mat(im)%T0+real(i-1,PR)*mat(im)%DT
     enddo
     allocate(mat(im)%rho(1:mat(im)%N_rho))
     do i=1,mat(im)%N_rho
       mat(im)%rho(i)=mat(im)%rho0*(mat(im)%rhof/mat(im)%rho0)**(real(i-1,PR)/real(mat(im)%N_rho-1,PR))
     enddo
    
    endif
 
    print*, 'Tmin=', mat(im)%T0, 'Tmax=', mat(im)%Tf
    print*, 'Emin=', mat(im)%Emin, 'Emax=', mat(im)%Emax
    print*, 'Pmin=', mat(im)%Pmin, 'Pmax=', mat(im)%Pmax
    print*, 'Pminloc=', minloc(mat(im)%P(:,:))
    print*, 'Pmaxloc=', maxloc(mat(im)%P(:,:))


end subroutine read_table

subroutine get_index_rT(im,rho,T,ir,iT)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T
   integer, intent(out) :: ir, iT
   real(PR) :: x,y
   !!! floor: integer(4) < 2e9

   call get_ir(im,rho,ir)

   call get_iT_T(im,T,iT)


end subroutine get_index_rT


subroutine get_index_rE(im,rho,u,ir,iT)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, u
   integer, intent(out) :: ir, iT
   real(PR) :: x, x1, x2, cr1, cr2, u1, u2
   integer :: i, iTmin, iTmax
   logical :: found 

   call get_ir(im,rho,ir)

   if(u.ge.mat(im)%Emax)then

     iT=mat(im)%N_T
     outtabE=1

   elseif(u.le.mat(im)%Emin)then

     iT=1
     outtabE=-1

   else

     x=log(rho)
     x1=mat(im)%logr0+real(ir-1,PR)*mat(im)%Dlogr 
     !x1=log(mat(im)%rho(ir)) 
     x2=x1+mat(im)%Dlogr   
     x=max(x,x1) ; x=min(x,x2)
  
     cr1=(x2-x)/(x2-x1)
     cr2=(x-x1)/(x2-x1)
  
     iTmin=1
     iTmax=mat(im)%N_T
     found=.false.
     i=0
  
     do while(.not.found)
  
        iT=(iTmin+iTmax)/2
  
        u1=cr1*mat(im)%E(ir1,iT)+cr2*mat(im)%E(ir2,iT)
        u2=cr1*mat(im)%E(ir1,iT+1)+cr2*mat(im)%E(ir2,iT+1)
  
        if(u1.gt.u)then
          iTmax=iT 
        elseif(u2.lt.u)then
          iTmin=iT
        else
          found=.true.
        endif
  
        i=i+1
  
        if(i.gt.990)then
            print*, iTmin, iTmax, u1, u2, u
        endif
        if(i.gt.1000)then
            print*, 'PROBLEM get_index_rE'
            stop
        endif
  
     enddo
  
     !print*, 'get_index_rE it=', i

   endif

   iT1=iT
   iT2=min(iT+1,mat(im)%N_T)

end subroutine get_index_rE

subroutine get_index_rP(im,rho,P,ir,iT)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, P
   integer, intent(out) :: ir, iT
   real(PR) :: x, x1, x2, cr1, cr2, P1, P2
   integer :: i, iTmin, iTmax
   logical :: found 

   call get_ir(im,rho,ir)

   if(P.gt.mat(im)%Pmax)then

     iT=mat(im)%N_T
     outtabP=1

   elseif(P.lt.mat(im)%Pmin)then

     iT=1
     outtabP=-1

   else

     x=log(rho)
     x1=mat(im)%logr0+real(ir-1,PR)*mat(im)%Dlogr 
     !x1=log(mat(im)%rho(ir)) 
     x2=x1+mat(im)%Dlogr   
     x=max(x,x1) ; x=min(x,x2)
  
     cr1=(x2-x)/(x2-x1)
     cr2=(x-x1)/(x2-x1)
  
     iTmin=1
     iTmax=mat(im)%N_T
     found=.false.
     i=0

     do while(.not.found)
  
        iT=(iTmin+iTmax)/2
  
        P1=cr1*mat(im)%P(ir1,iT)+cr2*mat(im)%P(ir2,iT)
        P2=cr1*mat(im)%P(ir1,iT+1)+cr2*mat(im)%P(ir2,iT+1)
  
        if(P1.gt.P)then
          iTmax=iT 
        elseif(P2.lt.P)then
          iTmin=iT
        else
          found=.true.
        endif
  
        i=i+1
  
        if(i.gt.990)then
            print*, iTmin, iTmax, P1, P2, P
        endif
        if(i.gt.1000)then
            print*, 'PROBLEM get_index_rP'
            stop
        endif
  
     enddo
  
     !print*, 'get_index_rP it=', i

   endif

   iT1=iT
   iT2=min(iT+1,mat(im)%N_T)

end subroutine get_index_rP

subroutine get_index_TP(im,T,P,ir,iT)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: T, P
   integer, intent(out) :: ir, iT
   real(PR) :: x, x1, x2, cr1, cr2, P1, P2
   integer :: i, irmin, irmax
   logical :: found 

   !print*, 'get_iT_T'
   call get_iT_T(im,T,iT)

   if(P.gt.mat(im)%Pmax)then

     ir=mat(im)%N_rho
     outtabP=1

   elseif(P.lt.mat(im)%Pmin)then

     ir=1
     outtabP=-1

   else

     x=T
     x1=mat(im)%T0+real(iT-1,PR)*mat(im)%DT
     x2=x1+mat(im)%DT 
     x=max(x,x1) ; x=min(x,x2)
  
     cr1=(x2-x)/(x2-x1)
     cr2=(x-x1)/(x2-x1)
  
     irmin=1
     irmax=mat(im)%N_rho
     found=.false.
     i=0

     !print*, 'start dowhile'

     do while(.not.found)
  
        ir=(irmin+irmax)/2
  
        P1=cr1*mat(im)%P(ir,iT1)+cr2*mat(im)%P(ir,iT2)
        P2=cr1*mat(im)%P(ir+1,iT1)+cr2*mat(im)%P(ir+1,iT2)
  
        if(P1.gt.P)then
          irmax=ir 
        elseif(P2.lt.P)then
          irmin=ir
        else
          found=.true.
        endif
  
        i=i+1
        !print*, 'i=',i
 
        if(i.gt.990)then
            print*, irmin, irmax, P1, P2, P
        endif
        if(i.gt.1000)then
            print*, 'PROBLEM get_index_PT'
            stop
        endif
  
     enddo

   endif

   ir1=ir
   ir2=min(ir+1,mat(im)%N_rho)

end subroutine get_index_TP




subroutine get_ir(im,rho,ir)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: rho
   integer, intent(out) :: ir
   real(PR) :: x
   !!! floor: integer(4) < 2e9

   if(rho.gt.mat(im)%rhof)then
    ir=mat(im)%N_rho
    outtabr=1
   elseif(rho.lt.mat(im)%rho0)then
    ir=1
    outtabr=-1
   else  
    x=log(rho) 
    ir=floor(min((x-mat(im)%logr0)/mat(im)%Dlogr,1.0e9_PR))+1
    ir=max(ir,1) ; ir=min(ir,mat(im)%N_rho)
   endif

   ir1=ir
   ir2=min(ir+1,mat(im)%N_rho)

end subroutine get_ir

subroutine get_iT_T(im,T,iT)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: T
   integer, intent(out) :: iT
   real(PR) :: y
   !!! floor: integer(4) < 2e9
   if(T.gt.mat(im)%Tf)then
    iT=mat(im)%N_T
    outtabT=1
   elseif(T.lt.mat(im)%T0)then
    iT=1
    outtabT=-1
   else  
    y=T 
    iT=floor(min((y-mat(im)%T0)/mat(im)%DT,1.0e9_PR))+1
    iT=max(iT,1) ; iT=min(iT,mat(im)%N_T)
   endif
   iT1=iT
   iT2=min(iT+1,mat(im)%N_T)

end subroutine get_iT_T


subroutine get_indexsat_T(im,T,iT)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: T
   integer, intent(out) :: iT
   real(PR) :: x, x1, x2, T1, T2
   integer :: i, iTmin, iTmax
   logical :: found 

   if(T.ge.mat(im)%Tsatmax)then

     iT=mat(im)%N_sat

   elseif(T.le.mat(im)%Tsatmin)then

     iT=1

   else

     iTmin=1
     iTmax=mat(im)%N_sat
     found=.false.
     i=0
  
     do while(.not.found)
  
        iT=(iTmin+iTmax)/2
 
         
        T1=mat(im)%T_sat(iT)   
        T2=mat(im)%T_sat(iT+1) 
  
        if(T1.gt.T)then
          iTmax=iT 
        elseif(T2.lt.T)then
          iTmin=iT
        else
          found=.true.
        endif
  
        i=i+1
  
        if(i.gt.990)then
            print*, iTmin, iTmax, T1, T2, T
        endif
        if(i.gt.1000)then
            print*, 'PROBLEM get_indexsat_E'
            stop
        endif
  
     enddo
  
     !print*, 'get_index_rE it=', i

   endif

   iTsat1=iT
   iTsat2=min(iT+1,mat(im)%N_T)

end subroutine get_indexsat_T



!!!----fonctions utilisateurs

!!! Grandeurs from rT : nécessite un appel à get_index_rT

real(PR) function P_from_rT(im,rho,T)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T

   zNE=mat(im)%P(ir2,iT2)
   zNW=mat(im)%P(ir1,iT2)
   zSE=mat(im)%P(ir2,iT1)
   zSW=mat(im)%P(ir1,iT1)

   P_from_rT=get_from_rT_lin(im,rho,T)
   !P_from_rT_log=get_from_rT_log(im,rho,T)

end function P_from_rT

real(PR) function P_from_rT_log(im,rho,T)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T

   zNE=mat(im)%P(ir2,iT2)
   zNW=mat(im)%P(ir1,iT2)
   zSE=mat(im)%P(ir2,iT1)
   zSW=mat(im)%P(ir1,iT1)

   P_from_rT_log=get_from_rT_log(im,rho,T)

end function P_from_rT_log

real(PR) function E_from_rT(im,rho,T)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T

   zNE=mat(im)%E(ir2,iT2)
   zNW=mat(im)%E(ir1,iT2)
   zSE=mat(im)%E(ir2,iT1)
   zSW=mat(im)%E(ir1,iT1)

   E_from_rT=get_from_rT_lin(im,rho,T)
   !E_from_rT=get_from_rT_log(im,rho,T)

end function E_from_rT

real(PR) function c_from_rT(im,rho,T)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T

   zNE=mat(im)%c(ir2,iT2)
   zNW=mat(im)%c(ir1,iT2)
   zSE=mat(im)%c(ir2,iT1)
   zSW=mat(im)%c(ir1,iT1)

   c_from_rT=get_from_rT_lin(im,rho,T)
   !P_from_rT_log=get_from_rT_log(im,rho,T)

end function c_from_rT

real(PR) function Cv_from_rT(im,rho,T)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T

   zNE=mat(im)%Cv(ir2,iT2)
   zNW=mat(im)%Cv(ir1,iT2)
   zSE=mat(im)%Cv(ir2,iT1)
   zSW=mat(im)%Cv(ir1,iT1)

   Cv_from_rT=get_from_rT_lin(im,rho,T)
   !Cv_from_rT=get_from_rT_log(im,rho,T)

end function Cv_from_rT

real(PR) function g_from_rT(im,rho,T)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T

   zNE=mat(im)%g(ir2,iT2)
   zNW=mat(im)%g(ir1,iT2)
   zSE=mat(im)%g(ir2,iT1)
   zSW=mat(im)%g(ir1,iT1)

   g_from_rT=get_from_rT_lin(im,rho,T)
   !P_from_rT_log=get_from_rT_log(im,rho,T)

end function g_from_rT

real(PR) function pinf_from_rT(im,rho,T)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T

   zNE=mat(im)%pinf(ir2,iT2)
   zNW=mat(im)%pinf(ir1,iT2)
   zSE=mat(im)%pinf(ir2,iT1)
   zSW=mat(im)%pinf(ir1,iT1)

   pinf_from_rT=get_from_rT_lin(im,rho,T)

end function pinf_from_rT

real(PR) function Z_from_rT(im,rho,T)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T

   zNE=mat(im)%Zb(ir2,iT2)
   zNW=mat(im)%Zb(ir1,iT2)
   zSE=mat(im)%Zb(ir2,iT1)
   zSW=mat(im)%Zb(ir1,iT1)

   Z_from_rT=get_from_rT_lin(im,rho,T)
   !Z_from_rT=get_from_rT_log(im,rho,T)

end function Z_from_rT

!!! Grandeurs from ru : nécessite un appel à get_index_ru

!real(PR) function P_from_ru(im,rho,u)
!
!   implicit none
!   integer, intent(in) :: im
!   real(PR), intent(in) :: rho, u
!
!   zNE=mat(im)%P(ir2,iT2)
!   zNW=mat(im)%P(ir1,iT2)
!   zSE=mat(im)%P(ir2,iT1)
!   zSW=mat(im)%P(ir1,iT1)
!
!   P_from_ru=get_from_ru_lin(im,rho,u)
!
!end function P_from_ru
!
!real(PR) function P_from_ru_log(im,rho,u)
!
!   implicit none
!   integer, intent(in) :: im
!   real(PR), intent(in) :: rho, u
!
!   zNE=mat(im)%P(ir2,iT2)
!   zNW=mat(im)%P(ir1,iT2)
!   zSE=mat(im)%P(ir2,iT1)
!   zSW=mat(im)%P(ir1,iT1)
!
!   P_from_ru_log=get_from_ru_log(im,rho,u)
!
!end function P_from_ru_log

real(PR) function T_from_ru(im,rho,u)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, u
   real(PR) :: u_NE,u_NW,u_SE,u_SW
   real(PR) :: u1,u2,T1,T2,lrho1,lrho2,lrho
   integer :: ir, iT

   call get_index_rE(im,rho,u,ir,iT)

   u_NE=mat(im)%E(ir2,iT2)
   u_NW=mat(im)%E(ir1,iT2)
   u_SE=mat(im)%E(ir2,iT1)
   u_SW=mat(im)%E(ir1,iT1)

   lrho1=log(mat(im)%rho(ir1))
   lrho2=log(mat(im)%rho(ir2))
   lrho=log(rho)

   T1=mat(im)%T(iT1)
   T2=mat(im)%T(iT2)

   u1=( (lrho2-lrho)*u_SW+(lrho-lrho1)*u_SE )*mat(im)%Dlogrm1
   u2=( (lrho2-lrho)*u_NW+(lrho-lrho1)*u_NE )*mat(im)%Dlogrm1

   T_from_ru=( T2*(u-u1) + T1*(u2-u) )/(u2-u1)

end function T_from_ru

real(PR) function T_from_rP(im,rho,P)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, P
   real(PR) :: P_NE,P_NW,P_SE,P_SW
   real(PR) :: P1,P2,T1,T2,lrho1,lrho2,lrho
   integer :: ir, iT

   call get_index_rP(im,rho,P,ir,iT)

   P_NE=mat(im)%P(ir2,iT2)
   P_NW=mat(im)%P(ir1,iT2)
   P_SE=mat(im)%P(ir2,iT1)
   P_SW=mat(im)%P(ir1,iT1)

   lrho1=log(mat(im)%rho(ir1))
   lrho2=log(mat(im)%rho(ir2))
   lrho=log(rho)

   T1=mat(im)%T(iT1)
   T2=mat(im)%T(iT2)

   P1=( (lrho2-lrho)*P_SW+(lrho-lrho1)*P_SE )*mat(im)%Dlogrm1
   P2=( (lrho2-lrho)*P_NW+(lrho-lrho1)*P_NE )*mat(im)%Dlogrm1

   T_from_rP=( T2*(P-P1) + T1*(P2-P) )/(P2-P1)

end function T_from_rP

real(PR) function r_from_TP(im,T,P)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: T, P
   real(PR) :: P_NE,P_NW,P_SE,P_SW
   real(PR) :: P1,P2,T1,T2,lrho1,lrho2,lrho, rho1, rho2
   integer :: ir, iT

   call get_index_TP(im,T,P,ir,iT)

   P_NE=mat(im)%P(ir2,iT2)
   P_NW=mat(im)%P(ir1,iT2)
   P_SE=mat(im)%P(ir2,iT1)
   P_SW=mat(im)%P(ir1,iT1)

   T1=mat(im)%T(iT1)
   T2=mat(im)%T(iT2)

   rho1=mat(im)%rho(ir1)
   rho2=mat(im)%rho(ir2) 
 
   lrho1=log(rho1)
   lrho2=log(rho2)

   P1=( (T2-T)*P_SW+(T-T1)*P_NW )*mat(im)%DTm1
   P2=( (T2-T)*P_SE+(T-T1)*P_NE )*mat(im)%DTm1

   !r_from_TP=( rho2*(P-P1) + rho1*(P2-P) )/(P2-P1)

   lrho=( lrho2*(P-P1) + lrho1*(P2-P) )/(P2-P1)
   r_from_TP=exp(lrho)

end function r_from_TP



!real(PR) function c_from_ru(im,rho,u)
!
!   implicit none
!   integer, intent(in) :: im
!   real(PR), intent(in) :: rho, u
!
!   zNE=mat(im)%c(ir2,iT2)
!   zNW=mat(im)%c(ir1,iT2)
!   zSE=mat(im)%c(ir2,iT1)
!   zSW=mat(im)%c(ir1,iT1)
!
!   c_from_ru=get_from_ru_lin(im,rho,u)
!
!end function c_from_ru
!
!real(PR) function Cv_from_ru(im,rho,u)
!
!   implicit none
!   integer, intent(in) :: im
!   real(PR), intent(in) :: rho, u
!
!   zNE=mat(im)%Cv(ir2,iT2)
!   zNW=mat(im)%Cv(ir1,iT2)
!   zSE=mat(im)%Cv(ir2,iT1)
!   zSW=mat(im)%Cv(ir1,iT1)
!
!   Cv_from_ru=get_from_ru_lin(im,rho,u)
!
!end function Cv_from_ru
!
!real(PR) function g_from_ru(im,rho,u)
!
!   implicit none
!   integer, intent(in) :: im
!   real(PR), intent(in) :: rho, u
!
!   zNE=mat(im)%g(ir2,iT2)
!   zNW=mat(im)%g(ir1,iT2)
!   zSE=mat(im)%g(ir2,iT1)
!   zSW=mat(im)%g(ir1,iT1)
!
!   g_from_ru=get_from_ru_lin(im,rho,u)
!
!end function g_from_ru
!
!real(PR) function Z_from_ru(im,rho,u)
!
!   implicit none
!   integer, intent(in) :: im
!   real(PR), intent(in) :: rho, u
!
!   zNE=mat(im)%Zb(ir2,iT2)
!   zNW=mat(im)%Zb(ir1,iT2)
!   zSE=mat(im)%Zb(ir2,iT1)
!   zSW=mat(im)%Zb(ir1,iT1)
!
!   Z_from_ru=get_from_ru_lin(im,rho,u)
!
!end function Z_from_ru
!

!!!----temperature de fusion

real(PR) function Tm_from_r(im,rho)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: rho
   real(PR) :: x,x1,x2,y1,y2

   x1=mat(im)%nu(ir1)
   x2=mat(im)%nu(ir2)

   y1=mat(im)%Tm(ir1)
   y2=mat(im)%Tm(ir2)

   x=1.0_PR/rho

   Tm_from_r=interplin1D(x1,x2,y1,y2,x)

end function Tm_from_r

subroutine nusat_from_T(im,T,nusat1,nusat2)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: T
   real(PR), intent(out) :: nusat1, nusat2
   real(PR) :: x,x1,x2,y1,y2

   x1=mat(im)%T_sat(iTsat1)
   x2=mat(im)%T_sat(iTsat2)
   x=T

   y1=mat(im)%nu_sat(1,iTsat1)
   y2=mat(im)%nu_sat(1,iTsat2)

   nusat1=interplin1D(x1,x2,y1,y2,x)

   y1=mat(im)%nu_sat(2,iTsat1)
   y2=mat(im)%nu_sat(2,iTsat2)

   nusat2=interplin1D(x1,x2,y1,y2,x)

end subroutine nusat_from_T

subroutine pinfsat_from_T(im,T,pinf1,pinf2)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: T
   real(PR), intent(out) :: pinf1, pinf2
   real(PR) :: x,x1,x2,y1,y2

   x1=mat(im)%T_sat(iTsat1)
   x2=mat(im)%T_sat(iTsat2)
   x=T

   y1=mat(im)%pinf_sat(1,iTsat1)
   y2=mat(im)%pinf_sat(1,iTsat2)

   pinf1=interplin1D(x1,x2,y1,y2,x)

   y1=mat(im)%pinf_sat(2,iTsat1)
   y2=mat(im)%pinf_sat(2,iTsat2)

   pinf2=interplin1D(x1,x2,y1,y2,x)

end subroutine pinfsat_from_T

subroutine gsat_from_T(im,T,g1,g2)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: T
   real(PR), intent(out) :: g1, g2
   real(PR) :: x,x1,x2,y1,y2

   x1=mat(im)%T_sat(iTsat1)
   x2=mat(im)%T_sat(iTsat2)
   x=T

   y1=mat(im)%g_sat(1,iTsat1)
   y2=mat(im)%g_sat(1,iTsat2)

   g1=interplin1D(x1,x2,y1,y2,x)

   y1=mat(im)%g_sat(2,iTsat1)
   y2=mat(im)%g_sat(2,iTsat2)

   g2=interplin1D(x1,x2,y1,y2,x)

end subroutine gsat_from_T

subroutine Esat_from_T(im,T,e1,e2)

   implicit none
   !type(mat_tab), intent(in) :: mat
   integer, intent(in) :: im
   real(PR), intent(in) :: T
   real(PR), intent(out) :: e1, e2
   real(PR) :: x,x1,x2,y1,y2

   x1=mat(im)%T_sat(iTsat1)
   x2=mat(im)%T_sat(iTsat2)
   x=T

   y1=mat(im)%E_sat(1,iTsat1)
   y2=mat(im)%E_sat(1,iTsat2)

   e1=interplin1D(x1,x2,y1,y2,x)

   y1=mat(im)%E_sat(2,iTsat1)
   y2=mat(im)%E_sat(2,iTsat2)

   e2=interplin1D(x1,x2,y1,y2,x)

end subroutine Esat_from_T




!!!---------  FONCTION D'INTERPOLATION


!!!---from rT

real(PR) function get_from_rT_lin(im,rho,T)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T
   real(PR) :: x,y,x1,x2,y1,y2

   x1=log(mat(im)%rho(ir1)) 
   x2=x1+mat(im)%Dlogr   

   y1=mat(im)%T(iT1)
   y2=y1+mat(im)%DT

   x=log(rho) ; x=max(x,x1) ; x=min(x,x2)
   y=T ; y=max(y,y1) ; y=min(y,y2)

   get_from_rT_lin=interplin2D(x1,x2,y1,y2,zNE,zNW,zSE,zSW,x,y)

end function get_from_rT_lin

real(PR) function get_from_rT_log(im,rho,T)

   implicit none
   integer, intent(in) :: im
   real(PR), intent(in) :: rho, T
   real(PR) :: x,y,x1,x2,y1,y2,lzNE,lzNW,lzSE,lzSW

   x=log(rho) ; x=max(x,mat(im)%logr0) ; x=min(x,mat(im)%logrf)
   x1=log(mat(im)%rho(ir1)) 
   x2=x1+mat(im)%Dlogr   

   y1=mat(im)%T(iT1)
   y2=y1+mat(im)%DT

   y=T ; y=max(y,y1) ; y=min(y,y2)

   lzNE=log(zNE)
   lzNW=log(zNW)
   lzSE=log(zSE)
   lzSW=log(zSW)

   get_from_rT_log=exp(interplin2D(x1,x2,y1,y2,lzNE,lzNW,lzSE,lzSW,x,y))

end function get_from_rT_log


!!!---from ru
!
!real(PR) function get_from_ru_lin(im,rho,u)
!
!   implicit none
!   integer, intent(in) :: im
!   real(PR), intent(in) :: rho, u
!   real(PR) :: x,y,x1,x2,y1,y2
!
!   x1=log(mat(im)%rho(ir1)) 
!   x2=x1+mat(im)%Dlogr   
!
!   y1=mat(im)%E(iT1)
!   y2=mat(im)%E(iT2)
!
!   x=log(rho) ; x=max(x,x1) ; x=min(x,x2)
!   y=u ; y=max(y,y1) ; y=min(y,y2)
!
!   get_from_ru_lin=interplin2D(x1,x2,y1,y2,zNE,zNW,zSE,zSW,x,y)
!
!end function get_from_ru_lin
!
!real(PR) function get_from_ru_log(im,rho,u)
!
!   implicit none
!   integer, intent(in) :: im
!   real(PR), intent(in) :: rho, u
!   real(PR) :: x,y,x1,x2,y1,y2,lzNE,lzNW,lzSE,lzSW
!
!   x=log(rho) ; x=max(x,mat(im)%logr0) ; x=min(x,mat(im)%logrf)
!   x1=log(mat(im)%rho(ir1)) 
!   x2=x1+mat(im)%Dlogr   
!
!   y1=mat(im)%E(iT1) ; y2=mat(im)%E(iT2)
!   y=u ; y=max(y,y1) ; y=min(y,y2)
!
!   lzNE=log(zNE)
!   lzNW=log(zNW)
!   lzSE=log(zSE)
!   lzSW=log(zSW)
!
!   get_from_ru_log=exp(interplin2D(x1,x2,y1,y2,lzNE,lzNW,lzSE,lzSW,x,y))
!
!end function get_from_ru_log
!



real(PR) function interplin1D(x1,x2,y1,y2,x)

   implicit none
   real(PR), intent(in) :: x1,x2,y1,y2
   real(PR), intent(in) :: x
   real(PR) :: invs, s1, s2

   invs=1.0_PR/(x2-x1)

   s1=(x2-x)
   s2=(x-x1)
   
   interplin1D=invs*(s1*y1+s2*y2)

end function interplin1D

real(PR) function interplin2D(x1,x2,y1,y2,zNE,zNW,zSE,zSW,x,y)

   implicit none
   real(PR), intent(in) :: x1,x2,y1,y2
   real(PR), intent(in) :: zNE,zNW,zSE,zSW
   real(PR), intent(in) :: x,y
   real(PR) :: invs, sNE, sNW, sSE, sSW

   invs=1.0_PR/((x2-x1)*(y2-y1))

   sNE=(y-y1)*(x-x1)
   sNW=(y-y1)*(x2-x)
   sSE=(y2-y)*(x-x1)
   sSW=(y2-y)*(x2-x)
   
   interplin2D=invs*(sNE*zNE+sNW*zNW+sSE*zSE+sSW*zSW)

end function interplin2D



subroutine check_tab

   implicit none

   print*, '  check tab...'

   if(outtabr.eq.1)then

     print*, '   UPPER BOUND OF TAB IN DENSITY !!!  :('

   elseif(outtabr.eq.-1)then
  
     print*, '   LOWER BOUND OF TAB IN DENSITY !!!  :('

   endif

   if(outtabT.eq.1)then

     print*, '   UPPER BOUND OF TAB IN TEMPERATURE !!!  :('

   elseif(outtabT.eq.-1)then
  
     print*, '   LOWER BOUND OF TAB IN TEMPERATURE !!!  :('

   endif

   if(outtabE.eq.1)then

     print*, '   UPPER BOUND OF TAB IN ENERGY !!!  :('

   elseif(outtabE.eq.-1)then
  
     print*, '   LOWER BOUND OF TAB IN ENERGY !!!  :('

   endif

   if(outtabr.eq.0.and.outtabT.eq.0.and.outtabE.eq.0)then

     print*, '  -> everything is OK :) !!!'      

   endif

   if(outtabP.eq.1)then

     print*, '   UPPER BOUND OF TAB IN PRESSURE !!!  :('

   elseif(outtabP.eq.-1)then
  
     print*, '   LOWER BOUND OF TAB IN PRESSURE !!!  :('

   endif

end subroutine check_tab


subroutine test_lecture_tab

   implicit none
   integer :: iout=6
   real(PR) :: Temp, rho, Ener, P
   integer :: i, it
 
   write(iout,*) ""
   write(iout,*) '==========================================='
   write(iout,*) '            READ TABLE !!!                 '
   write(iout,*) '==========================================='
   
   allocate(mat(1:1))
   call read_table(1,'tab_Cu_401.dat')
   
   print*, 'READ TABLE FINI !!!'
   
   print*, '--------------- test 1 -------------------'
   rho=287.0_PR ; Temp=3310.0_PR
   call get_index_rT(1,rho,Temp,i,iT)
   print*, 'ir=', i, 'vs 662  iT=', iT, 'vs 125'
   print*, '    rho=',rho, 'rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '    T=',Temp, 'T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1) 

   print*, '--------------- test 2 -------------------'
   rho=8.93e-2_PR ; Temp=201.0_PR
   call get_index_rT(1,rho,Temp,i,iT)
   print*, 'ir=', i, 'vs 1  iT=', iT, 'vs 1'
   print*, '    rho=',rho, 'rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '    T=',Temp, 'T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1) 

   print*, '--------------- test 3 -------------------'
   rho=8.80e3_PR ; Temp=10170.0_PR
   call get_index_rT(1,rho,Temp,i,iT)
   print*, 'ir=', i, 'vs 942  iT=', iT, 'vs 399'
   print*, '    rho=',rho, 'rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '    T=',Temp, 'T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1) 
   print*, 'P1(iT)=', mat(1)%P(i,iT)  , 'P2(iT)=',mat(1)%P(i+1,iT)
   print*, 'P1(iT+1)=', mat(1)%P(i,iT+1), 'P2(iT+1)=', mat(1)%P(i+1,iT+1)
   !call get_from_rT_lin(1,rho,Temp,i,iT,mat(1)%P,P)
   P=P_from_rT(1,rho,Temp)
   print*, 'Plin=', P
   !call get_from_rT_log(1,rho,Temp,i,iT,mat(1)%P,P)
   P=P_from_rT_log(1,rho,Temp)
   print*, 'Plog=', P
   print*, '--------------- test 4 -------------------'
   rho=1.0e12_PR ; Temp=1e20_PR
   call get_index_rT(1,rho,Temp,i,iT)
   print*, 'ir=', i, 'vs 1000  iT=', iT, 'vs 400'
   print*, '    rho=',rho, 'rho1=', mat(1)%rho(i) 
   print*, '    T=',Temp, 'T1=', mat(1)%T(iT)
   !call get_from_rT_lin(1,rho,Temp,i,iT,mat(1)%P,P)
   P=P_from_rT(1,rho,Temp)
   print*, 'Plin=', P, 'vs ', mat(1)%P(i,iT)
   !call get_from_rT_log(1,rho,Temp,i,iT,mat(1)%P,P)
   P=P_from_rT_log(1,rho,Temp)
   print*, 'Plog=', P, 'vs ', mat(1)%P(i,iT)
   print*, '--------------- test 5 -------------------'
   rho=1.0e-12_PR ; Temp=1e-20_PR
   call get_index_rT(1,rho,Temp,i,iT)
   print*, 'ir=', i, 'vs 1  iT=', iT, 'vs 1'
   print*, '    rho=',rho, 'rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '    T=',Temp, 'T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1) 
   !call get_from_rT_lin(1,rho,Temp,i,iT,mat(1)%P,P)
   P=P_from_rT(1,rho,Temp)
   print*, 'Plin=', P, 'vs ', mat(1)%P(i,iT)
   !call get_from_rT_log(1,rho,Temp,i,iT,mat(1)%P,P)
   P=P_from_rT_log(1,rho,Temp)
   print*, 'Plog=', P, 'vs', mat(1)%P(i,iT)
   print*, '--------------- test 6 -------------------'
   rho=1.57e-1_PR ; Temp=4500.0_PR
   print*, '    rho=',rho, 'T=', Temp
   print*, '> get index_rT:'
   call get_index_rT(1,rho,Temp,i,iT)
   Ener=0.5_PR*(mat(1)%E(i,iT)+mat(1)%E(i,iT+1))
   print*, '   ir=', i, 'iT=', iT
   print*, '   E=', Ener
   print*, '> get index rE:'
   call get_index_rE(1,rho,Ener,i,iT)
   print*, 'ir=', i, ' iT=', iT
   print*, '... check if indexes are same !'
   print*, '    rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '    T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1), 'vs', 3880.0_PR  
   print*, '    E1(iT)=',mat(1)%E(i,iT), 'E2(iT)=',mat(1)%E(i+1,iT)
   print*, '    E1(iT+1)=',mat(1)%E(i,iT+1), 'E2(iT+1)=',mat(1)%E(i+1,iT+1)
   print*, '--------------- test 7 -------------------'
   rho=7.1e3_PR ; Ener=1.0e20_PR  !Temp=3880.0_PR
   print*, 'test Emax E=',Ener
   print*, 'Emax=', mat(1)%Emax
   call get_index_rE(1,rho,Ener,i,iT)
   print*, 'ir=', i, 'vs 924  iT=', iT, 'vs 400'

   print*, '--------------- test 8 -------------------'
   rho=1.57e-1_PR ; Temp=4501.0_PR
   print*, '    rho=',rho, 'T=', Temp
   print*, '> get index_rT:'
   call get_index_rT(1,rho,Temp,i,iT)
   print*, '    rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '    T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1) 

   !P=0.5_PR*(mat(1)%P(i,iT)+mat(1)%P(i,iT+1))
   print*, '   ir=', i, 'iT=', iT
   print*, '   >E_from_rT:'
   Ener=E_from_rT(1,rho,Temp)
   print*, '   E=', Ener
   print*, '   > get_index_rE:'
   call get_index_rE(1,rho,Ener,i,iT)
   print*, '   ir=', i, ' iT=', iT
   print*, '   > T_from_ru:'
   print*, '   T=', T_from_ru(1,rho,Ener)
   print*, '   ... check if indexes are same !'
   print*, '   rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '   T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1) 
   print*, '   E1(iT)=',mat(1)%E(i,iT), 'E2(iT)=',mat(1)%E(i+1,iT)
   print*, '   E1(iT+1)=',mat(1)%E(i,iT+1), 'E2(iT+1)=',mat(1)%E(i+1,iT+1)


   print*, '--------------- test 9 -------------------'
   rho=1.57e-1_PR ; Temp=4501.0_PR
   print*, '   rho=',rho, 'T=', Temp
   print*, '   > get index_rT:'
   call get_index_rT(1,rho,Temp,i,iT)
   print*, '   rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '   T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1) 

   !P=0.5_PR*(mat(1)%P(i,iT)+mat(1)%P(i,iT+1))
   print*, '   ir=', i, 'iT=', iT
   print*, '   P_from_rT:'
   P=P_from_rT(1,rho,Temp)
   print*, '   P=', P
   print*, '   > get index rP:'
   call get_index_rP(1,rho,P,i,iT)
   print*, '   ir=', i, ' iT=', iT
   print*, '   > T_from_rP:'
   print*, '   T=', T_from_rP(1,rho,P)
   print*, '   ... check if indexes are same !'
   print*, '   rho1=', mat(1)%rho(i),'rho2=',mat(1)%rho(i+1) 
   print*, '   T1=', mat(1)%T(iT),'T2=',mat(1)%T(iT+1) 
   print*, '   P1(iT)=',mat(1)%P(i,iT), 'P2(iT)=',mat(1)%P(i+1,iT)
   print*, '   P1(iT+1)=',mat(1)%P(i,iT+1), 'P2(iT+1)=',mat(1)%P(i+1,iT+1)

   print*, '--------------- test 10 -------------------'
   rho=1.0e-3_PR*mat(1)%rhos ; P=1.0e5_PR
   print*, '   rho=',rho, 'P=', P
   print*, '   > get index_rP:'
   call get_index_rP(1,rho,P,i,iT)
   print*, '   ir=', i, ' iT=', iT
   print*, '   > T_from_rP:'
   Temp=T_from_rP(1,rho,P)
   print*, '   T=', Temp 
   print*, '   > P_from_rT:'
   P=P_from_rT(1,rho,Temp)
   print*, '   P=', P
   print*, '   ... check if pressures are same !'

   print*, '--------------- test 11 -------------------'
   !rho=0.97e0_PR*mat(1)%rhos ; Temp=300.0_PR   !!!P=1.0e5_PR

   P=1.0e5_PR ; Temp=300.0_PR
   print*, '   P=',P, 'Temp=', Temp
   print*, '   > get index_TP:'
   call get_index_TP(1,Temp,P,i,iT)
   print*, '   ir=', i, ' iT=', iT
   print*, mat(1)%rho(i), mat(1)%rho(i+1)
   print*, mat(1)%P(i,iT), mat(1)%P(i+1,iT)
   print*, mat(1)%P(i,iT+1), mat(1)%P(i+1,iT+1)
   print*, '   > r_from_TP:'
   rho=r_from_TP(1,Temp,P)
   print*, '   rho=', rho
   print*, '   > get index_rT:'
   call get_index_rT(1,rho,Temp,i,iT)
   print*, '   ir=', i, ' iT=', iT
   print*, '   > P_from_rT:'
   print*, '   P=', P_from_rT(1,rho,Temp)


   call check_tab

end subroutine test_lecture_tab




end module mod_mat_tab
