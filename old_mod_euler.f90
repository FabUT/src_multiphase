module mod_euler


use mod_data, only : PR, iout !!!, Nt, Nl, Nx, Ny, Nz, t, dt

!use mod_snes, only : init_SNES, solve_SNES, restart_snes !SNESinterface

implicit none

integer :: i_sortie=-1
logical :: verb=.false.  

real(PR) :: dtmin=1.0e10_PR

integer :: ierr=0

!!!---routines
procedure(relax_p0), pointer :: relax_p => NULL() 
procedure(NOREC), pointer :: RECO => NULL()
procedure(MINMOD), pointer :: Limiteur => NULL()
procedure(src_euler1D), pointer :: SRC => NULL()
procedure(RK1), pointer :: Timescheme => NULL()

!!!---temps pour calcul des perfos
real(PR) :: Teuler, Trelax, t10,t1,t20,t2,t30,t3

!!!---nombre de variables conservatives: 4+Nl
integer :: Neq
real(PR), allocatable, dimension(:,:) :: U, dUdt

!!!---variables non conservatives : Nl*11
real(PR), allocatable, dimension(:,:,:) :: G, dGdt !!! Nl, 11, Nx

!!!---tableau et variables pour la relaxation de la pression
real(PR), allocatable :: flrhl0(:), el0(:), nul0(:), pl0(:), Zl(:), fl0(:), fl0_gl(:), coef1(:)
real(PR), allocatable :: xin(:), bin(:)
integer :: Nsnes=0
real(PR) :: pI0
integer :: itrelax=0
integer :: mode_relax=1 !!! 1: equation SGE (Ndanou 2015 eq 25 p 534)
                        !!! 2: système général tout EOS (Saurel 2009 p 1691)
integer :: mode_pI=1 !!! 1: pI=p ; 2: pI=sum(Zlpl)/sum(Zl) (cf Saurel 2009 p1691)
logical :: compute_c=.true. !!! desactive le calcul de c à l'entrée de relax_p

!!!---Propriétés matière

real(PR), allocatable :: mu(:)
real(PR), allocatable :: sigy(:)
real(PR), allocatable :: rho0(:)

!!!---EOS Stiffened gas Equation

!!!---eau:
!real(PR) :: g_sge_H2O=6.1_PR
!real(PR) :: p_sge_H2O=2.0e9_PR
!real(PR) :: mu_H2O=0.0_PR
!real(PR) :: sigy_H2O=0.0_PR
!real(PR) :: rho0_H2O=1000.0_PR
!!!---eau 
real(PR) :: g_sge_eau=4.4_PR
real(PR) :: p_sge_eau=6.0e8_PR
real(PR) :: mu_eau=0.0_PR
real(PR) :: sigy_eau=0.0_PR
real(PR) :: rho0_eau=1000.0_PR

!!!---air:
real(PR) :: g_sge_air=1.4_PR
real(PR) :: p_sge_air=0.0_PR
real(PR) :: mu_air=0.0_PR
real(PR) :: sigy_air=0.0_PR
real(PR) :: rho0_air=1.0_PR

!!!---Fe:
real(PR) :: g_sge_Fe=3.9_PR
real(PR) :: p_sge_Fe=43.6e9_PR
real(PR) :: mu_Fe=82.0e9_PR
real(PR) :: sigy_Fe=200.0e6_PR
real(PR) :: rho0_Fe=7860.0_PR

!!!---Al:
real(PR) :: g_sge_Al=3.5_PR
real(PR) :: p_sge_Al=32.0e9_PR
real(PR) :: mu_Al=26.0e9_PR
real(PR) :: sigy_Al=60.0e6_PR
real(PR) :: rho0_Al=2712.0_PR

!!!---Ti:
real(PR) :: g_sge_Ti=2.6_PR
real(PR) :: p_sge_Ti=44.0e9_PR
real(PR) :: mu_Ti=42.0e9_PR
real(PR) :: sigy_Ti=1030.0e6_PR
real(PR) :: rho0_Ti=4527.0_PR

!!!---Cu:
real(PR) :: g_sge_Cu=4.22_PR
real(PR) :: p_sge_Cu=34.2e9_PR
real(PR) :: mu_Cu=9.2e10_PR
real(PR) :: sigy_Cu=40.0e6_PR
real(PR) :: rho0_Cu=8900.0_PR



real(PR), allocatable :: p_sge(:), g_sge(:) 

contains

!!!===================== ROUTINE D'APPEL =====================

subroutine solve_euler(dt,cfl)

   use mod_data, only : Nl,Nx,Ny,Nz,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        dxm,dym,dzm
                        

   implicit none
   real(PR), intent(in) :: dt, cfl
   integer :: idir,i,j,k,iph,i2,im1,ip1,ip2
   real(PR) :: tloc, dtloc
  
   tloc=0.0_PR

   do while(tloc.lt.dt)

      dtloc=dtcfl(nx,ny,nz,dxm,c,vx,vy,vz,cfl)

      dtloc=min(dt-tloc,dtloc,dtmin) 

      !!!-------------------------  X  ----------------------------
 
      idir=1

      dUdt=0.0_PR
      dGdt=0.0_PR

      DO k=1,Nz ; DO j=1,Ny

      do i=1,Nx
      call primitive2conservative(Nl,U(:,i),G(:,:,i),&
                                  fl(:,i,j,k),Yl(:,i,j,k),rhl(:,i,j,k),&
                                  pl(:,i,j,k),elh(:,i,j,k),&
                                  al(:,:,i,j,k),bl(:,:,i,j,k),cl(:,:,i,j,k),&
                                  rh(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),&
                                  p(i,j,k),e(i,j,k),c(i,j,k))
      enddo
 
                  

      call Timescheme(Nl,Nx,dxm,U,G,dtloc)
 
      tloc=tloc+dtloc   

      do i=1,nx
      call conservative2primitive(Nl,U(:,i),G(:,:,i),&
                                  fl(:,i,j,k),Yl(:,i,j,k),rhl(:,i,j,k),&
                                  pl(:,i,j,k),elh(:,i,j,k),&
                                  al(:,:,i,j,k),bl(:,:,i,j,k),cl(:,:,i,j,k),&
                                  rh(i,j,k),vx(i,j,k),vy(i,j,k),vz(i,j,k),&
                                  p(i,j,k),e(i,j,k),c(i,j,k))
      enddo


      ENDDO ; ENDDO

   enddo !!! dowhile

end subroutine solve_euler


!!!======================== SUBROUTINE D'INITIALISATION =============================

subroutine init_euler(nx_i,ny_i,nz_i,Nl_i) 

   use mod_data, only : Nl,Nx,Ny,Nz,&
                        fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        Lx,Ly,Lz,&
                        dxm,dym,dzm,dx,dy,dz,&
                        xm,ym,zm,x,y,z

   implicit none
   integer, intent(in) :: nx_i,ny_i,nz_i,Nl_i
   integer :: i,j,k,iph

   write(iout,*) 'call init START !!!'


   Nx=nx_i
   Ny=ny_i
   Nz=Nz_i
   Nl=Nl_i
   
   !!!---MAILLAGE :
   
    Lx=1.0_PR ; Ly=1.0_PR ; Lz=1.0_PR 
   !Lx=0.017_PR ; Ly=1.0_PR ; Lz=1.0_PR

   write(iout,*) ""
   write(iout,fmt='(A)') "CREATION DU MAILLAGE"
   write(iout,*) ""
   write(iout,*) "Nx=",Nx, "Ny=", Ny, "Nz=",Nz 
   write(iout,*) "Lx=",Lx, "Ly=", Ly, "Lz=",Lz 
   write(iout,*) ""
   write(iout,*) '   allocation des grandeurs géométriques'
   allocate(x(1:nx))
   allocate(xm(0:nx))
   allocate(dx(1:nx))
   allocate(dxm(0:nx))
   allocate(y(1:ny))
   allocate(ym(0:ny))
   allocate(dy(1:ny))
   allocate(dym(0:ny))
   allocate(z(1:nz))
   allocate(zm(0:nz))
   allocate(dz(1:nz))
   allocate(dzm(0:nz))

   write(iout,*) '   valorisation'

   !!!---X:
   dx(:)=Lx/real(Nx,PR)
   xm(0)=0.0_PR
   do i=1,nx
     xm(i)=xm(i-1)+dx(i)
   enddo
   do i=1,nx
     x(i)=0.5_PR*(xm(i-1)+xm(i))
   enddo
   dxm(0)=2.0_PR*x(1)
   do i=1,nx-1
     dxm(i)=x(i+1)-x(i)
   enddo
   dxm(nx)=2.0_PR*(xm(nx)-x(nx))
 
  !!!---Y:
   dy(:)=Ly/real(Ny,PR)
   ym(0)=0.0_PR
   do i=1,ny
     ym(i)=ym(i-1)+dy(i)
   enddo
   do i=1,ny
     y(i)=0.5_PR*(ym(i-1)+ym(i))
   enddo
   dym(0)=2.0_PR*y(1)
   do i=1,ny-1
     dym(i)=y(i+1)-y(i)
   enddo
   dym(ny)=2.0_PR*(ym(ny)-y(ny))

   !!!---Z:
   dz(:)=Lz/real(Nz,PR)
   zm(0)=0.0_PR
   do i=1,nz
     zm(i)=zm(i-1)+dz(i)
   enddo
   do i=1,nz
     z(i)=0.5_PR*(zm(i-1)+zm(i))
   enddo
   dzm(0)=2.0_PR*z(1)
   do i=1,nz-1
     dzm(i)=z(i+1)-z(i)
   enddo
   dzm(nz)=2.0_PR*(zm(nz)-z(nz))


!!!---ALLOCATION DES TABLEAUX
   write(iout,*) ""
   write(iout,fmt='(A)') "ALLOCATION DES TABLEAUX"
   write(iout,*) ""
  
   !!!---grandeurs de mélange 
   !!!---vitesses hydrodynamiques:
   allocate(vx(1:nx,1:ny,1:nz))  ; vx(:,:,:)=0.0_PR
   allocate(vy(1:nx,1:ny,1:nz))  ; vy(:,:,:)=0.0_PR
   allocate(vz(1:nx,1:ny,1:nz))  ; vz(:,:,:)=0.0_PR
   !!!---tenseur des contraintes totales
   allocate(sig(1:3,1:3,1:nx,1:ny,1:nz)) ; sig(:,:,:,:,:)=0.0_PR

   !allocate(S11(1:nx,1:ny,1:nz)) ; S11(:,:,:)=0.0_PR
   !allocate(S22(1:nx,1:ny,1:nz)) ; S22(:,:,:)=0.0_PR
   !allocate(S33(1:nx,1:ny,1:nz)) ; S33(:,:,:)=0.0_PR
   !allocate(S12(1:nx,1:ny,1:nz)) ; S12(:,:,:)=0.0_PR
   !allocate(S13(1:nx,1:ny,1:nz)) ; S13(:,:,:)=0.0_PR
   !allocate(S23(1:nx,1:ny,1:nz)) ; S23(:,:,:)=0.0_PR
   !!!---densitée totale
   allocate(rh(1:nx,1:ny,1:nz))  ; rh(:,:,:)=0.0_PR
   !!!---énergie totale
   allocate(e(1:nx,1:ny,1:nz))   ; e(:,:,:)=0.0_PR
   !!!---vitesse du son totale
   allocate(c(1:nx,1:ny,1:nz))   ; c(:,:,:)=0.0_PR
!!!---pression totale
   allocate(p(1:nx,1:ny,1:nz))  ; p(:,:,:)=0.0_PR

   !!!---grandeurs phasiques

   !!!---propriétés physiques des matériaux (SG EOS ):
   allocate(p_sge(1:Nl)) ; allocate(g_sge(1:Nl))
   allocate(mu(1:Nl))    ; allocate(sigy(1:Nl))  
   allocate(rho0(1:Nl))
   !!!---champs :
   allocate(fl(1:Nl,1:nx,1:ny,1:nz))  ; fl(:,:,:,:)=1.0E-12_PR !!! volume fractions (alphal)
   allocate(rhl(1:Nl,1:nx,1:ny,1:nz)) ; rhl(:,:,:,:)=0.0_PR   !!! densité des phases
   allocate(Yl(1:Nl,1:nx,1:ny,1:nz))  ; Yl(:,:,:,:)=0.0_PR    !!! fraction massique
   allocate(pl(1:Nl,1:nx,1:ny,1:nz))  ; pl(:,:,:,:)=0.0_PR    !!! pression de phase
   allocate(elh(1:Nl,1:nx,1:ny,1:nz)) ; elh(:,:,:,:)=0.0_PR   !!! énergie hydro de phase
   allocate(ele(1:Nl,1:nx,1:ny,1:nz)) ; ele(:,:,:,:)=0.0_PR   !!! énergie elastique de phase
   allocate(al(1:Nl,1:3,1:nx,1:ny,1:nz)) ; al(:,:,:,:,:)=0.0_PR !!! composantes x des vecteurs de base
   allocate(bl(1:Nl,1:3,1:nx,1:ny,1:nz)) ; bl(:,:,:,:,:)=0.0_PR !!! composantes y des vecteurs de base
   allocate(cl(1:Nl,1:3,1:nx,1:ny,1:nz)) ; cl(:,:,:,:,:)=0.0_PR !!! composantes z des vecteurs de base

   !!!---tenseur des contraintes phasique (symmétrique) :
   allocate(sigl(1:Nl,1:3,1:3,1:nx,1:ny,1:nz)) ; sigl(:,:,:,:,:,:)=0.0_PR
   !!!--vecteurs L/R d'entrée pour la subroutine HLLC:
   !allocate(fl_L(1:Nl))     ; allocate(fl_R(1:Nl))   
   !allocate(rhl_L(1:Nl))    ; allocate(rhl_R(1:Nl))   
   !allocate(pl_L(1:Nl))     ; allocate(pl_R(1:Nl))   
   !allocate(elh_L(1:Nl))    ; allocate(elh_R(1:Nl))   
   !allocate(al_L(1:Nl,1:3)) ; allocate(al_R(1:Nl,1:3))   
   !allocate(bl_L(1:Nl,1:3)) ; allocate(bl_R(1:Nl,1:3))   
   !allocate(cl_L(1:Nl,1:3)) ; allocate(cl_R(1:Nl,1:3))   

   Neq=Nl+4 !!! nombre d'équations à résoudre partie conservative
   !!!!---vecteurs flux pour la partie conservative (eq 19 Ndanou 2015)
   !allocate(F1(1:Neq))     ; F1(:)=0.0_PR
   !!!!---vecteurs flux pour la partie non-conservative (eq 19 Ndanou 2015)
   !allocate(F2(1:Nl,1:11)) ; F2(:,:)=0.0_PR
   !allocate(H2(1:Nl,1:11)) ; H2(:,:)=0.0_PR
   !allocate(K2(1:Nl,1:11)) ; K2(:,:)=0.0_PR
   !allocate(M2(1:Nl,1:11)) ; M2(:,:)=0.0_PR

   !!!---vecteur des variavles conservatives :
   allocate(U(1:Neq,1:nx))    ; U(:,:)=0.0_PR
   allocate(dUdt(1:Neq,1:nx)) ; dUdt(:,:)=0.0_PR
   !!!---Partie non-conservative : fraction volumique + Finger tensor +
   !                               énergie hydro de chaque phase   
   allocate(G(1:Nl,1:11,1:nx))    ; G(:,:,:)=0.0_PR
   allocate(dGdt(1:Nl,1:11,1:nx)) ; dGdt(:,:,:)=0.0_PR

   !!!---initialisation du tenseur des contraintes

   write(iout,*) "   initialisation du tenseur des contraintes" 

   do k=1,Nz
     do j=1,Ny
       do i=1,Nx 
          !!!---composante x des vecteurs de base
          al(1:Nl,1,i,j,k)    = 1.0_PR
          al(1:Nl,2,i,j,k)    = 0.0_PR
          al(1:Nl,3,i,j,k)    = 0.0_PR
          !!!---composante y des vecteurs de base
          bl(1:Nl,1,i,j,k)    = 0.0_PR
          bl(1:Nl,2,i,j,k)    = 1.0_PR
          bl(1:Nl,3,i,j,k)    = 0.0_PR
          !!!---composante z des vecteurs de base
          cl(1:Nl,1,i,j,k)    = 0.0_PR
          cl(1:Nl,2,i,j,k)    = 0.0_PR
          cl(1:Nl,3,i,j,k)    = 1.0_PR
       enddo
      enddo
    enddo
 
   sigl(:,:,:,:,:,:)=0.0_PR

   !!!---initialisation du solveur SNES pour la relaxation de la pression

   call init_relaxp(Nl)


   !!!---vecteur 1D des variables primitives
   !Nprim=9+Nl*13

   !write(iout,*) '  Nprim=', Nprim  
   !allocate(W1(1:Nprim,1:Nx))
   !allocate(W_L(1:Nprim,1:Nx))
   !allocate(W_R(1:Nprim,1:Nx))

   !!!---Reconstruction

   RECO => MUSCL

   !!!---Limiteur:

   Limiteur => MINMOD

   !!!Limiteur => VANLEER

   !!!---Timescheme:

   Timescheme => RK2

   !!!---SRC:

   SRC => src_euler1D

   !!!---Temps:
   Teuler=0.0_PR
   Trelax=0.0_PR

   write(iout,*) 'call init FINI !!!'

end subroutine init_euler 

!!!========================== SOLVEUR DE RIEMANN ===============================

subroutine HLLC_flux(Nl,U_L,U_R,G_L,G_R,us,vs,ws,F1,F2)

    !!! input:
    !!! Nl : nombre de phases
    !!! fl      : fractions volumiques des phases
    !!! rl      : densité des phases
    !!! pl      : pression de phase
    !!! elh     : énergie hydrodynamique de phase
    !!! r       : densité totale (sum(fl rl))
    !!! u, v, w : vitesse x, y, z
    !!! c       : vitesse du son
    !!! E       : énergie spécifique totale (sum(Ylel + 1/2v^2)
    !!! sig(1:3)=s11, s12, s13 : composantes du tenseur de contraintes
    !!! al, bl, cl : composantes x,y,z des vecteurs de la base locale de la phase l
    !!! output :
    !!! F1 : flux des variables conservatives :( (fr)1*u, ... (fr)l*u ...,ru^2-s11,ruv-s12,ruw-s13,(rE-s11)u-s12v-s13w )   
    !!! F2,H2,K2,M2: flux des variables non conservatives 

    implicit none
  
    integer, intent(in) :: Nl
    real(PR), intent(in) :: U_L(1:Nl+4), U_R(1:Nl+4)
    real(PR), intent(in) :: G_L(1:Nl,1:11), G_R(1:Nl,1:11)  
    real(PR), intent(out) :: us,vs,ws
    real(PR), intent(out) :: F1(1:Nl+4)
    real(PR), intent(out) :: F2(1:Nl,1:11)
    !!!---Left/Right variables:   
    real(PR), dimension(Nl) :: fl_L, fl_R, Yl_L, Yl_R, pl_L, pl_R, rl_L, rl_R, elh_L, elh_R, ele_L, ele_R
    real(PR), dimension(Nl,3) :: al_L,bl_L,cl_L,al_R,bl_R,cl_R
    real(PR), dimension(Nl,3,3) :: sigl_L, sigl_R !!! s11, s12, s13
    real(PR), dimension(3,3) :: sig_L, sig_R !!! s11, s12, s13
    real(PR) :: r_L,r_R,vx_L,vx_R,vy_L,vy_R,vz_L,vz_R,e_L,e_R,c_L,c_R
    real(PR) :: S_R,S_L,C0,C1,C2,C3,p_L,p_R,s11_L,S12_L,s13_L,s11_R,S12_R,s13_R
    !!!---star variables :
    real(PR) :: frls(1:Nl)
    real(PR) :: rs,s11s,s12s,s13s,Es,ps
    real(PR) ::rls(1:Nl),pls(1:Nl),elhs(1:Nl),als(1:Nl,1:3),bls(1:Nl,1:3),cls(1:Nl,1:3)
    integer :: i, iph


    !!!---CONVERSION CONSERVATIVE -> PRIMITIVE

    call conservative2primitive(Nl,U_L,G_L,fl_L,Yl_L,rl_L,pl_L,elh_L,al_L,bl_L,cl_L,r_L,vx_L,vy_L,vz_L,p_L,e_L,c_L)
    call MAJ_meca(Nl,mu,rho0,rl_L,fl_L,al_L,bl_L,cl_L,pl_L,sigl_L,sig_L,ele_L)
    call MAJ_mixt(Nl,rl_L,fl_L,Yl_L,pl_L,sigl_L,elh_L,ele_L,vx_L,vy_L,vz_L,r_L,p_L,sig_L,e_L,c_L)

    call conservative2primitive(Nl,U_R,G_R,fl_R,Yl_R,rl_R,pl_R,elh_R,al_R,bl_R,cl_R,r_R,vx_R,vy_R,vz_R,p_R,e_R,c_R)
    call MAJ_meca(Nl,mu,rho0,rl_R,fl_R,al_R,bl_R,cl_R,pl_R,sigl_R,sig_R,ele_R)
    call MAJ_mixt(Nl,rl_R,fl_R,Yl_R,pl_R,sigl_R,elh_R,ele_R,vx_R,vy_R,vz_R,r_R,p_R,sig_R,e_R,c_R)

    s11_L=sig_L(1,1) ; s11_R=sig_R(1,1)
    s12_L=sig_L(1,2) ; s12_R=sig_R(1,2)
    s13_L=sig_L(1,3) ; s13_R=sig_R(1,3)

    !!!---FLUX
    S_R=max(vx_L+c_L,vx_R+c_R)
    S_L=min(vx_L-c_L,vx_R-c_R)

    IF(S_L.ge.0.0_PR)THEN

       us=vx_L ; vs=vy_L ; ws=vz_L

       F1(1:Nl)=fl_L(1:Nl)*rl_L(1:Nl)*vx_L
       F1(Nl+1)=r_L*vx_L**2-s11_L
       F1(Nl+2)=r_L*vx_L*vy_L-s12_L
       F1(Nl+3)=r_L*vx_L*vz_L-s13_L
       F1(Nl+4)=(r_L*E_L-s11_L)*vx_L-s12_L*vy_L-s13_L*vz_L
 
       do i=1,Nl+4
       if(F1(i).ne.F1(i))then
         print*, 'F1 1:', i          
         stop
       endif
       enddo
     
       F2(1:Nl,1)=vx_L*fl_L(1:Nl)
       do i=1,3
       F2(1:Nl,1+i)=vx_L*al_L(1:Nl,i)
       enddo
       do i=1,3
       F2(1:Nl,4+i)=vx_L*bl_L(1:Nl,i)
       enddo
       do i=1,3
       F2(1:Nl,7+i)=vx_L*cl_L(1:Nl,i)
       enddo
       F2(1:Nl,11)=fl_L(1:Nl)*rl_L(1:Nl)*vx_L*elh_L(1:Nl)

    ELSEIF(S_R.le.0.0_PR)THEN

       us=vx_R ; vs=vy_R ; ws=vz_R

       F1(1:Nl)=fl_R(1:Nl)*rl_R(1:Nl)*vx_R
       F1(Nl+1)=r_R*vx_R**2-s11_R
       F1(Nl+2)=r_R*vx_R*vy_R-s12_R
       F1(Nl+3)=r_R*vx_R*vz_R-s13_R
       F1(Nl+4)=(r_R*E_R-s11_R)*vx_R-s12_R*vy_R-s13_R*vz_R
       do i=1,Nl+4
       if(F1(i).ne.F1(i))then
         print*, 'F1 2:', i
         stop          
       endif
       enddo

       F2(1:Nl,1)=vx_R*fl_R(1:Nl)
       do i=1,3
       F2(1:Nl,1+i)=vx_R*al_R(1:Nl,i)
       enddo
       do i=1,3
       F2(1:Nl,4+i)=vx_R*bl_R(1:Nl,i)
       enddo
       do i=1,3
       F2(1:Nl,7+i)=vx_R*cl_R(1:Nl,i)
       enddo
       F2(1:Nl,11)=fl_R(1:Nl)*rl_R(1:Nl)*vx_R*elh_R(1:Nl)

    ELSE

       us=((r_L*vx_L**2-s11_L)-(r_R*vx_R**2-s11_R)-S_L*r_L*vx_L+S_R*r_R*vx_R) / &
          (r_L*vx_L-r_R*vx_R-S_L*r_L+S_R*r_R)


       IF(S_L.le.0.0_PR.and.us.ge.0.0_PR)THEN
        
          frls(1:Nl)=fl_L(1:Nl)*rl_L(1:Nl)*(S_L-vx_L)/(S_L-us)

          rs=sum(frls(1:Nl))

          C0=1.0_PR/( (vx_R-S_R)*r_R-(vx_L-S_L)*r_L )
          C1=(vx_R-S_R)*r_R
          C2=(vx_L-S_L)*r_L
          C3=(vx_L-S_L)*r_L*(vx_R-S_R)*r_R    
          s11s=C0*( C1*s11_L-C2*s11_R+C3*(vx_R-vx_L) )
          s12s=C0*( C1*s12_L-C2*s12_R+C3*(vy_R-vy_L) )
          s13s=C0*( C1*s13_L-C2*s13_R+C3*(vz_R-vz_L) )

          C0=1.0_PR/((vx_L-S_L)*r_L)
          vs=vy_L+(s12s-s12_L)*C0           
          ws=vz_L+(s13s-s13_L)*C0     
      
          Es=(r_L*E_L*(vx_L-S_L)-s11_L*vx_L-s12_L*vy_L-s13_L*vz_L+s11s*us+s12s*vs+s13s*ws)/&
           (rs*(us-S_L))
  
          F1(1:Nl)=frls(1:Nl)*us
          F1(Nl+1)=rs*us**2-s11s
          F1(Nl+2)=rs*us*vs-s12s
          F1(Nl+3)=rs*us*ws-s13s
          F1(Nl+4)=(rs*Es-s11s)*us-s12s*vs-s13s*ws
          do i=1,Nl+4
          if(F1(i).ne.F1(i))then
            print*, 'F1 3:', i
            stop          
          endif
          enddo

          !!!---F2:
          do i=1,3
          als(1:Nl,i)=(al_L(1:Nl,i)*(vx_L-S_L)+bl_L(1:Nl,i)*(vy_L-vs)+cl_L(1:Nl,i)*(vz_L-ws))/&
                     (us-S_L)
          enddo

          bls(1:Nl,1:3)=bl_L(1:Nl,1:3)
          cls(1:Nl,1:3)=cl_L(1:Nl,1:3)
 
          !!! Saurel 2009 p1689
          !!!---étant donné la densité initiale rl_L, la densité finale rls et la
          !pression initiale pl_L, la pression finale pls est imposée par la 
          !condition de Hugoniot

          rls(1:Nl)=rl_L(1:Nl)*(vx_L-S_L)/(us-S_L)

          do iph=1,Nl
             pls(iph)=hugoniot_pressure(iph,pl_L(iph),rl_L(iph),rls(iph))   
             elhs(iph)=energie_hydro(iph,rls(iph),pls(iph))
          enddo

          F2(1:Nl,1)=us*fl_L(1:Nl)
          do i=1,3
          F2(1:Nl,1+i)=us*als(1:Nl,i)
          enddo
          do i=1,3
          F2(1:Nl,4+i)=us*bls(1:Nl,i)
          enddo
          do i=1,3
          F2(1:Nl,7+i)=us*cls(1:Nl,i)
          enddo
          !!!---(alpha.rho)_l* ou "alpha_l*.rho_l*"
          !F2(1:Nl,11)=frls(1:Nl)*us*elhs(1:Nl)
           F2(1:Nl,11)=fl_L(1:Nl)*rls(1:Nl)*us*elhs(1:Nl)

       ELSEIF(S_R.ge.0.0_PR.and.us.le.0.0_PR)THEN

          frls(1:Nl)=fl_R(1:Nl)*rl_R(1:Nl)*(S_R-vx_R)/(S_R-us)

          rs=sum(frls(1:Nl))

          C0=1.0_PR/( (vx_R-S_R)*r_R-(vx_L-S_L)*r_L )
          C1=(vx_R-S_R)*r_R
          C2=(vx_L-S_L)*r_L
          C3=(vx_L-S_L)*r_L*(vx_R-S_R)*r_R    
          s11s=C0*( C1*s11_L-C2*s11_R+C3*(vx_R-vx_L) )
          s12s=C0*( C1*s12_L-C2*s12_R+C3*(vy_R-vy_L) )
          s13s=C0*( C1*s13_L-C2*s13_R+C3*(vz_R-vz_L) )

          C0=1.0_PR/((vx_R-S_R)*r_R)
          vs=vy_R+(s12s-s12_R)*C0           
          ws=vz_R+(s13s-s13_R)*C0     
      
          Es=(r_R*E_R*(vx_R-S_R)-s11_R*vx_R-s12_R*vy_R-s13_R*vz_R+s11s*us+s12s*vs+s13s*ws)/&
           (rs*(us-S_R))
  
          F1(1:Nl)=frls(1:Nl)*us
          F1(Nl+1)=rs*us**2-s11s
          F1(Nl+2)=rs*us*vs-s12s
          F1(Nl+3)=rs*us*ws-s13s
          F1(Nl+4)=(rs*Es-s11s)*us-s12s*vs-s13s*ws
       do i=1,Nl+4
       if(F1(i).ne.F1(i))then
         print*, 'F1 4:', i, Es, rs, us, vs, ws, s11s, s12s, s13s
         print*, ''
         print*, 'F1 4:', C0, vx_R-S_R
         stop          
       endif
       enddo

          !!!---F2:
          do i=1,3
          als(1:Nl,i)=(al_R(1:Nl,i)*(vx_R-S_R)+bl_R(1:Nl,i)*(vy_R-vs)+cl_R(1:Nl,i)*(vz_R-ws))/&
                     (us-S_R)
          enddo
          bls(1:Nl,1:3)=bl_R(1:Nl,1:3)
          cls(1:Nl,1:3)=cl_R(1:Nl,1:3)


          !!! Saurel 2009 p1689----------------------
          !!!---étant donné la densité initiale rl_R, la densité finale rls et la
          !pression initiale pl_R, la pression finale pls est imposée par la 
          !condition de Hugoniot

          rls(1:Nl)=rl_R(1:Nl)*(vx_R-S_R)/(us-S_R)
         
          do iph=1,Nl
             pls(iph)=hugoniot_pressure(iph,pl_R(iph),rl_R(iph),rls(iph))   
             elhs(iph)=energie_hydro(iph,rls(iph),pls(iph))
          enddo

          !!!----------------------------------------

          F2(1:Nl,1)=us*fl_R(1:Nl)
          do i=1,3
          F2(1:Nl,1+i)=us*als(1:Nl,i)
          enddo
          do i=1,3
          F2(1:Nl,4+i)=us*bls(1:Nl,i)
          enddo
          do i=1,3
          F2(1:Nl,7+i)=us*cls(1:Nl,i)
          enddo
          !!!---(alpha.rho)_l* ou "alpha_l*.rho_l*"
          !F2(1:Nl,11)=frls(1:Nl)*us*elhs(1:Nl)
          F2(1:Nl,11)=fl_R(1:Nl)*rls(1:Nl)*us*elhs(1:Nl)

       ENDIF

    ENDIF

end subroutine HLLC_flux

!!!========================== SUBROUTINE DE CONVERSION =========================
!!!---1D

subroutine primitive2conservative(Nl,U,G,fl,Yl,rhl,pl,elh,al,bl,cl,rh,vx,vy,vz,p,e,c)

   implicit none

   integer, intent(in) :: Nl
   real(PR), intent(out) :: U(1:Nl+4), G(1:Nl,1:11)
   real(PR), intent(in) :: fl(1:Nl), Yl(1:Nl)
   real(PR), intent(in) :: rhl(1:Nl), pl(1:Nl)
   real(PR), intent(in) :: elh(1:Nl)
   real(PR), intent(in) :: al(1:Nl,1:3), bl(1:Nl,1:3), cl(1:Nl,1:3)
   real(PR), intent(in) :: rh
   real(PR), intent(in) :: vx, vy, vz
   real(PR), intent(in) :: p, e, c
   real(PR) :: c1

   G(1:Nl,1)=fl(1:Nl)
   G(1:Nl,2:4) =al(1:Nl,1:3)
   G(1:Nl,5:7) =bl(1:Nl,1:3)
   G(1:Nl,8:10)=cl(1:Nl,1:3)
   G(1:Nl,11)=fl(1:Nl)*rhl(1:Nl)*elh(1:Nl)
   U(1:Nl)=fl(1:Nl)*rhl(1:Nl)
   U(Nl+1)=rh*vx
   U(Nl+2)=rh*vy
   U(Nl+3)=rh*vz
   U(Nl+4)=rh*e

 
end subroutine primitive2conservative

subroutine conservative2primitive(Nl,U,G,fl,Yl,rhl,pl,elh,al,bl,cl,rh,vx,vy,vz,p,e,c)

   implicit none

   integer, intent(in) :: Nl
   real(PR), intent(in) :: U(1:Nl+4), G(1:Nl,1:11)
   real(PR), intent(out) :: fl(1:Nl), Yl(1:Nl)
   real(PR), intent(out) :: rhl(1:Nl), pl(1:Nl)
   real(PR), intent(out) :: elh(1:Nl)
   real(PR), intent(out) :: al(1:Nl,1:3), bl(1:Nl,1:3), cl(1:Nl,1:3)
   real(PR), intent(out) :: rh
   real(PR), intent(out) :: vx, vy, vz
   real(PR), intent(out) :: p, e, c
   real(PR) :: c1
   integer ::  iph, i,j,k,i2

     !!!---variables phasiques
     fl(1:Nl)=G(1:Nl,1)
     al(1:Nl,1:3)=G(1:Nl,2:4)  
     bl(1:Nl,1:3)=G(1:Nl,5:7)  
     cl(1:Nl,1:3)=G(1:Nl,8:10)
     do iph=1,Nl
      pl(iph)=(g_sge(iph)-1.0_PR)*G(iph,11)/fl(iph)-g_sge(iph)*p_sge(iph)
      !if(pl(iph).lt.0.0_PR)then
      !   print*, 'pl<0:', fl(:)
      !   stop
      !endif
     enddo

     rhl(1:Nl)=U(1:Nl)/fl(1:Nl)
     rh=sum(U(1:Nl))
     Yl(1:Nl)=fl(1:Nl)*rhl(1:Nl)/rh
     elh(1:Nl)=G(1:Nl,11)/(fl(1:Nl)*rhl(1:Nl))

     !if(rh.lt.0.0_PR)then
     !   print*, 'rh<0', U(1:Nl), fl(1:Nl)
     !   stop
     !endif

     !do iph=1,Nl
     !   pl(iph)=pressure(iph,rhl(iph),elh(iph))
     !enddo

     !!!---variables de mélange
     vx=U(Nl+1)/rh
     vy=U(Nl+2)/rh
     vz=U(Nl+3)/rh
     e=U(Nl+4)/rh

     !!!---p ne dépend que de G may be <0 !
     p=sum(fl(1:Nl)*pl(1:Nl))

     !print*, 'p2=', p

     !if(p.lt.0.0_PR)then
     !   print*, 'p2=', p
     !   stop
     !endif

     if(compute_c)then
      c=soundspeed_mixt(Nl,fl(1:Nl),Yl(1:Nl),rhl(1:Nl),pl(1:Nl),rh)
     else 
      c=0.0_PR
     endif
 
end subroutine conservative2primitive


!!!============================ SUBROUTINE EOS =================================
!!!---0D

subroutine MAJ_mixt(Nl,rhl,fl,Yl,pl,sigl,elh,ele,vx,vy,vz,rh,p,sig,e,c)

   implicit none
   integer, intent(in) :: Nl
   real(PR), intent(in) :: rhl(1:Nl),fl(1:Nl),Yl(1:Nl),pl(1:Nl)
   real(PR), intent(in) :: sigl(1:Nl,1:3,1:3), elh(1:Nl), ele(1:Nl)
   real(PR), intent(in) :: vx,vy,vz
   real(PR), intent(out) :: rh, p, sig(1:3,1:3), e, c
   integer :: i,j,k,iph

   !!!write(iout,*) '   calcul des variables de mélange : rh, Yl, E, Sij' 

   rh=sum(fl(1:Nl)*rhl(1:Nl)) 
   e=sum(Yl(1:Nl)*(elh(1:Nl)+ele(1:Nl)))+0.5_PR*(vx**2+vy**2+vz**2)
   p=sum( fl(1:Nl)*pl(1:Nl) )
   forall(i=1:3,j=1:3) sig(i,j)=sum(fl(1:Nl)*sigl(1:Nl,i,j))
   c=soundspeed_mixt(Nl,fl(1:Nl),Yl(1:Nl),rhl(1:Nl),pl(1:Nl),rh)

end subroutine MAJ_mixt

real(PR) function energie_hydro(iph,rhl,pl)

   implicit none
   !!! p=rho(Gam-1)e-GammP0  !!! Pa
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! rh=(p+GammP0)/(e*(Gam-1)) !!! kg/m^3
   !!! c^2=Gam(p+p0)/rho !!! m/s
   !!! eau: Gam=6.1 p0=2e9 ~ molecular attraction
   integer, intent(in) :: iph
   real(PR), intent(in) :: rhl, pl

   energie_hydro=( pl + g_sge(iph)*p_sge(iph) )/( (g_sge(iph) - 1.0_PR)*rhl )

end function energie_hydro

real(PR) function ddnu_energie_hydro(iph,rh,p)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! de/drho=-(p+Gam*P0)/( (Gam-1)*rho^2 )
   integer, intent(in) :: iph
   real(PR), intent(in) :: rh, p

   !ddrh_energie_hydro=-( p + g_sge(iph)*p_sge(iph) )/( (g_sge(iph) - 1.0_PR)*rh**2 )

   ddnu_energie_hydro= ( p + g_sge(iph)*p_sge(iph) )/( g_sge(iph) - 1.0_PR )

end function ddnu_energie_hydro

real(PR) function ddp_energie_hydro(iph,rh,p)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! de/dp=1.0/( (Gam-1)*rho )
   integer, intent(in) :: iph
   real(PR), intent(in) :: rh, p

   ddp_energie_hydro=1.0_PR/( (g_sge(iph) - 1.0_PR)*rh )

end function ddp_energie_hydro

real(PR) function density(iph,e,p)

   implicit none
   !!! stiffened gas equation
   !!! p=rho(Gam-1)e-GammP0
   !!! e=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(e*(Gam-1))
   !!! c^2=Gam(p+p0)/rho
   integer, intent(in) :: iph
   real(PR), intent(in) :: e, p

   density=( p + g_sge(iph)*p_sge(iph) )/( (g_sge(iph) - 1.0_PR)*e )

end function density

real(PR) function pressure(iph,rh,e)

   implicit none
   !!! stiffened gas equation
   !!! p=rho(Gam-1)e-GammP0
   !!! e=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(e*(Gam-1))
   !!! c^2=Gam(p+p0)/rho
   !real(PR) :: pressure
   integer, intent(in) :: iph
   real(PR), intent(in) :: e, rh

   !pressure=max(rh*(g_sge(iph) - 1.0_PR)*e - g_sge(iph)*p_sge(iph), 0.0_PR)

   pressure=rh*(g_sge(iph) - 1.0_PR)*e - g_sge(iph)*p_sge(iph)


   if(pressure.lt.0.0_PR)then
      
      print*, 'p<0:'
      print*, 'iph, rh, rh*(g_sge(iph)-1.0_PR)*e, e, g_sge(iph)*p_sge(iph)'
      print*, iph, rh, rh*(g_sge(iph)-1.0_PR)*e, e, g_sge(iph)*p_sge(iph)
      stop
   endif

end function pressure

real(PR) function soundspeed(iph,rh,p)

   implicit none
   !!! stiffened gas equation
   !!! p=rho(Gam-1)e-GammP0
   !!! e=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(e*(Gam-1))
   !!! c^2=Gam*(p+p0)/rho
   integer, intent(in) :: iph
   real(PR), intent(in) :: p, rh

   !soundspeed=sqrt(g_sge(iph)*(p+p_sge(iph))/rh)
   if((g_sge(iph)*(p+p_sge(iph)))/rh.lt.0.0_PR)then
    ierr=1
    print*, 'PROBLEM soundspeed:', iph, p, rh
    stop
   endif

   soundspeed=sqrt(g_sge(iph)*(p+p_sge(iph))/rh)

end function soundspeed

real(PR) function soundspeed_mixt(Nl,fl,Yl,rhl,pl,rh)

   implicit none 
   !real(PR) :: soundspeed_mixt
   integer, intent(in) :: Nl
   real(PR), intent(in) :: rh
   real(PR), dimension(1:Nl), intent(in) :: fl, Yl, rhl, pl
   real(PR) :: c1
   integer :: iph

   !!!---frozen speed of sound:
   c1=0.0_PR
   do iph=1,Nl
     c1=c1+Yl(iph)*soundspeed(iph,rhl(iph),pl(iph))**2
   enddo

   soundspeed_mixt=sqrt(c1)

   !!!---Wood speed of sound:
   !c1=0.0_PR
   !do iph=1,Nl
   !  c1=c1+fl(iph)/(rhl(iph)*soundspeed(iph,rhl(iph),pl(iph))**2)
   !enddo
   !c1=c1*rh
   !soundspeed_mixt=sqrt(1.0_PR/c1)

   !c1=0.0_PR
   !do iph=1,Nl
   !  c1=max(c1,soundspeed(iph,rhl(iph,i,j,k),pl(iph,i,j,k)))
   !enddo
   !c(i,j,k)=c1

end function soundspeed_mixt

real(PR) function pressure_mixt(Nl,rh,e,fl)

   !==================================
   ! pression d'un mélange de fluides suivant
   ! une EOS de type SGE (p1692 de Saurel 2009)
   !=================================
   implicit none 
   integer, intent(in) :: Nl
   real(PR), intent(in) :: rh, e
   real(PR), dimension(1:Nl), intent(in) :: fl
   real(PR) :: C1,denom
   integer :: iph

   pressure_mixt=rh*e
   denom=0.0_PR

   do iph=1,Nl
     C1=fl(iph)/(g_sge(iph)-1.0_PR)
     pressure_mixt=pressure_mixt-C1*g_sge(iph)*p_sge(iph)
     denom=denom+C1
   enddo

   pressure_mixt=pressure_mixt/denom

end function pressure_mixt

real(PR) function Hugoniot_pressure(iph,p0,r0,rf)

   implicit none
   !!! stiffened gas equation :
   !!! p=rho(Gam-1)e-GammP0
   !!! e=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(e*(Gam-1))
   !!! c^2=Gam(p+p0)/rho
   integer, intent(in) :: iph
   real(PR), intent(in) :: p0, r0, rf
   real(PR) :: Gm1, Gp1, pinf 


   Gm1=(g_sge(iph)-1.0_PR)
   Gp1=(g_sge(iph)+1.0_PR)
   pinf=p_sge(iph)

   Hugoniot_pressure=(p0+pinf)*(Gm1*r0-Gp1*rf)/(Gm1*rf-Gp1*r0)-pinf 

end function Hugoniot_pressure

!!!========================== SUBROUTINES DE MECANIQUE =========================
!!!---0D

subroutine MAJ_meca(Nl,mu,rho0,rhl,fl,al,bl,cl,pl,sigl,sig,ele)

   implicit none
   integer, intent(in) :: Nl
   real(PR), intent(in) :: mu(1:Nl)
   real(PR), intent(in) :: rho0(1:Nl)
   real(PR), intent(in) :: rhl(1:Nl)
   real(PR), intent(in) :: fl(1:Nl)
   real(PR), intent(in) :: al(1:Nl,1:3)
   real(PR), intent(in) :: bl(1:Nl,1:3)
   real(PR), intent(in) :: cl(1:Nl,1:3)
   real(PR), intent(in) :: pl(1:Nl)
   real(PR), intent(out) :: sigl(1:Nl,1:3,1:3)
   real(PR), intent(out) :: sig(1:3,1:3)
   real(PR), intent(out) :: ele(1:Nl)
   integer :: i,j,k,iph


   sigl=0.0_PR

   !!!---partie deviatotique
   !do iph=1,Nl
   ! sigl(iph,1:3,1:3)=& !!! perfo: mettre iph après 1:3;1:3 ?
   !      Deviateur(mu(iph),rhl(iph),rho0(iph),al(iph,1:3),bl(iph,1:3),cl(iph,1:3))
   !enddo

   !!!---partie hydrostatique
   forall(iph=1:Nl,i=1:3) sigl(iph,i,i)=sigl(iph,i,i)-pl(iph)

   !!!---tenseur des contraintes totales
   forall(i=1:3,j=1:3) sig(i,j)=sum(fl(1:Nl)*sigl(1:Nl,i,j))

   !!!---Energie élastique
   ele=0.0_PR
   !do iph=1,Nl
   ! ele(iph)=energie_elastique(mu(iph),rho0(iph),al(iph,1:3),bl(iph,1:3),cl(iph,1:3))
   !enddo

end subroutine MAJ_meca

function energie_elastique(mu,rho0,a,b,c)

   implicit none
   real(PR) :: energie_elastique
   real(PR), intent(in) :: mu, rho0
   real(PR), intent(in) :: a(1:3), b(1:3), c(1:3)
   real(PR) :: e(1:3,1:3), G(1:3,1:3), g2(1:3,1:3)
   real(PR) :: detG, C1, j2l
   integer :: i,j,k

   !!!---vecteurs de base :
   e(1,1:3)=(/a(1),b(1),c(1)/)
   e(2,1:3)=(/a(2),b(2),c(2)/)
   e(3,1:3)=(/a(3),b(3),c(3)/)

   !!!---tenseur de Finger :
   do j=1,3
    do i=1,3
      G(i,j)=DOT_PRODUCT(e(i,1:3),e(j,1:3))
    enddo
   enddo

   !!!---Division par le déterminant :
   detG=G(1,1)*G(2,2)*G(3,3)+G(2,1)*G(3,2)*G(1,3)+G(3,1)*G(1,2)*G(2,3)&
       -G(3,1)*G(2,2)*G(1,3)-G(1,1)*G(3,2)*G(2,3)-G(2,1)*G(1,2)*G(3,3)

   if(detG.ne.0)then

      C1=detG**(-1.0_PR/3.0_PR)
      
      g=G*C1

      !!!---tenseur g**2 : g2
      do j=1,3
        do i=1,3
           g2(i,j)=sum(G(i,1:3)*G(1:3,j))   
        enddo
      enddo

      j2l=g2(1,1)+g2(2,2)+g2(3,3)
     
      energie_elastique=mu/(8.0_PR*rho0)*(j2l-3.0_PR)

   else

      energie_elastique=0.0_PR
    
   endif

  ! if(energie_elastique.ne.energie_elastique) stop
   
end function energie_elastique

function Deviateur(mu,rho,rho0,a,b,c)

   implicit none
   real(PR), dimension(1:3,1:3) :: Deviateur
   real(PR), intent(in) :: mu, rho, rho0
   real(PR), intent(in) :: a(1:3), b(1:3), c(1:3)
   real(PR) :: e(1:3,1:3), G(1:3,1:3), g2(1:3,1:3)
   real(PR) :: detG, C1, C2 
   integer :: i,j,k

   IF(mu.gt.0.0_PR)THEN

      !!!---vecteurs de base :
      e(1,1:3)=(/a(1),b(1),c(1)/)
      e(2,1:3)=(/a(2),b(2),c(2)/)
      e(3,1:3)=(/a(3),b(3),c(3)/)

      !!!---tenseur de Finger :
      do j=1,3
       do i=1,3
         G(i,j)=DOT_PRODUCT(e(i,1:3),e(j,1:3))
       enddo
      enddo

      !!!---Division par le déterminant :
      detG=G(1,1)*G(2,2)*G(3,3)+G(2,1)*G(3,2)*G(1,3)+G(3,1)*G(1,2)*G(2,3)&
          -G(3,1)*G(2,2)*G(1,3)-G(1,1)*G(3,2)*G(2,3)-G(2,1)*G(1,2)*G(3,3)

      if(detG.ne.0.0_PR)then

         C1=detG**(-1.0_PR/3.0_PR)
         
         g=G*C1

         !!!---tenseur g**2 : g2
         do j=1,3
           do i=1,3
              g2(i,j)=sum(G(i,1:3)*G(1:3,j))   
           enddo
         enddo

         C1=-0.5_PR*rho/rho0*mu
         C2=C1*1.0_PR/3.0_PR*(g2(1,1)+g2(2,2)+g2(3,3))
         
         Deviateur(1:3,1:3)=C1*g2(1:3,1:3)

         do i=1,3
           Deviateur(i,i)=Deviateur(i,i)-C2
         enddo
        
      else
 
         Deviateur=0.0_PR

      endif

   ELSE

      Deviateur=0.0_PR

   ENDIF !!! mu>0

end function Deviateur




!!!================= RELAXTATION =========================


subroutine relax(Nl,N,U,G)

   implicit none
   integer, intent(in) :: Nl,N
   real(PR), intent(inout) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)

   !!!---champs 1D locaux
   real(PR), dimension(Nl,N) :: fl, Yl, pl, rhl, elh, ele
   real(PR), dimension(Nl,3,N) :: al,bl,cl
   real(PR), dimension(Nl,3,3,N) :: sigl
   real(PR), dimension(3,3,N) :: sig
   real(PR), dimension(N) :: rh, vx, vy,vz, z, c, e, p
       
   integer :: i,j,k

   call cpu_time(t30)

   compute_c=.false.
   do i=1,N
   call conservative2primitive(Nl,U(:,i),G(:,:,i),fl(:,i),Yl(:,i),rhl(:,i),pl(:,i),elh(:,i),&
                               al(:,:,i),bl(:,:,i),cl(:,:,i),rh(i),vx(i),vy(i),vz(i),p(i),e(i),c(i))
   enddo
   compute_c=.true.

   !print*, 'minval p=', minval(p(:)), minloc(p,1) 

   if(ierr.ne.0) call crash('relaxp1')

   do i=1,N

       !if(i.eq.i_sortie)then
       !  verb=.true.
       !else
       !  verb=.false.
       !endif

       call relax_p(Nl,fl(1:Nl,i),Yl(1:Nl,i),rhl(1:Nl,i),&
                    elh(1:Nl,i),pl(1:Nl,i),rh(i),p(i),&
                    vx(i),vy(i),vz(i),c(i),e(i))

    enddo

   !print*, 'minval p2=', minval(p(:)), minloc(p,1) 

   if(ierr.ne.0) call crash('relaxp2')

   do i=1,N

      call MAJ_meca(Nl,mu(1:Nl),rho0(1:Nl),rhl(1:Nl,i),fl(1:Nl,i),&
                       al(1:Nl,1:3,i),bl(1:Nl,1:3,i),cl(1:Nl,1:3,i),&
                       pl(1:Nl,i),&
                       sigl(1:Nl,1:3,1:3,i),sig(1:3,1:3,i),ele(1:Nl,i))!!!S

      call MAJ_mixt(Nl,rhl(1:Nl,i),fl(1:Nl,i),Yl(1:Nl,i),pl(1:Nl,i),&
                      sigl(1:Nl,1:3,1:3,i),elh(1:Nl,i),ele(1:Nl,i),&
                      vx(i), vy(i), vz(i),&
                      rh(i),p(i),sig(1:3,1:3,i),e(i),c(i))
   enddo

   do i=1,N 
    call primitive2conservative(Nl,U(:,i),G(:,:,i),fl(:,i),Yl(:,i),rhl(:,i),pl(:,i),elh(:,i),&
                                al(:,:,i),bl(:,:,i),cl(:,:,i),rh(i),vx(i),vy(i),vz(i),p(i),e(i),c(i))
   enddo
   if(ierr.ne.0) call crash('relaxp3')

   call cpu_time(t3)

   Trelax=Trelax+t3-t30
 
end subroutine relax


!!!================ SUBROUTINE DE RELAXATION DE LA PRESSION ====================

subroutine init_relaxp(Nl)

     implicit none
     integer, intent(in) :: Nl

     IF(Nl.eq.1)THEN

        relax_p => relax_p0
        write(iout,*) '   relax_p => relax_p0' 

     ELSE

     if(allocated(xin))    deallocate(xin)  
     if(allocated(bin))    deallocate(bin)    
     if(allocated(flrhl0)) deallocate(flrhl0)
     if(allocated(el0))    deallocate(el0)
     if(allocated(nul0))   deallocate(nul0)
     if(allocated(pl0))    deallocate(pl0)
     if(allocated(Zl))     deallocate(Zl)
     if(allocated(fl0_gl))  deallocate(fl0_gl)
     if(allocated(fl0))    deallocate(fl0)
     if(allocated(coef1))    deallocate(coef1)

     allocate(flrhl0(1:Nl)); flrhl0(:)=0.0_PR
     allocate(el0(1:Nl))   ; el0(:)=0.0_PR
     allocate(nul0(1:Nl))  ; nul0(:)=0.0_PR
     allocate(pl0(1:Nl))   ; pl0(:)=0.0_PR
     allocate(Zl(1:Nl))    ; Zl(:)=0.0_PR
     allocate(fl0_gl(1:Nl)); fl0_gl(:)=0.0_PR
     allocate(fl0(1:Nl))   ; fl0(:)=0.0_PR
     allocate(coef1(1:Nl))   ; coef1(:)=0.0_PR


     if(mode_relax.eq.1)then
        !!!---valable pour SGE : 1 equation sur p
        Nsnes=1
        allocate(xin(1:Nsnes))   ; xin(:)=0.0_PR
        allocate(bin(1:Nsnes))   ; bin(:)=0.0_PR 

        if(mode_pI.eq.1)then
           relax_p => relax_p_sge_1
          ! call init_SNES(Nsnes,xin,bin,F_function=F_relaxp_sge_1,J_function=J_relaxp_sge_1)
           write(iout,*) '   relax_p => relax_p_sge_1' 
           write(iout,*) '   F_function => F_relax_p_sge_1' 
        elseif(mode_pI.eq.2)then
           relax_p => relax_p_sge_2
          ! call init_SNES(Nsnes,xin,bin,F_function=F_relaxp_sge_2,J_function=J_relaxp_sge_2)
           write(iout,*) '   relax_p => relax_p_sge_2' 
           write(iout,*) '   F_function => F_relax_p_sge_2' 
        endif

     elseif(mode_relax.eq.2)then

        !!!---equations : consevation de e+p/rho par phase + sum fl=1
        !!!---variables : nuk=1/rhk et p
        Nsnes=Nl+1
        allocate(xin(1:Nsnes)) ; xin(:)=0.0_PR
        allocate(bin(1:Nsnes)) ; bin(:)=0.0_PR 

        if(mode_pI.eq.1)then
           relax_p => relax_p_eos
           !call init_SNES(Nsnes,xin,bin,F_function=F_relaxp_1,J_function=J_relaxp_1)
           write(iout,*) '   relax_p => relax_p_eos' 
           write(iout,*) '   F_function => F_relax_p_1' 
        elseif(mode_pI.eq.2)then
           relax_p => relax_p_eos
           !call init_SNES(Nsnes,xin,bin,F_function=F_relaxp_2,J_function=J_relaxp_2)
           write(iout,*) '   relax_p => relax_p_eos' 
           write(iout,*) '   F_function => F_relax_p_2' 
        endif

     endif

     ENDIF

end subroutine init_relaxp

subroutine relax_p0(Nl,fl,Yl,rhl,elh,pl,rh,p,vx,vy,vz,c,e)

     implicit none
     integer, intent(in) :: Nl
     real(PR), intent(inout) :: fl(1:Nl)
     real(PR), intent(inout) :: Yl(1:Nl)
     real(PR), intent(inout) :: rhl(1:Nl)
     real(PR), intent(inout) :: elh(1:Nl)
     real(PR), intent(inout) :: pl(1:Nl)
     real(PR), intent(inout) :: p, rh, c
     real(PR), intent(in) :: vx,vy,vz, e
     real(PR) :: pmin, pmax
     integer :: iph

     p=pressure_mixt(Nl,rh,e-0.5_PR*(vx**2+vy**2+vz**2),fl(1:Nl)) 

     do iph=1,Nl
        pl(iph)=p
     enddo

     do iph=1,Nl
       elh(iph)=energie_hydro(iph,rhl(iph),p)
     enddo

     !print*, 'coucouz1' 
     !print*, fl(1:), Yl(1:), rhl(1:), pl(1:), rh, e-0.5_PR*(vx**2+vy**2+vz**2)
     c=soundspeed_mixt(Nl,fl(1:Nl),Yl(1:Nl),rhl(1:Nl),pl(1:Nl),rh)
     !print*, 'coucouz2' 

end subroutine relax_p0

subroutine relax_p_sge_1(Nl,fl,Yl,rhl,elh,pl,rh,p,vx,vy,vz,c,e)

     implicit none
     integer, intent(in) :: Nl
     real(PR), intent(inout) :: fl(1:Nl)
     real(PR), intent(inout) :: Yl(1:Nl)
     real(PR), intent(inout) :: rhl(1:Nl)
     real(PR), intent(inout) :: elh(1:Nl)
     real(PR), intent(inout) :: pl(1:Nl)
     real(PR), intent(inout) :: p, rh, c
     real(PR), intent(in) :: vx,vy,vz, e
     real(PR) :: pmin, pmax, ehtot
     !real(PR), allocatable :: xin(:), bin(:)
     integer :: iph

        !pmin=minval(pl(1:Nl))
        !pmax=maxval(pl(1:Nl))
        
        !IF(0.5_PR*abs(pmax-pmin)/(pmax+pmin).gt.1.0e-3_PR)THEN

           do iph=1,Nl
             flrhl0(iph)=fl(iph)*rhl(iph)
             fl0_gl(iph)=fl(iph)/g_sge(iph)
             pl0(iph)=pl(iph)  
             coef1(iph)=fl0_gl(iph)*(p_sge(iph)+pl0(iph))
             fl0(iph)=fl(iph)
             el0(iph)=elh(iph)
           enddo 
 
           xin(1)=maxval(pl(1:Nl))
           bin(:)=0.0_PR

           if(verb)then
              write(iout,*) '-----------------------------------'
              write(iout,*) '   x=', xin(1:Nsnes)
              write(iout,*) '   rhl=', rhl(1:Nl)
              write(iout,*) '   sumY=',sum(Yl(1:Nl)) 
              write(iout,*) '   sumYelh=',sum(Yl(1:Nl)*elh(1:Nl)) 
              write(iout,*) '   sumfl=',sum(fl(1:Nl)) 
              write(iout,*) '   elh=', elh(1:Nl)
              write(iout,*) '   p=', p
              write(iout,*) '   pl=', pl(1:Nl)
              write(iout,*) '   fl=', fl(1:Nl)
              write(iout,*) '   rho=', rh
           endif     

           !call solve_SNES(Nsnes,xin,bin)

           !call My_Newton(1,xin(1:1),minval(pl),maxval(pl),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)

           !call My_Newton(1,xin(1:1),0.0_PR,maxval(pl),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)

           call My_Newton(1,xin(1:1),1.0e-30_PR,maxval(pl),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)

           p=xin(1)

           !if(p.lt.0.0_PR)then
           !    verb=.true.
           !    print*, 'coucou1'
           !endif

           if(verb)then
              write(iout,*) '   e+p/rh=', elh(1:Nl)+p/rhl(1:Nl)
              write(iout,*) '-----------------------------------'     
              write(iout,*) '   p=', p, ':<'
           endif

           do iph=1,Nl
             fl(iph)=fl0_gl(iph)*(pl0(iph)+g_sge(iph)*p_sge(iph)+(g_sge(iph)-1.0_PR)*p)/(p_sge(iph)+p)
             rhl(iph)=flrhl0(iph)/fl(iph)
             Yl(iph)=fl(iph)*rhl(iph)/rh
           enddo

           !rh=sum(fl(1:Nl)*rhl(1:Nl))
 
           if(verb)then     
              do iph=1,Nl
               elh(iph)=energie_hydro(iph,rhl(iph),p)
               pl(iph)=pressure(iph,rhl(iph),elh(iph))
              enddo
              write(iout,*) '-----------------------------------'
              write(iout,*) '   x=', xin(1:Nsnes)
              write(iout,*) '   rhl=', rhl(1:Nl)
              write(iout,*) '   sumY=',sum(Yl(1:Nl)) 
              write(iout,*) '   sumfl=',sum(fl(1:Nl)) 
              write(iout,*) '   elh=', elh(1:Nl)
              write(iout,*) '   pl=', pl(1:Nl)
              write(iout,*) '   fl=', fl(1:Nl)
              write(iout,*) '   rho=', rh
              write(iout,*) '   e+p/rh=', elh(1:Nl)+p/rhl(1:Nl)
              write(iout,*) '-----------------------------------'
              write(iout,*) '---Set pressure from e---'
           endif    

        !ENDIF !!! pmin-pmax/(pmin+pmax)>1e-3

     ehtot=e-0.5_PR*(vx**2+vy**2+vz**2)

     p=pressure_mixt(Nl,rh,ehtot,fl(1:Nl))
 
     !if(p.lt.0.0_PR)then
     !    verb=.true.
     !    print*, 'coucou2'
     !endif

     do iph=1,Nl
        pl(iph)=p
     enddo

     do iph=1,Nl
       elh(iph)=energie_hydro(iph,rhl(iph),p)
     enddo

     if(abs(sum(Yl(1:Nl)*elh(1:Nl))-ehtot)/ehtot.gt.1.0e-6_PR)then
        write(iout,*) 'PROBLEM ehtot ne sum(Ylehl):'
        write(iout,*) sum(Yl(1:Nl)*elh(1:Nl)), ehtot
        stop
     endif

 
     c=soundspeed_mixt(Nl,fl(1:Nl),Yl(1:Nl),rhl(1:Nl),pl(1:Nl),rh)
 
     if(verb)then 
        write(iout,*) 'P=', p
        write(iout,*) 'elh1=', elh(1)
        write(iout,*) 'elh2=', elh(2)
     endif     

end subroutine relax_p_sge_1

subroutine relax_p_sge_2(Nl,fl,Yl,rhl,elh,pl,rh,p,vx,vy,vz,c,e)

     implicit none
     integer, intent(in) :: Nl
     real(PR), intent(inout) :: fl(1:Nl)
     real(PR), intent(inout) :: Yl(1:Nl)
     real(PR), intent(inout) :: rhl(1:Nl)
     real(PR), intent(inout) :: elh(1:Nl)
     real(PR), intent(inout) :: pl(1:Nl)
     real(PR), intent(inout) :: p, rh, c
     real(PR), intent(in) :: vx,vy,vz, e
     real(PR) :: pmin, pmax, ehtot
     !real(PR), allocatable :: xin(:), bin(:)
     integer :: iph

        pmin=minval(pl(1:Nl))
        pmax=maxval(pl(1:Nl))
        
        !IF(0.5_PR*abs(pmax-pmin)/(pmax+pmin).gt.1.0e-3_PR)THEN

           do iph=1,Nl
             flrhl0(iph)=fl(iph)*rhl(iph)
             fl0(iph)=fl(iph)
             pl0(iph)=pl(iph)  
           enddo 
 
           do iph=1,Nl
            if(pl(iph).gt.0.0_PR)then
            Zl(iph)=rhl(iph)*soundspeed(iph,rhl(iph),pl(iph))
            else
            Zl(iph)=0.0_PR
            endif
           enddo
           pI0=sum(Zl(1:Nl)*pl(1:Nl))/sum(Zl(1:Nl))
           !print*, 'pI0=', pI0

           xin(1)=maxval(pl(1:Nl))
           bin(:)=0.0_PR

           if(verb)then
              write(iout,*) 'pI0=', pI0
              write(iout,*) '-----------------------------------'
              write(iout,*) '   x=', xin(1:Nsnes)
              write(iout,*) '   rhl=', rhl(1:Nl)
              write(iout,*) '   sumY=',sum(Yl(1:Nl)) 
              write(iout,*) '   sumfl=',sum(fl(1:Nl)) 
              write(iout,*) '   elh=', elh(1:Nl)
              write(iout,*) '   p=', p
              write(iout,*) '   pl=', pl(1:Nl)
              write(iout,*) '   fl=', fl(1:Nl)
              write(iout,*) '   rho=', rh
           endif     

          ! call solve_SNES(Nsnes,xin,bin)
           call My_Newton(1,xin(1:1),1.0e-10_PR,maxval(pl),F=F_relaxp_sge_2,DF=J_relaxp_sge_2)

           p=xin(1)

           !if(p.lt.0.0_PR)then
           !    verb=.true.
           !    print*, 'coucou1'
           !endif

           if(verb)then
              write(iout,*) '   e+p/rh=', elh(1:Nl)+p/rhl(1:Nl)
              write(iout,*) '-----------------------------------'     
              write(iout,*) '   p=', p
           endif

           do iph=1,Nl
              fl(iph)=fl0(iph)*(pl0(iph)+g_sge(iph)*p_sge(iph)+(g_sge(iph)-1.0_PR)*pI0)/(p+g_sge(iph)*p_sge(iph)+(g_sge(iph)-1.0_PR)*pI0)
              rhl(iph)=flrhl0(iph)/fl(iph)
              Yl(iph)=fl(iph)*rhl(iph)/rh
           enddo

           do iph=1,Nl
            elh(iph)=energie_hydro(iph,rhl(iph),p)
            pl(iph)=pressure(iph,rhl(iph),elh(iph))
           enddo

           !!! rh=sum(fl(1:Nl)*rhl(1:Nl))
 
           if(verb)then     
              write(iout,*) '-----------------------------------'
              write(iout,*) '   x=', xin(1:Nsnes)
              write(iout,*) '   rhl=', rhl(1:Nl)
              write(iout,*) '   sumY=',sum(Yl(1:Nl)) 
              write(iout,*) '   sumfl=',sum(fl(1:Nl)) 
              write(iout,*) '   elh=', elh(1:Nl)
              write(iout,*) '   pl=', pl(1:Nl)
              write(iout,*) '   fl=', fl(1:Nl)
              write(iout,*) '   rho=', rh
              write(iout,*) '   e+p/rh=', elh(1:Nl)+p/rhl(1:Nl)
              write(iout,*) '-----------------------------------'
              write(iout,*) '---Set pressure from e---'
           endif    

        !ENDIF !!! pmin-pmax/(pmin+pmax)>1e-3

     ehtot=e-0.5_PR*(vx**2+vy**2+vz**2)

     p=pressure_mixt(Nl,rh,ehtot,fl(1:Nl))
 
     !if(p.lt.0.0_PR)then
     !    verb=.true.
     !    print*, 'coucou2'
     !endif

     do iph=1,Nl
        pl(iph)=p
     enddo

     do iph=1,Nl
       elh(iph)=energie_hydro(iph,rhl(iph),p)
     enddo

     if(abs(sum(Yl(1:Nl)*elh(1:Nl))-ehtot)/ehtot.gt.1.0e-6_PR)then
        write(iout,*) 'PROBLEM ehtot ne sum(Ylehl):'
        write(iout,*) sum(Yl(1:Nl)*elh(1:Nl)), ehtot
        stop
     endif
 
     c=soundspeed_mixt(Nl,fl(1:Nl),Yl(1:Nl),rhl(1:Nl),pl(1:Nl),rh)
 
     if(verb)then 
        write(iout,*) 'P=', p
        write(iout,*) 'elh1=', elh(1)
        write(iout,*) 'elh2=', elh(2)
     endif     

end subroutine relax_p_sge_2

subroutine relax_p_eos(Nl,fl,Yl,rhl,elh,pl,rh,p,vx,vy,vz,c,e)


     implicit none
     integer, intent(in) :: Nl
     real(PR), intent(inout) :: fl(1:Nl)
     real(PR), intent(inout) :: Yl(1:Nl)
     real(PR), intent(inout) :: rhl(1:Nl)
     real(PR), intent(inout) :: elh(1:Nl)
     real(PR), intent(inout) :: pl(1:Nl)
     real(PR), intent(inout) :: p, rh, c
     real(PR), intent(in) :: vx,vy,vz, e
     real(PR) :: pmin, pmax
     !real(PR), allocatable :: xin(:), bin(:)
     integer :: iph
 
        pmin=minval(pl(1:Nl))
        pmax=maxval(pl(1:Nl))
        
        !IF(0.5_PR*abs(pmax-pmin)/(pmax+pmin).gt.1.0e-3_PR)THEN

           do iph=1,Nl
             flrhl0(iph)=fl(iph)*rhl(iph)
             el0(iph)=elh(iph)
             nul0(iph)=1.0_PR/rhl(iph)
             pl0(iph)=pl(iph)  
           enddo 
 
           xin(1:Nl)=nul0(1:Nl)
           xin(Nl+1)=maxval(pl(1:Nl))
           bin(:)=0.0_PR

           do iph=1,Nl
            Zl(iph)=rhl(iph)*soundspeed(iph,rhl(iph),pl(iph))
           enddo
           pI0=sum(Zl(1:Nl)*pl(1:Nl))/sum(Zl(1:Nl))

           if(verb)then
              write(iout,*) 'pI0=', pI0
              write(iout,*) '-----------------------------------'
              write(iout,*) '   x=', xin(1:Nsnes)
              write(iout,*) '   rhl=', rhl(1:Nl)
              write(iout,*) '   sumY=',sum(Yl(1:Nl)) 
              write(iout,*) '   sumfl=',sum(fl(1:Nl)) 
              write(iout,*) '   elh=', elh(1:Nl)
              write(iout,*) '   p=', p
              write(iout,*) '   pl=', pl(1:Nl)
              write(iout,*) '   fl=', fl(1:Nl)
              write(iout,*) '   rho=', rh
           endif     

           !call solve_SNES(Nsnes,xin,bin)
           call My_Newton(1,xin(1:1),minval(pl),maxval(pl),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)

           p=xin(Nl+1)

           if(p.lt.0.0_PR)then
               verb=.true.
               print*, 'coucou1'
           endif

           if(verb)then
              if(mode_pI.eq.2)then
                 write(iout,*) '   e+pI0/rh=', elh(1:Nl)+pI0/rhl(1:Nl)
              else
                 write(iout,*) '   e+p/rh=', elh(1:Nl)+p/rhl(1:Nl)
              endif
              write(iout,*) '-----------------------------------'     
              write(iout,*) '   x=', xin(1:Nl+1)
           endif

           rhl(1:Nl)=1.0_PR/xin(1:Nl)
           fl(1:Nl)=flrhl0(1:Nl)*xin(1:Nl)
           Yl(1:Nl)=fl(1:Nl)*rhl(1:Nl)/rh

           do iph=1,Nl
            elh(iph)=energie_hydro(iph,rhl(iph),p)
            pl(iph)=pressure(iph,rhl(iph),elh(iph))
           enddo

           rh=sum(fl(1:Nl)*rhl(1:Nl))
 
           if(verb)then     
              write(iout,*) '-----------------------------------'
              write(iout,*) '   x=', xin(1:Nl+1)
              write(iout,*) '   rhl=', rhl(1:Nl)
              write(iout,*) '   sumY=',sum(Yl(1:Nl)) 
              write(iout,*) '   sumfl=',sum(fl(1:Nl)) 
              write(iout,*) '   elh=', elh(1:Nl)
              write(iout,*) '   p=', p
              write(iout,*) '   pl=', pl(1:Nl)
              write(iout,*) '   fl=', fl(1:Nl)
              write(iout,*) '   rho=', rh
              if(mode_pI.eq.2)then
                 write(iout,*) '   e+pI0/rh=', elh(1:Nl)+pI0/rhl(1:Nl)
              else
                 write(iout,*) '   e+p/rh=', elh(1:Nl)+p/rhl(1:Nl)
              endif
              write(iout,*) '-----------------------------------'
              write(iout,*) '---Set pressure from e---'
           endif    

        !ENDIF !!! pmin-pmax/(pmin+pmax)>1e-3

     p=pressure_mixt(Nl,rh,e-0.5_PR*(vx**2+vy**2+vz**2),fl(1:Nl)) 
     if(p.lt.0.0_PR)then
         verb=.true.
         print*, 'coucou2'
     endif

     do iph=1,Nl
        pl(iph)=p
     enddo

     do iph=1,Nl
       elh(iph)=energie_hydro(iph,rhl(iph),p)
     enddo
 
     c=soundspeed_mixt(Nl,fl(1:Nl),Yl(1:Nl),rhl(1:Nl),pl(1:Nl),rh)
 
     if(verb)then 
        write(iout,*) 'P=', p
        write(iout,*) 'elh1=', elh(1)
        write(iout,*) 'elh2=', elh(2)
     endif     

end subroutine relax_p_eos

subroutine F_relaxp_sge_1(N,x,f)

     !!!==== RELAXATION DE LA PRESSION SGE:
     !!! on cherche p telle que:
     !!! sum(fl_k*rho_k)*nu_k=1
     !!! avec :
     !!! nu_k(p)=nu_0*(p0+gkpl0+(gk-1)pI0)/(p+gkpl0+(gk-1)pI0)
     !!!    N=Nl+1
     !!!    x(1:Nl)=nu_k(1:Nl)
     !!!    x(Nl+1)=p

     use mod_data, only : Nl

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     real(PR) :: p

     f(1)=sum( fl0_gl(1:Nl)*(pl0(1:Nl)-x(1))/(p_sge(1:Nl)+x(1)) )

end subroutine F_relaxp_sge_1

subroutine J_relaxp_sge_1(N,x,J)

     use mod_data, only : Nl

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)
     real(PR) :: p
     integer :: k

     !!!---dérivée par rapport à p
     !p=x(1)
     J(1,1)=-sum(fl0_gl(1:Nl)*(p_sge(1:Nl)+pl0(1:Nl))/(p_sge(1:Nl)+x(1))**2)

end subroutine J_relaxp_sge_1


subroutine F_relaxp_sge_2(N,x,f)

     !!!==== RELAXATION DE LA PRESSION SGE:
     !!! on cherche p telle que:
     !!! sum(fl_k*rho_k)*nu_k=1
     !!! avec :
     !!! nu_k(p)=nu_0*(p0+gkpl0+(gk-1)pI0)/(p+gkpl0+(gk-1)pI0)
     !!!    N=Nl+1
     !!!    x(1:Nl)=nu_k(1:Nl)
     !!!    x(Nl+1)=p

     use mod_data, only : Nl

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     real(PR) :: p, c1
     integer :: k

     f(1)=sum( fl0(1:Nl)*(pl0(1:Nl)-x(1))/&
             (g_sge(1:Nl)*p_sge(1:Nl) + (g_sge(1:Nl)-1.0_PR)*PI0 + x(1)) )


end subroutine F_relaxp_sge_2

subroutine J_relaxp_sge_2(N,x,J)

     use mod_data, only : Nl

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)
     real(PR) :: p
     integer :: k

     p=x(1)

     !!!---dérivée par rapport à p
     J(1,1)=-sum(fl0(1:Nl)*(g_sge(1:Nl)*p_sge(1:Nl)+(g_sge(1:Nl)-1.0_PR)*pI0+pl0(1:Nl))/&
                (g_sge(1:Nl)*p_sge(1:Nl)+(g_sge(1:Nl)-1.0_PR)*pI0+p)**2 )
  
end subroutine J_relaxp_sge_2

subroutine F_relaxp_1(N,x,f)

     !!!==== RELAXATION DE LA PRESSION :
     !!! on cherche p, nu_k telles que:
     !!! conservation des énbergies phasiques:
     !!!    ek(p,nu_k)-el0+p*(nu_k-nu_0)=0
     !!! conservation de la masse:
     !!!    Sum Yk=1 => Sum fl_k*rho_k*nu_k =1
     !!!    N=Nl+1
     !!!    x(1:Nl)=nu_k(1:Nl)
     !!!    x(Nl+1)=p

     use mod_data, only : Nl

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     real(PR) :: p, rhk
     integer :: k

     p=x(N)

     do k=1,N-1
        rhk=1.0_PR/x(k)
        f(k)=energie_hydro(k,rhk,p)-el0(k)+p*(x(k)-nul0(k))
     enddo   

     f(N)=-1.0_PR
 
     do k=1,N-1
       f(N)=f(N)+flrhl0(k)*x(k)
     enddo

end subroutine F_relaxp_1

subroutine F_relaxp_2(N,x,f)

     !!!==== RELAXATION DE LA PRESSION :
     !!! on cherche p, nu_k telles que:
     !!! conservation des énbergies phasiques:
     !!!    ek(p,nu_k)-el0+p*(nu_k-nu_0)=0
     !!! conservation de la masse:
     !!!    Sum Yk=1 => Sum fl_k*rho_k*nu_k =1
     !!!    N=Nl+1
     !!!    x(1:Nl)=nu_k(1:Nl)
     !!!    x(Nl+1)=p

     use mod_data, only : Nl

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     real(PR) :: p, rhk
     integer :: k

     p=x(N)

     do k=1,N-1
        rhk=1.0_PR/x(k)
        f(k)=energie_hydro(k,rhk,p)-el0(k)+pI0*(x(k)-nul0(k))
     enddo   

     f(N)=-1.0_PR
 
     do k=1,N-1
       f(N)=f(N)+flrhl0(k)*x(k)
     enddo

end subroutine F_relaxp_2

subroutine J_relaxp_1(N,x,J)

     use mod_data, only : Nl

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)
     real(PR) :: p, rhk
     integer :: k

     J(:,:)=0.0_PR

     p=x(N)

     !!!---dérivée par rapport à nu
     do k=1,N-1
        rhk=1.0_PR/x(k)
        J(k,k)=ddnu_energie_hydro(k,rhk,p)+p
     enddo
     !!!---dérivée par rapport à p
     do k=1,N-1
        rhk=1.0_PR/x(k)
        J(k,N)=ddp_energie_hydro(k,rhk,p)+(x(k)-nul0(k))
     enddo

     do k=1,N-1
        J(N,k)=flrhl0(k)
     enddo
 
end subroutine J_relaxp_1

subroutine J_relaxp_2(N,x,J)

     use mod_data, only : Nl

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)
     real(PR) :: p, rhk
     integer :: k

     J(:,:)=0.0_PR

     p=x(N)

     !!!---dérivée par rapport à nu
     do k=1,N-1
        rhk=1.0_PR/x(k)
        J(k,k)=ddnu_energie_hydro(k,rhk,p)+pI0
     enddo
     !!!---dérivée par rapport à p 
     do k=1,N-1
        rhk=1.0_PR/x(k)
        J(k,N)=ddp_energie_hydro(k,rhk,p)
     enddo

     do k=1,N-1
        J(N,k)=flrhl0(k)
     enddo
 
end subroutine J_relaxp_2

!!!========================== SUBROUTINE DE TEST ===============================

subroutine test_fil

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z
 
   implicit none
   integer :: i,j,k,iph

   write(iout,*) ' Cas test tube à choc de SOD air pure'

   Nt=100
   dt=1.0e-5_PR
   t=0.0_PR

   !call init_euler(1000,1,1,1)
   call init_euler(200,1,1,1)

   p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
 
   write(iout,*) '   initialisation de fv,p,v,r' 
   do i=1,nx
      fl(1,i,:,:)=1.0_PR
   enddo

   pl(:,:,:,:)=1.0e5_PR

   do k=1,Nz ; do j=1,ny ; do i=1,nx 
      do iph=1,Nl
        rhl(iph,i,j,k)=rho0(iph)
      enddo
   enddo; enddo; enddo


   rhl(1,1:nx,:,:)=1.0_PR*rho0(1)
   vx(1:nx,:,:)=0.0_PR
   rhl(1,1:nx/2,:,:)=1.0_PR*rho0(1)  ; rhl(1,nx/2+1:nx,:,:)=0.125_PR*rho0(1)
   pl(1,1:nx/2,:,:)=1.0e5_PR ; pl(1,nx/2+1:nx,:,:)=1.0e4_PR

   forall(i=1:nx,j=1:ny,k=1:nz) rh(i,j,k)=sum(fl(1:Nl,i,j,k)*rhl(1:Nl,i,j,k))

   forall(iph=1:Nl,i=1:nx,j=1:ny,k=1:nz) Yl(iph,i,j,k)=fl(iph,i,j,k)*rhl(iph,i,j,k)/rh(i,j,k)

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx ; do iph=1,Nl
     elh(iph,i,j,k)=energie_hydro(iph,rhl(iph,i,j,k),pl(iph,i,j,k))
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx

    call MAJ_meca(Nl,mu(1:Nl),rho0(1:Nl),rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),&
                  al(1:Nl,1:3,i,j,k),bl(1:Nl,1:3,i,j,k),cl(1:Nl,1:3,i,j,k),&
                  pl(1:Nl,i,j,k),&
                  sigl(1:Nl,1:3,1:3,i,j,k),sig(1:3,1:3,i,j,k),ele(1:Nl,i,j,k))!!!S

   enddo ; enddo ; enddo

   write(iout,*) ' MAJ mixt'

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx

   call MAJ_mixt(Nl,rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),Yl(1:Nl,i,j,k),pl(1:Nl,i,j,k),&
                 sigl(1:Nl,1:3,1:3,i,j,k),elh(1:Nl,i,j,k),ele(1:Nl,i,j,k),&
                 vx(i,j,k), vy(i,j,k), vz(i,j,k),&
                 rh(i,j,k),p(i,j,k),sig(1:3,1:3,i,j,k),e(i,j,k),c(i,j,k))

   enddo ; enddo ; enddo

end subroutine test_fil





subroutine test_Ndanou9

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z
   implicit none
   integer :: i,j,k,iph
  
   write(iout,*) '  cas de test Ndanou 9'

   !!!------------------     cas test Fig 9 Ndanou 2015   ---------------------------------------
    call init_euler(2000,1,1,3)
   !! phase 1 : Air 2 : Aluminium 3 : Titane
   p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
   p_sge(2)=p_sge_Al  ; g_sge(2)=g_sge_Al  ; mu(2)=mu_Al  ; sigy(2)=sigy_Al  ; rho0(2)=rho0_Al 
   p_sge(3)=p_sge_Ti  ; g_sge(3)=g_sge_Ti  ; mu(3)=mu_Ti  ; sigy(3)=sigy_Ti  ; rho0(3)=rho0_Ti 

   mu(1:3)=0.0_PR
 
   write(iout,*) '   initialisation de fv,p,v,r' 
   do i=1,nx
     if(x(i).lt.0.001_PR)then
      fl(1,i,:,:)=1.0_PR
     elseif(x(i).lt.0.003_PR)then
      fl(2,i,:,:)=1.0_PR ; vx(i,:,:)=700.0_PR
     elseif(x(i).lt.0.0128_PR)then
      fl(3,i,:,:)=1.0_PR
     else
      fl(1,i,:,:)=1.0_PR
     endif
   enddo
   
   pl(:,:,:,:)=1.0e5_PR

   do k=1,Nz ; do j=1,ny ; do i=1,nx 
      do iph=1,Nl
        rhl(iph,i,j,k)=rho0(iph)
      enddo
   enddo; enddo; enddo

end subroutine test_Ndanou9

subroutine test_SOD_air

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z
 
   implicit none
   integer :: i,j,k,iph

   write(iout,*) ' Cas test tube à choc de SOD air pure'

   Nt=100
   dt=1.0e-5_PR
   t=0.0_PR

   !call init_euler(1000,1,1,1)
   call init_euler(200,1,1,1)

   p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
 
   write(iout,*) '   initialisation de fv,p,v,r' 
   do i=1,nx
      fl(1,i,:,:)=1.0_PR
   enddo

   pl(:,:,:,:)=1.0e5_PR

   do k=1,Nz ; do j=1,ny ; do i=1,nx 
      do iph=1,Nl
        rhl(iph,i,j,k)=rho0(iph)
      enddo
   enddo; enddo; enddo


   rhl(1,1:nx,:,:)=1.0_PR*rho0(1)
   vx(1:nx,:,:)=0.0_PR
   rhl(1,1:nx/2,:,:)=1.0_PR*rho0(1)  ; rhl(1,nx/2+1:nx,:,:)=0.125_PR*rho0(1)
   pl(1,1:nx/2,:,:)=1.0e5_PR ; pl(1,nx/2+1:nx,:,:)=1.0e4_PR

   forall(i=1:nx,j=1:ny,k=1:nz) rh(i,j,k)=sum(fl(1:Nl,i,j,k)*rhl(1:Nl,i,j,k))

   forall(iph=1:Nl,i=1:nx,j=1:ny,k=1:nz) Yl(iph,i,j,k)=fl(iph,i,j,k)*rhl(iph,i,j,k)/rh(i,j,k)

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx ; do iph=1,Nl
     elh(iph,i,j,k)=energie_hydro(iph,rhl(iph,i,j,k),pl(iph,i,j,k))
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx

    call MAJ_meca(Nl,mu(1:Nl),rho0(1:Nl),rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),&
                  al(1:Nl,1:3,i,j,k),bl(1:Nl,1:3,i,j,k),cl(1:Nl,1:3,i,j,k),&
                  pl(1:Nl,i,j,k),&
                  sigl(1:Nl,1:3,1:3,i,j,k),sig(1:3,1:3,i,j,k),ele(1:Nl,i,j,k))!!!S

   enddo ; enddo ; enddo

   write(iout,*) ' MAJ mixt'

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx

   call MAJ_mixt(Nl,rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),Yl(1:Nl,i,j,k),pl(1:Nl,i,j,k),&
                 sigl(1:Nl,1:3,1:3,i,j,k),elh(1:Nl,i,j,k),ele(1:Nl,i,j,k),&
                 vx(i,j,k), vy(i,j,k), vz(i,j,k),&
                 rh(i,j,k),p(i,j,k),sig(1:3,1:3,i,j,k),e(i,j,k),c(i,j,k))

   enddo ; enddo ; enddo

end subroutine test_SOD_air

subroutine test_convection

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z
 
   implicit none
   integer :: i,j,k,iph
 
   write(iout,*) ' cas test convection'  

   call init_euler(200,1,1,2)
   Nt=2000
   dt=1.0e-6_PR
   t=0.0_PR

   !!!                      Nl=2, Lx=1, Nx=200, Nit=2000 dt=1e-6
   p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
   p_sge(2)=p_sge_eau ; g_sge(2)=g_sge_eau ; mu(2)=mu_eau ; sigy(2)=sigy_eau ; rho0(2)=rho0_eau 

   write(iout,*) '   initialisation de fv,p,v,r'
 
   do i=1,nx/2.0_PR
      fl(1,i,:,:)=1.0e-8_PR
      fl(2,i,:,:)=1.0_PR-1.0e-8_PR
   enddo
   do i=nx/2+1,nx
      fl(1,i,:,:)=1.0_PR-1.0e-8_PR
      fl(2,i,:,:)=1.0e-8_PR
   enddo

   pl(:,:,:,:)=1.0e5_PR

   do k=1,Nz ; do j=1,ny ; do i=1,nx 
      rhl(1,i,j,k)=1.0_PR
      rhl(2,i,j,k)=1000.0_PR
   enddo; enddo; enddo

   vx(1:nx,:,:)=100.0_PR

end subroutine test_convection

subroutine test_eau_air_strong

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z

   implicit none
   integer :: i,j,k,iph

    write(iout,*) ' Cas test choc eau-air'

    call init_euler(1000,1,1,2)
    Nt=90
    dt=1.0e-7_PR
    t=0.0_PR

    !!!                      Nl=2 Lx=1 Nx=100 Nit=240 dt=1e-7
     p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
     p_sge(2)=p_sge_eau ; g_sge(2)=g_sge_eau ; mu(2)=mu_eau ; sigy(2)=sigy_eau ; rho0(2)=rho0_eau 

     write(iout,*) '   initialisation de fv,p,v,r'
 
     do i=1,nx
        if(x(i).le.0.75_PR)then
           fl(1,i,:,:)=1.0e-6_PR
           fl(2,i,:,:)=1.0_PR-1.0e-6_PR
           pl(:,i,:,:)=1.0e12_PR
        endif
     enddo
     do i=1,nx
        if(x(i).gt.0.75_PR)then
           fl(1,i,:,:)=1.0_PR-1.0e-6_PR
           fl(2,i,:,:)=1.0e-6_PR
           pl(:,i,:,:)=1.0e5_PR
        endif
     enddo

     do k=1,Nz ; do j=1,ny ; do i=1,nx 
        rhl(1,i,j,k)=10.0_PR
        rhl(2,i,j,k)=1000.0_PR
     enddo; enddo; enddo

     forall(i=1:nx,j=1:ny,k=1:nz) rh(i,j,k)=sum(fl(1:Nl,i,j,k)*rhl(1:Nl,i,j,k))

     forall(iph=1:Nl,i=1:nx,j=1:ny,k=1:nz) Yl(iph,i,j,k)=fl(iph,i,j,k)*rhl(iph,i,j,k)/rh(i,j,k)

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx ; do iph=1,Nl
       elh(iph,i,j,k)=energie_hydro(iph,rhl(iph,i,j,k),pl(iph,i,j,k))
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx

      call MAJ_meca(Nl,mu(1:Nl),rho0(1:Nl),rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),&
                    al(1:Nl,1:3,i,j,k),bl(1:Nl,1:3,i,j,k),cl(1:Nl,1:3,i,j,k),&
                    pl(1:Nl,i,j,k),&
                    sigl(1:Nl,1:3,1:3,i,j,k),sig(1:3,1:3,i,j,k),ele(1:Nl,i,j,k))!!!S

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx

     call MAJ_mixt(Nl,rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),Yl(1:Nl,i,j,k),pl(1:Nl,i,j,k),&
                   sigl(1:Nl,1:3,1:3,i,j,k),elh(1:Nl,i,j,k),ele(1:Nl,i,j,k),&
                   vx(i,j,k), vy(i,j,k), vz(i,j,k),&
                   rh(i,j,k),p(i,j,k),sig(1:3,1:3,i,j,k),e(i,j,k),c(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_eau_air_strong


subroutine test_eau_air

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z

   implicit none
   integer :: i,j,k,iph

    write(iout,*) ' Cas test choc eau-air'

    call init_euler(1000,1,1,2)
    Nt=240
    dt=1.0e-6_PR
    t=0.0_PR

    !!!                      Nl=2 Lx=1 Nx=100 Nit=240 dt=1e-7
     p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
     p_sge(2)=p_sge_eau ; g_sge(2)=g_sge_eau ; mu(2)=mu_eau ; sigy(2)=sigy_eau ; rho0(2)=rho0_eau 

     write(iout,*) '   initialisation de fv,p,v,r'
 
     do i=1,nx
        if(x(i).le.0.75_PR)then
           fl(1,i,:,:)=1.0e-6_PR
           fl(2,i,:,:)=1.0_PR-1.0e-6_PR
           pl(:,i,:,:)=1.0e9_PR
        endif
     enddo
     do i=1,nx
        if(x(i).gt.0.75_PR)then
           fl(1,i,:,:)=1.0_PR-1.0e-6_PR
           fl(2,i,:,:)=1.0e-6_PR
           pl(:,i,:,:)=1.0e5_PR
        endif
     enddo

     do k=1,Nz ; do j=1,ny ; do i=1,nx 
        rhl(1,i,j,k)=1.0_PR
        rhl(2,i,j,k)=1000.0_PR
     enddo; enddo; enddo

     forall(i=1:nx,j=1:ny,k=1:nz) rh(i,j,k)=sum(fl(1:Nl,i,j,k)*rhl(1:Nl,i,j,k))

     forall(iph=1:Nl,i=1:nx,j=1:ny,k=1:nz) Yl(iph,i,j,k)=fl(iph,i,j,k)*rhl(iph,i,j,k)/rh(i,j,k)

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx ; do iph=1,Nl
       elh(iph,i,j,k)=energie_hydro(iph,rhl(iph,i,j,k),pl(iph,i,j,k))
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx

      call MAJ_meca(Nl,mu(1:Nl),rho0(1:Nl),rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),&
                    al(1:Nl,1:3,i,j,k),bl(1:Nl,1:3,i,j,k),cl(1:Nl,1:3,i,j,k),&
                    pl(1:Nl,i,j,k),&
                    sigl(1:Nl,1:3,1:3,i,j,k),sig(1:3,1:3,i,j,k),ele(1:Nl,i,j,k))!!!S

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx

     call MAJ_mixt(Nl,rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),Yl(1:Nl,i,j,k),pl(1:Nl,i,j,k),&
                   sigl(1:Nl,1:3,1:3,i,j,k),elh(1:Nl,i,j,k),ele(1:Nl,i,j,k),&
                   vx(i,j,k), vy(i,j,k), vz(i,j,k),&
                   rh(i,j,k),p(i,j,k),sig(1:3,1:3,i,j,k),e(i,j,k),c(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_eau_air

subroutine test_rar_rar

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z

   implicit none
   integer :: i,j,k,iph

   write(iout,*) ' Cas test rar-rar'

   Nt=100
   dt=1.0e-5_PR
   t=0.0_PR

   call init_euler(200,1,1,1)

   p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
   !p_sge(1)=p_sge_eau ; g_sge(1)=g_sge_eau ; mu(1)=mu_eau ; sigy(1)=sigy_eau ; rho0(1)=rho0_eau 

   write(iout,*) '   initialisation de fv,p,v,r' 
   do i=1,nx
      fl(1,i,:,:)=1.0_PR
   enddo

   pl(:,:,:,:)=1.0e5_PR

   do k=1,Nz ; do j=1,ny ; do i=1,nx 
      do iph=1,Nl
        rhl(iph,i,j,k)=rho0(iph)
      enddo
   enddo; enddo; enddo


   rhl(1,1:nx,:,:)=1.0_PR*rho0(1)
   vx(1:nx/2,:,:)=-100.0_PR
   vx(nx/2+1:nx,:,:)=100.0_PR

   forall(i=1:nx,j=1:ny,k=1:nz) rh(i,j,k)=sum(fl(1:Nl,i,j,k)*rhl(1:Nl,i,j,k))

   forall(iph=1:Nl,i=1:nx,j=1:ny,k=1:nz) Yl(iph,i,j,k)=fl(iph,i,j,k)*rhl(iph,i,j,k)/rh(i,j,k)

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx ; do iph=1,Nl
     elh(iph,i,j,k)=energie_hydro(iph,rhl(iph,i,j,k),pl(iph,i,j,k))
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx

    call MAJ_meca(Nl,mu(1:Nl),rho0(1:Nl),rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),&
                  al(1:Nl,1:3,i,j,k),bl(1:Nl,1:3,i,j,k),cl(1:Nl,1:3,i,j,k),&
                  pl(1:Nl,i,j,k),&
                  sigl(1:Nl,1:3,1:3,i,j,k),sig(1:3,1:3,i,j,k),ele(1:Nl,i,j,k))!!!S

   enddo ; enddo ; enddo

   write(iout,*) ' MAJ mixt'

   do k=1,Nz ; do j=1,Ny ; do i=1,Nx

   call MAJ_mixt(Nl,rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),Yl(1:Nl,i,j,k),pl(1:Nl,i,j,k),&
                 sigl(1:Nl,1:3,1:3,i,j,k),elh(1:Nl,i,j,k),ele(1:Nl,i,j,k),&
                 vx(i,j,k), vy(i,j,k), vz(i,j,k),&
                 rh(i,j,k),p(i,j,k),sig(1:3,1:3,i,j,k),e(i,j,k),c(i,j,k))

   enddo ; enddo ; enddo

end subroutine test_rar_rar

subroutine test_cavitation

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z

   implicit none
   integer :: i,j,k,iph

    write(iout,*) ' Cas test cavitation'

    call init_euler(1000,1,1,2)
    !call init_euler(2000,1,1,2)

    Nt=200
    dt=1.0e-5_PR
    t=0.0_PR
    !dtmin=1.0e-7_PR

    !!!                      Nl=2 Lx=1 Nx=100 Nit=240 dt=1e-7
     p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
     p_sge(2)=p_sge_eau ; g_sge(2)=g_sge_eau ; mu(2)=mu_eau ; sigy(2)=sigy_eau ; rho0(2)=rho0_eau 
     !p_sge(1)=p_sge_eau ; g_sge(1)=g_sge_eau ; mu(1)=mu_eau ; sigy(1)=sigy_eau ; rho0(1)=rho0_eau 
     !p_sge(2)=p_sge_air ; g_sge(2)=g_sge_air ; mu(2)=mu_air ; sigy(2)=sigy_air ; rho0(2)=rho0_air 

     write(iout,*) '   initialisation de fv,p,v,r'

     do i=1,nx
        fl(1,i,:,:)=1.0e-2_PR
        fl(2,i,:,:)=1.0_PR-1.0e-2_PR
        pl(:,i,:,:)=1.0e5_PR
     enddo
 
     do i=1,nx
        if(x(i).le.0.5_PR)then
           vx(i,:,:)=-100.0_PR
        elseif(x(i).gt.0.5_PR)then
           vx(i,:,:)=100.0_PR
        endif

        !if(x(i).le.0.45_PR)then
        !   vx(i,:,:)=-100.0_PR
        !elseif(x(i).gt.0.55_PR)then
        !   vx(i,:,:)=100.0_PR
        !else
        !   vx(i,:,:)=-100.0_PR+200.0_PR/0.10_PR*(x(i)-0.45_PR)
        !endif
     enddo

     do k=1,Nz ; do j=1,ny ; do i=1,nx 
        rhl(1,i,j,k)=1.0_PR*rho0(1)
        rhl(2,i,j,k)=rho0(2)
     enddo; enddo; enddo

     forall(i=1:nx,j=1:ny,k=1:nz) rh(i,j,k)=sum(fl(1:Nl,i,j,k)*rhl(1:Nl,i,j,k))

     forall(iph=1:Nl,i=1:nx,j=1:ny,k=1:nz) Yl(iph,i,j,k)=fl(iph,i,j,k)*rhl(iph,i,j,k)/rh(i,j,k)

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx ; do iph=1,Nl
       elh(iph,i,j,k)=energie_hydro(iph,rhl(iph,i,j,k),pl(iph,i,j,k))
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx

      call MAJ_meca(Nl,mu(1:Nl),rho0(1:Nl),rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),&
                    al(1:Nl,1:3,i,j,k),bl(1:Nl,1:3,i,j,k),cl(1:Nl,1:3,i,j,k),&
                    pl(1:Nl,i,j,k),&
                    sigl(1:Nl,1:3,1:3,i,j,k),sig(1:3,1:3,i,j,k),ele(1:Nl,i,j,k))!!!S

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx

     call MAJ_mixt(Nl,rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),Yl(1:Nl,i,j,k),pl(1:Nl,i,j,k),&
                   sigl(1:Nl,1:3,1:3,i,j,k),elh(1:Nl,i,j,k),ele(1:Nl,i,j,k),&
                   vx(i,j,k), vy(i,j,k), vz(i,j,k),&
                   rh(i,j,k),p(i,j,k),sig(1:3,1:3,i,j,k),e(i,j,k),c(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_cavitation

subroutine test_cavitation2

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z

   implicit none
   integer :: i,j,k,iph

    write(iout,*) ' Cas test cavitation 2'

    call init_euler(2000,1,1,2)

    Nt=100
    dt=1.0e-6_PR
    t=0.0_PR
    dtmin=5.0e-8_PR

    !!!                      Nl=2 Lx=1 Nx=100 Nit=240 dt=1e-7
     p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
     p_sge(2)=p_sge_Cu ; g_sge(2)=g_sge_Cu ; mu(2)=mu_Cu ; sigy(2)=sigy_Cu ; rho0(2)=rho0_Cu 

     write(iout,*) '   initialisation de fv,p,v,r'

     do i=1,nx
        fl(1,i,:,:)=1.0e-4_PR
        fl(2,i,:,:)=1.0_PR-1.0e-4_PR
        pl(:,i,:,:)=1.0e5_PR
     enddo
 
     do i=1,nx
        if(x(i).le.0.5_PR)then
           vx(i,:,:)=-3000.0_PR
        endif
     enddo
     do i=1,nx
        if(x(i).gt.0.5_PR)then
           vx(i,:,:)=3000.0_PR
        endif
     enddo

     do k=1,Nz ; do j=1,ny ; do i=1,nx 
        rhl(1,i,j,k)=1.0_PR
        rhl(2,i,j,k)=8900.0_PR
     enddo; enddo; enddo

     forall(i=1:nx,j=1:ny,k=1:nz) rh(i,j,k)=sum(fl(1:Nl,i,j,k)*rhl(1:Nl,i,j,k))

     forall(iph=1:Nl,i=1:nx,j=1:ny,k=1:nz) Yl(iph,i,j,k)=fl(iph,i,j,k)*rhl(iph,i,j,k)/rh(i,j,k)

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx ; do iph=1,Nl
       elh(iph,i,j,k)=energie_hydro(iph,rhl(iph,i,j,k),pl(iph,i,j,k))
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx

      call MAJ_meca(Nl,mu(1:Nl),rho0(1:Nl),rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),&
                    al(1:Nl,1:3,i,j,k),bl(1:Nl,1:3,i,j,k),cl(1:Nl,1:3,i,j,k),&
                    pl(1:Nl,i,j,k),&
                    sigl(1:Nl,1:3,1:3,i,j,k),sig(1:3,1:3,i,j,k),ele(1:Nl,i,j,k))!!!S

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,Nz ; do j=1,Ny ; do i=1,Nx

     call MAJ_mixt(Nl,rhl(1:Nl,i,j,k),fl(1:Nl,i,j,k),Yl(1:Nl,i,j,k),pl(1:Nl,i,j,k),&
                   sigl(1:Nl,1:3,1:3,i,j,k),elh(1:Nl,i,j,k),ele(1:Nl,i,j,k),&
                   vx(i,j,k), vy(i,j,k), vz(i,j,k),&
                   rh(i,j,k),p(i,j,k),sig(1:3,1:3,i,j,k),e(i,j,k),c(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_cavitation2




subroutine test_compression

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z

   implicit none
   integer :: i,j,k,iph

   write(iout,*) ' Cas_ test cylindre en compression'

   call init_euler(10,1,1,3)

   !!!------------------     cas test cylindre compression  ---------------------------------------
   !! phase 1 : Air 2 : Aluminium 3 : Titane
   p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
   p_sge(2)=p_sge_Al  ; g_sge(2)=g_sge_Al  ; mu(2)=mu_Al  ; sigy(2)=sigy_Al  ; rho0(2)=rho0_Al 
   p_sge(3)=p_sge_Ti  ; g_sge(3)=g_sge_Ti  ; mu(3)=mu_Ti  ; sigy(3)=sigy_Ti  ; rho0(3)=rho0_Ti 
 
   write(iout,*) '   initialisation de fv,p,v,r' 
   do i=1,nx
      fl(2,i,:,:)=1.0_PR
   enddo   
   pl(:,:,:,:)=1.0e5_PR
   do k=1,Nz ; do j=1,ny ; do i=1,nx 
      do iph=1,Nl
        rhl(iph,i,j,k)=rho0(iph)
      enddo
   enddo; enddo; enddo

   !!!---compression selon -z
   al(:,1,:,:,:)=1.1_PR


end subroutine test_compression

!subroutine test_snes
!
!   implicit none
!
!   !!!---test petsc SNES
!    call init_euler(1,1,1,1)
!    Nt=0
!
!   print*, 'SNESinterface start !'
!   allocate(xin(1:2)) ; xin(:)=0.0_PR
!   allocate(bin(1:2)) ; bin(:)=0.0_PR
!   
!   xin(1)=1.0_PR
!   xin(2)=5.0_PR
!   
!   write(iout,*) 'TEST 1 solution: x = 1 2'
!
!   call init_SNES(2,xin,bin,F_function=F_test1,J_function=J_test1)
!   call solve_SNES(2,xin,bin)
!  
!   write(iout,*) '   x=', xin(1:2)
!   
!   write(iout,*) 'TEST 2 solution: x = 1.07716 0.0897632'
!
!   restart_snes=.true.
!   call init_SNES(2,xin,bin,F_function=F_test2,J_function=J_test2)
!   call solve_SNES(2,xin,bin)
!
!   write(iout,*) '   x=', xin(1:2)
!
!   stop
!
!end subroutine test_snes


subroutine test_relax

   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
                        al,bl,cl,sigl,sig,&
                        rh,vx,vy,vz,c,p,e,&
                        x,y,z

   implicit none
   integer :: iph

   write(iout,*) ' Cas test relaxation pression'

   call init_euler(1,1,1,2)

   p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
   p_sge(2)=p_sge_eau ; g_sge(2)=g_sge_eau ; mu(2)=mu_eau ; sigy(2)=sigy_eau ; rho0(2)=rho0_eau 


   rhl(1,1,1,1)=1.0053_PR
   rhl(2,1,1,1)=874.11_PR
   pl(1,1,1,1)=100749.0_PR
   pl(2,1,1,1)=270044072.0_PR
   fl(1,1,1,1)=0.99236_PR
   fl(2,1,1,1)=1.0_PR-fl(1,1,1,1)
   p(1,1,1)=2.162638e6_PR 
   vx(1,1,1)=334.91569_PR
   e(1,1,1)=966143.18905909313_PR

 
   rh(1,1,1)=sum(fl(1:Nl,1,1,1)*rhl(1:Nl,1,1,1))
 
   do iph=1,Nl
     elh(iph,1,1,1)=energie_hydro(iph,rhl(iph,1,1,1),pl(iph,1,1,1))
     Yl(iph,1,1,1)=fl(iph,1,1,1)*rhl(iph,1,1,1)/rh(1,1,1)
   enddo 
   c(1,1,1)=soundspeed_mixt(Nl,fl(1:Nl,1,1,1),Yl(1:Nl,1,1,1),rhl(1:Nl,1,1,1),pl(1:Nl,1,1,1),rh(1,1,1))

   verb=.true.

   call relax_p(Nl,fl(1:Nl,1,1,1),Yl(1:Nl,1,1,1),rhl(1:Nl,1,1,1),&
                elh(1:Nl,1,1,1),pl(1:Nl,1,1,1),rh(1,1,1),p(1,1,1),&
                vx(1,1,1),vy(1,1,1),vz(1,1,1),c(1,1,1),e(1,1,1))

   stop
 
end subroutine test_relax

!!!==== test SNES : problem 1
!!! solution x(1)=1, x(2)=2
subroutine F_test1(N,x,f)

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)

     !!---exemple 1 : x(1)=1 ; x(2)=2
     f(1)=x(1)**2+x(1)*x(2)-3.0_PR
     f(2)=x(1)*x(2)+x(2)**2-6.0_PR
    
end subroutine F_test1

subroutine J_test1(N,x,J)

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)

     !!!---exemple 1
     J(1,1)=2.0_PR*x(1)+x(2) ; J(1,2)=x(1)
     J(2,1)=x(2)             ; J(2,2)=x(1)+2.0_PR*x(2)

end subroutine J_test1

!!!==== test SNES : problem 1
!!! solution x(1)=1.07716, x(2)=0.0897632
subroutine F_test2(N,x,f)

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)

     !!---exemple 2: result: x=1.07716 ; y=0.0897632
     f(1)=sin(3.0_PR*x(1))+x(2)
     f(2)=x(1)*x(2)-12.0_PR*x(2)**2
    
end subroutine F_test2

subroutine J_test2(N,x,J)

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)

     !!!---exemple 2:
     J(1,1)=3.0_PR*cos(3.0_PR*x(1)) ; J(1,2)=1.0_PR
     J(2,1)=x(2)                    ; J(2,2)=x(1)-24.0_PR*x(2)

end subroutine J_test2

!!!=============================== RECONSTRUCTION ===================================

subroutine NOREC(Wim1,Wi,Wip1,Wip2,WL,WR)

   implicit none
   real(PR), intent(in) :: Wim1, Wi,Wip1,Wip2
   real(PR), intent(out) :: WL, WR

   WL=Wi  
   WR=Wip1

end subroutine NOREC

!SUBROUTINE MUSCL(v,v_W,v_E,nx,nr,nvar,ngh,limit,dir) !!maillage uniforme
!
!implicit none
!INTEGER, intent(in) :: nx,nr,nvar,ngh
!real(dp),dimension(-ngh+1:nx+ngh,-ngh+1:nr+ngh,1:nvar),intent(in) ::  v
!!real(dp),intent(in) :: dx(-ngh+1:nx+ngh), dxm(-ngh+1:nx+ngh)
!!real(dp),intent(in) :: dr(-ngh+1:nr+ngh), drm(-ngh+1:nr+ngh)
!real(dp),dimension(-ngh+1:nx+ngh,-ngh+1:nr+ngh,1:nvar),intent(out) :: v_E,v_W
!real(dp) :: PHI !!! 1:centre, -1:decentre, 1/3 : ordre 3
!real(dp) :: psiL1,psiL2,psiR1,psiR2,rL,rR
!logical,intent(in) :: limit 
!character(len=1), intent(in) :: dir
!integer :: i,j,k
!PHI=0.33333333d0 !!!!1.d0/3.d0
!
!if(dir.eq.'X')then
! do k=1,nvar
!  do j=1,nr
!   do i=0,nx+1
!  
!   if((v(i+1,j,k)-v(i,j,k)).eq.0.d0.and.(v(i,j,k)-v(i-1,j,k)).eq.0.d0)then
!    rL=1.d0 ; rR=1.d0
!   else
!    rL=(v(i+1,j,k)-v(i,j,k))/(v(i,j,k)-v(i-1,j,k))
!    rR=(v(i,j,k)-v(i-1,j,k))/(v(i+1,j,k)-v(i,j,k))
!   endif
!
!call limiteur(rL,psiL1); call limiteur(1.d0/rL,psiL2)
!call limiteur(rR,psiR1); call limiteur(1.d0/rR,psiR2)
!
!v_E(i,j,k)=v(i,j,k)+0.25d0*((1.d0-PHI)*psiL1*(v(i,j,k)-v(i-1,j,k))+(1.d0+PHI)*psiL2*(v(i+1,j,k)-v(i,j,k)))
!v_W(i,j,k)=v(i,j,k)+0.25d0*(-(1.d0+PHI)*psiR1*(v(i,j,k)-v(i-1,j,k))-(1.d0-PHI)*psiR2*(v(i+1,j,k)-v(i,j,k)))
!
!!v_E(i,j,k)=v(i,j,k)+capa(i-1,j)*0.25d0*((1.d0-PHI)*psiL1*(v(i,j,k)-v(i-1,j,k))+(1.d0+PHI)*psiL2*(v(i+1,j,k)-v(i,j,k)))
!!v_W(i,j,k)=v(i,j,k)+capa(i+1,j)*0.25d0*(-(1.d0+PHI)*psiR1*(v(i,j,k)-v(i-1,j,k))-(1.d0-PHI)*psiR2*(v(i+1,j,k)-v(i,j,k)))
!   enddo
! enddo
!enddo
!

subroutine MUSCL(Wim1,Wi,Wip1,Wip2,WL,WR)

   implicit none
   real(PR), intent(in) :: Wim1, Wi,Wip1,Wip2
   real(PR), intent(out) :: WL, WR
   real(PR) :: ri,rip1
   real(PR) :: PHI=1.0_PR/3.0_PR

   !!!---http://chimeracfd.com/programming/gryphon/muscl.html
   if((Wip1-Wi).eq.0.0_PR.and.(Wi-Wim1).eq.0)then
      WL=Wi 
   elseif((Wip1-Wi).eq.0.0_PR)then
      WL=Wi+0.25_PR*Limiteur(0.0_PR)*(1.0_PR-PHI)*(Wi-Wim1)
   elseif((Wi-Wim1).eq.0.0_PR)then
      WL=Wi+0.25_PR*Limiteur(0.0_PR)*(1.0_PR+PHI)*(Wip1-Wi)
   else
      ri=(Wi-Wim1)/(Wip1-Wi)
      WL=Wi+0.25_PR*( Limiteur(1.0_PR/ri)*(1.0_PR-PHI)*(Wi-Wim1) + Limiteur(ri)*(1.0_PR+PHI)*(Wip1-Wi) )
   endif

   if((Wip2-Wip1).eq.0.0_PR.and.(Wip1-Wi).eq.0.0_PR)then
      WR=Wip1
   elseif((Wip2-Wip1).eq.0.0_PR)then
      WR=Wip1-0.25_PR*Limiteur(0.0_PR)*(1.0_PR+PHI)*(Wip1-Wi)
   elseif((Wip1-Wi).eq.0.0_PR)then
      WR=Wip1-0.25_PR*Limiteur(0.0_PR)*(1.0_PR-PHI)*(Wip2-Wip1)
   else
      rip1=(Wip1-Wi)/(Wip2-Wip1)
      WR=Wip1-0.25_PR*( Limiteur(1.0_PR/rip1)*(1.0_PR+PHI)*(Wip1-Wi) + Limiteur(rip1)*(1.0_PR-PHI)*(Wip2-Wip1) )
   endif

   !!!---wiki :
   !if((Wip1-Wi).eq.0.0_PR.and.(Wi-Wim1).eq.0)then
   !   WL=Wi 
   !elseif((Wip1-Wi).eq.0.0_PR)then
   !elseif((Wi-Wim1).eq.0.0_PR)then
   !else
   !   ri=(Wi-Wim1)/(Wip1-Wi)
   !   WL=Wi+0.25_PR*Limiteur(ri)*( (1.0_PR-PHI)*(Wi-Wim1) + (1.0_PR+PHI)*(Wip1-Wi) )
   !endif

   !if((Wip2-Wip1).eq.0.0_PR.and.(Wip1-Wi).eq.0.0_PR)then
   !   WR=Wip1
   !elseif((Wip2-Wip1).eq.0.0_PR)then
   !   rip1=sign(1.0_PR,Wip1-Wi)*1.0e30_PR
   !   WR=Wip1-0.25_PR*Limiteur(rip1)*(1.0_PR+PHI)*(Wip1-Wi)
   !elseif((Wip1-Wi).eq.0.0_PR)then
   !   WR=Wip1-0.25_PR*Limiteur(0.0_PR)*(1.0_PR-PHI)*(Wip2-Wip1)
   !else
   !   rip1=(Wip1-Wi)/(Wip2-Wip1)
   !   WR=Wip1-0.25_PR*Limiteur(rip1)*( (1.0_PR+PHI)*(Wip1-Wi) + (1.0_PR-PHI)*(Wip2-Wip1) )
   !endif

end subroutine MUSCL

!!!=============================== LIMITEURS ===================================

real(PR) function MINMOD(r)

   implicit none
   real(PR),intent(in) :: r

   minmod=max(0.d0,min(1.d0,r))

end function MINMOD

real(PR) function VANLEER(r)

   implicit none
   real(PR),intent(in) :: r
   vanleer=(r+abs(r))/(1.d0+abs(r))

end function VANLEER

real(PR) function VANALBADA(r)

   implicit none
   real(PR),intent(in) :: r

   vanalbada=max(0.d0,(r+r*r)/(1.d0+r*r))

end function VANALBADA

real(PR) function SUPERBEE(r)

   implicit none
   real(PR),intent(in) :: r

   superbee=max(0.d0,min(1.d0,2.d0*r),min(2.d0,r))

end function SUPERBEE

!!!================================== 1D ================================

subroutine src_euler1D(Nl,N,U,G,dxm,dUdt,dGdt)

   implicit none
   integer, intent(in) :: Nl, N
   real(PR), intent(in) :: dxm(1:N)
   real(PR), intent(in) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
   real(PR), intent(out) :: dUdt(1:Nl+4,1:N), dGdt(1:Nl,1:11,1:N)
   !!!---Left/right
   real(PR) :: UL(1:Nl+4), GL(1:Nl,1:11)
   real(PR) :: UR(1:Nl+4), GR(1:Nl,1:11)

   !!!---champs 1D locaux
   !real(PR), dimension(Nl,N) :: fl, Yl, pl, rhl, elh, ele
   !real(PR), dimension(Nl,3,N) :: al,bl,cl
   !real(PR), dimension(Nl,3,3,N) :: sigl
   !real(PR), dimension(3,3,N) :: sig
   !real(PR), dimension(N) :: rh, vx, vy,vz, z, c, e, p
   ! !!!---Left/Right variables:   
   ! real(PR), dimension(Nl) :: fl_L, fl_R, pl_L, pl_R, rhl_L, rhl_R, elh_L, elh_R
   ! real(PR), dimension(Nl,3) :: al_L,bl_L,cl_L,al_R,bl_R,cl_R
   ! real(PR), dimension(3) :: sig_L, sig_R !!! s11, s12, s13
   ! real(PR) :: rh_L,rh_R,vx_L,vx_R,vy_L,vy_R,vz_L,vz_R,e_L,e_R,c_L,c_R
   !!!---Grandeurs étoile
   real(PR) :: vx_s, vy_s, vz_s
   !!!---vecteurs flux pour la partie conservative (eq 19 Ndanou 2015)
   real(PR) :: F1(1:Nl+4)
   !!!---vecteurs flux pour la partie non-conservative (eq 19 Ndanou 2015)
   real(PR) :: F2(1:Nl,1:11),H2(1:Nl,1:11),K2(1:Nl,1:11),M2(1:Nl,1:11)
   integer :: i, ip1, im1, ip2, iph, i2

   call cpu_time(t10) 

   dUdt=0.0_PR
   dGdt=0.0_PR

   do i=1,N-1

      im1=max(i-1,1) ; ip1=i+1 ; ip2=min(i+2,N)

      do i2=1,Nl+4
        call Reco(U(i2,im1),U(i2,i),U(i2,ip1),U(i2,ip2),UL(i2),UR(i2))
      enddo  

      GL=G(1:Nl,1:11,i) ; GR=G(1:Nl,1:11,ip1)

      call HLLC_flux(Nl,UL,UR,GL,GR,wL,wR,vx_s,vy_s,vz_s,F1,F2)

      !!!-----------------  PARTIE CONSERVATIVE --------------------
      !if(2Daxi)then

      !dUdt(:,i)=dUdt(:,i)-S(i)/vol(i)*F1(:)
      !dUdt(:,i+1)=dUdt(:,i+1)+S(i)/vol(i+1)*F1(:)
      !dUdt(Nl+1,i)=dUdt(Nl+1,i)+sum(H2(1:Nl,11))/x(i)
      !dUdt(Nl+1,i+1)=dUdt(Nl+1,i+1)+sum(H2(1:Nl,11))/x(i)


      !else

      dUdt(:,i)=dUdt(:,i)-1.0_PR/dxm(i)*F1(:)
      dUdt(:,i+1)=dUdt(:,i+1)+1.0_PR/dxm(i)*F1(:)

      !endif

      !!!---------------  PARTIE NON CONSERVATIVE ------------------

      !GL=G(1:Nl,1:11,i) ; GR=G(1:Nl,1:11,ip1)

      !!!   terme source i :

      call get_HKM(Nl,GL,H2,K2,M2)

      dGdt(1:Nl,1:11,i)  =dGdt(1:Nl,1:11,i)-1.0_PR/dxm(i)*(F2(1:Nl,1:11)&
                                                       +vx_s*H2(1:Nl,1:11)&
                                                       +vy_s*K2(1:Nl,1:11)&
                                                       +vz_s*M2(1:Nl,1:11))
      !!!   terme source i+1 :

      call get_HKM(Nl,GR,H2,K2,M2)

      dGdt(1:Nl,1:11,ip1)=dGdt(1:Nl,1:11,ip1)+1.0_PR/dxm(i)*(F2(1:Nl,1:11)&
                                                                      +vx_s*H2(1:Nl,1:11)&
                                                                      +vy_s*K2(1:Nl,1:11)&
                                                                      +vz_s*M2(1:Nl,1:11))
 

  enddo

   !!!---CL:
   dUdt(:,1)=0.0_PR
   dUdt(:,N)=0.0_PR
   dGdt(:,:,1)=0.0_PR
   dGdt(:,:,N)=0.0_PR

   call cpu_time(t1) 

   Teuler=Teuler+t1-t10

end subroutine src_euler1D

subroutine get_HKM(Nl,G,H2,K2,M2)

   implicit none

   integer, intent(in) :: Nl
   real(PR), intent(in) :: G(1:Nl,1:11)
   real(PR), intent(out) :: H2(1:Nl,1:11)
   real(PR), intent(out) :: K2(1:Nl,1:11)
   real(PR), intent(out) :: M2(1:Nl,1:11)
   real(PR) :: fl(1:Nl), al(1:Nl,1:3), bl(1:Nl,1:3), cl(1:Nl,1:3), pl(1:Nl)
   integer :: iph, i2

   fl(1:Nl)=G(1:Nl,1)
   al(1:Nl,1:3)=G(1:Nl,2:4)  
   bl(1:Nl,1:3)=G(1:Nl,5:7)  
   cl(1:Nl,1:3)=G(1:Nl,8:10) 
   do iph=1,Nl
      pl(iph)=(g_sge(iph)-1.0_PR)*G(iph,11)/fl(iph)-g_sge(iph)*p_sge(iph)
   enddo

   H2(1:Nl,1)=-1.0_PR*fl(1:Nl)
   H2(1:Nl,2:4)=0.0_PR
   do i2=1,3
    H2(1:Nl,4+i2)=-bl(1:Nl,i2)        
   enddo
   do i2=1,3
    H2(1:Nl,7+i2)=-cl(1:Nl,i2)        
   enddo
   H2(1:Nl,11)=fl(1:Nl)*pl(1:Nl)

   K2(1:Nl,1)=0.0_PR
   do i2=1,3
    K2(1:Nl,1+i2)=bl(1:Nl,i2)
   enddo
   K2(1:Nl,5:11)=0.0_PR

   M2(1:Nl,1)=0.0_PR
   do i2=1,3
    M2(1:Nl,1+i2)=cl(1:Nl,i2)
   enddo
   M2(1:Nl,5:11)=0.0_PR

end subroutine get_HKM

!!!========================= CALCUL DU PAS DE TEMPS =============================

function dtcfl(nx,ny,nz,dx,c,vx,vy,vz,cfl)

   implicit none
   real(PR) :: dtcfl
   integer, intent(in) :: nx,ny,nz
   real(PR), intent(in) :: cfl
   real(PR), dimension(1:nx) :: dx
   real(PR), dimension(1:nx,1:ny,1:nz), intent(in) :: c,vx,vy,vz
   real(PR) :: vmax
   integer :: i,j,k

   vmax=0.0_PR
   dtcfl=1.0e30_PR

   do k=1,nz
    do j=1,ny
     do i=1,nx
         vmax=max(abs(vx(i,j,k)),abs(vx(i,j,k)+c(i,j,k)),abs(vx(i,j,k)-c(i,j,k)))
         dtcfl=min(dtcfl,dx(i)/vmax)
         !dtcfl=min(dtcfl,vx,vx+c,vx-c)
         !dtcfl=min(dtcfl,vx,vx+c,vx-c)
     enddo
    enddo
   enddo

   dtcfl=cfl*dtcfl
   
   if(dtcfl.lt.1.0e-12_PR)then
       print*, 'PROBLEM CFL:', dtcfl
       stop
   endif

end function dtcfl

!!!===================== ROUTINES D'INTEGRATION TEMPORELLE =====================

subroutine RK1(Nl,N,dxm,U,G,dt)

   implicit none
   integer, intent(in) :: Nl, N
   real(PR), intent(in) :: dxm(1:N)
   real(PR), intent(in) :: dt
   real(PR), intent(inout) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
   real(PR) :: dUdt(1:Nl+4,1:N), dGdt(1:Nl,1:11,1:N)

       call SRC(Nl,N,U(1:Nl+4,1:N),G(1:Nl,1:11,1:N),dxm(1:N),&
                     dUdt(1:Nl+4,1:N),dGdt(1:Nl,1:11,1:N))
     
       U=U+dt*dUdt
       G=G+dt*dGdt

       call relax(Nl,N,U,G)

end subroutine RK1

subroutine RK2(Nl,N,dxm,U,G,dt)

   implicit none
   integer, intent(in) :: Nl, N
   real(PR), intent(in) :: dxm(1:N)
   real(PR), intent(in) :: dt
   real(PR), intent(inout) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
   real(PR) :: U2(1:Nl+4,1:N), G2(1:Nl,1:11,1:N)
   real(PR) :: dUdt(1:Nl+4,1:N), dGdt(1:Nl,1:11,1:N)

       call SRC(Nl,N,U,G,dxm,dUdt,dGdt)

       U2=U+0.5_PR*dt*dUdt
       G2=G+0.5_PR*dt*dGdt

       call relax(Nl,N,U2,G2)

       call SRC(Nl,N,U2,G2,dxm,dUdt,dGdt)
      
       U=U+dt*dUdt
       G=G+dt*dGdt

       call relax(Nl,N,U,G)

end subroutine RK2

subroutine RK3(Nl,N,dxm,U,G,dt)

   implicit none
   integer, intent(in) :: Nl, N
   real(PR), intent(in) :: dxm(1:N)
   real(PR), intent(in) :: dt
   real(PR), intent(inout) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
   real(PR) :: U1(1:Nl+4,1:N), G1(1:Nl,1:11,1:N)
   real(PR) :: U2(1:Nl+4,1:N), G2(1:Nl,1:11,1:N)
   real(PR) :: dUdt(1:Nl+4,1:N), dGdt(1:Nl,1:11,1:N)

       call SRC(Nl,N,U,G,dxm,dUdt,dGdt)

       U1=U+dt*dUdt
       G1=G+dt*dGdt

       call relax(Nl,N,U1,G1)

       call SRC(Nl,N,U1,G1,dxm,dUdt,dGdt)
      
       U2=0.75_PR*U+0.25_PR*(U1+dt*dUdt)
       G2=0.75_PR*G+0.25_PR*(G1+dt*dGdt)

       call relax(Nl,N,U2,G2)

       call SRC(Nl,N,U2,G2,dxm,dUdt,dGdt)
 
       U=1.0_PR/3.0_PR*(U+2.0_PR*(U2+dt*dUdt))
       G=1.0_PR/3.0_PR*(G+2.0_PR*(G2+dt*dGdt))
     
       call relax(Nl,N,U,G)

end subroutine RK3

!SUBROUTINE RK4(ww,nx,dx,dt)
!implicit none
!INTEGER,INTENT(IN)                      ::  nx
!real(dp), INTENT(IN)        ::  dx,dt
!real(dp),dimension(1:nx,3),intent(inout) ::  ww
!real(dp),dimension(1:nx,3) ::  rw_L,rw_R,rw_C1,rw_C2,w_C,w_R
!
!call Scheme_2D(ww(1:nx,1:3),rw_L(1:nx,1:3),nx,dx)
!
!w_C   = ww + 0.5d0*dt*rw_L
!call Scheme_2D(w_C(1:nx,1:3),rw_C1(1:nx,1:3),nx,dx)
!
!w_C  = ww + 0.5d0*dt*rw_C1
!call Scheme_2D(w_C(1:nx,1:3),rw_C2(1:nx,1:3),nx,dx)
!
!w_R   = ww + dt*rw_C2
!call Scheme_2D(w_R(1:nx,1:3),rw_R(1:nx,1:3),nx,dx)
!
!ww   =  ww + dt*(rw_L+rw_R+2.d0*(rw_C1+rw_C2))/6.d0
!
!END SUBROUTINE RK4
!

!!!===================== ROUTINES MATHEMATIQUES =====================

subroutine My_Newton(N,x,xmin,xmax,F,DF)
   !!! Newton une seule variable
   implicit none
   integer, intent(in) :: N
   real(PR), intent(inout) :: x(1:N)
   real(PR), intent(in) :: xmin,xmax 
   procedure(F_test1) :: F !!! fonction
   procedure(J_test1) :: DF !!! jacobienne
   real(PR) :: tol, fx(1:1), dfdx(1:1), dx
   integer :: it    

   call F(1,x,fx)
   tol=abs(fx(1))
   it=0

   do while(tol.gt.1.0e-9_PR.and.it.lt.10000)
 
      call F(1,x(1),fx(1:1))
      call DF(1,x(1),dfdx(1:1))

      dx=-fx(1)/dfdx(1)

      x(1)=x(1)+dx
      x(1)=max(x(1),xmin)
      x(1)=min(x(1),xmax)

      it=it+1

      tol=abs(fx(1))

      !if(it.ge.9900)then
      ! print*, 'it=', it, x(1), fx(1), dx, dfdx(1)
      !endif

   enddo

   if(it.ge.10000)then
     print*, 'PROBLEM NEWTON'
     print*, 'fl0:', fl0(:)
     print*, 'flrhl0(:)', flrhl0(:)
     print*, 'rhl0(:)', flrhl0(:)/fl0(:)
     print*, 'el0(:)', el0(:)

     print*, 'pl0:', pl0(:)
     print*, 'xmin=',xmin
     print*, 'xmax=',xmax
     stop 
   endif

end subroutine My_newton

!!!===================== ROUTINES DE CRASH =====================

subroutine crash(message)

   implicit none
   character(len=*), intent(in) :: message
 
   write(iout,*) '-----  PROBLEM solve_euler :c --------'
   write(iout,*) message
   write(iout,*) '-----------------------------------'

   !stop

   !call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)

end subroutine crash




end module mod_euler
