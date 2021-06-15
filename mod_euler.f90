module mod_euler


use mod_data, only : it, PR, iout, mesh, multifluide, fluide, radial, Input, initT

use mod_input, only : init_fields_from_input

use air_ETL, only : specific_heat 

use mod_mat_tab, only : mat_tab, mat, read_table,&
                       get_index_rT,P_from_rT, E_from_rT,&
                       Cv_from_rT, c_from_rT, g_from_rT,&
                       pinf_from_rT, Z_from_rT,&
                       get_index_rE,T_from_ru,& 
                       get_index_rP, T_from_rP,&
                       get_index_TP, r_from_TP

!!!, Nt, Nl, Nx, Ny, Nz, t, dt

!use mod_snes, only : init_SNES, solve_SNES, restart_snes !SNESinterface

implicit none

integer :: i_sortie=-1
logical :: verb=.false. 
logical :: just_thermo=.false.
logical :: plneg=.false.
logical :: phasechange=.false.
logical :: Reconstruct=.false.
logical :: verbnewton=.false.
integer :: modeG=1 !!! 1: conv de elh 2: convection de elh+ele

real(PR) :: dtmin=1.0e10_PR
real(PR) :: rtemp(1:10)

integer :: ierr=0

character(len=10) :: tscheme='RK1'

!!!---taille du maillage :

!!!---routines
procedure(relax_p0), pointer :: relax_p => NULL() 
procedure(NOREC), pointer :: RECO => NULL()
procedure(MINMOD), pointer :: Limiteur => NULL()
procedure(src_euler1D), pointer :: SRC => NULL()
procedure(RK1), pointer :: Timescheme => NULL()
procedure(soundspeed_mixt_frozen), pointer :: soundspeed_mixt => NULL()
procedure(pressure_mixt_sge), pointer :: pressure_mixt => NULL()

type :: phase
   integer :: typ=1 !!! 1: sge, 2: tab
   integer ::imat=0 !!! si typ=3
   character(len=50) :: filetab=''
   character(len=50) :: nom=''
   real(PR) :: p_sge, g_sge, mu, sigy, rho0, cv, sig, ener0, q, qp
   real(PR) :: Tvap=1.0e9_PR
   procedure(soundspeed_sge), nopass, pointer :: soundspeed
   procedure(pressure_sge), nopass, pointer :: pressure
   procedure(energie_hydro_sge), nopass, pointer :: energie_hydro
   procedure(pressure_star_sge), nopass, pointer :: pressure_star
   procedure(T_from_rP_sge), nopass, pointer :: T_from_rP
   procedure(rh_from_PT_sgeT), nopass, pointer :: rh_from_PT
end type phase

type(phase), allocatable :: ph(:)

!!!---taille du maillage
integer :: N1D,Nx,Ny,Nz

!!!---nombre de phases
integer :: Nl


!!!---temps pour calcul des perfos
real(PR) :: t1,t2,t3,t10,t20,t30
real(PR) :: t_euler=0.0_PR
real(PR) :: t_EM=0.0_PR
real(PR) :: t_src=0.0_PR
real(PR) :: t_init_Finger=0.0_PR
real(PR) :: t_relaxp=0.0_PR
real(PR) :: t_reset=0.0_PR
real(PR) :: t_reco=0.0_PR
real(PR) :: t_p2c=0.0_PR
real(PR) :: t_c2p=0.0_PR
real(PR) :: t_nc2p=0.0_PR
real(PR) :: t_tot=0.0_PR


!!!---nombre de variables conservatives: 4+Nl
integer :: Neq

!!!---tableaux de travail 1D

!!!---grandeurs 1D geometriques:
real(PR), allocatable :: SVL(:)
real(PR), allocatable :: SVR(:)
real(PR), allocatable :: Sgeo(:)

!!!---vecteur 1D primitif:
type(multifluide), allocatable :: MF1D(:)

!!!---vecteur 1D conservatif:
real(PR), allocatable :: U(:,:), dUdt(:,:) !!! 1:Nl+4 , 1:Nx

!!!---vecteur 1D non-conservatif :
!!!variables non conservatives : Nl*11
real(PR), allocatable :: G(:,:,:), dGdt(:,:,:) !!! Nl, 11, Nx

!!!---etats intermédiaires RK :
real(PR), allocatable :: U0(:,:), U1(:,:), U2(:,:) !!! 1:Nl+4 , 1:Nx
real(PR), allocatable :: G0(:,:,:), G1(:,:,:), G2(:,:,:) !!! Nl, 11, Nx

!!!---Vecteurs d'etats primitifs pour le schéma COMPACT
real(PR), allocatable :: AdW(:), WL(:), WR(:), Wi(:), Wim1(:), Wip1(:), Deltai(:)

!!!---variables pour relaxation
type(multifluide) :: MF0, MF01
real(PR) :: PNewton
integer :: err_relax=0
integer :: Nit_Newton=0

!!!---variable pour la meca : FT: Finger tensor FT2: FT^2
real(PR), dimension(1:3,1:3) :: Id, FT, FT2
real(PR) :: J1, J2, murr0, mu, alpham1, alpharhom1, rho
real(PR) :: detFT, detFT12, detFT13, detFT23, detFTm16, detFTm13, detFTm23, detFTm43, detFTm53
real(PR) :: a(1:3), b(1:3), c(1:3)
real(PR) :: aa, bb, cc, ab, bc, ac, ab2, ac2, bc2
real(PR) :: r13=1.0_PR/3.0_PR
real(PR) :: r23=2.0_PR/3.0_PR
real(PR) :: r43=4.0_PR/3.0_PR
logical  :: a2A=.false.  !!!! dans p2c, c2p, nonc2p, get_HKM et HLLC :
                        !!!  al=Al=al*fl**1/3 pour le transport !

!!!---tableaux de travail 1D pour la reconstruction
!real(PR), allocatable :: Urec_L(:,:), Urec_R(:,:)
!real(PR), allocatable :: Grec_L(:,:,:),Grec_R(:,:,:) !!! Nl, 11, Nx
type(multifluide), allocatable :: MFrec_L(:), MFrec_R(:) 

!!!---tableau et variables pour la relaxation de la pression
integer :: Nl0
real(PR), allocatable :: g_sge0(:), p_sge0(:), flrhl0(:), el0(:), nul0(:), pl0(:), Zl(:), fl0(:), fl0_gl(:), coef1(:), rhl0(:)
real(PR), allocatable :: xin(:), bin(:)
integer :: Nsnes=0
real(PR) :: pI0
integer :: itrelax=0
integer :: mode_relax=1 !!! 1: equation SGE (Ndanou 2015 eq 25 p 534)
                        !!! 2: système général tout EOS (Saurel 2009 p 1691)
integer :: mode_pI=1 !!! 1: pI=p ; 2: pI=sum(Zlpl)/sum(Zl) (cf Saurel 2009 p1691)
logical :: compute_c=.true. !!! desactive le calcul de c à l'entrée de relax_p

!!!--- Tolérances:
real(PR) :: toldet=1.0e-20_PR
real(PR) :: tolmu=1.0e-20_PR
real(PR) :: tolmuscl=1.0e-20_PR

!!!---Propriétés matière

!real(PR), allocatable :: mu(:)
!real(PR), allocatable :: sigy(:)
!real(PR), allocatable :: rho0(:)

!!!---EOS Stiffened gas Equation

!!!---eau:
!real(PR) :: g_sge_H2O=6.1_PR
!real(PR) :: p_sge_H2O=2.0e9_PR
!real(PR) :: mu_H2O=0.0_PR
!real(PR) :: sigy_H2O=0.0_PR
!real(PR) :: rho0_H2O=1000.0_PR

!!!---eau 1: Saurel 1999
real(PR) :: g_sge_eau=4.4_PR
real(PR) :: p_sge_eau=6.0e8_PR
real(PR) :: mu_eau=0.0_PR
real(PR) :: sigy_eau=0.0_PR
real(PR) :: rho0_eau=1000.0_PR
real(PR) :: cv_eau=4185.0_PR
real(PR) :: sig_eau=1.0e2_PR

!!!---eau 2: Pepitas 2009
real(PR) :: g_sge_eau2=2.35_PR
real(PR) :: p_sge_eau2=1.0e9_PR
real(PR) :: mu_eau2=0.0_PR
real(PR) :: sigy_eau2=0.0_PR
real(PR) :: rho0_eau2=1000.0_PR
real(PR) :: cv_eau2=1816.0_PR
real(PR) :: sig_eau2=1.0e2_PR
real(PR) :: q_eau2=-1167.0_PR
real(PR) :: qp_eau2=-23000.0_PR

!!!---air:
real(PR) :: g_sge_air=1.4_PR
real(PR) :: p_sge_air=0.0_PR
real(PR) :: mu_air=0.0_PR
real(PR) :: sigy_air=0.0_PR
real(PR) :: rho0_air=1.0_PR
real(PR) :: cv_air=786.0_PR !!! J/kg/K
real(PR) :: sig_air=1.0e-20_PR

!!!---Fe:
real(PR) :: g_sge_Fe=3.9_PR
real(PR) :: p_sge_Fe=43.6e9_PR
real(PR) :: mu_Fe=82.0e9_PR
real(PR) :: sigy_Fe=200.0e6_PR
real(PR) :: rho0_Fe=7860.0_PR
real(PR) :: cv_Fe=440.0_PR !!! J/kg/K
real(PR) :: sig_Fe=9.93e6_PR

!!!---Al:
real(PR) :: g_sge_Al=3.4_PR
real(PR) :: p_sge_Al=21.5e9_PR
real(PR) :: rho0_Al=2700.0_PR
!real(PR) :: g_sge_Al=3.5_PR
!real(PR) :: p_sge_Al=32.0e9_PR
!real(PR) :: rho0_Al=2712.0_PR
real(PR) :: mu_Al=26.0e9_PR
real(PR) :: sigy_Al=60.0e6_PR
real(PR) :: cv_Al=897.0_PR !!! J/kg/K
real(PR) :: sig_Al=37.7e6_PR

!!!---Ti:
real(PR) :: g_sge_Ti=2.6_PR
real(PR) :: p_sge_Ti=44.0e9_PR
real(PR) :: mu_Ti=42.0e9_PR
real(PR) :: sigy_Ti=1030.0e6_PR
real(PR) :: rho0_Ti=4527.0_PR
real(PR) :: cv_Ti=520.0_PR !!! J/kg/K
real(PR) :: sig_Ti=2.34e6_PR

!!!---Cu:
real(PR) :: g_sge_Cu=4.22_PR    !50.0_PR            !!! 4.22_PR
real(PR) :: p_sge_Cu=34.2e9_PR  !1.243361795e0_PR   !!! 34.2e9_PR
real(PR) :: mu_Cu=9.2e10_PR
real(PR) :: sigy_Cu=40.0e6_PR
real(PR) :: rho0_Cu=8900.0_PR
real(PR) :: cv_Cu=380.0_PR !!! J/kg/K
real(PR) :: sig_Cu=59.6e6_PR
real(PR) :: Tvap_Cu=2835.0_PR

!!!---Cug:
real(PR) :: g_sge_Cug=4.22_PR
real(PR) :: p_sge_Cug=0.0_PR
real(PR) :: mu_Cug=0.0_PR
real(PR) :: sigy_Cug=0.0_PR
real(PR) :: rho0_Cug=0.001_PR*8900.0_PR
real(PR) :: cv_Cug=380.0_PR !!! J/kg/K
real(PR) :: sig_Cug=1.0e-10_PR

!!!---epoxy:
real(PR) :: g_sge_epoxy=2.43_PR
real(PR) :: p_sge_epoxy=5.3e9_PR
real(PR) :: mu_epoxy=0.0_PR
real(PR) :: sigy_epoxy=0.0_PR
real(PR) :: rho0_epoxy=1185.0_PR
real(PR) :: cv_epoxy=380.0_PR !!! J/kg/K
real(PR) :: sig_epoxy=1.0e-10_PR

!!!---spinel:
real(PR) :: g_sge_spinel=1.62_PR
real(PR) :: p_sge_spinel=141.0e9_PR
real(PR) :: mu_spinel=0.0_PR
real(PR) :: sigy_spinel=0.0_PR
real(PR) :: rho0_spinel=3622.0_PR
real(PR) :: cv_spinel=380.0_PR !!! J/kg/K
real(PR) :: sig_spinel=1.0e-10_PR

!real(PR), allocatable :: p_sge(:), g_sge(:) 

contains

!!!===================== ROUTINE D'APPEL =====================

subroutine solve_euler(M,dt,cfl) !!!!(M,dt,cfl)

   !use mod_data, only : Nl,Nx,Ny,Nz,fl,Yl,pl,rhl,elh,ele,&
   !                     al,bl,cl,sigl,sig,&
   !                     rh,vx,vy,vz,c,p,e,&
   !                     dxm,dym,dzm
 
   !use mod_data, only : dxm,dym,dzm
                        

   implicit none
   type(mesh), intent(inout) :: M
   real(PR), intent(in) :: dt, cfl
   integer :: idir,i,j,k,iph,i2,im1,ip1,ip2
   real(PR) :: tloc, dtloc

   tloc=0.0_PR

   t_c2p=0.0_PR ; t_src=0.0_PR ; t_nc2p=0.0_PR ; t_relaxp=0.0_PR
   t_reset=0.0_PR ; t_init_Finger=0.0_PR

   do while(tloc.lt.dt)

      dtloc=dtcfl(M,cfl)

      dtloc=min(dt-tloc,dtloc,M%dtmin) 
      !print*, 'dtloc=', dtloc, 'dt=', dt 

      !!!-------------------------  X  ----------------------------
 
      idir=1
      N1D=Nx 
      MF1D(1:Nx)=M%MF(1:Nx,1,1)     
      SVL(0:Nx+1)=M%SVLx(0:Nx+1)
      SVR(0:Nx+1)=M%SVRx(0:Nx+1)
      Sgeo(0:Nx+1)=M%Sgeox(0:Nx+1)   

      DO k=1,Nz ; DO j=1,Ny

      !print*, 'PP1:', MF1D(101)%p,MF1D(101)%pl(1),MF1D(101)%pl(2) 


      call Timescheme(dtloc)
 

      !if(.not.plneg)then

         M%MF(1:Nx,j,k)=MF1D(1:Nx) 

         tloc=tloc+dtloc  

      !else

      !   dtloc=dtloc/2.0_PR
      !   print*, 'dtloc/2 !!!'
      !   stop   
 
      !   if(dtloc.eq.M%dtmin)then
      !    dtloc=M%dtmin
      !    print*, 'dtloc/2=', dtmin
      !    stop
      !   endif

      !endif
 
      ENDDO ; ENDDO


      !!!-------------------------  Y  ----------------------------
      !!!-------------------------  Z  ----------------------------

   enddo !!! dowhile

   t_tot=t_c2p+t_src+t_relaxp+t_nc2p+t_reset

         
   !print*, 'c2p:',t_c2p/t_tot, 'src:', t_src/t_tot, 'nc2p:', t_nc2p/t_tot, &
   !             'relax:',t_relaxp/t_tot,'reset:',t_reset/t_tot
!!!'Fing:', t_init_Finger/t_tot 
 
end subroutine solve_euler


!!!======================== SUBROUTINE D'INITIALISATION =============================

subroutine init_euler_from_input(M) 

    implicit none
    type(mesh), intent(inout) :: M
    integer :: i,j,k,iph

    write(iout,*) 'call init START !!!'

    Nx=Input%Nx
    Ny=Input%Ny
    Nz=Input%Nz
    Nl=Input%Nl

    M%Nx=Nx
    M%Ny=Ny
    M%Nz=Nz
    M%Nl=Nl
    M%Lx=Input%Lx
    M%Ly=Input%Ly
    M%Lz=Input%Lz
    M%Nt=Input%Noutput
    M%dt=Input%dtoutput

    !!!---MAILLAGE :

    call compute_mesh(M)

    !!!---ALLOCATION DES TABLEAUX
    write(iout,*) ""
    write(iout,fmt='(A)') "ALLOCATION DES TABLEAUX"
    write(iout,*) ""

    !!!---allocation de l'objet Fluide:

    allocate(M%MF(1:Nx,1:Ny,1:Nz))
  
    !!!---initialization de l'objet fluide
    DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
     call init_multifluide(M%MF(i,j,k)) 
    ENDDO ; ENDDO ; ENDDO

    !!!---nombre d'équations à résoudre partie conservative
    Neq=Nl+4 

    !!!!---vecteur 1D geometrique
    allocate(SVL(0:Nx+1))  ; SVL=0.0_PR
    allocate(SVR(0:Nx+1))  ; SVR=0.0_PR
    allocate(Sgeo(0:Nx+1)) ; Sgeo=0.0_PR

    !!!---vecteur 1D pour spltting directionnel:
    allocate(MF1D(0:Nx+1))
    do i=0,Nx+1
       call init_multifluide(MF1D(i)) 
    enddo

    !!!!---vecteur 1D des variavles conservatives :
    allocate(U(1:Neq,0:nx+1))    ; U(:,:)=0.0_PR
    allocate(dUdt(1:Neq,0:nx+1)) ; dUdt(:,:)=0.0_PR
    !!!!---vecteur 1D non-conservatif : fraction volumique + Finger tensor +
    !!                               énergie hydro de chaque phase   
    allocate(G(1:Nl,1:11,0:nx+1))    ; G(:,:,:)=0.0_PR
    allocate(dGdt(1:Nl,1:11,0:nx+1)) ; dGdt(:,:,:)=0.0_PR

    !!!---tableaux 1D pour la reconstruction

    if(trim(adjustl(Input%reco)).eq.'NOREC')then
         Reconstruct=.false.
    elseif(trim(adjustl(Input%reco)).eq.'MUSCL')then
         Reconstruct=.true.
    endif

    !!!---Limiteur:
    if(trim(adjustl(Input%reco)).eq.'MINMOD')then
       Limiteur => MINMOD
    elseif(trim(adjustl(Input%reco)).eq.'VANLEER')then
       Limiteur => VANLEER
    endif


    if(Reconstruct)then
       allocate(MFrec_L(0:Nx+1))
       allocate(MFrec_R(0:Nx+1))
       DO i=1,Nx
        call init_multifluide(MFrec_L(i)) 
        call init_multifluide(MFrec_R(i)) 
       ENDDO
    endif

    !!!---initialisation du solveur SNES pour la relaxation de la pression

    call init_relaxp

    !!!---vecteur 1D des variables primitives
    !Nprim=9+Nl*13

    !write(iout,*) '  Nprim=', Nprim  
    !allocate(W1(1:Nprim,1:Nx))
    !allocate(W_L(1:Nprim,1:Nx))
    !allocate(W_R(1:Nprim,1:Nx))

    !!!================= AFFECTATION PHASE -> MATERIAUX =================

    allocate(ph(1:Nl))
    do iph=1,Nl
        if(Input%EOS.eq.1)then
            ph(iph)%typ=3
        else
            ph(iph)%typ=4
        endif
        ph(iph)%imat=Input%f2m(iph)
    enddo

    call init_phase(M)

    !!!---Soundspeed of the mixture:
    if(Input%soundspeed_mixt.eq.2)then
        soundspeed_mixt => soundspeed_mixt_Wood
    else
        soundspeed_mixt => soundspeed_mixt_frozen
    endif

    !!!---Pressure of the mixture:
    if(Input%EOS.eq.1)then
        pressure_mixt => pressure_mixt_sge
    else
        pressure_mixt => pressure_mixt_sgeT
    endif

    !!!---Reconstruction

    if(trim(adjustl(Input%reco)).eq.'NOREC')then
        Reconstruct=.false.
    elseif(trim(adjustl(Input%reco)).eq.'MUSCL')then
        Reconstruct=.true.
    endif

    !!!---Limiteur:
    if(trim(adjustl(Input%reco)).eq.'MINMOD')then
        Limiteur => MINMOD
    elseif(trim(adjustl(Input%reco)).eq.'VANLEER')then
        Limiteur => VANLEER
    endif

    !!!---Timescheme:
    if(trim(adjustl(Input%tscheme)).eq.'RK1')then

        Timescheme => RK1

    elseif(trim(adjustl(Input%tscheme)).eq.'RK2')then

        !!!---tableaux de travail RK2
        allocate(U0(1:Neq,0:nx+1))      ; U0(:,:)=0.0_PR
        allocate(G0(1:Nl,1:11,0:nx+1))  ; G0(:,:,:)=0.0_PR

        Timescheme => RK2 !!!_SSP

    elseif(trim(adjustl(Input%tscheme)).eq.'RK2_SSP')then

        !!!---tableaux de travail RK2
        allocate(U0(1:Neq,0:nx+1))      ; U0(:,:)=0.0_PR
        allocate(G0(1:Nl,1:11,0:nx+1))  ; G0(:,:,:)=0.0_PR

        Timescheme => RK2_SSP

    elseif(trim(adjustl(Input%tscheme)).eq.'RK3')then

        !!!---tableaux de travail RK3
        allocate(U0(1:Neq,0:nx+1))      ; U0(:,:)=0.0_PR
        allocate(G0(1:Nl,1:11,0:nx+1))  ; G0(:,:,:)=0.0_PR
        allocate(U1(1:Neq,0:nx+1))      ; U1(:,:)=0.0_PR
        allocate(G1(1:Nl,1:11,0:nx+1))  ; G1(:,:,:)=0.0_PR

        Timescheme => RK3

    elseif(trim(adjustl(Input%tscheme)).eq.'RK4')then

        !!!---tableaux de travail RK3
        allocate(U0(1:Neq,0:nx+1))      ; U0(:,:)=0.0_PR
        allocate(G0(1:Nl,1:11,0:nx+1))  ; G0(:,:,:)=0.0_PR
        allocate(U1(1:Neq,0:nx+1))      ; U1(:,:)=0.0_PR
        allocate(G1(1:Nl,1:11,0:nx+1))  ; G1(:,:,:)=0.0_PR
        allocate(U2(1:Neq,0:nx+1))      ; U2(:,:)=0.0_PR
        allocate(G2(1:Nl,1:11,0:nx+1))  ; G2(:,:,:)=0.0_PR

        Timescheme => RK4

    elseif(trim(adjustl(Input%tscheme)).eq.'COMPACT')then

        allocate(WL(1:12*Nl+3))   ; WL=0.0_PR
        allocate(WR(1:12*Nl+3))   ; WR=0.0_PR
        allocate(Wi(1:12*Nl+3))   ; Wi=0.0_PR
        allocate(Wip1(1:12*Nl+3)) ; Wip1=0.0_PR
        allocate(Wim1(1:12*Nl+3)) ; Wim1=0.0_PR
        allocate(Deltai(1:12*Nl+3)) ; Deltai=0.0_PR
        allocate(AdW(1:12*Nl+3))  ; AdW=0.0_PR

        Timescheme => COMPACT
 
    else

        print*, ' ö please select a time integrator' ; stop

    endif

    !!!---SRC:

    if(Input%just_thermo.eq.1)then
        SRC => src_EM
    else
        SRC => src_euler1D
    endif

    !!!--- Initialisation des champs:

    write(iout,*) '   initialisation de fv,p,v,r' 

    call init_fields_from_input(M)

    if(Input%isoT.eq.1)then
        if(initT)then
           forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%rh=ph(iph)%rh_from_PT(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%p,M%MF(i,j,k)%F(iph)%T)
        else
           write(iout,*) "      > ERREUR init_euler_from_input: si isoT=1, specifier la temperature dans le input.dat !!!"
           stop
        endif
    endif

    forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

    forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

    forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)

    write(iout,*) ' MAJ meca'

    do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

     call MAJ_meca(M%MF(i,j,k))

    enddo ; enddo ; enddo

    write(iout,*) ' MAJ mixt'

    do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

    call MAJ_mixt(M%MF(i,j,k))

    enddo ; enddo ; enddo

    write(iout,*) 'call init from input FINI !!!'

end subroutine init_euler_from_input

subroutine compute_mesh(M)

    implicit none
    type(mesh), intent(inout) :: M
    real(PR) :: Lx,Ly,Lz
    integer :: i

    write(iout,*) ""
    write(iout,fmt='(A)') " call to compute mesh..."
    write(iout,*) ""
    Lx=M%Lx ; Ly=M%Ly ; Lz=M%Lz
    write(iout,*) "Nx=",Nx, "Ny=", Ny, "Nz=",Nz 
    write(iout,*) "Lx=",Lx, "Ly=", Ly, "Lz=",Lz 
    write(iout,*) ""
    write(iout,*) '   allocation des grandeurs géométriques'
    allocate(M%x(1:nx))
    allocate(M%xm(0:nx))
    allocate(M%dx(1:nx))
    allocate(M%dxm(0:nx))
    allocate(M%y(1:ny))
    allocate(M%ym(0:ny))
    allocate(M%dy(1:ny))
    allocate(M%dym(0:ny))
    allocate(M%z(1:nz))
    allocate(M%zm(0:nz))
    allocate(M%dz(1:nz))
    allocate(M%dzm(0:nz))
    allocate(M%SVLx(0:nx+1))
    allocate(M%SVRx(0:nx+1))
    allocate(M%Sgeox(0:nx+1))
    write(iout,*) '   valorisation'

    !!!---X:
    M%dx(:)=Lx/real(Nx,PR)
    M%xm(0)=0.0_PR
    M%xm(0)=2.0_PR*M%dx(1) 
    do i=1,nx
      M%xm(i)=M%xm(i-1)+M%dx(i)
    enddo
    do i=1,nx
      M%x(i)=0.5_PR*(M%xm(i-1)+M%xm(i))
    enddo
    M%dxm(0)=2.0_PR*M%x(1)
    do i=1,nx-1
      M%dxm(i)=M%x(i+1)-M%x(i)
    enddo
    M%dxm(nx)=2.0_PR*(M%xm(nx)-M%x(nx))
    if(radial)then
       do i=1,Nx
          M%SVLx(i)=M%xm(i)/(M%x(i)*M%dx(i))
          M%SVRx(i)=M%xm(i)/((M%x(i)+M%dx(i))*M%dx(min(i+1,Nx)))
          M%Sgeox(i)=1.0_PR/M%x(i)
       enddo
       M%Sgeox(1)=0.0_PR
       !M%SVLx(0)=0.0_PR ; M%SVRx(0)=0.0_PR
       M%SVLx(0)=M%SVLx(1) ; M%SVRx(0)=M%SVLx(1)
       M%SVLx(Nx+1)=(M%xm(Nx)+M%dxm(Nx))/((M%x(Nx)+M%dx(Nx))*M%dx(Nx))
       M%SVRx(Nx+1)=(M%xm(Nx)+M%dxm(Nx))/((M%x(Nx)+2.0_PR*M%dx(Nx))*M%dx(Nx))
       M%Sgeox(0)=0.0_PR ; M%Sgeox(Nx+1)=1.0_PR/(M%x(Nx)+M%dx(Nx))
    else
      do i=1,Nx
          M%SVLx(i)=1.0_PR/M%dx(i)
          M%SVRx(i)=1.0_PR/M%dx(min(i+1,Nx))
          M%Sgeox(i)=0.0_PR
      enddo
      M%SVLx(0)=M%SVLx(1) ; M%SVRx(0)=M%SVRx(1)
      M%SVLx(Nx+1)=M%SVLx(Nx)
      M%SVRx(Nx+1)=M%SVRx(Nx)
      M%Sgeox(0)=0.0_PR ; M%Sgeox(Nx+1)=0.0_PR
    endif 

    !!---Y:
    M%dy(:)=Ly/real(Ny,PR)
    M%ym(0)=0.0_PR
    do i=1,ny
      M%ym(i)=M%ym(i-1)+M%dy(i)
    enddo
    do i=1,ny
      M%y(i)=0.5_PR*(M%ym(i-1)+M%ym(i))
    enddo
    M%dym(0)=2.0_PR*M%y(1)
    do i=1,ny-1
      M%dym(i)=M%y(i+1)-M%y(i)
    enddo
    M%dym(ny)=2.0_PR*(M%ym(ny)-M%y(ny))
    
    !!!---Z:
    M%dz(:)=Lz/real(Nz,PR)
    M%zm(0)=0.0_PR
    do i=1,nz
      M%zm(i)=M%zm(i-1)+M%dz(i)
    enddo
    do i=1,nz
      M%z(i)=0.5_PR*(M%zm(i-1)+M%zm(i))
    enddo
    M%dzm(0)=2.0_PR*M%z(1)
    do i=1,nz-1
      M%dzm(i)=M%z(i+1)-M%z(i)
    enddo
    M%dzm(nz)=2.0_PR*(M%zm(nz)-M%z(nz))

end subroutine compute_mesh


subroutine init_euler(nx_i,ny_i,nz_i,nl_i,Lx,Ly,Lz,M) 

   !use mod_data, only : Nl,Nx,Ny,Nz,&
   !                     fl,Yl,pl,rhl,elh,ele,&
   !                     al,bl,cl,sigl,sig,&
   !                     rh,vx,vy,vz,c,p,e,&
   !                     Lx,Ly,Lz,&
   !                     dxm,dym,dzm,dx,dy,dz,&
   !                     xm,ym,zm,x,y,z

   implicit none
   integer, intent(in) :: nx_i, ny_i, nz_i, Nl_i
   real(PR), intent(in) :: Lx,Ly,Lz
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph

   write(iout,*) 'call init START !!!'

   Nx=Nx_i
   Ny=Ny_i
   Nz=Nz_i
   Nl=Nl_i

   M%Nx=nx
   M%Ny=ny
   M%Nz=Nz
   M%Nl=Nl
   M%Lx=Lx
   M%Ly=Ly
   M%Lz=Lz

   !!!---MAILLAGE :

   write(iout,*) ""
   write(iout,fmt='(A)') "CREATION DU MAILLAGE"
   write(iout,*) ""
   write(iout,*) "Nx=",Nx, "Ny=", Ny, "Nz=",Nz 
   write(iout,*) "Lx=",Lx, "Ly=", Ly, "Lz=",Lz 
   write(iout,*) ""
   write(iout,*) '   allocation des grandeurs géométriques'
   allocate(M%x(1:nx))
   allocate(M%xm(0:nx))
   allocate(M%dx(1:nx))
   allocate(M%dxm(0:nx))
   allocate(M%y(1:ny))
   allocate(M%ym(0:ny))
   allocate(M%dy(1:ny))
   allocate(M%dym(0:ny))
   allocate(M%z(1:nz))
   allocate(M%zm(0:nz))
   allocate(M%dz(1:nz))
   allocate(M%dzm(0:nz))
   allocate(M%SVLx(0:nx+1))
   allocate(M%SVRx(0:nx+1))
   allocate(M%Sgeox(0:nx+1))
   write(iout,*) '   valorisation'

   !!!---X:
   M%dx(:)=Lx/real(Nx,PR)
   M%xm(0)=0.0_PR
   M%xm(0)=2.0_PR*M%dx(1) 
   do i=1,nx
     M%xm(i)=M%xm(i-1)+M%dx(i)
   enddo
   do i=1,nx
     M%x(i)=0.5_PR*(M%xm(i-1)+M%xm(i))
   enddo
   M%dxm(0)=2.0_PR*M%x(1)
   do i=1,nx-1
     M%dxm(i)=M%x(i+1)-M%x(i)
   enddo
   M%dxm(nx)=2.0_PR*(M%xm(nx)-M%x(nx))
   if(radial)then
      do i=1,Nx
         M%SVLx(i)=M%xm(i)/(M%x(i)*M%dx(i))
         M%SVRx(i)=M%xm(i)/((M%x(i)+M%dx(i))*M%dx(min(i+1,Nx)))
         M%Sgeox(i)=1.0_PR/M%x(i)
      enddo
      M%Sgeox(1)=0.0_PR
      !M%SVLx(0)=0.0_PR ; M%SVRx(0)=0.0_PR
      M%SVLx(0)=M%SVLx(1) ; M%SVRx(0)=M%SVLx(1)
      M%SVLx(Nx+1)=(M%xm(Nx)+M%dxm(Nx))/((M%x(Nx)+M%dx(Nx))*M%dx(Nx))
      M%SVRx(Nx+1)=(M%xm(Nx)+M%dxm(Nx))/((M%x(Nx)+2.0_PR*M%dx(Nx))*M%dx(Nx))
      M%Sgeox(0)=0.0_PR ; M%Sgeox(Nx+1)=1.0_PR/(M%x(Nx)+M%dx(Nx))
   else
     do i=1,Nx
         M%SVLx(i)=1.0_PR/M%dx(i)
         M%SVRx(i)=1.0_PR/M%dx(min(i+1,Nx))
         M%Sgeox(i)=0.0_PR
     enddo
     M%SVLx(0)=M%SVLx(1) ; M%SVRx(0)=M%SVRx(1)
     M%SVLx(Nx+1)=M%SVLx(Nx)
     M%SVRx(Nx+1)=M%SVRx(Nx)
     M%Sgeox(0)=0.0_PR ; M%Sgeox(Nx+1)=0.0_PR
   endif 

  !!!---Y:
   M%dy(:)=Ly/real(Ny,PR)
   M%ym(0)=0.0_PR
   do i=1,ny
     M%ym(i)=M%ym(i-1)+M%dy(i)
   enddo
   do i=1,ny
     M%y(i)=0.5_PR*(M%ym(i-1)+M%ym(i))
   enddo
   M%dym(0)=2.0_PR*M%y(1)
   do i=1,ny-1
     M%dym(i)=M%y(i+1)-M%y(i)
   enddo
   M%dym(ny)=2.0_PR*(M%ym(ny)-M%y(ny))

   !!!---Z:
   M%dz(:)=Lz/real(Nz,PR)
   M%zm(0)=0.0_PR
   do i=1,nz
     M%zm(i)=M%zm(i-1)+M%dz(i)
   enddo
   do i=1,nz
     M%z(i)=0.5_PR*(M%zm(i-1)+M%zm(i))
   enddo
   M%dzm(0)=2.0_PR*M%z(1)
   do i=1,nz-1
     M%dzm(i)=M%z(i+1)-M%z(i)
   enddo
   M%dzm(nz)=2.0_PR*(M%zm(nz)-M%z(nz))


!!!---ALLOCATION DES TABLEAUX
   write(iout,*) ""
   write(iout,fmt='(A)') "ALLOCATION DES TABLEAUX"
   write(iout,*) ""

   !!!---allocation de l'objet Fluide:

   allocate(M%MF(1:Nx,1:Ny,1:Nz))
  
   !!!---initialization de l'objet fluide
   DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
    call init_multifluide(M%MF(i,j,k)) 
   ENDDO ; ENDDO ; ENDDO

   !!!---nombre d'équations à résoudre partie conservative
   Neq=Nl+4 

   !!!!---vecteur 1D geometrique
   allocate(SVL(0:Nx+1))  ; SVL=0.0_PR
   allocate(SVR(0:Nx+1))  ; SVR=0.0_PR
   allocate(Sgeo(0:Nx+1)) ; Sgeo=0.0_PR

   !!!---vecteur 1D primitif:
   allocate(MF1D(0:Nx+1))
   do i=0,Nx+1
      call init_multifluide(MF1D(i)) 
   enddo

   !!!!---vecteur 1D des variavles conservatives :
   allocate(U(1:Neq,0:nx+1))    ; U(:,:)=0.0_PR
   allocate(dUdt(1:Neq,0:nx+1)) ; dUdt(:,:)=0.0_PR
   !!!!---vecteur 1D non-conservatif : fraction volumique + Finger tensor +
   !!                               énergie hydro de chaque phase   
   allocate(G(1:Nl,1:11,0:nx+1))    ; G(:,:,:)=0.0_PR
   allocate(dGdt(1:Nl,1:11,0:nx+1)) ; dGdt(:,:,:)=0.0_PR

   !!!---tableaux 1D pour la reconstruction

   if(Reconstruct)then
      allocate(MFrec_L(0:Nx+1))
      allocate(MFrec_R(0:Nx+1))
      DO i=1,Nx
       call init_multifluide(MFrec_L(i)) 
       call init_multifluide(MFrec_R(i)) 
      ENDDO
   endif

   !!!---initialisation du solveur SNES pour la relaxation de la pression

   call init_relaxp

   !!!---vecteur 1D des variables primitives
   !Nprim=9+Nl*13

   !write(iout,*) '  Nprim=', Nprim  
   !allocate(W1(1:Nprim,1:Nx))
   !allocate(W_L(1:Nprim,1:Nx))
   !allocate(W_R(1:Nprim,1:Nx))

   !!!================= POINTEURS DE PROCEDURES =================

   allocate(ph(1:Nl))
   do iph=1, Nl
    ph(iph)%soundspeed => soundspeed_sge
   enddo

   !!!---Reconstruction

   RECO => NOREC

   !!!---Limiteur:

   Limiteur => MINMOD

   !!!Limiteur => VANLEER

   !!!---Timescheme:
   if(trim(tscheme).eq.'RK3')then

      !!!---tableaux de travail RK3
      allocate(U0(1:Neq,0:nx+1))      ; U0(:,:)=0.0_PR
      allocate(G0(1:Nl,1:11,0:nx+1))  ; G0(:,:,:)=0.0_PR
      !allocate(F0(0:Nx+1))
      allocate(U1(1:Neq,0:nx+1))      ; U1(:,:)=0.0_PR
      allocate(G1(1:Nl,1:11,0:nx+1))  ; G1(:,:,:)=0.0_PR

      Timescheme => RK3

   elseif(trim(tscheme).eq.'RK2')then

      !!!---tableaux de travail RK2
      allocate(U0(1:Neq,0:nx+1))      ; U0(:,:)=0.0_PR
      allocate(G0(1:Nl,1:11,0:nx+1))  ; G0(:,:,:)=0.0_PR
      !allocate(F0(0:Nx+1))

      Timescheme => RK2 !!!_SSP

   elseif(trim(tscheme).eq.'RK2_SSP')then

      !!!---tableaux de travail RK2
      allocate(U0(1:Neq,0:nx+1))      ; U0(:,:)=0.0_PR
      allocate(G0(1:Nl,1:11,0:nx+1))  ; G0(:,:,:)=0.0_PR
      !allocate(F0(0:Nx+1))

      Timescheme => RK2_SSP

   elseif(trim(tscheme).eq.'RK1')then

      Timescheme => RK1
 
   elseif(trim(tscheme).eq.'COMPACT')then

      Timescheme => COMPACT
 
   else
      print*, ' ö please select a time integrator' ; stop
   endif

   !!!---SRC:

   if(just_thermo)then
      SRC => src_EM
   else
      SRC => src_euler1D
   endif

   write(iout,*) 'call init FINI !!!'

end subroutine init_euler 

!!!========================== SOLVEUR DE RIEMANN ===============================

subroutine HLLC_flux(MF_L,MF_R,us,vs,ws,F1,F2)
   
!!!(F_L,F_R,U_L,U_R,G_L,G_R,us,vs,ws,F1,F2)

    !!! input:
    !!! Nl : nombre de phases
    !!! MF_L, MF_R   : conteneurs des variables primitives left et right
    !!!   MF%F%f       : fractions volumiques des phases
    !!!   MF%F%rh      : densité des phases
    !!!   MF%F%p       : pression de phase
    !!!   MF%F%eh      : énergie hydrodynamique de phase
    !!!   MF%rh       : densité totale (sum(F%f%F%rh))
    !!!   MF%u, v, w : vitesse x, y, z
    !!!   MF%c       : vitesse du son
    !!!   MF%E       : énergie spécifique totale (sum(Ylel) + 1/2v^2)
    !!!   MF%F%sig(1:3)=s11, s12, s13 : composantes du tenseur de contraintes
    !!!   MF%F%a, b, c : composantes x,y,z des vecteurs de la base locale de la phase l
    !!! output :
    !!! F1 : flux des variables conservatives :( (fr)1*u, ... (fr)l*u ...,ru^2-s11,ruv-s12,ruw-s13,(rE-s11)u-s12v-s13w )   
    !!! F2 : flux des variables non conservatives 

    implicit none
  
    type(multifluide), intent(in) :: MF_L, MF_R
    !real(PR), intent(in) :: U_L(1:Nl+4), U_R(1:Nl+4)
    !real(PR), intent(in) :: G_L(1:Nl,1:11), G_R(1:Nl,1:11)  
    real(PR), intent(out) :: us,vs,ws
    real(PR), intent(out) :: F1(1:Nl+4)
    real(PR), intent(out) :: F2(1:Nl,1:11)
    !!!---Left/Right variables:   
    real(PR) :: fl_L(1:Nl), fl_R(1:Nl), rl_L(1:Nl), rl_R(1:Nl)
    real(PR) :: Yl_L(1:Nl), Yl_R(1:Nl), pl_L(1:Nl), pl_R(1:Nl)
    real(PR) :: elh_L(1:Nl), elh_R(1:Nl), ele_L(1:Nl), ele_R(1:Nl)
    real(PR) :: al_L(1:Nl,1:3),bl_L(1:Nl,1:3),cl_L(1:Nl,1:3)
    real(PR) :: al_R(1:Nl,1:3),bl_R(1:Nl,1:3),cl_R(1:Nl,1:3)
    real(PR) :: r_L,r_R,vx_L,vx_R,vy_L,vy_R,vz_L,vz_R,e_L,e_R,c_L,c_R
    !real(PR) :: p_L,p_R,
    real(PR) :: s11_L,S12_L,s13_L,s11_R,S12_R,s13_R
    real(PR) :: S_R,S_L,C0,C1,C2,C3
    !!!---star variables :
    type(multifluide) :: MFs
    real(PR) :: frls(1:Nl)
    real(PR) :: rs,s11s,s12s,s13s,Es,ps
    real(PR) ::rls(1:Nl),pls(1:Nl),elhs(1:Nl),als(1:Nl,1:3),bls(1:Nl,1:3),cls(1:Nl,1:3),eles(1:Nl)
    integer :: i, iph
    real(PR) :: coeff


    fl_L=MF_L%F(1:Nl)%f   ; fl_R=MF_R%F(1:Nl)%f
    rl_L=MF_L%F(1:Nl)%rh  ; rl_R=MF_R%F(1:Nl)%rh
    pl_L=MF_L%F(1:Nl)%p   ; pl_R=MF_R%F(1:Nl)%p
    elh_L=MF_L%F(1:Nl)%eh ; elh_R=MF_R%F(1:Nl)%eh
    do iph=1,Nl
       al_L(iph,:)=MF_L%F(iph)%a   ; al_R(iph,:)=MF_R%F(iph)%a
       bl_L(iph,:)=MF_L%F(iph)%b   ; bl_R(iph,:)=MF_R%F(iph)%b
       cl_L(iph,:)=MF_L%F(iph)%c   ; cl_R(iph,:)=MF_R%F(iph)%c
    enddo

    if(a2A)then
      do iph=1,Nl
       coeff=fl_L(iph)**(r13)
       al_L=al_L*coeff          
       bl_L=bl_L*coeff          
       cl_L=cl_L*coeff
       coeff=fl_R(iph)**(r13)
       al_R=al_R*coeff
       bl_R=bl_R*coeff
       cl_R=cl_R*coeff          
      enddo
    endif
 

    r_L  =MF_L%rh          ; r_R  =MF_R%rh
    vx_L =MF_L%vx          ; vx_R =MF_R%vx
    vy_L =MF_L%vy          ; vy_R =MF_R%vy
    vz_L =MF_L%vz          ; vz_R =MF_R%vz
    E_L  =MF_L%E           ; E_R  =MF_R%E
    c_L  =MF_L%c           ; c_R  =MF_R%c
    s11_L=MF_L%sig(1,1)    ; s11_R=MF_R%sig(1,1)
    s12_L=MF_L%sig(1,2)    ; s12_R=MF_R%sig(1,2)
    s13_L=MF_L%sig(1,3)    ; s13_R=MF_R%sig(1,3)
    ele_L =MF_L%ee          ; ele_R =MF_R%ee
 
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
 
       !!!do i=1,Nl+4
       !!!if(F1(i).ne.F1(i))then
       !!!  print*, 'F1 1:', i          
       !!!  stop
       !!!endif
       !!!enddo
     
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
       if(modeG.eq.1)then
          F2(1:Nl,11)=fl_L(1:Nl)*rl_L(1:Nl)*vx_L*elh_L(1:Nl)
       else
          F2(1:Nl,11)=fl_L(1:Nl)*rl_L(1:Nl)*vx_L*(elh_L(1:Nl)+ele_L(1:Nl))
       endif

    ELSEIF(S_R.le.0.0_PR)THEN

       us=vx_R ; vs=vy_R ; ws=vz_R

       F1(1:Nl)=fl_R(1:Nl)*rl_R(1:Nl)*vx_R
       F1(Nl+1)=r_R*vx_R**2-s11_R
       F1(Nl+2)=r_R*vx_R*vy_R-s12_R
       F1(Nl+3)=r_R*vx_R*vz_R-s13_R
       F1(Nl+4)=(r_R*E_R-s11_R)*vx_R-s12_R*vy_R-s13_R*vz_R
       !do i=1,Nl+4
       !if(F1(i).ne.F1(i))then
       !  print*, 'F1 2:', i
       !  stop          
       !endif
       !enddo

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
       if(modeG.eq.1)then
          F2(1:Nl,11)=fl_R(1:Nl)*rl_R(1:Nl)*vx_R*elh_R(1:Nl)
       else
          F2(1:Nl,11)=fl_R(1:Nl)*rl_R(1:Nl)*vx_R*(elh_R(1:Nl)+ele_R(1:Nl))
       endif

    ELSE

       us=((r_L*vx_L**2-s11_L)-(r_R*vx_R**2-s11_R)-S_L*r_L*vx_L+S_R*r_R*vx_R) / &
          (r_L*vx_L-r_R*vx_R-S_L*r_L+S_R*r_R)

       C0=1.0_PR/( (vx_R-S_R)*r_R-(vx_L-S_L)*r_L )
       C1=(vx_R-S_R)*r_R
       C2=(vx_L-S_L)*r_L
       C3=(vx_L-S_L)*r_L*(vx_R-S_R)*r_R    
       s11s=C0*( C1*s11_L-C2*s11_R+C3*(vx_R-vx_L) )
       s12s=C0*( C1*s12_L-C2*s12_R+C3*(vy_R-vy_L) )
       s13s=C0*( C1*s13_L-C2*s13_R+C3*(vz_R-vz_L) )


       IF(S_L.le.0.0_PR.and.us.ge.0.0_PR)THEN
          !print*, 'coucou S_L'

          frls(1:Nl)=fl_L(1:Nl)*rl_L(1:Nl)*(S_L-vx_L)/(S_L-us)

          rs=sum(frls(1:Nl))

          C0=1.0_PR/((vx_L-S_L)*r_L)
          vs=vy_L+(s12s-s12_L)*C0           
          ws=vz_L+(s13s-s13_L)*C0     
      
          Es=(r_L*E_L*(vx_L-S_L)-s11_L*vx_L-s12_L*vy_L-s13_L*vz_L+s11s*us+s12s*vs+s13s*ws)/&
           (rs*(us-S_L))
 

          F1(1:Nl)=fl_L(1:Nl)*rl_L(1:Nl)*vx_L
          F1(Nl+1)=r_L*vx_L**2-s11_L
          F1(Nl+2)=r_L*vx_L*vy_L-s12_L
          F1(Nl+3)=r_L*vx_L*vz_L-s13_L
          F1(Nl+4)=(r_L*E_L-s11_L)*vx_L-s12_L*vy_L-s13_L*vz_L

          F1(1:Nl)=F1(1:Nl)+S_L*(frls(1:Nl)-fl_L(1:Nl)*rl_L(1:Nl))
          F1(Nl+1)=F1(Nl+1)+S_L*(rs*us-r_L*vx_L)
          F1(Nl+2)=F1(Nl+2)+S_L*(rs*vs-r_L*vy_L)
          F1(Nl+3)=F1(Nl+3)+S_L*(rs*ws-r_L*vz_L)
          F1(Nl+4)=F1(Nl+4)+S_L*(rs*Es-r_L*E_L)

          !!!---old version : 
          !F1(1:Nl)=frls(1:Nl)*us
          !F1(Nl+1)=rs*us**2-s11s
          !F1(Nl+2)=rs*us*vs-s12s
          !F1(Nl+3)=rs*us*ws-s13s
          !F1(Nl+4)=(rs*Es-s11s)*us-s12s*vs-s13s*ws

          !do i=1,Nl+4
          !if(F1(i).ne.F1(i))then
          !  print*, 'F1 3:', i
          !  stop          
          !endif
          !enddo

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

          if(modeG.eq.2)then
            MFs=MF_L
            do iph=1,Nl
                 MFs%F(iph)%a=als(iph,:) ; MFs%F(iph)%b=bls(iph,:) ; MFs%F(iph)%c=cls(iph,:)
                 MFs%F(iph)%rh=rls(iph)
                 call init_Finger(MFs%F(iph))
                 !!!---uniquement la partie deviatotique du tenseur des contraintes
                 call Deviateur(MFs%F(iph)%sig(1:3,1:3))
                 !!!---Energie élastique
                 eles(iph)=energie_elastique(MFs%F(iph)%f,MFs%F(iph)%mu,MFs%F(iph)%rh0)
                 MFs%F(iph)%ee=eles(iph)
            enddo
          endif


          do iph=1,Nl

             if(modeG.eq.1)then
                pls(iph) = ph(iph)%pressure_star(MF_L%F(iph),rls(iph))  
             else
                pls(iph) = pressure_star_sge2(MF_L%F(iph),MFs%F(iph))
             endif  

             elhs(iph)=ph(iph)%energie_hydro(MF_L%F(iph),rls(iph),pls(iph))

          enddo

          !F2(1:Nl,1)=vx_L*fl_L(1:Nl)
          !do i=1,3
          !F2(1:Nl,1+i)=vx_L*al_L(1:Nl,i)
          !enddo
          !do i=1,3
          !F2(1:Nl,4+i)=vx_L*bl_L(1:Nl,i)
          !enddo
          !do i=1,3
          !F2(1:Nl,7+i)=vx_L*cl_L(1:Nl,i)
          !enddo
          !F2(1:Nl,11)=fl_L(1:Nl)*rl_L(1:Nl)*vx_L*elh_L(1:Nl)

          !do i=1,3
          !F2(1:Nl,1+i)=F2(1:Nl,1+i)+S_L*( als(1:Nl,i) - al_L(1:Nl,i))
          !enddo
          !F2(1:Nl,11)=F2(1:Nl,11)+S_L*( fl_L(1:Nl)*rls(1:Nl)*elhs(1:Nl) - fl_L(1:Nl)*rl_L(1:Nl)*elh_L(1:Nl) )


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
          !!!---(alpha.rho)_l* ou "alpha_l*.rho_l*" identique !
          !F2(1:Nl,11)=frls(1:Nl)*us*elhs(1:Nl)
          if(modeG.eq.1)then
             F2(1:Nl,11)=fl_L(1:Nl)*rls(1:Nl)*us*elhs(1:Nl)
          else
             F2(1:Nl,11)=fl_L(1:Nl)*rls(1:Nl)*us*(elhs(1:Nl)+eles(1:Nl))
          endif
          !!print*, 'F2 2',fl_L(2),rls(2),us,elhs(2)

       ELSEIF(S_R.ge.0.0_PR.and.us.le.0.0_PR)THEN

          !print*, 'coucou S_R'

          frls(1:Nl)=fl_R(1:Nl)*rl_R(1:Nl)*(S_R-vx_R)/(S_R-us)

          rs=sum(frls(1:Nl))

          C0=1.0_PR/((vx_R-S_R)*r_R)
          vs=vy_R+(s12s-s12_R)*C0           
          ws=vz_R+(s13s-s13_R)*C0     
      
          Es=(r_R*E_R*(vx_R-S_R)-s11_R*vx_R-s12_R*vy_R-s13_R*vz_R+s11s*us+s12s*vs+s13s*ws)/&
           (rs*(us-S_R))

          F1(1:Nl)=fl_R(1:Nl)*rl_R(1:Nl)*vx_R
          F1(Nl+1)=r_R*vx_R**2-s11_R
          F1(Nl+2)=r_R*vx_R*vy_R-s12_R
          F1(Nl+3)=r_R*vx_R*vz_R-s13_R
          F1(Nl+4)=(r_R*E_R-s11_R)*vx_R-s12_R*vy_R-s13_R*vz_R

          F1(1:Nl)=F1(1:Nl)+S_R*(frls(1:Nl)-fl_R(1:Nl)*rl_R(1:Nl))
          F1(Nl+1)=F1(Nl+1)+S_R*(rs*us-r_R*vx_R)
          F1(Nl+2)=F1(Nl+2)+S_R*(rs*vs-r_R*vy_R)
          F1(Nl+3)=F1(Nl+3)+S_R*(rs*ws-r_R*vz_R)
          F1(Nl+4)=F1(Nl+4)+S_R*(rs*Es-r_R*E_R)

          !---old version  
          !F1(1:Nl)=frls(1:Nl)*us
          !F1(Nl+1)=rs*us**2-s11s
          !F1(Nl+2)=rs*us*vs-s12s
          !F1(Nl+3)=rs*us*ws-s13s
          !F1(Nl+4)=(rs*Es-s11s)*us-s12s*vs-s13s*ws

       !!! do i=1,Nl+4
       !!! if(F1(i).ne.F1(i))then
       !!!   print*, 'F1 4:', i, Es, rs, us, vs, ws, s11s, s12s, s13s
       !!!   print*, C1,C2,C3
       !!!   print*, s11_L,s11_R,(vx_R-vx_L)
       !!!   print*, s12_L,s12_R,(vy_R-vy_L)
       !!!   print*, s13_L,s13_R,(vz_R-vz_L)

       !!!   print*, ''
       !!!   print*, 'F1 4:', C0, vx_R-S_R
       !!!   stop          
       !!! endif
       !!! enddo

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

          if(modeG.eq.2)then
            MFs=MF_R
            do iph=1,Nl
                 MFs%F(iph)%a=als(iph,:) ; MFs%F(iph)%b=bls(iph,:) ; MFs%F(iph)%c=cls(iph,:)
                 MFs%F(iph)%rh=rls(iph)
                 call init_Finger(MFs%F(iph))
                 !!!---uniquement la partie deviatotique du tenseur des contraintes
                 call Deviateur(MFs%F(iph)%sig(1:3,1:3))
                 !!!---Energie élastique
                 eles(iph)=energie_elastique(MFs%F(iph)%f,MFs%F(iph)%mu,MFs%F(iph)%rh0)
                 MFs%F(iph)%ee=eles(iph)
            enddo
          endif


          do iph=1,Nl
             !pls(iph)=hugoniot_pressure(iph,pl_R(iph),rl_R(iph),rls(iph),F_R)   
             !elhs(iph)=ph(iph)%energie_hydro(iph,rls(iph),pls(iph),F_R)
             !pls(iph)=hugoniot_pressure_sge(F_L%g_sge(iph),F_L%p_sge(iph),pl_L(iph),rl_L(iph),rls(iph))   
             !elhs(iph)=energie_hydro_sge(F_L%g_sge(iph),F_L%p_sge(iph),rls(iph),pls(iph))

             if(modeG.eq.1)then
                pls(iph)=ph(iph)%pressure_star(MF_R%F(iph),rls(iph))   
             else
                pls(iph)=pressure_star_sge2(MF_R%F(iph),MFs%F(iph))
             endif  

             !pls(iph) =pressure_star_sge2(iph,F_R,rls(iph),&
             !          als(iph,1:3),bls(iph,1:3),cls(iph,1:3))
  
             elhs(iph)=ph(iph)%energie_hydro(MF_R%F(iph),rls(iph),pls(iph))

          enddo

          !!!----------------------------------------
          !F2(1:Nl,1)=vx_R*fl_R(1:Nl)
          !do i=1,3
          !F2(1:Nl,1+i)=vx_R*al_R(1:Nl,i)
          !enddo
          !do i=1,3
          !F2(1:Nl,4+i)=vx_R*bl_R(1:Nl,i)
          !enddo
          !do i=1,3
          !F2(1:Nl,7+i)=vx_R*cl_R(1:Nl,i)
          !enddo
          !F2(1:Nl,11)=fl_R(1:Nl)*rl_R(1:Nl)*vx_R*elh_R(1:Nl)

          !do i=1,3
          !F2(1:Nl,1+i)=F2(1:Nl,1+i)+S_R*( als(1:Nl,i) - al_R(1:Nl,i) )
          !enddo
          !F2(1:Nl,11)=F2(1:Nl,11)+S_R*( fl_R(1:Nl)*rls(1:Nl)*elhs(1:Nl) - fl_R(1:Nl)*rl_R(1:Nl)*elh_R(1:Nl) )

          !!!---old
          F2(1:Nl,1)=us*MF_R%F(1:Nl)%f
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
          if(modeG.eq.1)then 
             F2(1:Nl,11)=MF_R%F(1:Nl)%f*rls(1:Nl)*us*elhs(1:Nl)
          else
             F2(1:Nl,11)=MF_R%F(1:Nl)%f*rls(1:Nl)*us*(elhs(1:Nl)+eles(1:Nl))
          endif
          !print*, 'F2 3',fl_R(2),rls(2),us,elhs(2)

       ENDIF

    ENDIF

end subroutine HLLC_flux

!!!========================== SUBROUTINE DE CONVERSION =========================
!!!---1D
subroutine primitive2conservative_1D

    implicit none
    integer :: i

    call cpu_time(t1)

    do i=1,N1D
        call primitive2conservative(MF1D(i),U(1:Nl+4,i),G(1:Nl,1:11,i))
        MF1D(i)%sonde(1)=U(Nl+1,i)
    enddo 

    call cpu_time(t2) ; t_p2c=t_p2c+t2-t1

end subroutine primitive2conservative_1D

subroutine primitive2conservative(MF,U,G)

   implicit none

   type(multifluide), intent(in) :: MF
   real(PR), intent(out) :: U(1:Nl+4), G(1:Nl,1:11)
   real(PR) :: c1, coeff
   integer :: iph

   G(1:Nl,1)=MF%F(1:Nl)%f

   do iph=1,Nl
      if(a2A)then
           coeff=MF%F(iph)%f**(r13)
      else
           coeff=1.0_PR
      endif
      G(iph,2:4) =MF%F(iph)%a(1:3)*coeff
      G(iph,5:7) =MF%F(iph)%b(1:3)*coeff
      G(iph,8:10)=MF%F(iph)%c(1:3)*coeff
   enddo

   if(modeG.eq.1)then
      G(1:Nl,11)=MF%F(1:Nl)%f*MF%F(1:Nl)%rh*MF%F(1:Nl)%eh
   else
      G(1:Nl,11)=MF%F(1:Nl)%f*MF%F(1:Nl)%rh*(MF%F(1:Nl)%eh+MF%F(1:Nl)%ee)
   endif

   U(1:Nl)=MF%F(1:Nl)%f*MF%F(1:Nl)%rh
   U(Nl+1)=MF%rh*MF%vx
   U(Nl+2)=MF%rh*MF%vy
   U(Nl+3)=MF%rh*MF%vz
   U(Nl+4)=MF%rh*MF%E
 
end subroutine primitive2conservative

subroutine conservative2primitive_1D

    implicit none
    integer :: i

    call cpu_time(t1)

    do i=1,N1D
       call conservative2primitive(U(1:Nl+4,i),G(1:Nl,1:11,i),MF1D(i)) 
    enddo 

    call cpu_time(t2) ; t_c2p=t_c2p+t2-t1

end subroutine conservative2primitive_1D

subroutine conservative2primitive(U,G,MF)

    !!!----------
    !!! Dérive les propriétés primitives a partir
    !!! de vecteur G des grandeurs non-conservatives:
    !!!
    !!! F(iph)%f=G(iph,1)
    !!! F(iph)%a,b,c: G(iph,2:10)
    !!! V1: F(iph)%rh=U(iph)/F(iph)%f
    !!! V2: F(iph)%rh=detFT12*F%rh0*F(iph)%fl0/F(iph)%fl
    !!! modeG=1:
    !!!  F(iph)%u=G(iph,11)/U(iph) (=f*rh*eh/f*rh)
    !!! modeG=2:  F(iph)%eh=G(iph,11)/U(iph)-F%ee ???
    !!!  F(iph)%p=ph(iph)%pressure(F,rh,u)
    !!!  F(iph)%ee=energie_elastique...
    !!!  MF%E=U(Nl+4)/MF%rh
    !!!  MF%ee=sum(Y*ee) MF%p=sum(Y*p) MF%u=MF%E-ee-&/év**2
    !!!----------

    implicit none
    
    real(PR), intent(in) :: U(1:Nl+4), G(1:Nl,1:11)
    type(multifluide), intent(inout) :: MF
    real(PR) :: c1, Ener
    integer ::  iph, i,j,k,i2
    real(PR) :: sumfl, coeff
    real(PR) :: rtemp(1:10)
    integer :: itemp(1:10)

    !!!---variables phasiques
    !call cpu_time(time_1)

    if(Nl.gt.1)then

        rtemp(1:Nl)=G(1:Nl,1)   
        forall(iph=1:Nl) itemp(iph)=iph 

        call hpsort(Nl,rtemp(1:Nl),itemp(1:Nl))

        sumfl=0.0_PR
 
        do j=1,Nl-1
          iph=itemp(j)  
          MF%F(iph)%f=G(iph,1)
          sumfl=sumfl+MF%F(iph)%f
        enddo

    endif

    iph=itemp(Nl)
    MF%F(iph)%f=1.0_PR-sumfl

    !MF%F(1:Nl)%f=G(1:Nl,1)
    !write(iout,*) ' al1 1:', MF%F(1)%a(1)

    do iph=1,Nl
       if(a2A)then
          coeff=MF%F(iph)%f**(-r13)
       else
          coeff=1.0_PR
       endif
       MF%F(iph)%a(1:3)=G(iph,2:4)*coeff 
       MF%F(iph)%b(1:3)=G(iph,5:7)*coeff
       MF%F(iph)%c(1:3)=G(iph,8:10)*coeff
    enddo
    !write(iout,*) ' al1 2:', MF%F(1)%a(1)

    !!!---V1

    MF%F(1:Nl)%rh=U(1:Nl)/MF%F(1:Nl)%f
    MF%rh=sum(U(1:Nl))
    MF%F(1:Nl)%Y=U(1:Nl)/MF%rh 

    !!!---V2
    !do iph=1,Nl

    !  call init_Finger(iph,F)

    !  F%rhl(1:Nl)=detFT12*F%rho0(1:Nl)*F%fl0(1:Nl)/F%fl(1:Nl)     !ph(iph)%rho0

    !  !!!---partie deviatotique du tenseur des contraintes
    !  call Deviateur(F%sigl(iph,1:3,1:3))

    !!  !!!---relaxation plastique ?

    !!  !!!---Energie élastique
    !  F%ele(iph)=energie_elastique(F%fl(iph),F%mu(iph),ph(iph)%rho0) 

    !enddo

    !F%rh=sum(F%fl(1:Nl)*F%rhl(1:Nl))
    !F%Yl(1:Nl)=F%fl(1:Nl)*F%rhl(1:Nl)/F%rh 
    !F%flrhl(1:Nl)=F%fl(1:Nl)*F%rhl(1:Nl)/F%rh 


    !!!------------------------------------------

    if(modeG.eq.1)then

       do iph=1,Nl
       
        !Ener=G(iph,11)/(F%fl(iph)*F%rhl(iph))
       
        Ener=G(iph,11)/U(iph)

        MF%F(iph)%eh0=MF%F(iph)%eh
        MF%F(iph)%eh=Ener

        MF%F(iph)%p=ph(iph)%pressure(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%eh) 

        !MF%F(iph)%p=(MF%F(iph)%g_sge-1.0_PR)*MF%F(iph)%rh*Ener-MF%F(iph)%g_sge*MF%F(iph)%p_sge

        if(MF%F(iph)%p.gt.1.e20_PR)then
           print*, 'conserv2prim: pl>1e20:',MF%F(iph)%p,G(iph,11),MF%F(iph)%f
           stop
        !elseif(F%pl(iph).le.0.e0_PR)then
        !   print*, 'conserv2prim: pl<0:',F%pl(iph),G(iph,11),F%fl(iph)
        !   stop 
        endif
       enddo

    elseif(modeG.eq.2)then !!! bizarre...

       do iph=1,Nl

        call init_Finger(MF%F(iph))

        !!!---Energie élastique
        MF%F(iph)%ee=energie_elastique(MF%F(iph)%f,MF%F(iph)%mu,MF%F(iph)%rh0) 
 
        Ener=G(iph,11)/U(iph)

        MF%F(iph)%eh0=MF%F(iph)%eh
        MF%F(iph)%eh=Ener-MF%F(iph)%ee

        !MF%F(iph)%p=(MF%F(iph)%g_sge-1.0_PR)*MF%F(iph)%rh*MF%F(iph)%eh-MF%F(iph)%g_sge*MF%F(iph)%p_sge

        MF%F(iph)%p=ph(iph)%pressure(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%eh) 

        if(MF%F(iph)%p.gt.1.e20_PR)then
           print*, 'conserv2prim: pl>1e20:',MF%F(iph)%p,G(iph,11),MF%F(iph)%f
           stop
        !elseif(F%pl(iph).le.0.e0_PR)then
        !   print*, 'conserv2prim: pl<0:',F%pl(iph),G(iph,11),F%fl(iph)
        !   stop 
        endif
       enddo


    endif

    !if(rh.lt.0.0_PR)then
    !   print*, 'rh<0', U(1:Nl), fl(1:Nl)
    !   stop
    !endif

    !do iph=1,Nl
    !   pl(iph)=pressure(iph,rhl(iph),elh(iph))
    !enddo

    !!!---variables de mélange
    MF%vx=U(Nl+1)/MF%rh
    MF%vy=U(Nl+2)/MF%rh
    MF%vz=U(Nl+3)/MF%rh
    MF%E=U(Nl+4)/MF%rh

    !!! ATTENTION sens de calculer energie meca 
    !!! même si p hors equilibre ?
    !!!---p ne dépend que de G may be <0 !
    MF%p=sum(MF%F(1:Nl)%f*MF%F(1:Nl)%p)

    call MAJ_meca(MF)

    MF%ee=sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%ee)

    MF%u=MF%E-MF%ee-0.5_PR*(MF%vx**2+MF%vy**2+MF%vz**2)

    !!!---Temperature:
    !call MAJ_thermo
    !do iph=1,Nl 
    !   F%Tl(iph)=F%Tl(iph)+(F%elh(iph)-F%elh0(iph))/F%cvl(iph)
    !   F%cvl(iph)=specific_heat(F%Tl(iph),F%pl(iph))/F%g_sge(iph)   
    !enddo
    !F%T=sum(F%fl(:)*F%Tl(:))

    !call cpu_time(time_2)

    !t_c2p=t_c2p+time_2-time_1

end subroutine conservative2primitive

subroutine nonconservative2primitive_1D

    implicit none
    integer :: i

    call cpu_time(t1)   

    do i=1,N1D
        call nonconservative2primitive(G(1:Nl,1:11,i),MF1D(i)) 
        !call conservative2primitive(U(1:Nl+4,i),G(1:Nl,1:11,i),MF1D(i)) 
        if(abs(MF1D(i)%ee).gt.1.0e15_PR)then
            print*, 'PROBLEM u 2 i=', i
            stop
        endif
    enddo 

    call cpu_time(t2) ; t_nc2p=t_nc2p+t2-t1

end subroutine nonconservative2primitive_1D

subroutine nonconservative2primitive(G,MF) 

    !!!----------
    !!! Dérive les propriétés primitives a partir
    !!! de vecteur G des grandeurs non-conservatives:
    !!!
    !!! F(iph)%f=G(iph,1)
    !!! F(iph)%a,b,c: G(iph,2:10)
    !!! rhoeh(iph) = G(iph,11)/F(iph)%fl
    !!! p=rhoe*(g-1)-gamp*inf !!! (sge)
    !!!----------

    implicit none
    
    real(PR), intent(in) :: G(1:Nl,1:11)
    type(multifluide), intent(inout) :: MF
    real(PR) :: c1, Ener
    integer ::  iph, i,j,k,i2
    real(PR) :: sumfl, coeff
    real(PR) :: rtemp(1:10)
    integer :: itemp(1:10)

    !!!---variables phasiques

    if(Nl.gt.1)then

        rtemp(1:Nl)=G(1:Nl,1)   
        forall(iph=1:Nl) itemp(iph)=iph 

        call hpsort(Nl,rtemp(1:Nl),itemp(1:Nl))

        sumfl=0.0_PR
 
        do j=1,Nl-1
          iph=itemp(j)  
          MF%F(iph)%f=G(iph,1)
          sumfl=sumfl+MF%F(iph)%f
        enddo

        iph=itemp(Nl)
        MF%F(iph)%f=1.0_PR-sumfl

    else

       MF%F(1)%f=1.0_PR

    endif

    if(a2A)then
       do iph=1,Nl
          coeff=MF%F(iph)%f**(-r13)
          MF%F(iph)%a(1:3)=G(iph,2:4)*coeff 
          MF%F(iph)%b(1:3)=G(iph,5:7)*coeff
          MF%F(iph)%c(1:3)=G(iph,8:10)*coeff
       enddo
    else
       do iph=1,Nl
          MF%F(iph)%a(1:3)=G(iph,2:4) 
          MF%F(iph)%b(1:3)=G(iph,5:7)
          MF%F(iph)%c(1:3)=G(iph,8:10)
       enddo
    endif

    call MAJ_meca(MF)

    !do iph=1,Nl

    !  call init_Finger(iph,F)

    !  !F%rhl(iph)=detFT12*F%rho0(iph)*F%fl0(iph)/F%fl(iph)     !ph(iph)%rho0

    !  !F%rhl(iph)=detFT12*F%rho0(iph)    !!!!!*F%fl0(iph)/F%fl(iph)     !ph(iph)%rho0

    !  !!!---partie deviatotique du tenseur des contraintes
    !  call Deviateur(F%sigl(iph,1:3,1:3))

    !!  !!!---relaxation plastique ?

    !!  !!!---Energie élastique
    !  F%ele(iph)=energie_elastique(F%fl(iph),F%mu(iph),ph(iph)%rho0) 

    !enddo

    !---V1
    !F%rh=sum(F%fl(1:Nl)*F%rhl(1:Nl))
    !F%Yl(1:Nl)=F%fl(1:Nl)*F%rhl(1:Nl)/F%rh 
    !F%flrhl(1:Nl)=F%fl(1:Nl)*F%rhl(1:Nl)/F%rh 
    !!F%ee=sum(F%Yl(1:Nl)*F%ele(1:Nl))

    !do iph=1,Nl

    !    Ener=G(iph,11)/(F%fl(iph)*F%rhl(iph))

    !    F%elh0(iph)=F%elh(iph)
    !    F%elh(iph)=Ener

    !    F%pl(iph)=(F%g_sge(iph)-1.0_PR)*F%rhl(iph)*Ener-F%g_sge(iph)*F%p_sge(iph)

    !enddo

    !F%u=sum(F%Yl(1:Nl)*F%elh(1:Nl))

    !---V2 : aucune hypothèse sur rhl

    do iph=1,Nl

        Ener=G(iph,11)/(MF%F(iph)%f)

        MF%F(iph)%p=(MF%F(iph)%g_sge-1.0_PR)*Ener-MF%F(iph)%g_sge*MF%F(iph)%p_sge

    enddo

end subroutine nonconservative2primitive





!!!============================ SUBROUTINE EOS =================================
!!!---0D

subroutine MAJ_mixt(MF)   !!!!rhl,fl,Yl,pl,sigl,elh,ele,vx,vy,vz,rh,p,sig,e,c)

   implicit none
   type(multifluide), intent(inout) :: MF
   integer :: i,j,k,iph

   !!!write(iout,*) '   calcul des variables de mélange : rh, Yl, E, Sij' 

   if(MF%F(1)%rh.lt.0.0_PR)then
   print*, 'MAJ_mixt:'
   print*, 'rhl:', MF%F(:)%rh
   print*, 'fl:' , MF%F(:)%f
   print*, 'Yl:' , MF%F(:)%Y
   print*, 'pl:' , MF%F(:)%p
   print*, 'elh:', MF%F(:)%eh
   print*, 'P:'  , MF%F(:)%p
   endif

   MF%rh=sum(MF%F(1:Nl)%f*MF%F(1:Nl)%rh) 
   MF%u =sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%eh)
   MF%ee=sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%ee)
   MF%E =MF%u+MF%ee+0.5_PR*(MF%vx**2+MF%vy**2+MF%vz**2)
   MF%p =sum( MF%F(1:Nl)%f*MF%F(1:Nl)%p )
   forall(i=1:3,j=1:3) MF%sig(i,j)=sum(MF%F(1:Nl)%f*MF%F(1:Nl)%sig(i,j))
   MF%c=soundspeed_mixt(MF)

end subroutine MAJ_mixt

subroutine MAJ_THERMO(MF)
   implicit none
   type(multifluide), intent(inout) :: MF
   integer :: iph


   !!!---Cv
   MF%cv=sum(MF%F(:)%Y*MF%F(:)%cv)

   !!!---Temperature:
   do iph=1,Nl 

     !if(ph(iph)%typ.eq.1)then

      !!!---elh0 mise a jour dans conservtive2primitive
      
      !MF%F(iph)%T=MF%F(iph)%T+(MF%F(iph)%eh-MF%F(iph)%eh0)/MF%F(iph)%cv

      MF%F(iph)%T=ph(iph)%T_from_rP(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%p)

      !F%Tl(iph)=300.0_PR+(F%elh(iph)-F%elh0(iph))/F%cvl(iph)

     ! endif

      !F%cvl(iph)=specific_heat(min(F%Tl(iph),70000.0_PR),min(F%pl(iph),1e8_PR))/F%g_sge(iph)   
   enddo

   MF%T=sum(MF%F(:)%Y*MF%F(:)%cv*MF%F(:)%T)/MF%cv

end subroutine MAJ_THERMO

!!!============== EOS SGET ==================

pure real(PR) function pressure_sgeT(F,rh,u)

   implicit none
   !!! stiffened gas equation with temperature
   !!! p=rho(Gam-1)(u-q) - GammP0
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: u, rh

   pressure_sgeT=rh*(F%g_sge - 1.0_PR)*(u - F%q ) - F%g_sge*F%p_sge

end function pressure_sgeT

pure real(PR) function energie_hydro_sgeT(F,rh,p)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) + q !!! J/kg
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, p

   energie_hydro_sgeT=( p + F%g_sge*F%p_sge )/( (F%g_sge - 1.0_PR)*rh ) + F%q

end function energie_hydro_sgeT

real(PR) function soundspeed_sgeT(F,rh,p)

   implicit none
   !!! stiffened gas equation
   !!! p=rho(Gam-1)e-GammP0
   !!! e=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(e*(Gam-1))
   !!! c^2=Gam*(p+p0)/rho
   !!! integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, p

   if((p+F%p_sge)/rh.lt.0.0_PR)then
    ierr=1
    print*, 'PROBLEM soundspeed:', F%iph, p, rh, F%p_sge
    !stop
   endif

   !!! c² = ( p/rh² - (de/drh)p ) / (de/dp)rh
   !!! same as soundspeed_sge

   soundspeed_sgeT=sqrt(F%g_sge*(p+F%p_sge)/rh+r43*F%mu/F%rh0)

end function soundspeed_sgeT

pure real(PR) function pressure_star_sgeT(F,rs)

   implicit none
   real(PR), intent(in) :: rs
   type(fluide), intent(in) :: F
   real(PR) :: G, Gm1, Gp1, nus, nu0, p0, pinf, e0
   real(PR) :: sig(1:3,1:3)

   !!!---formule de Hugoniot Saurel 2009 p1689 : idem sge ?
   G=F%g_sge
   Gm1=(G-1.0_PR)
   Gp1=(G+1.0_PR)
   pressure_star_sgeT=(F%p+F%p_sge)*(Gm1*F%rh-Gp1*rs)/(Gm1*rs-Gp1*F%rh)-F%p_sge 

end function pressure_star_sgeT

pure real(PR) function T_from_rP_sgeT(F,rh,p)

   implicit none
   real(PR), intent(in) :: rh, p
   type(fluide), intent(in) :: F

   !!!---formule 16 de Pepitas 2009 p 752
   T_from_rP_sgeT=(p+F%p_sge)/((F%g_sge-1.0_PR)*rh*F%Cv)

end function T_from_rP_sgeT

pure real(PR) function T_from_pe_sgeT(F,p,e)

   implicit none
   !!!integer, intent(in) :: iph
   real(PR), intent(in) :: p, e
   type(fluide), intent(in) :: F

   !!!--- deduit de formules 15-16 de Pepitas 2009 p 752
   T_from_pe_sgeT=(e-F%q)*(p+F%p_sge)/((p+F%g_sge*F%p_sge)*F%Cv)

end function T_from_pe_sgeT

pure real(PR) function T_from_re_sgeT(F,rh,e)

   implicit none
   real(PR), intent(in) :: rh, e
   type(fluide), intent(in) :: F

   !!!--- deduit de formules 15-16 de Pepitas 2009 p 752
   T_from_re_sgeT=(rh*(e-F%q)-F%p_sge)/(rh*F%Cv)

end function T_from_re_sgeT

pure real(PR) function rh_from_PT_sgeT(F,P,T)

   implicit none
   real(PR), intent(in) :: P, T
   type(fluide), intent(in) :: F

   !!!--- formule 16 de Pepitas 2009 p 752
   rh_from_PT_sgeT=(P+F%p_sge)/((F%g_sge-1.0_PR)*F%Cv*T)

end function rh_from_PT_sgeT

!!!--- derivees pour Newton

pure real(PR) function dedT_rho_sgeT(F,T,nu)

   implicit none
   real(PR), intent(in) :: T, nu
   type(fluide), intent(in) :: F

   dedT_rho_sgeT=F%Cv

end function dedT_rho_sgeT

pure real(PR) function dednu_T_sgeT(F,T,nu)

   implicit none
   real(PR), intent(in) :: T, nu
   type(fluide), intent(in) :: F

   dednu_T_sgeT=F%p_sge

end function dednu_T_sgeT

pure real(PR) function dPdT_rho_sgeT(F,T,nu)

   implicit none
   real(PR), intent(in) :: T, nu
   type(fluide), intent(in) :: F

   dPdT_rho_sgeT=(F%g_sge-1.0_PR)*F%Cv/nu

end function dPdT_rho_sgeT

pure real(PR) function dPdnu_T_sgeT(F,T,nu)

   implicit none
   real(PR), intent(in) :: T, nu
   type(fluide), intent(in) :: F

   dPdnu_T_sgeT=-(F%g_sge-1.0_PR)*F%Cv*T/nu**2

end function dPdnu_T_sgeT

!!!============== EOS SGE ==================

!!!---V2: iph et F en argument !!!! 

pure real(PR) function pressure_sge(F,rh,u)

   implicit none
   !!! stiffened gas equation
   !!! p=rho(Gam-1)u-GammP0
   !!! u=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(u*(Gam-1))
   !!! c^2=Gam(p+p0)/rho
   !real(PR) :: pressure
   !integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: u, rh

   !pressure=max(rh*(g_sge(iph) - 1.0_PR)*e - g_sge(iph)*p_sge(iph), 0.0_PR)

   !pressure=F%rhl(iph)*(F%g_sge(iph) - 1.0_PR)*F%elh(iph) - F%g_sge(iph)*F%p_sge(iph)

   pressure_sge=rh*(F%g_sge - 1.0_PR)*u - F%g_sge*F%p_sge

   !if(pressure_sge.lt.0.0_PR)then
   !   
   !   print*, 'p<0:'
   !   print*, 'iph, rh, rh*(g_sge(iph)-1.0_PR)*u, u, g_sge(iph)*p_sge(iph)'
   !   print*, iph, rh, rh*(F%g_sge(iph)-1.0_PR)*F%elh(iph), F%elh(iph), F%g_sge(iph)*F%p_sge(iph)
   !   stop
   !endif

end function pressure_sge

pure real(PR) function energie_hydro_sge(F,rh,p)

   implicit none
   !!! p=rho(Gam-1)e-GammP0  !!! Pa
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! rh=(p+GammP0)/(e*(Gam-1)) !!! kg/m^3
   !!! c^2=Gam(p+p0)/rho !!! m/s
   !!! eau: Gam=6.1 p0=2e9 ~ molecular attraction
   !!! integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, p

   energie_hydro_sge=( p + F%g_sge*F%p_sge )/( (F%g_sge - 1.0_PR)*rh )

end function energie_hydro_sge

real(PR) function soundspeed_sge(F,rh,p)

   implicit none
   !!! stiffened gas equation
   !!! p=rho(Gam-1)e-GammP0
   !!! e=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(e*(Gam-1))
   !!! c^2=Gam*(p+p0)/rho
   !!! integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, p

   if((p+F%p_sge)/rh.lt.0.0_PR)then
    ierr=1
    print*, 'PROBLEM soundspeed:', F%iph, p, rh, F%p_sge
    !stop
   endif

   !soundspeed_sge=sqrt(F%g_sge(iph)*(p+F%p_sge(iph))/rh)

   soundspeed_sge=sqrt(F%g_sge*(p+F%p_sge)/rh+r43*F%mu/F%rh0)

   !if((F%pl(iph)+F%p_sge(iph))/F%rhl(iph).lt.0.0_PR)then
   ! ierr=1
   ! print*, 'PROBLEM soundspeed:', iph, F%pl(iph),F%rhl(iph)
   ! stop
   !endif

   !soundspeed_sge=sqrt(F%g_sge(iph)*(F%pl(iph)+F%p_sge(iph))/F%rhl(iph))

end function soundspeed_sge

pure real(PR) function pressure_star_sge(F,rs)

   implicit none
   !!!integer, intent(in) :: iph
   real(PR), intent(in) :: rs
   type(fluide), intent(in) :: F
   real(PR) :: G, Gm1, Gp1, nus, nu0, p0, pinf, e0
   real(PR) :: sig(1:3,1:3)

   !!!---formule de Hugoniot Saurel 2009 p1689 :
   G=F%g_sge
   Gm1=(G-1.0_PR)
   Gp1=(G+1.0_PR)
   pressure_star_sge=(F%p+F%p_sge)*(Gm1*F%rh-Gp1*rs)/(Gm1*rs-Gp1*F%rh)-F%p_sge 

   !!!---formule de Hugoniot Saurel 2009 p1689 modif Favrie 2015 ? :
   !G=F%g_sge(iph)
   !Gm1=(G-1.0_PR)
   !Gp1=(G+1.0_PR)
   !pressure_star_sge=(-F%sigl(iph,1,1)+F%p_sge(iph))*(Gm1*F%rhl(iph)-Gp1*rs)/(Gm1*rs-Gp1*F%rhl(iph))-F%p_sge(iph) 

   !!!---Equation (57) p 6065 Favrie 2009 + eos sge
   !G=F%g_sge(iph)
   !Gm1=(G-1.0_PR)
   !p0=F%pl(iph)
   !pinf=F%p_sge(iph)
   !e0=F%elh(iph)
   !nus=1.0_PR/rs
   !nu0=1.0_PR/F%rhl(iph)

   !pressure_star_sge=-(p0*Gm1*(nus-nu0)+2.0_PR*(nus*G*pinf-Gm1*e0))/&
   !                   ((nus-nu0)*Gm1+2.0_PR*nus)

end function pressure_star_sge

pure real(PR) function pressure_star_sge2(F0,Fs)

   implicit none
   !integer, intent(in) :: iph
   type(fluide), intent(in) :: F0, Fs
   real(PR) :: G, Gm1, Gp1, nus, nu0, p0, pinf, e0, ee
   real(PR) :: p, s11, s12, s13, s11_0, s12_0, s13_0

   !!!---Equation (57) p 6065 Favrie 2009 + eos sge total

   G=F0%g_sge
   Gm1=(G-1.0_PR)
   p0=F0%p
   pinf=F0%p_sge
   e0=F0%eh+F0%ee
   nus=1.0_PR/Fs%rh
   nu0=1.0_PR/F0%rh

   s11_0    = F0%sig(1,1)
   s12_0    = F0%sig(1,2)
   s13_0    = F0%sig(1,3)
   s11      = Fs%sig(1,1)
   s12      = Fs%sig(1,2)
   s13      = Fs%sig(1,3)

   ee=Fs%ee

   p=e0-ee-nus/Gm1*G*pinf+0.5_PR*((nus-nu0)*(s11+s11_0)+&
                                   (nu0*F0%a(2)-nus*Fs%a(2))*(s12+s12_0)+&
                                   (nu0*F0%a(3)-nus*Fs%a(3))*(s13+s13_0))


   p=p/(nus/Gm1+0.5_PR*(nus-nu0))

   pressure_star_sge2=p 

end function pressure_star_sge2

!real(PR) function pressure_star_sge2(iph,F0,rs,als,bls,cls)
!
!   implicit none
!   integer, intent(in) :: iph
!   type(fluide), intent(in) :: F0
!   real(PR), intent(in) :: rs, als(1:3), bls(1:3), cls(1:3)
!   type(fluide) :: Fs
!   real(PR) :: G, Gm1, Gp1, nus, nu0, p0, pinf, e0, ele
!   real(PR) :: p, s11, s12, s13, s11_0, s12_0, s13_0
!   real(PR) :: al0(1:3), bl0(1:3), cl0(1:3)
!
!   !!!---Equation (57) p 6065 Favrie 2009 + eos sge total
!
!   G=F0%g_sge(iph)
!   Gm1=(G-1.0_PR)
!   p0=F0%pl(iph)
!   pinf=F0%p_sge(iph)
!   e0=F0%elh(iph)+F0%ele(iph)
!   nus=1.0_PR/rs
!   nu0=1.0_PR/F0%rhl(iph)
!
!   Fs=F0
!   Fs%rhl(iph)=rs
!   Fs%al(iph,1:3)=als(1:3)
!   Fs%bl(iph,1:3)=bls(1:3)
!   Fs%cl(iph,1:3)=cls(1:3)
!
!   call init_Finger(iph,Fs) 
!
!   !!!---uniquement la partie deviatotique du tenseur des contraintes
!   call Deviateur(Fs%sigl(iph,1:3,1:3))
!
!   !!!---Energie élastique
!   ele=energie_elastique(Fs%fl(iph),Fs%mu(iph),ph(iph)%rho0)
!
!   al0(1:3) = F0%al(iph,1:3)
!   s11_0    = F0%sigl(iph,1,1)
!   s12_0    = F0%sigl(iph,1,2)
!   s13_0    = F0%sigl(iph,1,3)
!   s11      = Fs%sigl(iph,1,1)
!   s12      = Fs%sigl(iph,1,2)
!   s13      = Fs%sigl(iph,1,3)
!
!   ele=Fs%ele(iph)
!
!   p=e0-ele-nus/Gm1*G*pinf+0.5_PR*((nus-nu0)*(s11+s11_0)+&
!                                   (nu0*al0(2)-nus*als(2))*(s12+s12_0)+&
!                                   (nu0*al0(3)-nus*als(3))*(s13+s13_0))
!
!
!   p=p/(nus/Gm1+0.5_PR*(nus-nu0))
!
!   pressure_star_sge2=p 
!
!end function pressure_star_sge2

pure real(PR) function ddnu_energie_hydro_sge(F,rh,p)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! de/drho=-(p+Gam*Pinf)/( (Gam-1)*rho^2 )
   !!! integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, p

   ddnu_energie_hydro_sge= ( p + F%g_sge*F%p_sge )/( F%g_sge - 1.0_PR )

end function ddnu_energie_hydro_sge

pure real(PR) function ddp_energie_hydro_sge(F,rh,p)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! de/dp=1.0/( (Gam-1)*rho )
   !!! integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, p

   ddp_energie_hydro_sge=1.0_PR/( (F%g_sge - 1.0_PR)*rh )

end function ddp_energie_hydro_sge

pure real(PR) function u_from_TP_sge(F,T,p)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! de/dp=1.0/( (Gam-1)*rho )
   !!!integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: T, p

   u_from_TP_sge=ph(F%iph)%ener0+F%cv*(T-300.0_PR)

end function u_from_TP_sge

pure real(PR) function r_from_uP_sge(F,u,P)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! de/dp=1.0/( (Gam-1)*rho )
   !!! integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: u, P

   r_from_uP_sge=(P+F%g_sge*F%p_sge)/&
                 (u*(F%g_sge-1.0_PR))

end function r_from_uP_sge

pure real(PR) function T_from_rP_sge(F,rh,p)

   implicit none
   !!!integer, intent(in) :: iph
   real(PR), intent(in) :: rh, p
   type(fluide), intent(in) :: F

   !!!---formule 16 de Pepitas 2009 p 752
   T_from_rP_sge=1.0_PR/rh*(p+F%p_sge)/((F%g_sge-1.0_PR)*F%Cv)

end function T_from_rP_sge

!real(PR) function energie_hydro_sge(iph,F,rh,p)
!
!   implicit none
!   !!! p=rho(Gam-1)e-GammP0  !!! Pa
!   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
!   !!! rh=(p+GammP0)/(e*(Gam-1)) !!! kg/m^3
!   !!! c^2=Gam(p+p0)/rho !!! m/s
!   !!! eau: Gam=6.1 p0=2e9 ~ molecular attraction
!   integer, intent(in) :: iph
!   type(multifluide), intent(in) :: F
!   real(PR), intent(in) :: rh, p
!
!   energie_hydro_sge=( p + F%g_sge(iph)*F%p_sge(iph) )/( (F%g_sge(iph) - 1.0_PR)*rh )
!
!end function energie_hydro_sge

!!!______________________________________

!!!---------V0 : g et pinf en argument
!!!______________________________________

real(PR) function pressure_sge0(rh,u,g_sge,p_sge)

   implicit none
   !!! stiffened gas equation
   !!! p=rho(Gam-1)u-GammP0
   !!! u=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(u*(Gam-1))
   !!! c^2=Gam(p+p0)/rho
   !real(PR) :: pressure
   real(PR), intent(in) :: g_sge,p_sge, u, rh

   !pressure=max(rh*(g_sge(iph) - 1.0_PR)*e - g_sge(iph)*p_sge(iph), 0.0_PR)

   pressure_sge0=rh*(g_sge - 1.0_PR)*u - g_sge*p_sge

   if(pressure_sge0.lt.0.0_PR)then
      
      print*, 'p<0:'
      print*, 'rh, rh*(g_sge-1.0)*u, u, g_sge*p_sge'
      print*, rh, rh*(g_sge-1.0_PR)*u, u, g_sge*p_sge
      stop
   endif

end function pressure_sge0

real(PR) function energie_hydro_sge0(rhl,pl,g,pinf)

   implicit none
   !!! p=rho(Gam-1)e-GammP0  !!! Pa
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! rh=(p+GammP0)/(e*(Gam-1)) !!! kg/m^3
   !!! c^2=Gam(p+p0)/rho !!! m/s
   !!! eau: Gam=6.1 p0=2e9 ~ molecular attraction
   real(PR), intent(in) :: g,pinf,rhl, pl

   energie_hydro_sge0=( pl + g*pinf )/( (g - 1.0_PR)*rhl )

end function energie_hydro_sge0


real(PR) function soundspeed_sge0(rhl,pl,g,pinf)

   implicit none
   !!! stiffened gas equation
   !!! p=rho(Gam-1)e-GammP0
   !!! e=(p+GammP0)/(rho*(Gam-1))
   !!! rh=(p+GammP0)/(e*(Gam-1))
   !!! c^2=Gam*(p+p0)/rho
   real(PR), intent(in) :: g,pinf,rhl, pl

   if((pl+pinf)/rhl.lt.0.0_PR)then
    ierr=1
    print*, 'PROBLEM soundspeed:', pl,rhl,pinf
    stop
   endif

   soundspeed_sge0=sqrt(g*(pl+pinf)/rhl)

end function soundspeed_sge0

real(PR) function pressure_star_sge0(p0,r0,rf,g,pinf)

   implicit none
   real(PR), intent(in) :: g,pinf,p0, r0, rf
   real(PR) :: Gm1, Gp1

   !!!---formule de Hugoniot SGE :
   Gm1=(g-1.0_PR)
   Gp1=(g+1.0_PR)

   pressure_star_sge0=(p0+pinf)*(Gm1*r0-Gp1*rf)/(Gm1*rf-Gp1*r0)-pinf 

end function pressure_star_sge0

real(PR) function ddnu_energie_hydro_sge0(rh,p,g,pinf)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! de/drho=-(p+Gam*P0)/( (Gam-1)*rho^2 )
   real(PR), intent(in) :: g,pinf,rh, p

   ddnu_energie_hydro_sge0= ( p + g*pinf )/( g - 1.0_PR )

end function ddnu_energie_hydro_sge0

real(PR) function ddp_energie_hydro_sge0(rh,p,g,pinf)

   implicit none
   !!! e=(p+GammP0)/(rho*(Gam-1)) !!! J/kg
   !!! de/dp=1.0/( (Gam-1)*rho )
   real(PR), intent(in) :: g, pinf, rh, p

   ddp_energie_hydro_sge0=1.0_PR/( (g - 1.0_PR)*rh )

end function ddp_energie_hydro_sge0


!!!============== EOS TABULEE ==================

real(PR) function pressure_tab(F,rh,u)

   implicit none
   !!!integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, u
   real(PR) :: T

   T=T_from_ru(F%iph,rh,u)
   !T=F%Tl(iph) 
   pressure_tab=P_from_rT(F%iph,rh,T)

end function pressure_tab

real(PR) function energie_hydro_tab(F,rh,p)

   implicit none
   !!!integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, p
   real(PR) :: T

   T=T_from_rP(F%iph,rh,P)
   !T=F%Tl(iph) 

   energie_hydro_tab=E_from_rT(F%iph,rh,T)

end function energie_hydro_tab

real(PR) function soundspeed_tab(F,rh,p)

   implicit none
   !!!integer, intent(in) :: iph
   type(fluide), intent(in) :: F
   real(PR), intent(in) :: rh, p
   real(PR) :: T

   T=T_from_rP(F%iph,rh,p)
   !T=F%Tl(iph) 

   soundspeed_tab=c_from_rT(F%iph,F%rh,T)

end function soundspeed_tab


!!!============== EOS MIXT ==================
real(PR) function soundspeed_mixt_Wood(MF) 

   implicit none 
   type(multifluide), intent(in) :: MF 
   real(PR) :: c1
   integer :: iph

   !!!---Wood speed of sound:
   c1=0.0_PR
   do iph=1,Nl
     c1=c1+MF%F(iph)%f/(MF%F(iph)%rh*ph(iph)%soundspeed(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%p)**2)
   enddo
   c1=c1*MF%rh
   soundspeed_mixt_Wood=sqrt(1.0_PR/c1)

end function soundspeed_mixt_Wood

real(PR) function soundspeed_mixt_frozen(MF) 

   implicit none 
   type(multifluide), intent(in) :: MF 
   real(PR) :: c1
   integer :: iph

   !!!---frozen speed of sound:
   c1=0.0_PR
   do iph=1,Nl
     c1=c1+MF%F(iph)%Y*(ph(iph)%soundspeed(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%p))**2
   enddo
   soundspeed_mixt_frozen=sqrt(c1)

   !!!---Wood speed of sound:
   !c1=0.0_PR
   !do iph=1,Nl
   !  c1=c1+MF%F(iph)%f/(MF%F(iph)%rh*ph(iph)%soundspeed(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%p)**2)
   !enddo
   !c1=c1*MF%rh
   !soundspeed_mixt=sqrt(1.0_PR/c1)

   !c1=0.0_PR
   !do iph=1,Nl
   !  c1=max(c1,soundspeed(iph,rhl(iph,i,j,k),pl(iph,i,j,k)))
   !enddo
   !c(i,j,k)=c1

end function soundspeed_mixt_frozen

real(PR) function pressure_mixt_sge(MF)

   !==================================
   ! pression d'un mélange de fluides suivant
   ! une EOS de type SGE (p1692 de Saurel 2009)
   ! p=( rho*e - sum(alpha*gam*pinf/(gam-1)) ) / sum( alpha/(gam-1) )
   !=================================
   implicit none 
   type(multifluide), intent(in) :: MF
   real(PR) :: C1,denom
   integer :: iph

   pressure_mixt_sge=MF%rh*MF%u
   denom=0.0_PR

   do iph=1,Nl
     C1=MF%F(iph)%f/(MF%F(iph)%g_sge-1.0_PR)
     pressure_mixt_sge=pressure_mixt_sge-C1*MF%F(iph)%g_sge*MF%F(iph)%p_sge
     denom=denom+C1
   enddo

   pressure_mixt_sge=pressure_mixt_sge/denom

   if(pressure_mixt_sge.le.0)then

     print*, '********************************************'
     print*, 'PROBLEM PRESSURE MIXT SGE <0'
     if(Nl.ge.2)then
       print*, 'rho*u=',MF%rh*MF%u, 'rho*ele=', MF%rh*MF%ee, 'rhoE=', MF%rh*MF%E
       print*, 'u=',MF%u, 'ele=', MF%ee, 'E=', MF%E
       print*, 'sum(flelh)=', sum(MF%F(1:Nl)%f*MF%F(1:Nl)%eh)
       do iph=1,Nl
        print*, 'phase', iph, '--------------------------' 
        print*, ' fl/(g-1):',MF%F(iph)%f/(MF%F(iph)%g_sge-1.0_PR)
        print*, ' fl/(g-1)*gp:', MF%F(iph)%f/(MF%F(iph)%g_sge-1.0_PR)*MF%F(iph)%g_sge*MF%F(iph)%p_sge
        print*, ' F%rhl,fl,pl=', MF%F(iph)%rh, MF%F(iph)%f, MF%F(iph)%p
        print*, ' F%ele=', MF%F(iph)%ee, 'F%elh=', MF%F(iph)%eh
       enddo
       stop
     endif
     print*, '********************************************'

     err_relax=1

   endif

end function pressure_mixt_sge

real(PR) function pressure_mixt_sgeT(MF)

   !==================================
   ! pression d'un mélange de fluides suivant
   ! une EOS de type SGET (cf. Pepitas 2009 p 753 eq 19)
   ! p=( rho*(e-sum(Yq)) - sum(alpha*gam*pinf/(gam-1)) ) / sum( alpha/(gam-1) )
   !=================================
   implicit none 
   type(multifluide), intent(in) :: MF
   real(PR) :: C1,denom
   integer :: iph

   pressure_mixt_sgeT=MF%rh*(MF%u-sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%q))
   denom=0.0_PR

   do iph=1,Nl
     C1=MF%F(iph)%f/(MF%F(iph)%g_sge-1.0_PR)
     pressure_mixt_sgeT=pressure_mixt_sgeT-C1*MF%F(iph)%g_sge*MF%F(iph)%p_sge
     denom=denom+C1
   enddo

   pressure_mixt_sgeT=pressure_mixt_sgeT/denom

   if(pressure_mixt_sgeT.le.0)then

     print*, '********************************************'
     print*, 'PROBLEM PRESSURE MIXT SGET <0'
     if(Nl.ge.2)then
       print*, 'rho*u=',MF%rh*MF%u, 'rho*ele=', MF%rh*MF%ee, 'rhoE=', MF%rh*MF%E
       print*, 'u=',MF%u, 'ele=', MF%ee, 'E=', MF%E
       print*, 'rho*u-sumYq=', MF%rh*(MF%u-sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%q))
       print*, 'sum(flelh)=', sum(MF%F(1:Nl)%f*MF%F(1:Nl)%eh)
       do iph=1,Nl
        print*, 'phase', iph, '--------------------------' 
        print*, ' fl/(g-1):',MF%F(iph)%f/(MF%F(iph)%g_sge-1.0_PR)
        print*, ' fl/(g-1)*gp:', MF%F(iph)%f/(MF%F(iph)%g_sge-1.0_PR)*MF%F(iph)%g_sge*MF%F(iph)%p_sge
        print*, ' F%rhl,fl,pl=', MF%F(iph)%rh, MF%F(iph)%f, MF%F(iph)%p
        print*, ' F%ele=', MF%F(iph)%ee, 'F%elh=', MF%F(iph)%eh
       enddo
       stop
     endif
     print*, '********************************************'

     err_relax=1

   endif

end function pressure_mixt_sgeT

subroutine sortie_relax(MF)


   implicit none
   type(multifluide), intent(in) :: MF
   integer :: i
   real(PR) :: p(1:1), p0, pfin, dp, Jaco(1:1,1:1), func(1:1)

    print*, '' 
    print*, '============================================'
    print*, ''
    print*, '            PROBLEM RELAXATION              '
    print*, ''
    print*, ' err_relax=', err_relax
    if(err_relax.eq.1)then
    print*, ' Pmixt < 0 :(' 
    elseif(err_relax.eq.2)then
    print*, ' sum(Ylelh) .ne. ehtot x(' 
    elseif(err_relax.eq.3)then
    print*, ' it Newton > 10000 :p'
    endif
    print*, '============================================'
    print*, ''
    print*, 'Nl=',Nl
    do i=1,Nl
       print*, 'phase:',i, trim(ph(i)%nom)  
    enddo
    print*, ''
    print*, 'itérations Newton:', Nit_Newton
    print*, ''
    print*, '------------AVANT NEWTON:--------------'
    print*, ''
    print*, '---variables phasiques:'
    print*, 'fl='   , MF0%F(1:Nl)%f
    print*, 'pl='   , MF0%F(1:Nl)%p
    print*, 'rhl='  , MF0%F(1:Nl)%rh
    print*, 'Yl='   , MF0%F(1:Nl)%Y
    print*, 'pinf=' , MF0%F(1:Nl)%p_sge
    print*, 'gamma=', MF0%F(1:Nl)%g_sge
    print*, 'elh='  , MF0%F(1:Nl)%eh
    print*, 'ele='  , MF0%F(1:Nl)%ee
    print*, 'sum(alelh)=',sum(MF0%F(1:Nl)%Y*MF0%F(1:Nl)%eh)
    print*, 'sum(alele)=',sum(MF0%F(1:Nl)%Y*MF0%F(1:Nl)%ee)

    print*, '---variables de mélange:'
    print*, 'rh=',MF0%rh
    print*, 'vx=',MF0%vx
    print*, 'u=' ,MF0%u
    print*, 'ee=',MF0%ee
    print*, 'e=' ,MF0%e
    print*, 'P=' ,MF0%p

    print*, ''
    print*, '------------APRES NEWTON:--------------'
    print*, ''
    print*, '---variables phasiques modifiées:'
    print*, 'fl=' ,MF01%F(1:Nl)%f
    print*, 'pl=' ,MF01%F(1:Nl)%p
    print*, 'rhl=',MF01%F(1:Nl)%rh
    print*, 'Yl=' ,MF01%F(1:Nl)%Y
    print*, 'elh=',MF01%F(1:Nl)%eh
    print*, 'ele=',MF01%F(1:Nl)%ee
    print*, ''
    print*, '---nouvelle pression:'
    print*, 'P=',PNewton
    print*, ''
    print*, '---vérif conservation masse:'
    print*, 'sum(fl*rhl)0=', sum(MF0%F(1:Nl)%f*MF0%F(1:Nl)%rh)
    print*, 'sum(fl*rhl) =', sum(MF01%F(1:Nl)%f*MF01%F(1:Nl)%rh)
    print*, ''
    print*, '---vérif Yl:'
    print*, 'sum(Yl)0=', sum(MF0%F(1:Nl)%Y)
    print*, 'sum(Yl) =', sum(MF01%F(1:Nl)%Y)
    print*, ''
    print*, '---vérif fl:'
    print*, 'sum(Yl)0=', sum(MF0%F(1:Nl)%Y)
    print*, 'sum(Yl) =', sum(MF01%F(1:Nl)%Y)
    print*, ''
    print*, '---vérif pl:'
    do i=1,Nl
    print*, '  pression_sge(',i,')=', ph(i)%pressure(MF01%F(i),MF01%F(i)%rh,MF01%F(i)%eh) 
    enddo
    print*, ''
    print*, '---vérif conservation entalpie:'
    print*, 'elh0+newp/rh0=', MF0%F(1:Nl)%eh+PNewton/MF0%F(1:Nl)%rh
    print*, 'elh+newp/rh='  , MF01%F(1:Nl)%eh+PNewton/MF01%F(1:Nl)%rh
     
    print*, ''
    print*, '------------APRES PMIXT:--------------'
    print*, ''
    print*, 'MF%p=',MF%p
    print*, ''
    print*, '--------------------------------------'
    print*, ''
    print*, ' > sortie vers testrelax.dat'
    print*, ''
    print*, '============================================'


   open(unit=12, file='testrelax.dat', status='replace')

   p0=PNewton-10.0_PR*PNewton   !!-1.0e6_PR
   pfin=PNewton+10.0_PR*PNewton  !!!1.0e6_PR

   dp=(pfin-p0)/real(10000,PR)   

   do i=1,10000

      p(1)=p0+real(i-1,PR)*dp
      Jaco(1,1)=0.0_PR

      call F_relaxp_sge_1(1,p,func)
      call J_relaxp_sge_1(1,p,Jaco)

      write(12,*) p(1), func(1), Jaco(1,1)

   enddo

   close(12)

end subroutine sortie_relax


!!!========================== SUBROUTINES DE MECANIQUE =========================
!!!---0D

subroutine MAJ_meca(MF) !!!mu,rho0) !!!!rhl,fl,al,bl,cl,pl,sigl,sig,ele)

   implicit none
   type(multifluide), intent(inout) :: MF
   integer :: i,j,k,iph

   MF%sig=0.0_PR

   do iph=1,Nl


      MF%F(iph)%sig(1:3,1:3)=0.0_PR

      call init_Finger(MF%F(iph)) 

      MF%F(iph)%detG=detFT

      !!!---partie deviatotique du tenseur des contraintes
      call Deviateur(MF%F(iph)%sig(1:3,1:3))

      !!!---partie hydrostatique
      do i=1,3
         MF%F(iph)%sig(i,i)=MF%F(iph)%sig(i,i)-MF%F(iph)%p
      enddo

      !!!---Energie élastique
      MF%F(iph)%ee=energie_elastique(MF%F(iph)%f,MF%F(iph)%mu,MF%F(iph)%rh0)

   enddo

   !!!---tenseur des contraintes totales
   forall(i=1:3,j=1:3) MF%sig(i,j)=sum(MF%F(1:Nl)%f*MF%F(1:Nl)%sig(i,j))

   !!!---energie elastique totale
   MF%ee=sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%ee)

end subroutine MAJ_meca

!subroutine phase_change(F,ischange) !!!mu,rho0) !!!!rhl,fl,al,bl,cl,pl,sigl,sig,ele)
!
!   implicit none
!   type(fluide), intent(inout) :: F
!   logical :: ischange
!   integer :: i,j,k,iph,iph2,im
!   real(PR) :: nu,nu1,nu2,p12,Y1,Y2,a1,a2,nu1,nu2
!   real(PR) :: rh2,elh2,Y12,rh12,a12,elh12,ax12,ax2,Yx12,Yx2
!   real(PR) :: rho_exp, u_exp, p_exp
!   type(fluide) :: F2
!   integer :: transit(1:2) 
!   
!
!   ischange=.false.
!   transit=0
!   transit(2)=3
!
!   do iph=1,Nl
! 
!      iph2=transit(iph)
!
!
!      if(iph2.gt.0)then
!
!          !!!---transition iph -> iph2
!
!           rho=F%rhl(iph)
!           Ener=F%elh(iph)
!           im=F%im(iph)
!           call get_index_rE(im,rho,Ener,ir,iT)
!
!           !P=MF1D(i)%P 
!           !call get_index_rP(iph,rho,P,ir,iT)
!           !Temp=T_from_rP(iph,rho,P)
!
!           Temp=T_from_ru(im,rho,Ener)
!
!
!           if(Temp.lt.ph(im)%Tsatmax.and.Temp.gt.ph(im)%Tsatmin)then
!
!              nu=1.0_PR/F%rhl(iph)
!              
!              call get_indexsat_T(im,T,iT)
!
!              !call get_sat(iph,F%elh(iph),nu1,nu2)
!
!              call nusat_from_T(im,Temp,nu1,nu2)
!              call pinfsat_from_T(im,Temp,pinf1,pinf2)
!              call gsat_from_T(im,Temp,g1,g2)
!              call Esat_from_T(im,Temp,e1,e2)
!
!              if(nu.gt.nu1.and.nu.lt.nu2)then
!
!                 ischange=.true.
!           
!                 !!!--propriété de phase 1+2
!                 a12=F%fl(iph)+F%fl(iph2)
!                 Y12=F%Yl(iph)+F%Yl(iph2)
!                 rh12=F%rhl(iph)
!                 !!!---propriétés de la phase 2
!                 F%Tl(iph2)=Temp
!                 F%Yl(iph2)=(nu-nu1)/(nu2-nu1)*Y12
!                 F%rhl(iph2)=1.0_PR/nu2 
!                 F%fl(iph2)=rh12/rh2*Y2
!                 F%p_sge(iph2)=pinf2
!                 F%g_sge(iph2)=g2
!                 F%elh(iph2)=e2
!                 F%pl(iph2)=ph(im)%pressure(iph2,F,F%rhl(iph2),F%elh(iph2))
!
!                 !!!---propriétés de la phase 1
!                 F%Tl(iph)=Temp
!                 F%Yl(iph)=(nu2-nu)/(nu2-nu1)*Y12
!                 F%rhl(iph)=1.0_PR/nu1 
!                 F%fl(iph)=rh12/rh1*Y1
!                 F%p_sge(iph)=pinf1
!                 F%g_sge(iph)=g1
!                 F%elh(iph)=e1
!                 F%pl(iph)=ph(im)%pressure(iph,F,F%rhl(iph),F%elh(iph))
!
!              endif
!
!
!          endif
!
!      endif
!
!   enddo
!
!end subroutine phase_change

!subroutine get_sat(iph,u,nu1,nu2)
!   implicit none
!   integer, intent(in) :: iph
!   real(PR), intent(in) :: u
!   real(PR), intent(out) :: nu1,nu2
!   real(PR) :: T
!
!   !if(T.lt.2500.0_PR)then
!
!   !    nu1=0.0_PR
!   !    nu2=0.0_PR
!
!   !else
!
!   !nu1=1.0_PR/(1.01_PR*rho0_Cu)
!   !nu2=1.0_PR/(0.001_PR*rho0_Cu)
!   !p12=1.0_PR/nu2*8.314_PR/63.0e-3_PR*T
!
!   !endif
!
!   T=u/cv_Cu
!   nu1=1.0_PR/(0.8_PR*rho0_Cu)
!   nu2=1.0_PR/(0.001_PR*rho0_Cug)
!   !p12=pressure_sge(g_sge_Cug,p_sge_Cug,rho0_Cug,u)
!
!end subroutine get_sat





real(PR) function energie_elastique(fl, mu,rho0)

   !!!--------------------------------------
   !!! retourne l'énergie elastique
   !!! Nécessite d'appeler init_Finger
   !!!--------------------------------------

   implicit none
   real(PR), intent(in) :: fl, mu, rho0
   real(PR) :: G(1:3,1:3), G2(1:3,1:3)
   integer :: i,j,k

   if(abs(detFT).gt.toldet.and.mu.gt.tolmu)then

      !!!--- EOS (46) Favrie 2009 p 6053

      !energie_elastique=alpharhom1*0.25_PR*mu*detFTm16*(&
      !                  J2-2.0_PR*detFT13*J1+3.0_PR*detFT23)


      !!!--- EOS Favrie 2009 p 6041 
      !   do j=1,3 
      !     do i=1,3
      !       G(i,j)=FT(i,j)*detFTm13
      !     enddo
      !   enddo
 
      !  !!!!---tenseur (g-I)**2 : g2
      !   G(1:3,1:3)=G(1:3,1:3)-Id(1:3,1:3)
      !   do j=1,3
      !     do i=1,3
      !        G2(i,j)=sum(G(i,1:3)*G(1:3,j))   
      !     enddo
      !   enddo
    
      !   energie_elastique=0.25_PR*mu/rho0*(G2(1,1)+G2(2,2)+G2(3,3))


      !!!---Favrie 2009 p 6041 :
      !energie_elastique=0.25_PR*mu/rho0*(J2-2.0_PR*J1+3.0_PR)
 
      !G(1:3,1:3)=FT(1:3,1:3)-Id(1:3,1:3)    
      !do j=1,3
      !  do i=1,3
      !     G2(i,j)=sum(G(i,1:3)*G(1:3,j))   
      !  enddo
      !enddo 

      !energie_elastique=0.25_PR*mu/rho0*(G(1,1)+G(2,2)+G(3,3))
 
      !!!--- EOS (5) Favrie 2015 p 526
      !energie_elastique=0.125_PR*mu/rho0*(J2-3.0_PR)


       !!!---mine
      !energie_elastique=0.25_PR*mu/rho*detFT12*(J2-2.0_PR*J1+3.0_PR)


      !!!---EOS 2015 p539
      energie_elastique=0.125_PR*mu/rho0*( (a(1)**4+b(2)**4+c(3)**4)/(a(1)*b(2)*c(3))**(4.0_PR/3.0_PR) - 3.0_PR ) 



   else

      energie_elastique=0.0_PR
    
   endif

   !energie_elastique=0.0_PR

end function energie_elastique

subroutine Deviateur(D)

    !!!--------------------------------------
    !!! retourne le deviateur de contraintes
    !!! Nécessite d'appeler init_Finger
    !!!--------------------------------------
    
    implicit none
    real(PR), dimension(1:3,1:3), intent(out) :: D
    integer :: i,j,k
    real(PR) :: tr
    
    IF(murr0.gt.0.0_PR)THEN
    
        if(abs(detFT).le.toldet)then
        
            D=0.0_PR
        
        else
        
            do j=1,3
                do i=1,3
                
                    !!!---Fravrie 2015 p 6053:
                    
                     !D(i,j)=-alpham1*mu*detFTm16*( FT2(i,j) - r13*J2*Id(i,j)&
                     !              -detFT13*(FT(i,j)  - r13*J1*Id(i,j)) )
                    
                    !D(i,j)=-murr0*( detFTm23*(FT2(i,j) - r13*J2*Id(i,j))&
                    !               -detFTm13*(FT(i,j)  - r13*J1*Id(i,j)) )
                    
                    !!!---Favrie 2009 p 6041 :
                    !D(i,j)=-murr0*( (FT2(i,j) - r13*J2*Id(i,j))&
                    !              -(FT(i,j)  - r13*J1*Id(i,j)) )
                    
                    !!!---Fravrie 2015 p 527:
                    !D(i,j)=-0.5_PR*murr0*( FT2(i,j) - r13*J2*Id(i,j) )
                    
                    
                    !!!---mine
                    !D(i,j)=-mu*detFT12*( (FT2(i,j) - r13*J2*Id(i,j))&
                    !               -(FT(i,j)  - r13*J1*Id(i,j)) )
                    
                    if(i.eq.1.and.j.eq.1)then
                    
                       D(i,j)=-0.5_PR*mu/(a(1)*b(2)*c(3))**(1.0_PR/3.0_PR)*(&
                       a(1)**4-1.0_PR/3.0_PR*(a(1)**4+b(2)**4+c(3)**4))
                    
                    else
                    
                        D(i,j)=0.0_PR
                    
                    endif
                
                enddo
            enddo       
            
            !tr=D(1,1)+D(2,2)+D(3,3)
            
            !if(tr.ne.0.0_PR)then
            !   print*, 'PROBLEM DEVIATEUR: tr=', tr
            !   print*, '----------------------------'
            !   print*, D(1,1), D(1,2), D(1,3)
            !   print*, D(2,1), D(2,2), D(2,3)
            !   print*, D(3,1), D(3,2), D(3,3)
            !   print*, '----------------------------'
            !   stop
            !endif 
        
        endif
    
    ELSE
    
        D=0.0_PR
    
    ENDIF !!! mu>0

end subroutine Deviateur




!!!================= RELAXTATION =========================


subroutine relaxation !!!!!(Nl,N,F) !!!,U,G)

   implicit none
   !integer, intent(in) :: Nl,N
   !type(fluide), intent(inout) :: F(1:N)
   !!! real(PR), intent(inout) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)

   !!!---champs 1D locaux
   !real(PR), dimension(Nl,N) :: fl, Yl, pl, rhl, elh, ele
   !real(PR), dimension(Nl,3,N) :: al,bl,cl
   !real(PR), dimension(Nl,3,3,N) :: sigl
   !real(PR), dimension(3,3,N) :: sig
   !real(PR), dimension(N) :: rh, vx, vy,vz, z, c, e, p
       
   integer :: i,j,k,iph,ir,iT
   logical :: ischange
   real(PR) :: rho, Ener, Temp, p_sge, g_sge, new_p, P, P2

   call cpu_time(t1)

   !print*, 'EQUILIBRAGE DE LA PRESSION 1 '

   err_relax=0

   do i=1,N1D

       call relax_p(MF1D(i))

       if(err_relax.ne.0)then

          call sortie_relax(MF1D(i))

          stop

       endif

    enddo


   !!!!!!-------------CHANGEMENT DE PHASE :

   !call phase_change(MF1D(i),ischange)       


   !!!!---------------------- MODIFICATION DE p_inf et gamma !


   !!!   !print*, 'MODIFICATION DU FLUIDE :'
   !!!   !print*, '----Before:'
   !!!   !do i=1,N1D
   !!!   !  MF1D(i)%B=MF1D(i)%g_sge(2) 
   !!!   !enddo
   !!!   !print*, '   max gamma2=', maxval(MF1D(1:N1D)%B) !g_sge(2))
   !!!   !print*, '   min gamma2=', minval(MF1D(1:N1D)%B) !g_sge(2))
   !!!   !do i=1,N1D
   !!!   !  MF1D(i)%B=MF1D(i)%p_sge(2) 
   !!!   !enddo
   !!!   !print*, '   max pinf2=' , maxval(MF1D(1:N1D)%B) !p_sge(2))
   !!!   !print*, '   min pinf2=' , minval(MF1D(1:N1D)%B) !p_sge(2))

   !!!   plneg=.false.

   !!!   do i=1,N1D

   !!!     do iph=1,Nl

   !!!        if(ph(iph)%typ.eq.2)then

   !!!           rho=MF1D(i)%rhl(iph)

   !!!           Ener=MF1D(i)%elh(iph)
   !!!           call get_index_rE(iph,rho,Ener,ir,iT)

   !!!           !P=MF1D(i)%P 
   !!!           !call get_index_rP(iph,rho,P,ir,iT)
   !!!           !Temp=T_from_rP(iph,rho,P)

   !!!           Temp=T_from_ru(iph,rho,Ener)

   !!!           MF1D(i)%Tl(iph)=Temp

   !!!           p_sge=pinf_from_rT(iph,rho,Temp)

   !!!           g_sge=g_from_rT(iph,rho,Temp)

   !!!           P=P_from_rT(iph,rho,Temp)
   !!!           P2=rho*(g_sge - 1.0_PR)*Ener - g_sge*p_sge

   !!!           if( abs((P2-P)/(0.5_PR*(P+P2))) .gt.1.0e-1_PR )then
   !!!           print*, 'P ne P2:', P, P2
   !!!            print*, 'old pinf,g :',  MF1D(i)%p_sge(iph), MF1D(i)%g_sge(iph)
   !!!            print*, 'new pinf,g :',  p_sge, g_sge
   !!!            print*, 'old p:',MF1D(i)%pl(iph)
   !!!            print*, 'old p verif:', rho*(MF1D(i)%g_sge(iph) - 1.0_PR)*Ener - MF1D(i)%g_sge(iph)*MF1D(i)%p_sge(iph)
   !!!            !print*, 'new_p:', new_p
   !!!            print*, 'old E:', Ener
   !!!            print*, 'verif E:',E_from_rT(iph,rho,Temp) 



   !!!           stop
   !!!           endif

   !!!           MF1D(i)%pl(iph)=P

   !!!           MF1D(i)%p_sge(iph)=p_sge
   !!!           MF1D(i)%g_sge(iph)=g_sge
   !!!           MF1D(i)%cvl(iph)=Cv_from_rT(iph,rho,Temp)
   !!!           MF1D(i)%ZTF(iph)=Z_from_rT(iph,rho,Temp)


 
   !!!           !MF1D(i)%pl(iph)=ph(iph)%pressure(iph,MF1D(i),MF1D(i)%rhl(iph),MF1D(i)%elh(iph))
  
   !!!           !Ener=MF1D(i)%elh(iph)

   !!!           !new_p=rho*(g_sge - 1.0_PR)*Ener - g_sge*p_sge

   !!!           !Ener=(P+g_sge*p_sge)/(rho*(g_sge-1.0_PR)) 

   !!!           !if(new_p.lt.0.0_PR)then
   !!!           ! print*, 'neg pressure after pinf modif'
   !!!           ! print*, 'i=', i, 'iph=', iph
   !!!           
   !!!           ! plneg=.true.

   !!!           ! stop
   !!!           !endif 
 
   !!!           !MF1D(i)%pl(iph)=new_p

   !!!           !MF1D(i)%elh(iph)=Ener
   !!!           !MF1D(i)%p_sge(iph)=p_sge
   !!!           !MF1D(i)%g_sge(iph)=g_sge
   !!!           !MF1D(i)%cvl(iph)=Cv_from_rT(iph,rho,Temp)
   !!!           !MF1D(i)%ZTF(iph)=Z_from_rT(iph,rho,Temp)

   !!!           !MF1D(i)%u=sum(MF1D(i)%Yl(:)*MF1D(i)%elh(:))

   !!!           !MF1D(i)%e=MF1D(i)%u+0.5_PR*sqrt(MF1D(i)%vx**2+MF1D(i)%vy**2+MF1D(i)%vz**2)
   !!!   
   !!!        elseif(ph(iph)%typ.eq.1)then

   !!!           if(MF1D(i)%Tl(iph).gt.ph(iph)%Tvap)then

   !!!                phasechange=.true.

   !!!                MF1D(i)%p_sge(iph)=ph(iph)%p_sge*exp(-(MF1D(i)%Tl(iph)-ph(iph)%Tvap)/ph(iph)%Tvap)      

   !!!                MF1D(i)%pl(iph)=ph(iph)%pressure(iph,MF1D(i),MF1D(i)%rhl(iph),MF1D(i)%elh(iph))

   !!!                !MF1D(i)%elh(iph)=ph(iph)%energie_hydro(iph,MF1D(i),MF1D(i)%rhl(iph),MF1D(i)%pl(iph))

   !!!           endif



   !!!        endif

   !!!      enddo

   !!!   enddo

   !!!   if(phasechange) print*, 'CHANGEMENT DE PHASE !!!'
   !!!   phasechange=.false.


   !!!   !print*, '----After:'
   !!!   !do i=1,N1D
   !!!   !  MF1D(i)%B=MF1D(i)%g_sge(2) 
   !!!   !enddo
   !!!   !print*, '   max gamma2=', maxval(MF1D(1:N1D)%B) !g_sge(2))
   !!!   !print*, '   min gamma2=', minval(MF1D(1:N1D)%B) !g_sge(2))
   !!!   !do i=1,N1D
   !!!   !  MF1D(i)%B=MF1D(i)%p_sge(2) 
   !!!   !enddo
   !!!   !print*, '   max pinf2=' , maxval(MF1D(1:N1D)%B) !p_sge(2))
   !!!   !print*, '   min pinf2=' , minval(MF1D(1:N1D)%B) !p_sge(2))



   !!!   !!!!------------------- EQUILIBRAGE DE LA PRESSION 2 ------------------------

   !!!   !!print*, 'minval p=', minval(p(:)), minloc(p,1) 
   !!!   !!if(ierr.ne.0) call crash('relaxp1')

   !!!   !!print*, 'EQUILIBRAGE DE LA PRESSION 2 '

   !!!   do i=1,N1D

   !!!       call relax_p(MF1D(i))

   !!!       if(pmixtlt0)then
   !!!          print*, 'equ pression2 : p<0 at i=', i
   !!!          stop
   !!!       endif

   !!!    enddo


   !!!------------------ CHANGEMENTS DE PHASE !!! --------------------------

   !!! do i=1,N1D


   !!!       !if(i.eq.123.and.it.eq.2123)then

   !!!       ! print*, '=============================='
   !!!       ! print*, 'fl=', MF1D(i)%fl(1), MF1D(i)%fl(2) 
   !!!       ! print*, 'Yl=', MF1D(i)%Yl(1), MF1D(i)%Yl(2) 
   !!!       ! print*, 'rhl=', MF1D(i)%rhl(1), MF1D(i)%rhl(2) 
   !!!       ! print*, 'pl=', MF1D(i)%pl(1), MF1D(i)%pl(2) 
   !!!       ! print*, 'elh=', MF1D(i)%elh(1), MF1D(i)%elh(2) 
   !!!       ! print*, '=============================='

   !!!       !endif




   !!!    call phase_change(MF1D(i),ischange)       

   !!!    if(ischange)then
   !!!       !print*, 'change i=', i


   !!!       !if(i.eq.123)then

   !!!       ! print*, '=============================='
   !!!       ! print*, 'fl=', MF1D(i)%fl(1), MF1D(i)%fl(2) 
   !!!       ! print*, 'Yl=', MF1D(i)%Yl(1), MF1D(i)%Yl(2) 
   !!!       ! print*, 'rhl=', MF1D(i)%rhl(1), MF1D(i)%rhl(2) 
   !!!       ! print*, 'pl=', MF1D(i)%pl(1), MF1D(i)%pl(2) 
   !!!       ! print*, 'elh=', MF1D(i)%elh(1), MF1D(i)%elh(2) 
   !!!       ! print*, '=============================='

   !!!       !endif

   !!!       call relax_p(MF1D(i))
   !!!       if(pmixtlt0)then
   !!!          print*, 'pression2 : p<0 at i=', i
   !!!          stop
   !!!       endif

   !!!       MF1D(i)%c=soundspeed_mixt(MF1D(i))

   !!!    endif

   !!! enddo

   !!!------------------- EQUILIBRAGE DE LA PRESSION 2------------------------

   !print*, 'minval p=', minval(p(:)), minloc(p,1) 

   ! if(ierr.ne.0) call crash('relaxp1')

   !do i=1,N1D

   !    call relax_p(MF1D(i))
   !    if(pmixtlt0)then
   !       print*, 'pression2 : p<0 at i=', i
   !       stop
   !    endif

   ! enddo


   !!!----------------- MAJ ---------------------------------------

   !print*, 'minval p2=', minval(p(:)), minloc(p,1) 

   !if(ierr.ne.0) call crash('relaxp2')

   do i=1,N1D
      call MAJ_meca(MF1D(i)) 
      call MAJ_mixt(MF1D(i)) 
      call MAJ_thermo(MF1D(i))
   enddo

!   do i=1,N1D 
!    call primitive2conservative(MF1D(i),U(:,i),G(:,:,i)) 
!   enddo
!   if(ierr.ne.0) call crash('relaxp3')

   call cpu_time(t2)

   t_relaxp=t_relaxp+t2-t1
 
end subroutine relaxation


!!!================ SUBROUTINE DE RELAXATION DE LA PRESSION ====================

subroutine init_relaxp

     implicit none

     IF(Nl.eq.1)THEN

        relax_p => relax_p0
        write(iout,*) '   relax_p => relax_p0' 

     ELSE

     if(allocated(xin))    deallocate(xin)  
     if(allocated(bin))    deallocate(bin)   
     if(allocated(p_sge0)) deallocate(p_sge0)
     if(allocated(g_sge0)) deallocate(g_sge0)
     if(allocated(flrhl0)) deallocate(flrhl0)
     if(allocated(el0))    deallocate(el0)
     if(allocated(nul0))   deallocate(nul0)
     if(allocated(pl0))    deallocate(pl0)
     if(allocated(Zl))     deallocate(Zl)
     if(allocated(fl0_gl))  deallocate(fl0_gl)
     if(allocated(fl0))    deallocate(fl0)
     if(allocated(coef1))    deallocate(coef1)

     allocate(p_sge0(1:Nl)); p_sge0(:)=0.0_PR
     allocate(g_sge0(1:Nl)); g_sge0(:)=0.0_PR
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
           write(iout,*) '     > relax_p    => relax_p_sge_1' 
           write(iout,*) '     > F_function => F_relax_p_sge_1' 
        elseif(mode_pI.eq.2)then
           relax_p => relax_p_sge_2
          ! call init_SNES(Nsnes,xin,bin,F_function=F_relaxp_sge_2,J_function=J_relaxp_sge_2)
           write(iout,*) '     > relax_p    => relax_p_sge_2' 
           write(iout,*) '     > F_function => F_relax_p_sge_2' 
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

subroutine relax_p0(MF)

     implicit none
     type(multifluide), intent(inout) :: MF
     real(PR) :: pmin, pmax
     integer :: iph

     MF%p=pressure_mixt(MF) 

     do iph=1,Nl
        MF%F(iph)%p=MF%p
     enddo

     do iph=1,Nl
       MF%F(iph)%eh=ph(iph)%energie_hydro(MF%F(iph),MF%F(iph)%rh,MF%p)
     enddo

end subroutine relax_p0

!subroutine relax_p_sge_bi(F)
!
!     implicit none
!     type(fluide), intent(inout) :: F 
!     real(PR) :: pmin, pmax, ehtot
!     !real(PR), allocatable :: xin(:), bin(:)
!     integer :: iph
!
!     !pmin=minval(pl(1:Nl))
!     !pmax=maxval(pl(1:Nl))
!     
!     !IF(0.5_PR*abs(pmax-pmin)/(pmax+pmin).gt.1.0e-3_PR)THEN
!
!    Nl0=Nl
!
!    !!!---Nombre de phases
!     do iph=1,Nl
!       if(ph(iph)%typ.eq.2)then
!         
!         Esatmax=maxval(mat(iph)%E_sat(1,:))
!         if(F%elh(iph).lt.mat(iph)%Esatmax)then 
!            call get_indexsat_E(iph,F%elh(iph),iTsat)
!            call get_nusat_E(iph,F%elh(iph),nusat1,nusat2)
!            nu=1.0_PR/F%rhl(iph)            
!            if(nu.gt.nu_sat1.and.nu.lt.nusat2)then
!               Nl0=Nl0+1
!                              
!               F%state=2
!
!            endif
!         endif
!       endif
!     enddo 
!
!
!    !!!---
!
!
!     Nl0=Nl
!     do iph=1,Nl
!       flrhl0(iph)=F%fl(iph)*F%rhl(iph)
!       fl0_gl(iph)=F%fl(iph)/F%g_sge(iph)
!       pl0(iph)=F%pl(iph)  
!       coef1(iph)=fl0_gl(iph)*(F%p_sge(iph)+pl0(iph))
!       fl0(iph)=F%fl(iph)
!       el0(iph)=F%elh(iph)
!       p_sge0(iph)=F%p_sge(iph)
!       g_sge0(iph)=F%g_sge(iph)
!     enddo 
!
!     xin(1)=maxval(F%pl(1:Nl))
!
!     bin(:)=0.0_PR
!
!     if(verb)then
!        write(iout,*) '==================================='
!        write(iout,*) '   x=', xin(1:Nsnes)
!        write(iout,*) '   rhl=', F%rhl(1:Nl)
!        write(iout,*) '   sumY=',sum(F%Yl(1:Nl)) 
!        write(iout,*) '   sumYelh=',sum(F%Yl(1:Nl)*F%elh(1:Nl)) 
!        write(iout,*) '   sumfl=',sum(F%fl(1:Nl)) 
!        write(iout,*) '   elh=', F%elh(1:Nl)
!        write(iout,*) '   p=', F%p
!        write(iout,*) '   pl=', F%pl(1:Nl)
!        write(iout,*) '   fl=', F%fl(1:Nl)
!        write(iout,*) '   rho=', F%rh
!     endif     
!
!     !call solve_SNES(Nsnes,xin,bin)
!
!     !call My_Newton(1,xin(1:1),minval(pl),maxval(pl),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)
!
!     !call My_Newton(1,xin(1:1),0.0_PR,maxval(pl),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)
!
!     call My_Newton(1,xin(1:1),1.0e-30_PR,maxval(F%pl),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)
!
!     F%p=xin(1)
!
!     !if(p.lt.0.0_PR)then
!     !    verb=.true.
!     !    print*, 'coucou1'
!     !endif
!
!     if(verb)then
!        write(iout,*) '----------- After Newton -----------'
!        write(iout,*) '   p=', F%p, ':<'
!        write(iout,*) '   e+p/rh=', F%elh(1:Nl)+F%p/F%rhl(1:Nl)
!     endif
!
!     do iph=1,Nl
!       F%fl(iph)=fl0_gl(iph)*(pl0(iph)+F%g_sge(iph)*F%p_sge(iph)+(F%g_sge(iph)-1.0_PR)*F%p)/(F%p_sge(iph)+F%p)
!       F%rhl(iph)=flrhl0(iph)/F%fl(iph)
!       F%Yl(iph)=F%fl(iph)*F%rhl(iph)/F%rh
!     enddo
!
!     !rh=sum(fl(1:Nl)*rhl(1:Nl))
! 
!     if(verb)then     
!        do iph=1,Nl
!         !if(ph(iph)%typ.eq.2) F%Tl(iph)=T_from_rP(iph,F%rhl(iph),F%p)
!         F%elh(iph)=ph(iph)%energie_hydro(iph,F,F%rhl(iph),F%p)
!         !!!!F%elh(iph)=energie_hydro_sge(F%g_sge(iph),F%p_sge(iph),F%rhl(iph),F%p)
!         F%pl(iph)=ph(iph)%pressure(iph,F,F%rhl(iph),F%elh(iph))
!        enddo
!        write(iout,*) '------------------------'
!        write(iout,*) '   x=', xin(1:Nsnes)
!        write(iout,*) '   rhl=', F%rhl(1:Nl)
!        write(iout,*) '   sumY=',sum(F%Yl(1:Nl)) 
!        write(iout,*) '   sumfl=',sum(F%fl(1:Nl)) 
!        write(iout,*) '   elh=', F%elh(1:Nl)
!        write(iout,*) '   pl=', F%pl(1:Nl)
!        write(iout,*) '   fl=', F%fl(1:Nl)
!        write(iout,*) '   rho=', F%rh
!        write(iout,*) '   e+p/rh=', F%elh(1:Nl)+F%p/F%rhl(1:Nl)
!        write(iout,*) '-------------------------'
!        write(iout,*) '---Set pressure from e---'
!     endif    
!
!        !ENDIF !!! pmin-pmax/(pmin+pmax)>1e-3
!
!     ehtot=F%u !!!!F%e-0.5_PR*(F%vx**2+F%vy**2+F%vz**2)
!
!     F%p=pressure_mixt(F)
! 
!     if(verb)then 
!        write(iout,*) 'u=', F%u
!        write(iout,*) 'P=', F%p
!     endif     
!
!
!     !if(p.lt.0.0_PR)then
!     !    verb=.true.
!     !    print*, 'coucou2'
!     !endif
!
!     do iph=1,Nl
!        F%pl(iph)=F%p
!     enddo
!
!     do iph=1,Nl
!       F%elh(iph)=ph(iph)%energie_hydro(iph,F,F%rhl(iph),F%p)
!       !F%elh(iph)=energie_hydro_sge(F%g_sge(iph),F%p_sge(iph),F%rhl(iph),F%p)
!     enddo
!
!     if(verb)then 
!        write(iout,*) 'elh1=', F%elh(1)
!        write(iout,*) 'elh2=', F%elh(2)
!     endif     
!
!
!     if(abs(sum(F%Yl(1:Nl)*F%elh(1:Nl))-ehtot)/ehtot.gt.1.0e-6_PR)then
!        write(iout,*) 'PROBLEM ehtot ne sum(Ylehl):'
!        write(iout,*) sum(F%Yl(1:Nl)*F%elh(1:Nl)), ehtot
!        stop
!     endif
!
! 
!     !F%c=soundspeed_mixt(F)
! 
!end subroutine relax_p_sge_bi



subroutine relax_p_sge_1(MF)

     implicit none
     type(multifluide), intent(inout) :: MF 
     real(PR) :: pmin, pmax, ehtot,Pmixt,coeff,t2,t1
     !real(PR), allocatable :: xin(:), bin(:)
     integer :: iph

     !pmin=minval(pl(1:Nl))
     !pmax=maxval(pl(1:Nl))
     
     !IF(0.5_PR*abs(pmax-pmin)/(pmax+pmin).gt.1.0e-3_PR)THEN

     MF0=MF

     Nl0=Nl
     do iph=1,Nl
       !flrhl0(iph)=F%fl(iph)*F%rhl(iph)
       fl0_gl(iph)=MF%F(iph)%f/MF%F(iph)%g_sge
       pl0(iph)=MF%F(iph)%p  
       !coef1(iph)=fl0_gl(iph)*(F%p_sge(iph)+pl0(iph))
       fl0(iph)=MF%F(iph)%f
       !el0(iph)=F%elh(iph)
       p_sge0(iph)=MF%F(iph)%p_sge
       g_sge0(iph)=MF%F(iph)%g_sge
     enddo 
 
     xin(1)=maxval(abs(MF%F(1:Nl)%p))
 
     bin(:)=0.0_PR

     call My_Newton(1,xin(1:1),minval(MF%F(:)%p),maxval(MF%F(:)%p),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)

     MF%p=xin(1)

     Pnewton=MF%p

     do iph=1,Nl
       coeff=1.0_PR/MF%F(iph)%f
       MF%F(iph)%f0=MF%F(iph)%f
       MF%F(iph)%f=fl0_gl(iph)*(pl0(iph)+MF%F(iph)%g_sge*MF%F(iph)%p_sge+(MF%F(iph)%g_sge-1.0_PR)*MF%p)/(MF%F(iph)%p_sge+MF%p)
       coeff=coeff*MF%F(iph)%f

       MF%F(iph)%p=MF%p

       !!!   F%rhl(iph)=flrhl0(iph)/F%fl(iph)
       !!!   !F%rhl(iph)=F%flrhl(iph)/F%fl(iph)

       !!!   F%Yl(iph)=F%fl(iph)*F%rhl(iph)/F%rh

       !!!   F%elh(iph)=ph(iph)%energie_hydro(iph,F,F%rhl(iph),F%p)

       !coeff=coeff**(-r13)
       !coeff=coeff**(-1.0_PR)
       !F%al(iph,1:3)=F%al(iph,1:3)*coeff
       !F%bl(iph,1:3)=F%bl(iph,1:3)*coeff
       !F%cl(iph,1:3)=F%cl(iph,1:3)*coeff 

       !!!   !F%al(iph,1)=F%rhl(iph)/F%rho0(iph)

     enddo

     MF01=MF

     !!!---Reset énergie interne 

     !!! ATTENTION: recalculer u à partir des nouveaux Yl*ele 

     !call MAJ_meca(F)

     !F%u=F%E-F%ee-0.5_PR*(F%vx**2+F%vy**2+F%vz**2)

     !ehtot=F%u 

     !Pmixt=pressure_mixt(F)

     !F%p=Pmixt

     !do iph=1,Nl
     !   F%pl(iph)=F%p
     !enddo

     !do iph=1,Nl
     !   F%elh(iph)=ph(iph)%energie_hydro(iph,F,F%rhl(iph),F%p)
     !enddo

     !if(abs(sum(F%Yl(1:Nl)*F%elh(1:Nl))-ehtot)/ehtot.gt.1.0e-6_PR)then

     !   print*, '********************************************'
     !   write(iout,*) 'PROBLEM ehtot ne sum(Ylehl):'
     !   write(iout,*) sum(F%Yl(1:Nl)*F%elh(1:Nl)), ehtot
     !   print*, '********************************************'
     !   err_relax=2
     !endif

end subroutine relax_p_sge_1

subroutine reset_internal_energy_1D

    implicit none
    integer :: i

    call cpu_time(t1)  

    do i=1,N1D
        call reset_internal_energy(U(1:Nl+4,i),MF1D(i))
    enddo
    do i=1,N1D
        if(abs(MF1D(i)%ee).gt.1.0e15_PR)then
            print*, 'PROBLEM u 3 i=', i
            stop
        endif
    enddo 

    call cpu_time(t2) ; t_reset=t_reset+t2-t1

end subroutine reset_internal_energy_1D

subroutine reset_internal_energy(U,MF)

    implicit none
    real(PR), intent(in) :: U(1:Nl+4)
    type(multifluide), intent(inout) :: MF
    real(PR) :: ehtot, Pmixt, coeff
    integer :: iph

    !call init_Finger(1,F) 
    !call cpu_time(time_1)

    MF%F(1:Nl)%rh=U(1:Nl)/MF%F(1:Nl)%f
    MF%rh=sum(U(1:Nl))
    MF%F(1:Nl)%Y=U(1:Nl)/MF%rh 
    
    !do iph=1,Nl
    !   call init_Finger(iph,F)
    !   coeff=(F%fl(iph)/F%fl0(iph))**(-r13)
    !   F%al(iph,1:3)=F%al(iph,1:3)*coeff
    !   F%bl(iph,1:3)=F%bl(iph,1:3)*coeff
    !   F%cl(iph,1:3)=F%cl(iph,1:3)*coeff
    !enddo
    
    do iph=1,Nl
       call init_Finger(MF%F(iph))
       coeff=(MF%F(iph)%rh/(MF%F(iph)%rh0*detFT12))**(r13)
       MF%F(iph)%a(1:3)=MF%F(iph)%a(1:3)*coeff
       MF%F(iph)%b(1:3)=MF%F(iph)%b(1:3)*coeff
       MF%F(iph)%c(1:3)=MF%F(iph)%c(1:3)*coeff
    !   !coeff=(F%rhl(iph)/(F%rho0(iph)*detFT12))
    !   !F%al(iph,1:3)=F%al(iph,1:3)*coeff
    enddo
    
    !call MAJ_meca(F)
    
    MF%vx=U(Nl+1)/MF%rh
    MF%vy=U(Nl+2)/MF%rh
    MF%vz=U(Nl+3)/MF%rh
    MF%E=U(Nl+4)/MF%rh
    MF%p=sum(MF%F(1:Nl)%f*MF%F(1:Nl)%p)
    
    MF%u=MF%E-MF%ee-0.5_PR*(MF%vx**2+MF%vy**2+MF%vz**2)
    
    ehtot=MF%u 
    
    Pmixt=pressure_mixt(MF)
    
    MF%p=Pmixt
    
    do iph=1,Nl
       MF%F(iph)%p=MF%p
    enddo
    
    do iph=1,Nl
       MF%F(iph)%eh=ph(iph)%energie_hydro(MF%F(iph),MF%F(iph)%rh,MF%p)
    enddo
    
    ehtot=sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%eh)
    
    if(abs(ehtot-MF%u)/(ehtot+MF%u).gt.1.0e-3_PR)then
       print*, 'PROBLEM ehtot:'
       print*, 'ehtot=', ehtot, 'MF%u=', MF%u
       stop
    endif
    
    call MAJ_meca(MF)
    call MAJ_mixt(MF)

end subroutine reset_internal_energy

subroutine relax_p_sge_2(MF)

     implicit none
     type(multifluide), intent(inout) :: MF
     real(PR) :: pmin, pmax, ehtot
     !real(PR), allocatable :: xin(:), bin(:)
     integer :: iph

     pmin=minval(MF%F(1:Nl)%p)
     pmax=maxval(MF%F(1:Nl)%p)
     
     !IF(0.5_PR*abs(pmax-pmin)/(pmax+pmin).gt.1.0e-3_PR)THEN


     Nl0=Nl
     do iph=1,Nl
       flrhl0(iph)=MF%F(iph)%f*MF%F(iph)%rh
       fl0_gl(iph)=MF%F(iph)%f/MF%F(iph)%g_sge
       pl0(iph)   =MF%F(iph)%p  
       coef1(iph) =fl0_gl(iph)*(MF%F(iph)%p_sge+pl0(iph))
       fl0(iph)   =MF%F(iph)%f
       el0(iph)   =MF%F(iph)%eh
       p_sge0(iph)=MF%F(iph)%p_sge
       g_sge0(iph)=MF%F(iph)%g_sge
     enddo 
 
     do iph=1,Nl
      if(MF%F(iph)%p.gt.0.0_PR)then
      Zl(iph)=MF%F(iph)%rh*ph(iph)%soundspeed(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%p)
      else
      Zl(iph)=0.0_PR
      endif
     enddo
     pI0=sum(Zl(1:Nl)*MF%F(1:Nl)%p)/sum(Zl(1:Nl))
     !print*, 'pI0=', pI0

     xin(1)=maxval(MF%F(1:Nl)%p)
     bin(:)=0.0_PR

     if(verb)then
        write(iout,*) 'pI0=', pI0
        write(iout,*) '-----------------------------------'
        write(iout,*) '   x=', xin(1:Nsnes)
        write(iout,*) '   rhl=', MF%F(1:Nl)%rh
        write(iout,*) '   sumY=',sum(MF%F(1:Nl)%Y) 
        write(iout,*) '   sumfl=',sum(MF%F(1:Nl)%f) 
        write(iout,*) '   elh=', MF%F(1:Nl)%eh
        write(iout,*) '   p=', MF%p
        write(iout,*) '   pl=', MF%F(1:Nl)%p
        write(iout,*) '   fl=', MF%F(1:Nl)%f
        write(iout,*) '   rho=', MF%rh
     endif     

     ! call solve_SNES(Nsnes,xin,bin)
     call My_Newton(1,xin(1:1),1.0e-10_PR,maxval(MF%F(:)%p),F=F_relaxp_sge_2,DF=J_relaxp_sge_2)

     MF%p=xin(1)

     !if(p.lt.0.0_PR)then
     !    verb=.true.
     !    print*, 'coucou1'
     !endif

     if(verb)then
        write(iout,*) '   e+p/rh=', MF%F(1:Nl)%eh+MF%p/MF%F(1:Nl)%rh
        write(iout,*) '-----------------------------------'     
        write(iout,*) '   p=', MF%p
     endif

     do iph=1,Nl
        MF%F(iph)%f=fl0(iph)*(pl0(iph)+MF%F(iph)%g_sge*MF%F(iph)%p_sge+(MF%F(iph)%g_sge-1.0_PR)*pI0)/&
        (MF%p+MF%F(iph)%g_sge*MF%F(iph)%p_sge+(MF%F(iph)%g_sge-1.0_PR)*pI0)
        MF%F(iph)%rh=flrhl0(iph)/MF%F(iph)%f
        MF%F(iph)%Y=MF%F(iph)%f*MF%F(iph)%rh/MF%rh
     enddo

     do iph=1,Nl
      MF%F(iph)%eh=ph(iph)%energie_hydro(MF%F(iph),MF%F(iph)%rh,MF%p)
      !F%elh(iph)=energie_hydro_sge(F%g_sge(iph),F%p_sge(iph),F%rhl(iph),F%p)
      MF%F(iph)%p=ph(iph)%pressure(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%eh)
     enddo

     !!! rh=sum(fl(1:Nl)*rhl(1:Nl))
 
     if(verb)then     
        write(iout,*) '-----------------------------------'
        write(iout,*) '   x=', xin(1:Nsnes)
        write(iout,*) '   rhl=', MF%F(1:Nl)%rh
        write(iout,*) '   sumY=',sum(MF%F(1:Nl)%Y) 
        write(iout,*) '   sumfl=',sum(MF%F(1:Nl)%f) 
        write(iout,*) '   elh=', MF%F(1:Nl)%eh
        write(iout,*) '   pl=', MF%F(1:Nl)%p
        write(iout,*) '   fl=', MF%F(1:Nl)%f
        write(iout,*) '   rho=', MF%rh
        write(iout,*) '   e+p/rh=', MF%F(1:Nl)%eh+MF%p/MF%F(1:Nl)%rh
        write(iout,*) '-----------------------------------'
        write(iout,*) '---Set pressure from e---'
     endif    

     !ENDIF !!! pmin-pmax/(pmin+pmax)>1e-3

     ehtot=MF%u   !!!F%e-0.5_PR*(F%vx**2+F%vy**2+F%vz**2)

     MF%p=pressure_mixt(MF)
 
     !if(p.lt.0.0_PR)then
     !    verb=.true.
     !    print*, 'coucou2'
     !endif

     do iph=1,Nl
        MF%F(iph)%p=MF%p
     enddo

     do iph=1,Nl
       MF%F(iph)%eh=ph(iph)%energie_hydro(MF%F(iph),MF%F(iph)%rh,MF%p)
       !F%elh(iph)=energie_hydro_sge(F%g_sge(iph),F%p_sge(iph),F%rhl(iph),F%p)
     enddo

     if(abs(sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%eh)-ehtot)/ehtot.gt.1.0e-6_PR)then
        write(iout,*) 'PROBLEM ehtot ne sum(Ylehl):'
        write(iout,*) sum(MF%F(1:Nl)%Y*MF%F(1:Nl)%eh), ehtot
        stop
     endif
 
     !F%c=soundspeed_mixt(F)
 
     if(verb)then 
        write(iout,*) 'P=', MF%p
        write(iout,*) 'elh1=', MF%F(1)%eh
        write(iout,*) 'elh2=', MF%F(2)%eh
     endif     

end subroutine relax_p_sge_2

subroutine relax_p_eos(MF)


     implicit none
     type(multifluide), intent(inout) :: MF
     real(PR) :: pmin, pmax
     !real(PR), allocatable :: xin(:), bin(:)
     integer :: iph
 
     pmin=minval(MF%F(1:Nl)%p)
     pmax=maxval(MF%F(1:Nl)%p)
     
     !IF(0.5_PR*abs(pmax-pmin)/(pmax+pmin).gt.1.0e-3_PR)THEN
     Nl0=Nl
     do iph=1,Nl
       flrhl0(iph)=MF%F(iph)%f*MF%F(iph)%rh
       fl0_gl(iph)=MF%F(iph)%f/MF%F(iph)%g_sge
       pl0(iph)=MF%F(iph)%p  
       nul0(iph)=1.0_PR/MF%F(iph)%rh
       coef1(iph)=fl0_gl(iph)*(MF%F(iph)%p_sge+pl0(iph))
       fl0(iph)=MF%F(iph)%f
       el0(iph)=MF%F(iph)%eh
       p_sge0(iph)=MF%F(iph)%p_sge
       g_sge0(iph)=MF%F(iph)%g_sge
     enddo 

     xin(1:Nl)=nul0(1:Nl)
     xin(Nl+1)=maxval(MF%F(1:Nl)%p)
     bin(:)=0.0_PR

     do iph=1,Nl
      Zl(iph)=MF%F(iph)%rh*ph(iph)%soundspeed(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%p)
     enddo
     pI0=sum(Zl(1:Nl)*MF%F(1:Nl)%p)/sum(Zl(1:Nl))

     if(verb)then
        write(iout,*) 'pI0=', pI0
        write(iout,*) '-----------------------------------'
        write(iout,*) '   x=', xin(1:Nsnes)
        write(iout,*) '   rhl=', MF%F(1:Nl)%rh
        write(iout,*) '   sumY=',sum(MF%F(1:Nl)%Y) 
        write(iout,*) '   sumfl=',sum(MF%F(1:Nl)%f) 
        write(iout,*) '   elh=', MF%F(1:Nl)%eh
        write(iout,*) '   p=', MF%p
        write(iout,*) '   pl=', MF%F(1:Nl)%p
        write(iout,*) '   fl=', MF%F(1:Nl)%f
        write(iout,*) '   rho=', MF%rh
     endif     

     !call solve_SNES(Nsnes,xin,bin)
     call My_Newton(1,xin(1:1),minval(MF%F(:)%p),maxval(MF%F(:)%p),F=F_relaxp_sge_1,DF=J_relaxp_sge_1)

     MF%p=xin(Nl+1)

     if(MF%p.lt.0.0_PR)then
         verb=.true.
         print*, 'coucou1'
     endif

     if(verb)then
        if(mode_pI.eq.2)then
           write(iout,*) '   e+pI0/rh=', MF%F(1:Nl)%eh+pI0/MF%F(1:Nl)%rh
        else
           write(iout,*) '   e+p/rh=', MF%F(1:Nl)%eh+MF%p/MF%F(1:Nl)%rh
        endif
        write(iout,*) '-----------------------------------'     
        write(iout,*) '   x=', xin(1:Nl+1)
     endif

     MF%F(1:Nl)%rh=1.0_PR/xin(1:Nl)
     MF%F(1:Nl)%f=flrhl0(1:Nl)*xin(1:Nl)
     MF%F(1:Nl)%Y=MF%F(1:Nl)%f*MF%F(1:Nl)%rh/MF%rh

     do iph=1,Nl
      MF%F(iph)%eh=ph(iph)%energie_hydro(MF%F(iph),MF%F(iph)%rh,MF%p)
      !F%elh(iph)=energie_hydro_sge(F%g_sge(iph),F%p_sge(iph),F%rhl(iph),F%p)
      MF%F(iph)%p=ph(iph)%pressure(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%eh)
     enddo

     MF%rh=sum(MF%F(1:Nl)%f*MF%F(1:Nl)%rh)
 
     if(verb)then     
        write(iout,*) '-----------------------------------'
        write(iout,*) '   x=', xin(1:Nl+1)
        write(iout,*) '   rhl=', MF%F(1:Nl)%rh
        write(iout,*) '   sumY=',sum(MF%F(1:Nl)%Y) 
        write(iout,*) '   sumfl=',sum(MF%F(1:Nl)%f) 
        write(iout,*) '   elh=', MF%F(1:Nl)%eh
        write(iout,*) '   p=', MF%p
        write(iout,*) '   pl=', MF%F(1:Nl)%p
        write(iout,*) '   fl=', MF%F(1:Nl)%f
        write(iout,*) '   rho=', MF%rh
        if(mode_pI.eq.2)then
           write(iout,*) '   e+pI0/rh=', MF%F(1:Nl)%eh+pI0/MF%F(1:Nl)%rh
        else
           write(iout,*) '   e+p/rh=', MF%F(1:Nl)%eh+MF%p/MF%F(1:Nl)%rh
        endif
        write(iout,*) '-----------------------------------'
        write(iout,*) '---Set pressure from e---'
     endif    

     !ENDIF !!! pmin-pmax/(pmin+pmax)>1e-3

     MF%p=pressure_mixt(MF) 
     if(MF%p.lt.0.0_PR)then
         verb=.true.
         print*, 'coucou2'
     endif

     do iph=1,Nl
        MF%F(iph)%p=MF%p
     enddo

     do iph=1,Nl
       MF%F(iph)%eh=ph(iph)%energie_hydro(MF%F(iph),MF%F(iph)%rh,MF%p)
       !F%elh(iph)=energie_hydro_sge(F%g_sge(iph),F%p_sge(iph),F%rhl(iph),F%p)
     enddo
 
     !F%c=soundspeed_mixt(F)
 
     if(verb)then 
        write(iout,*) 'P=', MF%p
        write(iout,*) 'elh1=', MF%F(1)%eh
        write(iout,*) 'elh2=', MF%F(2)%eh
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

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     !type(multifluide), intent(in) :: flu 
     real(PR) :: p
     integer :: iph 

     !f(1)=sum( fl0_gl(1:flu%Nl)*(pl0(1:flu%Nl)-x(1))/(flu%p_sge(1:flu%Nl)+x(1)) )

     !f(1)=sum( fl0_gl(1:Nl0)*(pl0(1:Nl0)-x(1))/(p_sge0(1:Nl0)+x(1)) )
    
      f(1)=-1.0_PR
      do iph=1,Nl0
         f(1)=f(1)+fl0_gl(iph)*(pl0(iph)+g_sge0(iph)*p_sge0(iph)+(g_sge0(iph)-1.0_PR)*x(1))/(p_sge0(iph)+x(1))
      enddo

end subroutine F_relaxp_sge_1

subroutine J_relaxp_sge_1(N,x,J)

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)
     !type(multifluide), intent(in) :: flu
     real(PR) :: p
     integer :: k, iph

     !!!---dérivée par rapport à p
     !p=x(1)
     !J(1,1)=-sum(fl0_gl(1:flu%Nl)*(flu%p_sge(1:flu%Nl)+pl0(1:flu%Nl))/(flu%p_sge(1:flu%Nl)+x(1))**2)
     !J(1,1)=-sum(fl0_gl(1:Nl0)*(p_sge0(1:Nl0)+pl0(1:Nl0))/(p_sge0(1:Nl0)+x(1))**2)

      J(1,1)=0.0_PR
      do iph=1,Nl0
         J(1,1)=J(1,1)-fl0_gl(iph)*(p_sge0(iph)+pl0(iph))/(p_sge0(iph)+x(1))**2 
      enddo

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

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     !type(fluide), intent(in) :: flu
     real(PR) :: p, c1
     integer :: k

     !f(1)=sum( fl0(1:flu%Nl)*(pl0(1:flu%Nl)-x(1))/&
     !        (flu%g_sge(1:flu%Nl)*flu%p_sge(1:flu%Nl) + (flu%g_sge(1:flu%Nl)-1.0_PR)*PI0 + x(1)) )

     f(1)=sum( fl0(1:Nl0)*(pl0(1:Nl0)-x(1))/&
             (g_sge0(1:Nl0)*p_sge0(1:Nl0) + (g_sge0(1:Nl0)-1.0_PR)*PI0 + x(1)) )


end subroutine F_relaxp_sge_2

subroutine J_relaxp_sge_2(N,x,J)

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)
     !type(fluide), intent(in) :: flu
     real(PR) :: p
     integer :: k

     p=x(1)

     !!!---dérivée par rapport à p
     !J(1,1)=-sum(fl0(1:flu%Nl)*(flu%g_sge(1:flu%Nl)*flu%p_sge(1:flu%Nl)+(flu%g_sge(1:flu%Nl)-1.0_PR)*pI0+pl0(1:flu%Nl))/&
      !          (flu%g_sge(1:flu%Nl)*flu%p_sge(1:flu%Nl)+(flu%g_sge(1:flu%Nl)-1.0_PR)*pI0+p)**2 )
 
     J(1,1)=-sum(fl0(1:Nl0)*(g_sge0(1:Nl0)*p_sge0(1:Nl0)+(g_sge0(1:Nl0)-1.0_PR)*pI0+pl0(1:Nl0))/&
                (g_sge0(1:Nl0)*p_sge0(1:Nl0)+(g_sge0(1:Nl0)-1.0_PR)*pI0+p)**2 )
  
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

     implicit none
     integer, intent(in) :: N !!! =Nl0+1
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     !type(fluide), intent(in) :: flu
     real(PR) :: p, rhk
     integer :: k

     p=x(N)

     do k=1,N-1
        rhk=1.0_PR/x(k)
        !f(k)=ph(k)%energie_hydro(k,flu,rhk,p)-el0(k)+p*(x(k)-nul0(k))
        f(k)=energie_hydro_sge0(rhk,p,g_sge0(k),p_sge0(k))-el0(k)+p*(x(k)-nul0(k))
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

     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: f(1:N)
     !type(fluide), intent(in) :: flu
     real(PR) :: p, rhk
     integer :: k

     p=x(N)

     do k=1,N-1
        rhk=1.0_PR/x(k)
        !f(k)=ph(iph)%energie_hydro(k,flu,rhk,p)-el0(k)+pI0*(x(k)-nul0(k))
        f(k)=energie_hydro_sge0(rhk,p,g_sge0(k),p_sge0(k))-el0(k)+pI0*(x(k)-nul0(k))
     enddo   

     f(N)=-1.0_PR
 
     do k=1,N-1
       f(N)=f(N)+flrhl0(k)*x(k)
     enddo

end subroutine F_relaxp_2

subroutine J_relaxp_1(N,x,J)


     implicit none
     integer, intent(in) :: N
     real(PR), intent(in) :: x(1:N)
     real(PR), intent(out) :: J(1:N,1:N)
     !type(fluide), intent(in) :: flu
     real(PR) :: p, rhk
     integer :: k

     J(:,:)=0.0_PR

     p=x(N)

     !!!---dérivée par rapport à nu
     do k=1,N-1
        rhk=1.0_PR/x(k)
        J(k,k)=ddnu_energie_hydro_sge0(rhk,p,g_sge0(k),p_sge0(k))+p
        !!!J(k,k)=ddnu_energie_hydro_sge(k,F0,rhk,p)+p
     enddo
     !!!---dérivée par rapport à p
     do k=1,N-1
        rhk=1.0_PR/x(k)
        J(k,N)=ddp_energie_hydro_sge0(rhk,p,g_sge0(k),p_sge0(k))+(x(k)-nul0(k))
        !!!J(k,N)=ddp_energie_hydro_sge(k,F0,rhk,p)+(x(k)-nul0(k))
     enddo

     do k=1,N-1
        J(N,k)=flrhl0(k)
        !!!J(N,k)=F0%fl(k)*F0%rhl(k)
     enddo
 
end subroutine J_relaxp_1

!!!  subroutine J_relaxp_2(N,x,J,flu)
!!!  
!!!       implicit none
!!!       integer, intent(in) :: N
!!!       real(PR), intent(in) :: x(1:N)
!!!       real(PR), intent(out) :: J(1:N,1:N)
!!!       type(multifluide), intent(in) :: flu
!!!       real(PR) :: p, rhk
!!!       integer :: k
!!!  
!!!       J(:,:)=0.0_PR
!!!  
!!!       p=x(N)
!!!  
!!!       !!!---dérivée par rapport à nu
!!!       do k=1,N-1
!!!          rhk=1.0_PR/x(k)
!!!          J(k,k)=ddnu_energie_hydro_sge0(rhk,p,g_sge0(k),p_sge0(k))+pI0
!!!       enddo
!!!       !!!---dérivée par rapport à p 
!!!       do k=1,N-1
!!!          rhk=1.0_PR/x(k)
!!!          J(k,N)=ddp_energie_hydro_sge0(rhk,p,g_sge0(k),p_sge0(k))
!!!       enddo
!!!  
!!!       do k=1,N-1
!!!          J(N,k)=flrhl0(k)
!!!       enddo
!!!   
!!!  end subroutine J_relaxp_2

!!!========================== SUBROUTINE DE TEST ===============================
subroutine test_protection(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph,ir,iT,iE
   real(PR) :: rho,T,P,E,Cv,Temp,Temp2,P2,Ener,resid

   write(iout,*) '   test protection'

   radial=.false.

   M%Nt=500
   M%dt=1.0e-9_PR

   !!!-------------Nx,Ny,Nz,Nl,Lx,Ly,Lz,M
   call init_euler(2000,1,1,3,2.0e-3_PR,2.e-4_PR,1.0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='air'
   ph(2)%typ=1 ; ph(2)%nom='Cug'
   ph(3)%typ=1 ; ph(3)%nom='epoxy'

   !ph(2)%typ=2 ; ph(2)%filetab='tab_Cu_401.dat'
   call init_phase(M)
 
   Temp=300.0_PR
   P=1.0e5_PR
   Temp2=3000.0_PR
   Ener=ph(2)%ener0+ph(2)%cv*(Temp2-300.0_PR)
   P2=rho0_Cu*Ener*(ph(2)%g_sge-1.0_PR)-ph(2)%g_sge*ph(2)%p_sge

   print*, 'cas test protection: Pcu=', P2, 'Ener=', Ener 

   !!!---fractions volumiques

   resid=1.0e-6_PR

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      if(M%x(i).lt.0.050e-3_PR)then !!! cug
         M%MF(i,j,k)%F(1)%f=0.5_PR*resid
         M%MF(i,j,k)%F(1)%f=1.0_PR-resid
         M%MF(i,j,k)%F(1)%f=0.5_PR*resid
         call set_TP0(M%MF(i,j,k),Temp2,P2)
         !call set_TP0(M%MF(i,j,k),Temp,P)

      elseif(M%x(i).lt.0.550e-3_PR)then !!! epoxy
         M%MF(i,j,k)%F(1)%f=0.5_PR*resid
         M%MF(i,j,k)%F(2)%f=0.5_PR*resid
         M%MF(i,j,k)%F(3)%f=1.0_PR-resid
         call set_TP0(M%MF(i,j,k),Temp,P)
      else !!! Air
         M%MF(i,j,k)%F(1)%f=1.0_PR-resid
         M%MF(i,j,k)%F(2)%f=0.5_PR*resid
         M%MF(i,j,k)%F(3)%f=0.5_PR*resid
         call set_TP0(M%MF(i,j,k),Temp,P)
      endif

      M%MF(i,j,k)%vx=0.0_PR

   enddo; enddo ; enddo

end subroutine test_protection



subroutine test_advection_Cu(M)


   !!!***Favire 2009 p6059-6064
   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph,ir,iT,iE
   real(PR) :: rho,T,P,E,Cv,Temp,Temp2,P2,Ener,resid

   write(iout,*) '   test advection Cu'

   radial=.false.

   M%Nt=380
   M%dt=1.0e-6_PR

   !!!-------------Nx,Ny,Nz,Nl,Lx,Ly,Lz,M
   call init_euler(100,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='Cu'
   ph(2)%typ=1 ; ph(2)%nom='air'

   call init_phase(M)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%rh  = ph(iph)%rho0
      M%MF(i,j,k)%F(iph)%p  = 1.0e5_PR
      !M%MF(i,j,k)%mu(iph)  = 0.0_PR
   enddo; enddo; enddo ; enddo
 
   !!!---fractions volumiques

   resid=1.0e-6_PR

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
         M%MF(i,j,k)%vx=1000.0_PR
         if(M%x(i).le.0.4_PR)then 
           M%MF(i,j,k)%F(1)%f=1.0_PR-resid
           M%MF(i,j,k)%F(2)%f=resid
         else
           M%MF(i,j,k)%F(1)%f=resid
           M%MF(i,j,k)%F(2)%f=1.0_PR-resid
         endif
   enddo; enddo ; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:Nl)%f*M%MF(i,j,k)%F(1:Nl)%rh)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh
      M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
    call MAJ_meca(M%MF(i,j,k))
   enddo ; enddo ; enddo


   write(iout,*) ' MAJ mixt'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_mixt(M%MF(i,j,k))
   enddo ; enddo ; enddo

end subroutine test_advection_Cu

subroutine test_solid_solid_shock_tube(M)


   !!!***Favrie 2009 p6060-6064
   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph,ir,iT,iE
   real(PR) :: rho,T,P,E,Cv,Temp,Temp2,P2,Ener,resid

   write(iout,*) '   test solid-solid shock tube'

   radial=.false.

   M%Nt=80  !!!134
   M%dt=1.0e-6_PR


   !!!-------------Nx,Ny,Nz,Nl,Lx,Ly,Lz,M
   call init_euler(1000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='Al'
   ph(2)%typ=1 ; ph(2)%nom='Cu'

   call init_phase(M)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%rh  = ph(iph)%rho0
      !M%MF(i,j,k)%mu(iph)  = 0.0_PR
   enddo; enddo; enddo ; enddo
 
   !!!---fractions volumiques

   resid=1.0e-6_PR

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
         if(M%x(i).le.0.5_PR)then 
           M%MF(i,j,k)%F(1)%f=1.0_PR-resid
           M%MF(i,j,k)%F(2)%f=resid
           M%MF(i,j,k)%F(:)%p  = 1.0e8_PR
         else
           M%MF(i,j,k)%F(1)%f=resid
           M%MF(i,j,k)%F(2)%f=1.0_PR-resid
           M%MF(i,j,k)%F(:)%p  = 1.0e5_PR
         endif
   enddo; enddo ; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh
      M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
    call MAJ_meca(M%MF(i,j,k))
   enddo ; enddo ; enddo


   write(iout,*) ' MAJ mixt'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_mixt(M%MF(i,j,k))
   enddo ; enddo ; enddo

end subroutine test_solid_solid_shock_tube


subroutine test_solid_gas_shock_tube(M)


   !!!***Favrie 2009 p6061-6065
   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph,ir,iT,iE
   real(PR) :: rho,T,P,E,Cv,Temp,Temp2,P2,Ener,resid,coeff

   write(iout,*) '   test solid-gas shock tube'

   radial=.false.

   M%Nt=87
   M%dt=1.0e-6_PR


   !!!-------------Nx,Ny,Nz,Nl,Lx,Ly,Lz,M
   call init_euler(1000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='Cu'
   ph(2)%typ=1 ; ph(2)%nom='air'

   call init_phase(M)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
      M%MF(i,j,k)%F(1)%rh  = ph(1)%rho0
      M%MF(i,j,k)%F(2)%rh  = 50.0_PR*ph(2)%rho0
      !M%MF(i,j,k)%mu(1)  = 0.0_PR

      M%MF(i,j,k)%F(1)%rh0=M%MF(i,j,k)%F(1)%rh
      M%MF(i,j,k)%F(2)%rh0=M%MF(i,j,k)%F(2)%rh

   enddo; enddo; enddo
 
   !!!---fractions volumiques

   resid=1.0e-6_PR

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
         if(M%x(i).le.0.6_PR)then 
           M%MF(i,j,k)%F(1)%f=1.0_PR-resid
           M%MF(i,j,k)%F(2)%f=resid
           M%MF(i,j,k)%F(:)%p  = 5.0e9_PR
         else
           M%MF(i,j,k)%F(1)%f=resid
           M%MF(i,j,k)%F(2)%f=1.0_PR-resid
           M%MF(i,j,k)%F(:)%p  = 1.0e5_PR
         endif
         do iph=1,Nl
            M%MF(i,j,k)%F(iph)%f0=M%MF(i,j,k)%F(iph)%f
         enddo

         !do iph=1,Nl
         !   coeff=M%MF(i,j,k)%fl(iph)**(r13)
         !   M%MF(i,j,k)%al(iph,1:3)=M%MF(i,j,k)%al(iph,1:3)*coeff
         !   M%MF(i,j,k)%bl(iph,1:3)=M%MF(i,j,k)%bl(iph,1:3)*coeff
         !   M%MF(i,j,k)%cl(iph,1:3)=M%MF(i,j,k)%cl(iph,1:3)*coeff
         !enddo

   enddo; enddo ; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh
      M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
    call MAJ_meca(M%MF(i,j,k))
   enddo ; enddo ; enddo


   write(iout,*) ' MAJ mixt'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_mixt(M%MF(i,j,k))
   enddo ; enddo ; enddo

end subroutine test_solid_gas_shock_tube

subroutine test_cavitation_Cu(M)

   !!!***Favire 2009 p6062-6065
   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph,ir,iT,iE
   real(PR) :: rho,T,P,E,Cv,Temp,Temp2,P2,Ener,resid,coeff

   write(iout,*) '   test cavitation Cu'

   radial=.false.

   M%Nt=88
   M%dt=1.0e-6_PR

   !!!-------------Nx,Ny,Nz,Nl,Lx,Ly,Lz,M
   call init_euler(2000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='Cu'
   ph(2)%typ=1 ; ph(2)%nom='air'

   !ph(2)%typ=2 ; ph(2)%filetab='tab_Cu_401.dat'
   call init_phase(M)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%rh  = ph(iph)%rho0
      M%MF(i,j,k)%F(iph)%p  = 1.0e5_PR
      !M%MF(i,j,k)%mu(iph)=0.0_PR
   enddo; enddo; enddo ; enddo
 
   !!!---fractions volumiques

   resid=1.0e-4_PR

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
         M%MF(i,j,k)%F(1)%f=1.0_PR-resid
         M%MF(i,j,k)%F(2)%f=resid
         if(M%x(i).le.0.5_PR)then 
           M%MF(i,j,k)%vx=-3000.0_PR
         else
           M%MF(i,j,k)%vx=3000.0_PR
         endif
         !do iph=1,Nl
         !   coeff=M%MF(i,j,k)%fl(iph)**(r13)
         !   M%MF(i,j,k)%al(iph,1:3)=M%MF(i,j,k)%al(iph,1:3)*coeff
         !   M%MF(i,j,k)%bl(iph,1:3)=M%MF(i,j,k)%bl(iph,1:3)*coeff
         !   M%MF(i,j,k)%cl(iph,1:3)=M%MF(i,j,k)%cl(iph,1:3)*coeff
         !enddo
   enddo; enddo ; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh
      M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
    call MAJ_meca(M%MF(i,j,k))
   enddo ; enddo ; enddo

      do i=1,M%Nx
        if(abs(M%MF(i,1,1)%ee).gt.1.0e15_PR)then
          print*, 'PROBLEM u 1 i=', i
          stop
        endif
      enddo 

   write(iout,*) ' MAJ mixt'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_mixt(M%MF(i,j,k))
   enddo ; enddo ; enddo

end subroutine test_cavitation_Cu

subroutine test_Ndanou2015(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph,ir,iT,iE
   real(PR) :: rho,T,P,E,Cv,Temp,Temp2,P2,Ener,resid

   write(iout,*) '   test protection'

   radial=.false.

   M%Nt=200
   M%dt=1.0e-9_PR

   !!!-------------Nx,Ny,Nz,Nl,Lx,Ly,Lz,M
   call init_euler(4000,1,1,3,1.7e-2_PR,1.0_PR,1.0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='air'
   ph(2)%typ=1 ; ph(2)%nom='Al'
   ph(3)%typ=1 ; ph(3)%nom='Ti'

   !ph(2)%typ=2 ; ph(2)%filetab='tab_Cu_401.dat'
   call init_phase(M)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%rh  = ph(iph)%rho0
      M%MF(i,j,k)%F(iph)%p  = 1.0e5_PR
      !M%MF(i,j,k)%rho0(1)=M%MF(i,j,k)%rhl(1)
      !M%MF(i,j,k)%rho0(2)=M%MF(i,j,k)%rhl(2)
      !M%MF(i,j,k)%mu(2)=0.0e-2_PR
      !M%MF(i,j,k)%mu(3)=0.0e-2_PR
   enddo; enddo; enddo ; enddo
 
   !!!---fractions volumiques

   resid=1.0e-6_PR

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      M%MF(i,j,k)%vx=0.0_PR
      if(M%x(i).lt.1.0e-3_PR)then !!! air
         M%MF(i,j,k)%F(1)%f=1.0_PR-resid
         M%MF(i,j,k)%F(2)%f=0.5_PR*resid
         M%MF(i,j,k)%F(3)%f=0.5_PR*resid
      elseif(M%x(i).lt.3.0e-3_PR)then !!! Al
         M%MF(i,j,k)%F(1)%f=0.5_PR*resid
         M%MF(i,j,k)%F(2)%f=1.0_PR-resid
         M%MF(i,j,k)%F(3)%f=0.5_PR*resid
         M%MF(i,j,k)%vx=700.0_PR
      elseif(M%x(i).lt.12.8e-3_PR)then !!! Ti
         M%MF(i,j,k)%F(1)%f=0.5_PR*resid
         M%MF(i,j,k)%F(2)%f=0.5_PR*resid
         M%MF(i,j,k)%F(3)%f=1.0_PR-resid
      else !!!Air
         M%MF(i,j,k)%F(1)%f=1.0_PR-resid
         M%MF(i,j,k)%F(2)%f=0.5_PR*resid
         M%MF(i,j,k)%F(3)%f=0.5_PR*resid
      endif

   enddo; enddo ; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh
      M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
    call MAJ_meca(M%MF(i,j,k))
   enddo ; enddo ; enddo
   write(iout,*) ' MAJ mixt'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_mixt(M%MF(i,j,k))
   enddo ; enddo ; enddo

end subroutine test_Ndanou2015

subroutine test_meca(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph,ir,iT,iE
   real(PR) :: rho,T,P,E,Cv,Temp,Temp2,P2,Ener,resid

   write(iout,*) '   test meca'

   radial=.false.

   M%Nt=1
   M%dt=1.0e-9_PR

   !!!-------------Nx,Ny,Nz,Nl,Lx,Ly,Lz,M
   call init_euler(1000,1,1,2,2.0e-3_PR,2.e-4_PR,1.0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='Ti'
   ph(2)%typ=1 ; ph(2)%nom='Air'

   call init_phase(M)
 
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%rh  = ph(iph)%rho0
      M%MF(i,j,k)%F(iph)%p  = 1.0e5_PR
      M%MF(i,j,k)%vx       = 0.0_PR
   enddo; enddo; enddo ; enddo
 
   !!!---fractions volumiques

   resid=1.0e-6_PR

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      if(M%x(i).lt.1.0e-3_PR)then !!! Ti
         M%MF(i,j,k)%F(1)%f=1.0_PR-resid
         M%MF(i,j,k)%F(2)%f=0.5_PR*resid
         M%MF(i,j,k)%F(:)%a(1)=0.9_PR
      else !!!Air
         M%MF(i,j,k)%F(2)%f=1.0_PR-resid
         M%MF(i,j,k)%F(1)%f=0.5_PR*resid
      endif

   enddo; enddo ; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
      M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh
      M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
    call MAJ_meca(M%MF(i,j,k))
   enddo ; enddo ; enddo
   write(iout,*) ' MAJ mixt'
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_mixt(M%MF(i,j,k))
   enddo ; enddo ; enddo

end subroutine test_meca

subroutine test_fil(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph,ir,iT,iE
   real(PR) :: rho,T,P,E,Cv,Temp   

   write(iout,*) '   test fil cuivre'

   radial=.false.

   M%Nt=5000
   M%dt=0.1e-9_PR

   !!!-------------Nx,Ny,Nz,Nl,Lx,Ly,Lz,M
   call init_euler(1000,1,1,2,1.0e-3_PR,2.e-4_PR,1.0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='air'
   ph(2)%typ=1 ; ph(2)%nom='Cug'
   !ph(2)%typ=2 ; ph(2)%filetab='tab_Cu_401.dat'
   call init_phase(M)
   
   !!!---fractions volumiques

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

     !if(M%x(i).lt.5.0e-5_PR)then 
     !if(M%x(i).lt.0.550e-3_PR.and.M%x(i).gt.0.495e-3_PR)then 
     !    M%MF(i,j,k)%fl(1)=1.0e-8_PR
     !    M%MF(i,j,k)%fl(2)=1.0_PR-1.0e-8_PR
     ! else
     !    M%MF(i,j,k)%fl(1)=1.0_PR-1.0e-8_PR
     !    M%MF(i,j,k)%fl(2)=1.0e-8_PR
     ! endif
     if(M%x(i).lt.0.050e-3_PR)then 
         M%MF(i,j,k)%F(1)%f=1.0e-8_PR
         M%MF(i,j,k)%F(2)%f=1.0_PR-1.0e-8_PR
      elseif(M%x(i).lt.0.550e-3_PR)then
         M%MF(i,j,k)%F(1)%f=1.0_PR-1.0e-8_PR
         M%MF(i,j,k)%F(2)%f=1.0e-8_PR
      endif


      M%MF(i,j,k)%vx=0.0_PR

   enddo; enddo ; enddo

   !!!---set TP
   P=1.0e5_PR
   Temp=300.0_PR

   call set_TP(M,Temp,P) 

   !ph(2)%typ=2 ; ph(2)%filetab='tab_Cu_401.dat'


end subroutine test_fil

subroutine test_fil2(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph
   real(PR) :: trace  
 
   write(iout,*) ' Cas test fil 2'

   M%Nt=5000 !!! 00
   M%dt=0.1e-9_PR

   !call init_euler(1000,1,1,1)

   call init_euler(1000,1,1,2,1.0e-3_PR,2.e-4_PR,1.e0_PR,M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 

      M%MF(i,j,k)%F(1)%p_sge=p_sge_Cu ; M%MF(i,j,k)%F(1)%g_sge=g_sge_Cu
      M%MF(i,j,k)%F(1)%mu=mu_Cu       ; M%MF(i,j,k)%F(1)%sigy=sigy_Cu
      M%MF(i,j,k)%F(1)%rh0=rho0_Cu    ; M%MF(i,j,k)%F(1)%cv=cv_Cu

      M%MF(i,j,k)%F(2)%p_sge=p_sge_Cug ; M%MF(i,j,k)%F(2)%g_sge=g_sge_Cug
      M%MF(i,j,k)%F(2)%mu=mu_Cug       ; M%MF(i,j,k)%F(2)%sigy=sigy_Cug
      M%MF(i,j,k)%F(2)%rh0=rho0_Cug    ; M%MF(i,j,k)%F(2)%cv=cv_Cug  
 
   enddo; enddo; enddo

   trace=1.0e-20_PR
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      if(abs(M%x(i)-0.0e-3_PR).lt.100.0e-6_PR)then

         M%MF(i,j,k)%F(1)%f=1.0_PR-trace
         !M%MF(i,j,k)%Yl(1)=1.0_PR-trace
         M%MF(i,j,k)%F(1)%p=1.0e5_PR
         M%MF(i,j,k)%F(1)%rh=1.0_PR*M%MF(i,j,k)%F(1)%rh0

         M%MF(i,j,k)%F(2)%f=trace
         !M%MF(i,j,k)%Yl(2)=trace
         M%MF(i,j,k)%F(2)%p=1.0e5_PR
         M%MF(i,j,k)%F(2)%rh=M%MF(i,j,k)%F(2)%rh0

      else

         M%MF(i,j,k)%F(2)%f=1.0_PR-trace
         !M%MF(i,j,k)%Yl(2)=1.0_PR-trace
         M%MF(i,j,k)%F(2)%p=1.0e5_PR
         M%MF(i,j,k)%F(2)%rh=M%MF(i,j,k)%F(2)%rh0

         M%MF(i,j,k)%F(1)%f=trace
         !M%MF(i,j,k)%Yl(1)=trace
         M%MF(i,j,k)%F(1)%p=1.0e5_PR
         M%MF(i,j,k)%F(1)%rh=M%MF(i,j,k)%F(1)%rh0

      endif

      M%MF(i,j,k)%vx=0.0_PR

   enddo; enddo ; enddo

   !---Yl imposé :
   !forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=1.0_PR/sum(M%MF(i,j,k)%Yl(1:M%Nl)/M%MF(i,j,k)%rhl(1:M%Nl))
   !forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%fl(iph)=M%MF(i,j,k)%Yl(iph)*M%MF(i,j,k)%rh/M%MF(i,j,k)%rhl(iph)


   !---fl imposé:
   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)
   forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
     !M%MF(i,j,k)%elh(iph)=energie_hydro_sge(M%MF(i,j,k)%g_sge(iph),M%MF(i,j,k)%p_sge(iph),M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph))
     M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   !!!---set energie de reference
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,Nl
!     M%MF(i,j,k)%elh0(iph)=energie_hydro_sge(M%MF(i,j,k)%g_sge(iph),M%MF(i,j,k)%p_sge(iph),M%MF(i,j,k)%rho0(iph),1.e5_PR)
     M%MF(i,j,k)%F(iph)%eh0=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh0,1.e5_PR)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
    call MAJ_meca(M%MF(i,j,k))
   enddo ; enddo ; enddo

   write(iout,*) ' MAJ mixt'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_mixt(M%MF(i,j,k))
   enddo ; enddo ; enddo

   !!!---THERMO
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
      M%MF(i,j,k)%F(:)%T=300.0_PR
      M%MF(i,j,k)%T=300.0_PR
      !do iph=1, Nl
      !M%MF(i,j,k)%cvl(iph)=specific_heat(min(M%MF(i,j,k)%Tl(iph),70000.0_PR),min(M%MF(i,j,k)%pl(iph),1.0e8_PR))&
      !                    /M%MF(i,j,k)%g_sge(iph)   
      !enddo
   enddo ; enddo ; enddo


end subroutine test_fil2



!subroutine test_Ndanou9(F)
!
!   use mod_data, only : Nt,t,dt,x
!
!   !use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
!   !                     al,bl,cl,sigl,sig,&
!   !                     rh,vx,vy,vz,c,p,e,&
!   !                     x,y,z
!
!   implicit none
!   type(fluide), allocatable, intent(inout) :: F(:,:,:)
!   integer :: i,j,k,iph
!  
!   write(iout,*) '  cas de test Ndanou 9'
!
!   !!!------------------     cas test Fig 9 Ndanou 2015   ---------------------------------------
!    call init_euler(2000,1,1,3,F)
!   !! phase 1 : Air 2 : Aluminium 3 : Titane
!
!
!   do k=1,Nz ; do j=1,ny ; do i=1,nx 
!   F(i,j,k)%p_sge(1)=p_sge_air ; F(i,j,k)%g_sge(1)=g_sge_air ; F(i,j,k)%mu(1)=mu_air ; F(i,j,k)%sigy(1)=sigy_air ; F(i,j,k)%rho0(1)=rho0_air 
!   F(i,j,k)%p_sge(2)=p_sge_Al  ; F(i,j,k)%g_sge(2)=g_sge_Al  ; F(i,j,k)%mu(2)=mu_Al  ; F(i,j,k)%sigy(2)=sigy_Al  ; F(i,j,k)%rho0(2)=rho0_Al 
!   F(i,j,k)%p_sge(3)=p_sge_Ti  ; F(i,j,k)%g_sge(3)=g_sge_Ti  ; F(i,j,k)%mu(3)=mu_Ti  ; F(i,j,k)%sigy(3)=sigy_Ti  ; F(i,j,k)%rho0(3)=rho0_Ti 
!   enddo; enddo; enddo
!
!   forall(i=1:nx,j=1:ny,k=1:nz) F(i,j,k)%mu(1:3)=0.0_PR
! 
!   write(iout,*) '   initialisation de fv,p,v,r' 
!
!   do k=1,Nz ; do j=1,ny ; do i=1,nx 
!     if(x(i).lt.0.001_PR)then
!      F(i,j,k)%fl(1)=1.0_PR
!     elseif(x(i).lt.0.003_PR)then
!      F(i,j,k)%fl(2)=1.0_PR ; F(i,j,k)%vx=700.0_PR
!     elseif(x(i).lt.0.0128_PR)then
!      F(i,j,k)%fl(3)=1.0_PR
!     else
!      F(i,j,k)%fl(1)=1.0_PR
!     endif
!     F(i,j,k)%pl(:)=1.0e5_PR
!   enddo; enddo; enddo
!   
!
!   do k=1,Nz ; do j=1,ny ; do i=1,nx 
!      do iph=1,Nl
!        F(i,j,k)%rhl(iph)=F(i,j,k)%rho0(iph)
!      enddo
!   enddo; enddo; enddo
!
!end subroutine test_Ndanou9

subroutine test_SOD_air(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph

   write(iout,*) ' Cas test tube à choc de SOD air pure'

   M%Nt=100
   M%dt=1.0e-5_PR

   call init_euler(200,1,1,1,1.0_PR,1.0_PR,1.0_PR,M)

   ph(1)%typ=1 ; ph(1)%nom='air' ; call init_phase(M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
      M%MF(i,j,k)%F(1)%p=1.0e5_PR
   enddo; enddo; enddo

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
      M%MF(i,j,k)%F(1)%f=1.0_PR
      M%MF(i,j,k)%vx=0.0_PR
      if(i.le.M%Nx/2)then
         M%MF(i,j,k)%F(1)%rh=1.0_PR*M%MF(i,j,k)%F(1)%rh0
         M%MF(i,j,k)%F(1)%p=1.0e5_PR 
      else
         M%MF(i,j,k)%F(1)%rh=0.125_PR*M%MF(i,j,k)%F(1)%rh0
         M%MF(i,j,k)%F(1)%p=1.0e4_PR
      endif
   enddo; enddo; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

   forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

   print*, 'energie_hydro'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
     M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

    call MAJ_meca(M%MF(i,j,k))

   enddo ; enddo ; enddo

   write(iout,*) ' MAJ mixt'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

   call MAJ_mixt(M%MF(i,j,k))

   enddo ; enddo ; enddo

end subroutine test_SOD_air

subroutine test_convection(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph
 
   write(iout,*) ' cas test convection'  

   call init_euler(200,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)
   M%Nt=2000
   M%dt=1.0e-6_PR

   ph(1)%typ=1 ; ph(1)%nom='air' 
   ph(2)%typ=1 ; ph(2)%nom='eau'
   call init_phase(M)

  write(iout,*) '   initialisation de fv,p,v,r'
 
   do k=1,M%Nz ; do j=1,M%Ny
   do i=1,M%Nx/2
      M%MF(i,j,k)%F(1)%f=1.0e-8_PR
      M%MF(i,j,k)%F(2)%f=1.0_PR-1.0e-8_PR
   enddo
   do i=M%Nx/2+1,M%Nx
      M%MF(i,j,k)%F(1)%f=1.0_PR-1.0e-8_PR
      M%MF(i,j,k)%F(2)%f=1.0e-8_PR
   enddo
   enddo; enddo

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
      M%MF(i,j,k)%F(:)%p=1.0e5_PR
      M%MF(i,j,k)%F(1)%rh=1.0_PR
      M%MF(i,j,k)%F(2)%rh=1000.0_PR
      M%MF(i,j,k)%vx=100.0_PR
   enddo; enddo; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)
   forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
     !M%MF(i,j,k)%elh(iph)=ph(iph)%energie_hydro(iph,M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph),M%MF(i,j,k))
     !M%MF(i,j,k)%elh(iph)=energie_hydro_sge(M%MF(i,j,k)%g_sge(iph),M%MF(i,j,k)%p_sge(iph),M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph))
     M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

    call MAJ_meca(M%MF(i,j,k))

   enddo ; enddo ; enddo

   write(iout,*) ' MAJ mixt'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

   call MAJ_mixt(M%MF(i,j,k))

   enddo ; enddo ; enddo

end subroutine test_convection

subroutine test_eau_air(M)

   implicit none
    type(mesh), intent(inout) :: M 
   integer :: i,j,k,iph

    write(iout,*) ' Cas test choc eau-air'

    !call init_euler(1000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

    call init_euler(500,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

    M%Nt=240
    M%dt=1.0e-6_PR

    ph(1)%typ=1 ; ph(1)%nom='air' 
    ph(2)%typ=1 ; ph(2)%nom='eau'
    call init_phase(M)

     write(iout,*) '   initialisation de fv,p,v,r'

     do k=1,M%Nz ; do j=1,M%Ny
     do i=1,M%Nx
        if(M%x(i).le.0.75_PR)then
           M%MF(i,j,k)%F(1)%f=1.0e-6_PR
           M%MF(i,j,k)%F(2)%f=1.0_PR-1.0e-6_PR
           M%MF(i,j,k)%F(:)%p=1.0e9_PR
        endif
     enddo
     do i=1,M%Nx
        if(M%x(i).gt.0.75_PR)then
           M%MF(i,j,k)%F(1)%f=1.0_PR-1.0e-6_PR
           M%MF(i,j,k)%F(2)%f=1.0e-6_PR
           M%MF(i,j,k)%F(:)%p=1.0e5_PR
        endif
     enddo
     enddo ; enddo

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
        M%MF(i,j,k)%F(1)%rh=1.0_PR
        M%MF(i,j,k)%F(2)%rh=1000.0_PR
     enddo; enddo; enddo

     forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

     forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
       M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      call MAJ_meca(M%MF(i,j,k))

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

     call MAJ_mixt(M%MF(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_eau_air

subroutine test_eau_air_strong(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph

    write(iout,*) ' Cas test choc eau-air'

    call init_euler(2000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)
    M%Nt=90
    M%dt=1.0e-7_PR

   ph(1)%typ=1 ; ph(1)%nom='air' 
   ph(2)%typ=1 ; ph(2)%nom='eau'
   call init_phase(M)

    write(iout,*) '   initialisation de fv,p,v,r'
 
    do k=1,M%Nz ; do j=1,M%Ny
       do i=1,M%Nx
          if(M%x(i).le.0.75_PR)then
             M%MF(i,j,k)%F(1)%f=1.0e-6_PR
             M%MF(i,j,k)%F(2)%f=1.0_PR-1.0e-6_PR
             M%MF(i,j,k)%F(:)%p=1.0e12_PR
          endif
       enddo
       do i=1,M%Nx
          if(M%x(i).gt.0.75_PR)then
             M%MF(i,j,k)%F(1)%f=1.0_PR-1.0e-6_PR
             M%MF(i,j,k)%F(2)%f=1.0e-6_PR
             M%MF(i,j,k)%F(:)%p=1.0e6_PR
          endif
       enddo
     enddo; enddo

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
        M%MF(i,j,k)%F(1)%rh=10.0_PR
        M%MF(i,j,k)%F(2)%rh=1000.0_PR
     enddo; enddo; enddo

     forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

     forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
       !M%MF(i,j,k)%elh(iph)=ph(iph)%energie_hydro(iph,M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph),M%MF(i,j,k))
       !M%MF(i,j,k)%elh(iph)=energie_hydro_sge(M%MF(i,j,k)%g_sge(iph),M%MF(i,j,k)%p_sge(iph),M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph))
        M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)

     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      call MAJ_meca(M%MF(i,j,k))

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

     call MAJ_mixt(M%MF(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_eau_air_strong


subroutine test_epoxy_spinel(M)

   implicit none
    type(mesh), intent(inout) :: M 
   integer :: i,j,k,iph

    write(iout,*) ' Cas test choc epoxy-spinel'

    !call init_euler(1000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

    call init_euler(500,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

    M%Nt=80
    M%dt=1.0e-6_PR

    ph(1)%typ=1 ; ph(1)%nom='epoxy' 
    ph(2)%typ=1 ; ph(2)%nom='spinel'
    call init_phase(M)

     write(iout,*) '   initialisation de fv,p,v,r'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
        M%MF(i,j,k)%F(1)%f=0.5954_PR
        M%MF(i,j,k)%F(2)%f=1.0_PR-M%MF(i,j,k)%F(1)%f
        if(M%x(i).le.0.6_PR)then
           M%MF(i,j,k)%F(:)%p=1.0e10_PR
        else
           M%MF(i,j,k)%F(:)%p=1.0e5_PR
        endif
        M%MF(i,j,k)%F(:)%rh=M%MF(i,j,k)%F(:)%rh0

     enddo ; enddo ; enddo

     forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

     forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
       M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      call MAJ_meca(M%MF(i,j,k))

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

     call MAJ_mixt(M%MF(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_epoxy_spinel

subroutine test_epoxy_spinel_strong(M)

   implicit none
    type(mesh), intent(inout) :: M 
   integer :: i,j,k,iph

    write(iout,*) ' Cas test choc epoxy-spinel strong'

    !call init_euler(1000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

    call init_euler(500,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

    M%Nt=80
    M%dt=1.0e-6_PR

    ph(1)%typ=1 ; ph(1)%nom='epoxy' 
    ph(2)%typ=1 ; ph(2)%nom='spinel'
    call init_phase(M)

     write(iout,*) '   initialisation de fv,p,v,r'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
        M%MF(i,j,k)%F(1)%f=0.5954_PR
        M%MF(i,j,k)%F(2)%f=1.0_PR-M%MF(i,j,k)%F(1)%f
        if(M%x(i).le.0.6_PR)then
           M%MF(i,j,k)%F(:)%p=2.0e11_PR
        else
           M%MF(i,j,k)%F(:)%p=1.0e5_PR
        endif
        M%MF(i,j,k)%F(:)%rh=M%MF(i,j,k)%F(:)%rh0

     enddo ; enddo ; enddo

     forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

     forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
       M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      call MAJ_meca(M%MF(i,j,k))

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

     call MAJ_mixt(M%MF(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_epoxy_spinel_strong

subroutine test_rar_rar(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph

   write(iout,*) ' Cas test rar-rar'

   M%Nt=100
   M%dt=1.0e-5_PR

   call init_euler(200,1,1,1,1.0_PR,1.0_PR,1.0_PR,M)

   ph(1)%typ=1 ; ph(1)%nom='air' ; call init_phase(M)

   write(iout,*) '   initialisation de fv,p,v,r' 

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
      M%MF(i,j,k)%F(1)%f=1.0_PR
      M%MF(i,j,k)%F(:)%p=1.0e5_PR
      do iph=1,M%Nl
        M%MF(i,j,k)%F(iph)%rh=M%MF(i,j,k)%F(iph)%rh0
      enddo
   enddo; enddo; enddo

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   M%MF(i,j,k)%F(1)%rh=1.0_PR*M%MF(i,j,k)%F(1)%rh0
   if(i.le.M%Nx/2)then
   M%MF(i,j,k)%vx=-100.0_PR
   else
   M%MF(i,j,k)%vx=100.0_PR
   endif
   enddo; enddo ; enddo

   forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

   forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
     !M%MF(i,j,k)%elh(iph)=ph(iph)%energie_hydro(iph,M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph),M%MF(i,j,k))
     !M%MF(i,j,k)%elh(iph)=energie_hydro_sge(M%MF(i,j,k)%g_sge(iph),M%MF(i,j,k)%p_sge(iph),M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph))
     M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)

   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

    call MAJ_meca(M%MF(i,j,k))

   enddo ; enddo ; enddo

   write(iout,*) ' MAJ mixt'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

   call MAJ_mixt(M%MF(i,j,k))

   enddo ; enddo ; enddo

end subroutine test_rar_rar

subroutine test_cavitation(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph

    write(iout,*) ' Cas test cavitation'

    call init_euler(1000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

    M%Nt=200
    M%dt=1.0e-5_PR

    ph(1)%typ=1 ; ph(1)%nom='air' 
    ph(2)%typ=1 ; ph(2)%nom='eau'
    call init_phase(M)

     write(iout,*) '   initialisation de fv,p,v,r'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
        M%MF(i,j,k)%F(1)%f=1.0e-2_PR
        M%MF(i,j,k)%F(2)%f=1.0_PR-1.0e-2_PR
        M%MF(i,j,k)%F(:)%p=1.0e5_PR
     enddo; enddo; enddo

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
        if(M%x(i).le.0.5_PR)then
           M%MF(i,j,k)%vx=-100.0_PR
        elseif(M%x(i).gt.0.5_PR)then
           M%MF(i,j,k)%vx=100.0_PR
        endif
     enddo; enddo; enddo

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
        M%MF(i,j,k)%F(1)%rh=1.0_PR*M%MF(i,j,k)%F(1)%rh0
        M%MF(i,j,k)%F(2)%rh=M%MF(i,j,k)%F(2)%rh0
     enddo; enddo; enddo

     forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

     forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
        M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      call MAJ_meca(M%MF(i,j,k))

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

     call MAJ_mixt(M%MF(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_cavitation

subroutine test_cavitation2(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: i,j,k,iph

    write(iout,*) ' Cas test cavitation 2'

    call init_euler(2000,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

    M%Nt=100
    M%dt=1.0e-6_PR
    M%dtmin=5.0e-8_PR

    ph(1)%typ=1 ; ph(1)%nom='air' 
    ph(2)%typ=1 ; ph(2)%nom='Cu'
    call init_phase(M)

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
        M%MF(i,j,k)%F(1)%f=1.0e-4_PR
        M%MF(i,j,k)%F(2)%f=1.0_PR-1.0e-4_PR
        M%MF(i,j,k)%F(:)%p=1.0e5_PR
        if(M%x(i).le.0.5_PR)then
           M%MF(i,j,k)%vx=-3000.0_PR
        elseif(M%x(i).gt.0.5_PR)then
           M%MF(i,j,k)%vx=3000.0_PR
        endif
     enddo; enddo; enddo

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
        M%MF(i,j,k)%F(1)%rh=1.0_PR
        M%MF(i,j,k)%F(2)%rh=8900.0_PR
     enddo; enddo; enddo

     forall(i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:M%Nl)%f*M%MF(i,j,k)%F(1:M%Nl)%rh)

     forall(iph=1:M%Nl,i=1:M%Nx,j=1:M%Ny,k=1:M%Nz) M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh


     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
       !M%MF(i,j,k)%elh(iph)=ph(iph)%energie_hydro(iph,M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph),M%MF(i,j,k))
       !M%MF(i,j,k)%elh(iph)=energie_hydro_sge(M%MF(i,j,k)%g_sge(iph),M%MF(i,j,k)%p_sge(iph),M%MF(i,j,k)%rhl(iph),M%MF(i,j,k)%pl(iph))
        M%MF(i,j,k)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh,M%MF(i,j,k)%F(iph)%p)
     enddo ; enddo ; enddo ; enddo

     write(iout,*) ' MAJ meca'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

      call MAJ_meca(M%MF(i,j,k))

     enddo ; enddo ; enddo

     write(iout,*) ' MAJ mixt'

     do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx

     call MAJ_mixt(M%MF(i,j,k))

     enddo ; enddo ; enddo

end subroutine test_cavitation2




!subroutine test_compression
!
!   use mod_data, only : Nl,Nx,Ny,Nz,Nt,t,dt,fl,Yl,pl,rhl,elh,ele,&
!                        al,bl,cl,sigl,sig,&
!                        rh,vx,vy,vz,c,p,e,&
!                        x,y,z
!
!   implicit none
!   integer :: i,j,k,iph
!
!   write(iout,*) ' Cas_ test cylindre en compression'
!
!   call init_euler(10,1,1,3)
!
!   !!!------------------     cas test cylindre compression  ---------------------------------------
!   !! phase 1 : Air 2 : Aluminium 3 : Titane
!   p_sge(1)=p_sge_air ; g_sge(1)=g_sge_air ; mu(1)=mu_air ; sigy(1)=sigy_air ; rho0(1)=rho0_air 
!   p_sge(2)=p_sge_Al  ; g_sge(2)=g_sge_Al  ; mu(2)=mu_Al  ; sigy(2)=sigy_Al  ; rho0(2)=rho0_Al 
!   p_sge(3)=p_sge_Ti  ; g_sge(3)=g_sge_Ti  ; mu(3)=mu_Ti  ; sigy(3)=sigy_Ti  ; rho0(3)=rho0_Ti 
! 
!   write(iout,*) '   initialisation de fv,p,v,r' 
!   do i=1,nx
!      fl(2,i,:,:)=1.0_PR
!   enddo   
!   pl(:,:,:,:)=1.0e5_PR
!   do k=1,Nz ; do j=1,ny ; do i=1,nx 
!      do iph=1,Nl
!        rhl(iph,i,j,k)=rho0(iph)
!      enddo
!   enddo; enddo; enddo
!
!   !!!---compression selon -z
!   al(:,1,:,:,:)=1.1_PR
!
!
!end subroutine test_compression

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


subroutine test_relax(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: iph, i, j, k

   write(iout,*) ' Cas test relaxation pression'

   call init_euler(1,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
   M%MF(i,j,k)%F(1)%p_sge=p_sge_air ; M%MF(i,j,k)%F(1)%g_sge=g_sge_air ; M%MF(i,j,k)%F(1)%mu=mu_air 
   M%MF(i,j,k)%F(1)%sigy =sigy_air  ; M%MF(i,j,k)%F(1)%rh0=rho0_air    ;
   M%MF(i,j,k)%F(2)%p_sge=p_sge_eau ; M%MF(i,j,k)%F(2)%g_sge=g_sge_eau ; M%MF(i,j,k)%F(2)%mu=mu_eau
   M%MF(i,j,k)%F(2)%sigy =sigy_eau  ; M%MF(i,j,k)%F(2)%rh0=rho0_eau 
   enddo ; enddo ; enddo

   M%MF(1,1,1)%F(1)%rh=1.0053_PR
   M%MF(1,1,1)%F(2)%rh=874.11_PR
   M%MF(1,1,1)%F(1)%p=100749.0_PR
   M%MF(1,1,1)%F(2)%p=270044072.0_PR
   M%MF(1,1,1)%F(1)%f=0.99236_PR
   M%MF(1,1,1)%F(2)%f=1.0_PR-M%MF(1,1,1)%F(1)%f
   M%MF(1,1,1)%p=2.162638e6_PR 
   M%MF(1,1,1)%c=334.91569_PR
   M%MF(1,1,1)%vx=0.0_PR
   M%MF(1,1,1)%u=966143.18905909313_PR
   M%MF(1,1,1)%e=M%MF(1,1,1)%u

   M%MF(1,1,1)%rh=sum(M%MF(1,1,1)%F(1:M%Nl)%f*M%MF(1,1,1)%F(1:M%Nl)%rh)
  
   do iph=1,M%Nl
     !M%MF(1,1,1)%elh(iph)=ph(iph)%energie_hydro(iph,M%MF(1,1,1)%rhl(iph),M%MF(1,1,1)%pl(iph),M%MF(1,1,1))
     M%MF(1,1,1)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(1,1,1)%F(iph),M%MF(1,1,1)%F(iph)%rh,M%MF(1,1,1)%F(iph)%p) 
     M%MF(1,1,1)%F(iph)%Y=M%MF(1,1,1)%F(iph)%f*M%MF(1,1,1)%F(iph)%rh/M%MF(1,1,1)%rh
   enddo 
   M%MF(1,1,1)%c=soundspeed_mixt(M%MF(1,1,1))

   verb=.true.

   call relax_p(M%MF(1,1,1))

   stop
 
end subroutine test_relax

subroutine test_relax2(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: iph, i, j, k

   write(iout,*) ' Cas test relaxation pression'

   call init_euler(1,1,1,2,1.0_PR,1.0_PR,1.0_PR,M)

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx 
   M%MF(i,j,k)%F(1)%p_sge=p_sge_Cu  ; M%MF(i,j,k)%F(1)%g_sge=g_sge_Cu ; M%MF(i,j,k)%F(1)%mu=mu_Cu 
   M%MF(i,j,k)%F(1)%sigy=sigy_Cu    ; M%MF(i,j,k)%F(1)%rh0=rho0_Cu
   M%MF(i,j,k)%F(2)%p_sge=p_sge_air ; M%MF(i,j,k)%F(2)%g_sge=g_sge_air ; M%MF(i,j,k)%F(2)%mu=mu_air
   M%MF(i,j,k)%F(2)%sigy=sigy_air   ; M%MF(i,j,k)%F(2)%rh0=rho0_air 
   enddo ; enddo ; enddo

   M%MF(1,1,1)%F(1)%rh=8776.5334167879500_PR
   M%MF(1,1,1)%F(2)%rh=1.0482138728300139_PR
   M%MF(1,1,1)%F(1)%p=21501633.5864868_PR
   M%MF(1,1,1)%F(2)%p=857.76873357
   M%MF(1,1,1)%F(1)%f=0.7045895_PR
   M%MF(1,1,1)%F(2)%f=1.0_PR-M%MF(1,1,1)%F(1)%f
   M%MF(1,1,1)%p=0.0_PR 
   M%MF(1,1,1)%c=0.0_PR
   M%MF(1,1,1)%vx=55.053112917261807_PR
   M%MF(1,1,1)%u=5106671.5340925017_PR
   M%MF(1,1,1)%e=5108186.9567134418_PR
   M%MF(1,1,1)%F(1)%eh=5106927.3197702039_PR
   M%MF(1,1,1)%F(2)%eh=2046.5133529692748_PR
   M%MF(1,1,1)%F(1)%eh0=5106744.2590157716_PR
   M%MF(1,1,1)%F(2)%eh0=3506.4075082757072_PR

   M%MF(1,1,1)%rh=sum(M%MF(1,1,1)%F(1:M%Nl)%f*M%MF(1,1,1)%F(1:M%Nl)%rh)
  
   do iph=1,M%Nl
     M%MF(1,1,1)%F(iph)%Y=M%MF(1,1,1)%F(iph)%f*M%MF(1,1,1)%F(iph)%rh/M%MF(1,1,1)%rh
   enddo 

   verb=.true.
   print*, ':( fl1=', M%MF(1,1,1)%F(1)%f

   call relax_p(M%MF(1,1,1))

   stop
 
end subroutine test_relax2

subroutine test_relax3(M)

   implicit none
   type(mesh), intent(inout) :: M
   integer :: iph, i, j, k
   real(PR) :: p(1:1), p0, pfin, dp, Jaco(1:1,1:1), func(1:1)

   write(iout,*) ' Cas test relaxation pression'

   call init_euler(1,1,1,3,1.0_PR,1.0_PR,1.0_PR,M)

   !!!---type de phase :
   ph(1)%typ=1 ; ph(1)%nom='air'
   ph(2)%typ=1 ; ph(2)%nom='Cug'
   ph(3)%typ=1 ; ph(3)%nom='epoxy'

   call init_phase(M)

   M%MF(1,1,1)%F(1)%rh =4.8528e-2_PR
   M%MF(1,1,1)%F(2)%rh =11.03616_PR
   M%MF(1,1,1)%F(3)%rh =1169.4885_PR

   M%MF(1,1,1)%F(1)%p  =800674.48_PR
   M%MF(1,1,1)%F(2)%p  =826484.38_PR
   M%MF(1,1,1)%F(3)%p  =-289853.0_PR

   M%MF(1,1,1)%F(1)%f  =1.752828e-6_PR
   M%MF(1,1,1)%F(2)%f  =3.979411774e-9_PR
   M%MF(1,1,1)%F(3)%f  =0.999998_PR

   M%MF(1,1,1)%rh=sum(M%MF(1,1,1)%F(:)%f*M%MF(1,1,1)%F(:)%rh)
   do iph=1,M%Nl
     M%MF(1,1,1)%F(iph)%Y=M%MF(1,1,1)%F(iph)%f*M%MF(1,1,1)%F(iph)%rh/M%MF(1,1,1)%rh
   enddo 

   do iph=1,Nl
      M%MF(1,1,1)%F(iph)%eh=ph(iph)%energie_hydro(M%MF(1,1,1)%F(iph),M%MF(1,1,1)%F(iph)%rh,M%MF(1,1,1)%F(iph)%p)
   enddo 

   print*, 'elh1:', M%MF(1,1,1)%F(1)%eh, 'vs 41247835.23'
   print*, 'elh2:', M%MF(1,1,1)%F(2)%eh, 'vs 23257.377'
   print*, 'elh3:', M%MF(1,1,1)%F(3)%eh, 'vs 7700880.05'

   call MAJ_mixt(M%MF(1,1,1))
   M%MF(1,1,1)%vx=0.0_PR

   verb=.true.
   call relax_p(M%MF(1,1,1))

!   open(unit=12, file='testrelax.dat', status='replace')
!
!
!   p0=-1.0e6_PR
!   pfin=1.0e6_PR
!
!   dp=(pfin-p0)/real(1000,PR)   
!
!   do i=1,1000
!          
!
!      p(1)=p0+real(i-1,PR)*dp
!      Jaco(1,1)=0.0_PR
!
!      call F_relaxp_sge_1(1,p,func)
!      call J_relaxp_sge_1(1,p,Jaco)
!
!      write(12,*) p(1), func(1), Jaco(1,1)
!
!   enddo
!
!   close(12)

   stop
 
end subroutine test_relax3

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

subroutine NOREC(N,MF,U,G)

   !!---pas de reconstruction

   implicit none
   integer, intent(in) :: N
   type(multifluide), intent(in) :: MF(1:N)
   real(PR), intent(in) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
   integer :: i, ip1, iph

   do i=1,N

      ip1=min(i+1,N) 

      MFrec_L(i)=MF(i)
      MFrec_R(i)=MF(ip1)
      !Urec_L(1:Nl+4,i)=U(1:Nl+4,i)
      !Urec_R(1:Nl+4,i)=U(1:Nl+4,ip1)       
      !do iph=1,Nl 
      ! Grec_L(iph,1:11,i)=G(iph,1:11,i)
      ! Grec_R(iph,1:11,i)=G(iph,1:11,ip1)
      !enddo

      !!!---CONVERSION CONSERVATIVE -> PRIMITIVE
      !call conservative2primitive(Nl,UL(:,i),GL(:,:,i),FL(i))  
      !call MAJ_meca(Nl,FL(i)) ; call MAJ_mixt(Nl,FL(i))
      !call conservative2primitive(Nl,UR(:,i),GR(:,:,i),FR(i))   
      !call MAJ_meca(Nl,FR(i)) ; call MAJ_mixt(Nl,FR(i))

   enddo

end subroutine NOREC

!!! subroutine NOREC
!!! 
!!!    !!---pas de reconstruction
!!! 
!!!    implicit none
!!!    integer :: i, ip1, iph, N
!!! 
!!!    N=N1D
!!! 
!!!    do i=1,N
!!! 
!!!       ip1=min(i+1,N) 
!!! 
!!!       MFrec_L(i)=MF1D(i)
!!!       MFrec_R(i)=MF1D(ip1)
!!! 
!!!    enddo
!!! 
!!! end subroutine NOREC

subroutine MUSCL

   !!---reconstruction MUSCL sur variables conservatives

   implicit none

   integer :: i,j,k,ip1, im1, iph, i2, N
   real(PR) :: rtemp(1:10)
   integer :: itemp(1:10)
   real(PR) :: sumfl_L, sumfl_R

   call cpu_time(t1)   

   N=N1D

   do i=1,N

      im1=max(i-1,1) ; ip1=min(i+1,N)

      MFrec_L(i)=MF1D(i)
      MFrec_R(i)=MF1D(i)

      !!!!---variables primitives
 
      !!!---reconstruction de la pression
      call Reco_muscl3(MF1D(im1)%p,MF1D(i)%p,MF1D(ip1)%p,MFrec_L(i)%p,MFrec_R(i)%p)

      !!!---reconstruction fractions volumiques 
      if(Nl.gt.1)then

         rtemp(1:Nl)=MF1D(i)%F(1:Nl)%f     
         forall(iph=1:Nl) itemp(iph)=iph 

         call hpsort(Nl,rtemp(1:Nl),itemp(1:Nl))

         sumfl_L=0.0_PR
         sumfl_R=0.0_PR

         do j=1,Nl-1

           iph=itemp(j)       

           call Reco_muscl3(MF1D(im1)%F(iph)%f,MF1D(i)%F(iph)%f,MF1D(ip1)%F(iph)%f,MFrec_L(i)%F(iph)%f,MFrec_R(i)%F(iph)%f)

           sumfl_L=sumfl_L+MFrec_L(i)%F(iph)%f
           sumfl_R=sumfl_R+MFrec_R(i)%F(iph)%f

         enddo
          
         iph=itemp(Nl)
         MFrec_L(i)%F(iph)%f=1.0_PR-sumfl_L
         MFrec_R(i)%F(iph)%f=1.0_PR-sumfl_R

      endif !!! Nl > 1

      do iph=1,Nl

        !!!---reconstruction de la densité
        call Reco_muscl3(MF1D(im1)%F(iph)%rh,MF1D(i)%F(iph)%rh,MF1D(ip1)%F(iph)%rh,MFrec_L(i)%F(iph)%rh,MFrec_R(i)%F(iph)%rh)

        if(MFrec_L(i)%F(iph)%rh.lt.0.0_PR)then
          print*, 'PROBLEM MUSCL: FL%rh<0'
          print*, MF1D(im1)%F(iph)%rh,MF1D(i)%F(iph)%rh,MF1D(ip1)%F(iph)%rh
          stop
        elseif(MFrec_R(i)%F(iph)%rh.lt.0.0_PR)then
          print*, 'PROBLEM MUSCL: FR%rh<0'
          print*, MF1D(im1)%F(iph)%rh,MF1D(i)%F(iph)%rh,MF1D(ip1)%F(iph)%rh
          stop
        endif

        !!!---reconstruction de la pression
        do k=1,3
           MFrec_L(i)%F(iph)%sig(k,k)=MFrec_L(i)%F(iph)%sig(k,k)+MFrec_L(i)%F(iph)%p
           MFrec_R(i)%F(iph)%sig(k,k)=MFrec_R(i)%F(iph)%sig(k,k)+MFrec_R(i)%F(iph)%p
        enddo

        MFrec_L(i)%F(iph)%p=MFrec_L(i)%p
        MFrec_R(i)%F(iph)%p=MFrec_R(i)%p

        do k=1,3
           MFrec_L(i)%F(iph)%sig(k,k)=MFrec_L(i)%F(iph)%sig(k,k)-MFrec_L(i)%F(iph)%p
           MFrec_R(i)%F(iph)%sig(k,k)=MFrec_R(i)%F(iph)%sig(k,k)-MFrec_R(i)%F(iph)%p
        enddo

        !!!!---energie à partir de rho et p
        MFrec_R(i)%F(iph)%eh=ph(iph)%energie_hydro(MFrec_R(i)%F(iph),MFrec_R(i)%F(iph)%rh,MFrec_R(i)%F(iph)%p)
        MFrec_L(i)%F(iph)%eh=ph(iph)%energie_hydro(MFrec_L(i)%F(iph),MFrec_L(i)%F(iph)%rh,MFrec_L(i)%F(iph)%p)

      enddo  

      !!!!---reconstruction de la vitesse
      call Reco_muscl3(MF1D(im1)%vx,MF1D(i)%vx,MF1D(ip1)%vx,MFrec_L(i)%vx,MFrec_R(i)%vx)

      !!!---MAJ
      MFrec_R(i)%rh=sum(MFrec_R(i)%F(1:Nl)%f*MFrec_R(i)%F(1:Nl)%rh)
      do iph=1,Nl
      MFrec_R(i)%F(iph)%Y=MFrec_R(i)%F(iph)%f*MFrec_R(i)%F(iph)%rh/MFrec_R(i)%rh
      enddo
      MFrec_R(i)%u=sum(MFrec_R(i)%F(1:Nl)%Y*MFrec_R(i)%F(1:Nl)%eh)
      MFrec_R(i)%e=MFrec_R(i)%u+MFrec_R(i)%ee+0.5_PR*MFrec_R(i)%vx**2

      MFrec_L(i)%rh=sum(MFrec_L(i)%F(1:Nl)%f*MFrec_L(i)%F(1:Nl)%rh)
      do iph=1,Nl
      MFrec_L(i)%F(iph)%Y=MFrec_L(i)%F(iph)%f*MFrec_L(i)%F(iph)%rh/MFrec_L(i)%rh
      enddo
      MFrec_L(i)%u=sum(MFrec_L(i)%F(1:Nl)%Y*MFrec_L(i)%F(1:Nl)%eh)
      MFrec_L(i)%e=MFrec_L(i)%u+MFrec_L(i)%ee+0.5_PR*MFrec_L(i)%vx**2

  enddo

   do i=1,N
      !call relax_p(Frec_L(i))
      !call relax_p(Frec_R(i))
      !call MAJ_meca(Frec_L(i))
      !call MAJ_meca(Frec_R(i))
!
      !  do k=1,3
      !  Frec_L(i)%sigl(:,k,k)=-Frec_L(i)%p
      !  Frec_R(i)%sigl(:,k,k)=-Frec_R(i)%p
      !  enddo

      do j=1,3
        do k=1,3 
          MFrec_L(i)%sig(j,k)=sum(MFrec_L(i)%F(1:Nl)%f*MFrec_L(i)%F(1:Nl)%sig(j,k))
          MFrec_R(i)%sig(j,k)=sum(MFrec_R(i)%F(1:Nl)%f*MFrec_R(i)%F(1:Nl)%sig(j,k))
        enddo
      enddo
      MFrec_L(i)%c=soundspeed_mixt(MFrec_L(i))
      MFrec_R(i)%c=soundspeed_mixt(MFrec_R(i))
   enddo

   call cpu_time(t2) ;  t_reco=t_reco+t2-t1

end subroutine MUSCL

!!subroutine MUSCL
!!
!!   !!---reconstruction MUSCL sur variables conservatives
!!
!!   implicit none
!!
!!   integer :: i, ip1, im1, ip2, iph, i2, N
!!
!!   N=N1D
!!
!!   do i=1,N
!!      call conservative2primitive(U(1:Nl+4,i),G(1:Nl,1:11,i),MF1D(i)) 
!!   enddo 
!!
!!   do i=1,N
!!
!!      im1=max(i-1,1) ; ip1=min(i+1,N) ; ip2=min(i+2,N)
!!
!!      Frec_L(i)=MF1D(i)
!!      Frec_R(i)=MF1D(ip1)
!!
!!      !!!!---variables primitives
!!      do iph=1,Nl
!!
!!        !!!---reconstruction de la densité
!!        call Reco_muscl(MF1D(im1)%rhl(iph),MF1D(i)%rhl(iph),MF1D(ip1)%rhl(iph),MF1D(ip2)%rhl(iph),Frec_L(i)%rhl(iph),Frec_R(i)%rhl(iph))
!!        !!!---reconstruction de l'énergie
!!        !call Reco_muscl(MF1D(im1)%elh(iph),MF1D(i)%elh(iph),MF1D(ip1)%elh(iph),MF1D(ip2)%elh(iph),Frec_L(i)%elh(iph),Frec_R(i)%elh(iph))
!!        !!!---reconstruction de la pression
!!        call Reco_muscl(MF1D(im1)%pl(iph),MF1D(i)%pl(iph),MF1D(ip1)%pl(iph),MF1D(ip2)%pl(iph),Frec_L(i)%pl(iph),Frec_R(i)%pl(iph))
!!        !!!---energie à partir de rho et p
!!        Frec_R(i)%elh(iph)=ph(iph)%energie_hydro(iph,Frec_R(i),Frec_R(i)%rhl(iph),Frec_R(i)%pl(iph))
!!        Frec_L(i)%elh(iph)=ph(iph)%energie_hydro(iph,Frec_L(i),Frec_L(i)%rhl(iph),Frec_L(i)%pl(iph))
!!
!!      enddo  
!!
!!        call Reco_muscl(MF1D(im1)%vx,MF1D(i)%vx,MF1D(ip1)%vx,MF1D(ip2)%vx,Frec_L(i)%vx,Frec_R(i)%vx)
!!
!!        Frec_R(i)%rh=sum(Frec_R(i)%fl(1:Nl)*Frec_R(i)%rhl(1:Nl))
!!        Frec_R(i)%u=sum(Frec_R(i)%fl(1:Nl)*Frec_R(i)%elh(1:Nl))
!!        Frec_R(i)%e=Frec_R(i)%u+0.5_PR*Frec_R(i)%vx**2
!!        Frec_R(i)%P=sum(Frec_R(i)%fl(1:Nl)*Frec_R(i)%pl(1:Nl))
!!
!!        Frec_L(i)%rh=sum(Frec_L(i)%fl(1:Nl)*Frec_L(i)%rhl(1:Nl))
!!        Frec_L(i)%u=sum(Frec_L(i)%fl(1:Nl)*Frec_L(i)%elh(1:Nl))
!!        Frec_L(i)%e=Frec_L(i)%u+0.5_PR*Frec_L(i)%vx**2
!!        Frec_L(i)%P=sum(Frec_L(i)%fl(1:Nl)*Frec_L(i)%pl(1:Nl))
!!
!!   enddo
!!
!!   !do i=1,N
!!   !   call relax_p(Frec_L(i))
!!   !   call relax_p(Frec_R(i))
!!   !enddo
!!
!!
!!   !   do i=1,N-1
!!   !
!!   !      im1=max(i-1,1) ; ip1=i+1 ; ip2=min(i+2,N)
!!   !
!!   !
!!   !
!!   !      !!!---variables conservatives : 
!!   !      do i2=1,Nl+4
!!   !        call Reco_muscl(U(i2,im1),U(i2,i),U(i2,ip1),U(i2,ip2),Urec_L(i2,i),Urec_R(i2,i))
!!   !      enddo  
!!   !      GL(:,:,i)=G(:,:,i) ; GR(:,:,i)=G(:,:,ip1) !!! consistance sur 11 ?
!!   !      FL(i)=F(i)            ; FR(i)=F(ip1) 
!!   !
!!   !      !!!---CONVERSION CONSERVATIVE -> PRIMITIVE
!!   !      call conservative2primitive(UL(:,i),GL(:,:,i),FL(i))  
!!   !      call MAJ_meca(FL(i)) ; call MAJ_mixt(FL(i))
!!   !      call conservative2primitive(UR(:,i),GR(:,:,i),FR(i))   
!!   !      call MAJ_meca(FR(i)) ; call MAJ_mixt(FR(i))
!!   !
!!   !   enddo
!!
!!end subroutine MUSCL



!subroutine MUSCL(N,F,U,G,FL,FR,UL,UR,GL,GR)
!
!   !!---reconstruction MUSCL sur variables conservatives
!
!   implicit none
!
!   integer, intent(in) :: N
!   type(fluide), intent(in) :: F(1:N)
!   real(PR), intent(in) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
!   real(PR), intent(out) :: UL(1:Nl+4,1:N), GL(1:Nl,1:11,1:N)
!   real(PR), intent(out) :: UR(1:Nl+4,1:N), GR(1:Nl,1:11,1:N)
!   type(fluide), intent(out) :: FL(1:N), FR(1:N) !!! variables primitives
!   integer :: i, ip1, im1, ip2, iph, i2
!
!   do i=1,N-1
!
!      im1=max(i-1,1) ; ip1=i+1 ; ip2=min(i+2,N)
!
!      do i2=1,Nl+4
!        call Reco_muscl(U(i2,im1),U(i2,i),U(i2,ip1),U(i2,ip2),UL(i2,i),UR(i2,i))
!      enddo  
!      GL(:,:,i)=G(:,:,i) ; GR(:,:,i)=G(:,:,ip1) !!! consistance sur 11 ?
!      FL(i)=F(i)            ; FR(i)=F(ip1) 
!      !!!---CONVERSION CONSERVATIVE -> PRIMITIVE
!      call conservative2primitive(UL(:,i),GL(:,:,i),FL(i))  
!      call MAJ_meca(FL(i)) ; call MAJ_mixt(FL(i))
!      call conservative2primitive(UR(:,i),GR(:,:,i),FR(i))   
!      call MAJ_meca(FR(i)) ; call MAJ_mixt(FR(i))
!
!   enddo
!
!end subroutine MUSCL


!subroutine MUSCL2(Nl,N,F,U,G,FL,FR,UL,UR,GL,GR)
!
!   !!---reconstruction MUSCL sur variables primitives
!
!   implicit none
!
!   type(fluide), intent(in) :: F(1:N)
!   real(PR), intent(in) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
!   real(PR), intent(out) :: UL(1:Nl+4,1:N), GL(1:Nl,1:11,1:N)
!   real(PR), intent(out) :: UR(1:Nl+4,1:N), GR(1:Nl,1:11,1:N)
!   type(fluide), intent(out) :: FL(1:N), FR(1:N) !!! variables primitives
!   integer :: i, ip1, im1, ip2, iph, i2
!
!   do i=1,N-1
!
!      im1=max(i-1,1) ; ip1=i+1 ; ip2=min(i+2,N)
!
!      !!---Reconstruction variables primitives
!        FL=F(i) ; FR=F(ip1) 
!        !!!---reco : rho
!        call Reco_muscl(F(im1)%rh,F(i)%rh,F(ip1)%rh,F(ip2)%rh,FL%rh,FR%rh)
!        !!!---reco : vx, vy, vz
!        call Reco_muscl(F(im1)%vx,F(i)%vx,F(ip1)%vx,F(ip2)%vx,FL%vx,FR%vx)
!        call Reco_muscl(F(im1)%vy,F(i)%vy,F(ip1)%vy,F(ip2)%vy,FL%vy,FR%vy)
!        call Reco_muscl(F(im1)%vz,F(i)%vz,F(ip1)%vz,F(ip2)%vz,FL%vz,FR%vz)
!        !!!---p
!        call Reco_muscl(F(im1)%p,F(i)%p,F(ip1)%p,F(ip2)%p,FL%p,FR%p)
!        !!!---
!        !!!   ATTENTION : relaxp a rajouter :( !!!
!        !!!---pl
!        do i2=1,Nl
!          FL%pl(i2)=FL%p
!          FR%pl(i2)=FR%p
!        enddo 
!        !!!---rhl
!        do i2=1,Nl
!          FL%rhl(i2)=FL(i2)%fl(i2)*FL%rh
!          FR%rhl(i2)=FR(i2)%fl(i2)*FR%rh
!        enddo 
!        !!!---elh
!        do i2=1,Nl
!          FL%ehl(i2)=ph(iph)%energie_hydroi2,FL%rhl(i2),FL%pl(i2),FL)
!          FR%ehl(i2)=ph(iph)%energie_hydroi2,FR%rhl(i2),FR%pl(i2),FR)
!        enddo
!        !!!---Yl
!        FL%Yl(1:Nl)=FL%fl(1:Nl)*FL%rhl(1:Nl)/FL%rh
!        FR%Yl(1:Nl)=FR%fl(1:Nl)*FR%rhl(1:Nl)/FR%rh
!        !!!---u
!        FL%u=sum(FL%Yl(1:Nl)*FL%elh(1:Nl)) 
!        FR%u=sum(FR%Yl(1:Nl)*FR%elh(1:Nl)) 
!        !!!---e
!        FL%e=FL%u+0.5_PR*(FL%vx**2+FL%vy**2+FL%vz**2)
!        FR%e=FR%u+0.5_PR*(FR%vx**2+FR%vy**2+FR%vz**2)
!        !!!---c !!! abesoin de Yl, rhl, pl
!        FL%c=soundspeed_mixt(Nl,FL)
!        FR%c=soundspeed_mixt(Nl,FR)
!
!        !!!---reco uint
!        !call Reco_muscl(F(im1)%u,F(i)%u,F(ip1)%u,F(ip2)%u,FL%u,FR%u)
!        !FL%e=FL%rh*(FL%u+0.5_PR*(FL%vx**2+FL%vy**2+FL%vz**2))        
!        !FR%e=FR%rh*(FR%u+0.5_PR*(FR%vx**2+FR%vy**2+FR%vz**2))        
!        !!!---reco elh
!        !do i2=1,Nl
!        !   call Reco_muscl(F(im1)%elh(i2),F(i)%elh(i2),F(ip1)%elh(i2),F(ip2)%elh(i2),FL%elh(i2),FR%elh(i2))
!        !enddo
!        !!!---p !!! utilise rho et u
!        !FL%p=pressure_mixt(Nl,FL)
!        !FR%p=pressure_mixt(Nl,FR)
!        !!!---reco rlh
!        !FL%rhl(1:Nl)=FL%fl(1:Nl)*FL%rh
!        !FR%rhl(1:Nl)=FR%fl(1:Nl)*FR%rh
!        !!!!---reco Yl
!        !FL%Yl(1:Nl)=FL%fl(1:Nl)*FL%rhl(1:Nl)/FL%rh
!        !FR%Yl(1:Nl)=FR%fl(1:Nl)*FR%rhl(1:Nl)/FR%rh
!        !!!!---reco pl
!        !do i2=1,Nl
!        !FL%pl(i2)=pressure(i2,FL%rh,FL%u)
!        !FR%pl(i2)=pressure(i2,FR%rh,FR%u)
!        !enddo 
!
!   enddo
!
!end subroutine MUSCL2

subroutine Reco_NOREC(Wim1,Wi,Wip1,Wip2,WL,WR)

   implicit none
   real(PR), intent(in) :: Wim1, Wi,Wip1,Wip2
   real(PR), intent(out) :: WL, WR

   WL=Wi  
   WR=Wip1

end subroutine Reco_NOREC



subroutine Reco_muscl3(Wim1,Wi,Wip1,WL,WR)

   implicit none
   real(PR), intent(in) :: Wim1, Wi,Wip1
   real(PR), intent(out) :: WL, WR
   real(PR) :: GDO,GUP,sgDO,sgUP,Grad
   

   !!!---NOREC
   !WL=Wi 
   !WR=Wi

   !!!---MUSCL
   GDO=(Wi-Wim1)
   GUP=(Wip1-Wi)
   sgDO=sign(0.5_PR,GDO)
   sgUP=sign(0.5_PR,GUP)

   Grad=(sgDO+sgUP)*min(abs(GDO),abs(GUP))

   WL=Wi-0.5_PR*Grad
   WR=Wi+0.5_PR*Grad

end subroutine Reco_muscl3




subroutine Reco_muscl2(Wim1,Wi,Wip1,WL,WR)

   implicit none
   real(PR), intent(in) :: Wim1, Wi,Wip1
   real(PR), intent(out) :: WL, WR
   real(PR) :: Gradm1, Gradp1,phi,r

   !!!---NOREC
   WL=Wi 
   WR=Wi

   !!!---MUSCL
   Gradm1=(Wi-Wim1)
   Gradp1=(Wip1-Wi)

   if(abs(Gradm1).le.tolmuscl.and.abs(Gradp1).le.tolmuscl)then

      WL=Wi
      WR=Wi

   elseif(abs(Gradp1).le.tolmuscl)then

      WL=Wi
      phi=Limiteur(1.0e20_PR)
      WR=Wi+0.5_PR*phi*Gradm1

   else

      r=Gradm1/Gradp1
      phi=Limiteur(r)
      WL=Wi-0.5_PR*phi*Gradp1
      WR=Wi+0.5_PR*phi*Gradm1

   endif

end subroutine Reco_muscl2

subroutine Reco_muscl(Wim1,Wi,Wip1,Wip2,WL,WR)

   implicit none
   real(PR), intent(in) :: Wim1, Wi,Wip1,Wip2
   real(PR), intent(out) :: WL, WR
   real(PR) :: ri,rip1
   real(PR) :: PHI=1.0_PR/3.0_PR

   !!!---
   !WL=Wi 
   !WR=Wip1

   !!!---http://chimeracfd.com/programming/gryphon/muscl.html
   if(abs(Wip1-Wi).le.tolmuscl.and.abs(Wi-Wim1).le.tolmuscl)then
      WL=Wi 
   elseif(abs(Wip1-Wi).le.tolmuscl)then
      WL=Wi+0.25_PR*Limiteur(0.0_PR)*(1.0_PR-PHI)*(Wi-Wim1)
   elseif(abs(Wi-Wim1).le.tolmuscl)then
      WL=Wi+0.25_PR*Limiteur(0.0_PR)*(1.0_PR+PHI)*(Wip1-Wi)
   else
      ri=(Wi-Wim1)/(Wip1-Wi)
      WL=Wi+0.25_PR*( Limiteur(1.0_PR/ri)*(1.0_PR-PHI)*(Wi-Wim1) + Limiteur(ri)*(1.0_PR+PHI)*(Wip1-Wi) )
   endif

   if(abs(Wip2-Wip1).le.tolmuscl.and.abs(Wip1-Wi).le.tolmuscl)then
      WR=Wip1
   elseif(abs(Wip2-Wip1).le.tolmuscl)then
      WR=Wip1-0.25_PR*Limiteur(0.0_PR)*(1.0_PR+PHI)*(Wip1-Wi)
      !WR=Wip1-0.25_PR*( Limiteur(1.0e30_PR)*(1.0_PR+PHI)*(Wip1-Wi))
   elseif(abs(Wip1-Wi).le.tolmuscl)then
      WR=Wip1-0.25_PR*Limiteur(0.0_PR)*(1.0_PR-PHI)*(Wip2-Wip1)
      !WR=Wip1+0.25_PR*(Limiteur(1.0e30_PR)*(1.0_PR-PHI)*(Wip2-Wip1) )
   else
      rip1=(Wip1-Wi)/(Wip2-Wip1)
      WR=Wip1-0.25_PR*( Limiteur(1.0_PR/rip1)*(1.0_PR+PHI)*(Wip1-Wi) + Limiteur(rip1)*(1.0_PR-PHI)*(Wip2-Wip1) )
      !WR=Wip1-0.25_PR*( Limiteur(rip1)*(1.0_PR+PHI)*(Wip1-Wi) - Limiteur(1.0_PR/rip1)*(1.0_PR-PHI)*(Wip2-Wip1) )
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

end subroutine Reco_muscl

!!!=============================== LIMITEURS ===================================

real(PR) function MINMOD(r)

   implicit none
   real(PR),intent(in) :: r
   real(PR) :: PHI=1.0_PR/3.0_PR

   minmod=max(0.d0,min(1.d0,r))

   !minmod=min(1.d0,(3.0_PR-PHI)/(1.0_PR-PHI)*r)

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

subroutine src_euler1D !!!(Nl,N,SVL,SVR,ir,F,U,G,dUdt,dGdt)

   implicit none

   !integer, intent(in) :: Nl, N
   !real(PR), intent(in) :: SVL(1:N), SVR(1:N), ir(1:N) 
   !type(fluide), intent(in) :: F(1:N)
   !real(PR), intent(in) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
   !real(PR), intent(out) :: dUdt(1:Nl+4,1:N), dGdt(1:Nl,1:11,1:N)

   !!!---Left/right
   !real(PR) :: UL(1:Nl+4), GL(1:Nl,1:11)
   !real(PR) :: UR(1:Nl+4), GR(1:Nl,1:11)
   type(multifluide) :: MFL, MFR !!! variables primitives
   !!!---Grandeurs étoile
   real(PR) :: vx_s, vy_s, vz_s, vx_smax, F2max(1:11)
   !!!---vecteurs flux pour la partie conservative (eq 19 Ndanou 2015)
   real(PR) :: F1(1:Nl+4)
   !!!---vecteurs flux pour la partie non-conservative (eq 19 Ndanou 2015)
   real(PR) :: F2(1:Nl,1:11),H2(1:Nl,1:11),K2(1:Nl,1:11),M2(1:Nl,1:11)
   integer :: i, ip1, im1, ip2, iph, i2, N
   real(PR) :: t1,t2

   call cpu_time(t1)

   !call cpu_time(t10) 

   dUdt=0.0_PR
   dGdt=0.0_PR

   N=N1D

   if(Reconstruct)then

      !!!---Cl:-----------------  
      MFrec_R(0)=MFrec_L(1) 
      if(radial) MFrec_R(0)%vx=-MFrec_L(1)%vx
      MFrec_L(Nx+1)=MFrec_R(Nx)
      MFrec_L(0)=MFrec_R(0)
      MFrec_R(Nx+1)=MFrec_L(Nx+1)
      !!!---Cl:-----------------  

   else

      !!!---Cl:-----------------  
      MF1D(0)=MF1D(1)
      if(radial) MF1D(0)%vx=-MF1D(1)%vx
      MF1D(Nx+1)=MF1D(Nx)
      !!!-----------------------

   endif

   vx_smax=0.0_PR
   F2max=0.0_PR
    

   do i=0,N !!!-1

       ip1=i+1
      
       if(Reconstruct)then
          MFL=MFrec_R(i)
          MFR=MFrec_L(ip1)
       else 
          MFL=MF1D(i)
          MFR=MF1D(ip1)
       endif

      !!!---Solveur de Riemann :

      call HLLC_flux(MFL,MFR,vx_s,vy_s,vz_s,F1(1:Nl+4),F2(1:Nl,1:11))

      !vx_smax=max(vx_s,vx_smax)

      !!!-----------------  PARTIE CONSERVATIVE --------------------

      dUdt(1:Nl+4,i)  =dUdt(1:Nl+4,i)   - SVL(i)*F1(1:Nl+4)
      dUdt(1:Nl+4,ip1)=dUdt(1:Nl+4,ip1) + SVR(i)*F1(1:Nl+4)
      dUdt(Nl+1,i)    =dUdt(Nl+1,i)     + MF1D(i)%p*Sgeo(i)

      !!!---------------  PARTIE NON CONSERVATIVE ------------------

      !!!   terme source i :

      !call get_HKM(FL,H2,K2,M2)

      call get_HKM(MF1D(i),H2,K2,M2)

      dGdt(1:Nl,1:11,i)  =dGdt(1:Nl,1:11,i)-SVL(i)*(F2(1:Nl,1:11)&
                                                    +vx_s*H2(1:Nl,1:11)&
                                                    +vy_s*K2(1:Nl,1:11)&
                                                    +vz_s*M2(1:Nl,1:11))

      !!!   terme source i+1 :

      !call get_HKM(FR,H2,K2,M2)

      call get_HKM(MF1D(ip1),H2,K2,M2)

      dGdt(1:Nl,1:11,ip1)=dGdt(1:Nl,1:11,ip1)+SVR(i)*(F2(1:Nl,1:11)&
                                                      +vx_s*H2(1:Nl,1:11)&
                                                      +vy_s*K2(1:Nl,1:11)&
                                                      +vz_s*M2(1:Nl,1:11))

      !!!---sources externes
      !dUdt(Nl+1,i)=dUdt(Nl+1,i)+MF1D(i)%Qpx
      !dUdt(Nl+2,i)=dUdt(Nl+2,i)+MF1D(i)%Qpy
      !dUdt(Nl+3,i)=dUdt(Nl+3,i)+MF1D(i)%Qpz
      !dUdt(Nl+4,i)=dUdt(Nl+4,i)+MF1D(i)%Qe
      !do i2=1,Nl
      !   dGdt(i2,11,i)=dGdt(i2,11,i)+MF1D(i)%Qel(i2)*MF1D(i)%fl(i2)
      !enddo


  enddo

   !call cpu_time(t1) 

   call cpu_time(t2)
   t_euler=t_euler+t2-t1

   !i=900
   !print*, '   ---------------------------------------'
   !print*, '   dUdt=', dUdt(1:Nl+4,i)
   !print*, '   dGdt(1)=', dGdt(1,1:11,i)
   !print*, '   dGdt(2)=', dGdt(2,1:11,i)
   !print*, '   ---------------------------------------'

end subroutine src_euler1D


subroutine src_EM !!!(Nl,N,SVL,SVR,ir,F,U,G,dUdt,dGdt)

   implicit none

   integer :: i, ip1, im1, ip2, iph, i2, N

   call cpu_time(t10) 

   dUdt=0.0_PR
   dGdt=0.0_PR

   N=N1D

   !!!---Cl:-----------------
   MF1D(0)=MF1D(1)
   MF1D(0)%vx=-MF1D(1)%vx
   MF1D(Nx+1)=MF1D(Nx)
   !!!-----------------------

   do i=0,N !!!-1

      !!!---sources externes
      dUdt(Nl+1,i)=MF1D(i)%Qpx
      dUdt(Nl+2,i)=MF1D(i)%Qpy
      dUdt(Nl+3,i)=MF1D(i)%Qpz
      dUdt(Nl+4,i)=MF1D(i)%Qe
      do i2=1,Nl
         dGdt(i2,11,i)=MF1D(i)%F(i2)%Qe*MF1D(i)%F(i2)%f
      enddo

  enddo

   call cpu_time(t1) 

   t_EM=t_EM+t1-t10

end subroutine src_EM


subroutine get_HKM(MF,H2,K2,M2) !!! (F,G,H2,K2,M2)

   implicit none

   type(multifluide), intent(in) :: MF
   !real(PR), intent(in) :: G(1:Nl,1:11)
   real(PR), intent(out) :: H2(1:Nl,1:11)
   real(PR), intent(out) :: K2(1:Nl,1:11)
   real(PR), intent(out) :: M2(1:Nl,1:11)
   real(PR) :: fl(1:Nl), al(1:Nl,1:3), bl(1:Nl,1:3), cl(1:Nl,1:3), pl(1:Nl)
   real(PR) :: coeff
   integer :: iph, i2

   fl(1:Nl)=MF%F(1:Nl)%f    !G(1:Nl,1)
   !do iph=1,Nl
   !   !print*, iph, G(:,11), fl(:)
   !   !pl(iph)=(F%g_sge(iph)-1.0_PR)*G(iph,11)/fl(iph)-F%g_sge(iph)*F%p_sge(iph)
   !   !pl(iph)=(F%g_sge(iph)-1.0_PR)*F%rhl(iph)*F%ehl(iph)-F%g_sge(iph)*F%p_sge(iph)
   !enddo

   if(a2A)then
     do iph=1,Nl
      coeff=fl(iph)**(r13)
      al(iph,1:3)=MF%F(iph)%a(1:3)*coeff  
      bl(iph,1:3)=MF%F(iph)%b(1:3)*coeff  
      cl(iph,1:3)=MF%F(iph)%c(1:3)*coeff
     enddo
   else
     do iph=1,Nl
      al(iph,1:3)=MF%F(iph)%a(1:3)  
      bl(iph,1:3)=MF%F(iph)%b(1:3)  
      cl(iph,1:3)=MF%F(iph)%c(1:3)
     enddo
   endif

   pl(1:Nl)=MF%F(1:Nl)%p

   H2(1:Nl,1)=-1.0_PR*fl(1:Nl)
   H2(1:Nl,2:4)=0.0_PR
   do i2=1,3
    H2(1:Nl,4+i2)=-bl(1:Nl,i2)        
   enddo
   do i2=1,3
    H2(1:Nl,7+i2)=-cl(1:Nl,i2)        
   enddo
   if(modeG.eq.1)then
      H2(1:Nl,11)=fl(1:Nl)*pl(1:Nl)
   else
      H2(1:Nl,11)=-fl(1:Nl)*MF%F(1:Nl)%sig(1,1)
  endif

   K2(1:Nl,1)=0.0_PR
   do i2=1,3
    K2(1:Nl,1+i2)=bl(1:Nl,i2)
   enddo
   K2(1:Nl,5:11)=0.0_PR
   if(modeG.eq.2)then
      K2(1:Nl,11)=-fl(1:Nl)*MF%F(1:Nl)%sig(1,2)
   endif

   M2(1:Nl,1)=0.0_PR
   do i2=1,3
    M2(1:Nl,1+i2)=cl(1:Nl,i2)
   enddo
   M2(1:Nl,5:11)=0.0_PR
   if(modeG.eq.2)then
      M2(1:Nl,11)=-fl(1:Nl)*MF%F(1:Nl)%sig(1,3)
   endif


end subroutine get_HKM

!!!========================= CALCUL DU PAS DE TEMPS =============================

real(PR) function dtcfl(M,cfl)

   implicit none
   type(mesh), intent(in) :: M
   real(PR), intent(in) :: cfl
   real(PR) :: vmax, dttemp, vx, c, vmaxcfl, dxm
   integer :: i,j,k,icfl

   vmax=0.0_PR
   dtcfl=1.0e30_PR

   !print*, 'max c1', maxval(M%MF(:,:,:)%F(1)%c)
   !print*, 'max c2', maxval(M%MF(:,:,:)%F(2)%c)

   do k=1,M%Nz
    do j=1,M%Ny
     do i=1,M%Nx
         vmax=max(abs(M%MF(i,j,k)%vx),abs(M%MF(i,j,k)%vx+M%MF(i,j,k)%c),abs(M%MF(i,j,k)%vx-M%MF(i,j,k)%c))
         
         !dtcfl=min(dtcfl,M%dxm(i)/vmax)

         dttemp=M%dxm(i)/vmax

         if(dttemp.lt.dtcfl)then
            dtcfl=dttemp
            icfl=i
            vx=M%MF(i,j,k)%vx
            c=M%MF(i,j,k)%c
            vmaxcfl=vmax
            dxm=M%dxm(i)
         endif

         !dtcfl=min(dtcfl,vx,vx+c,vx-c)
         !dtcfl=min(dtcfl,vx,vx+c,vx-c)
     enddo
    enddo
   enddo

   dtcfl=cfl*dtcfl
   

   !print*, 'dtcfl:', dtcfl, icfl, vx, c, vmaxcfl, dxm, cfl 

   if(dtcfl.lt.1.0e-12_PR)then
       print*, 'PROBLEM CFL:', dtcfl
       stop
   endif

end function dtcfl

!!!===================== ROUTINES D'INTEGRATION TEMPORELLE =====================

subroutine test_p(message)
    implicit none
    character(len=*), intent(in)  :: message
    real(PR) :: pmax, p1max, p2max
    integer :: imax

    
    !pmax=maxval(MF1D(1:Nx)%p)
    !p1max=maxval(MF1D(1:Nx)%pl(1))
    !p2max=maxval(MF1D(1:Nx)%pl(2))
    !imax=maxloc(MF1D(1:Nx)%p,1)

    !if(pmax.gt.1.1e5_PR)then
    !write(6,"(A,A1,3(ES14.7,' '))"), trim(adjustl(message)),' ', pmax, p1max, p2max 
    !write(6,"(A,A6,I4)"), trim(adjustl(message)),' imax=', imax 

    !stop
    !endif

end subroutine test_p

subroutine COMPACT(dt) !!!(Nl,N,SVL,SVR,ir,F,dt) !!!U,G,dt)

    implicit none
    !!!-----------------------------
    !!! On suppose que les variables primitives varient doucement
    !!! Alors on utilise la "forme compacte" :  
    !!! dW/dt = W + A(W)dW/dx
    !!! avec:
    !!! W=(alphai, rhi, vx, vy, vz, pi, a(1:3), b(1:3), c(1:3) )
    !!! cf. Saurel 2009 Appendix C.
    !!! cf/ Favrie et al. 2009 p 6053 section 4.1.1
    !!! W(1:Nw) Nw=3+12*Nl
    !!!1  alphai: W(i1:i2-1)
    !!!2  rhi   : W(i2:i3-1)
    !!!3  vx,y,z: W(i3:i4-1)
    !!!4  pi    : W(i4:i5-1)
    !!!5  a(1)i : W(i5:i6-1)
    !!!6  a(2)i : W(i6:i7-1)
    !!!7  a(3)i : W(i7:i8-1)
    !!!8  b(1)i : W(i8:i9-1)
    !!!9  b(2)i : W(i9:i10-1)
    !!!10 b(3)i : W(i10:i11-1)
    !!!11 c(1)i : W(i11:i12-1)
    !!!12 c(2)i : W(i12:i13-1)
    !!!13 c(3)i : W(i13:Nw)
    !!!-----------------------------
    real(PR), intent(in) :: dt
    real(PR) :: c0, r
    integer :: i, j, k, iph, ip1
    logical :: test1, test2

    dUdt=0.0_PR ; dGdt=0.0_PR

    !!!----V2
    !!!   i1=1 ; i2=i1+Nl ; i3=i2+Nl ; i4=i3+3 ; i5=i4+Nl ; i6=i5+Nl
    !!!   i7=i6+Nl ; i8=i7+Nl ; i9=i8+Nl ; i10=i9+Nl, i11=i10+Nl, i12=i11+Nl
    !!!   i13=i12+Nl

    !!!   !!!---Wi
    !!!   Wi(i1:i2-1)   = MF1D(1)%F(1:Nl)%f
    !!!   Wi(i2:i3-1)   = MF1D(1)%F(1:Nl)%rh
    !!!   Wi(i3)        = MF1D(1)%vx
    !!!   Wi(i3+1)      = MF1D(1)%vy
    !!!   Wi(i3+2)      = MF1D(1)%vz
    !!!   Wi(i4:i5-1)   = MF1D(1)%F(1:Nl)%p
    !!!   Wi(i5:i6-1)   = MF1D(1)%F(1:Nl)%a(1)
    !!!   Wi(i6:i7-1)   = MF1D(1)%F(1:Nl)%a(2)
    !!!   Wi(i7:i8-1)   = MF1D(1)%F(1:Nl)%a(3)
    !!!   Wi(i8:i9-1)   = MF1D(1)%F(1:Nl)%b(1)
    !!!   Wi(i9:i10-1)  = MF1D(1)%F(1:Nl)%b(2)
    !!!   Wi(i10:i11-1) = MF1D(1)%F(1:Nl)%b(3)
    !!!   Wi(i11:i12-1) = MF1D(1)%F(1:Nl)%c(1)
    !!!   Wi(i12:i13-1) = MF1D(1)%F(1:Nl)%c(2)
    !!!   Wi(i13:Nw)    = MF1D(1)%F(1:Nl)%c(3)

    !!!   !!!---Wim1
    !!!   Wim1=Wi

    !!!   do ix=1,N1D
    !!!       !!!---Wip1
    !!!       ixp1=min(ix+1,N1D)
    !!!       !!!---Wip1
    !!!       Wip1(i1:i2-1)   = MF1D(1)%F(1:Nl)%f
    !!!       Wip1(i2:i3-1)   = MF1D(1)%F(1:Nl)%rh
    !!!       Wip1(i3)        = MF1D(1)%vx
    !!!       Wip1(i3+1)      = MF1D(1)%vy
    !!!       Wip1(i3+2)      = MF1D(1)%vz
    !!!       Wip1(i4:i5-1)   = MF1D(1)%F(1:Nl)%p
    !!!       Wip1(i5:i6-1)   = MF1D(1)%F(1:Nl)%a(1)
    !!!       Wip1(i6:i7-1)   = MF1D(1)%F(1:Nl)%a(2)
    !!!       Wip1(i7:i8-1)   = MF1D(1)%F(1:Nl)%a(3)
    !!!       Wip1(i8:i9-1)   = MF1D(1)%F(1:Nl)%b(1)
    !!!       Wip1(i9:i10-1)  = MF1D(1)%F(1:Nl)%b(2)
    !!!       Wip1(i10:i11-1) = MF1D(1)%F(1:Nl)%b(3)
    !!!       Wip1(i11:i12-1) = MF1D(1)%F(1:Nl)%c(1)
    !!!       Wip1(i12:i13-1) = MF1D(1)%F(1:Nl)%c(2)
    !!!       Wip1(i13:Nw)    = MF1D(1)%F(1:Nl)%c(3)
    !!!       
    !!!       !!!---Deltaim
    !!!       Dxm=M%dxm(ix-1)
    !!!       Deltaim=(Wi-Wim1)/Dxm

    !!!       !!!---Deltaip
    !!!       Dxp=M%dxm(ix)
    !!!       Deltaim=(Wip1-Wi)/Dxp

    !!!       !!!--Deltai
    !!!       !!! r=(Wi-Wim1)/(Wip1-Wi)
    !!!       test1=abs(Wip1-Wi).le.tolmuscl
    !!!       test2=abs(Wi-Wim1).le.tolmuscl
    !!!       if(test1.and.test2)then
    !!!          Deltai=0.0_PR
    !!!       else
    !!!          if(test1)then
    !!!             r=1.0e20_PR
    !!!          elseif(test2)then
    !!!             r=0.0_PR
    !!!          else
    !!!             r=(Wi-Wim1)/(Wip1-Wi)
    !!!          endif
    !!!          Deltai=LIMITEUR(r)
    !!!       endif
    !!!       !!!--- WL
    !!!       WL=Wi-0.5_PR*Dxm*Deltai
    !!!       WR=Wi+0.5_PR*Dxp*Deltai

    !!!       !!!--- A
    !!!       A=0.0_PR
    !!!       ! vx sur diagonale
    !!!       do i=1,Nw
    !!!           A(i,i)=Wi(i3)
    !!!       enddo
    !!!       ! dependance de rhi à vx
    !!!       A(i2:i3-1,i3)=W(i2:i3-1)
    !!!       ! dependance de vx à alpha et pi
    !!!       rhm1=1.0_PR/MF1D(ix)%rh
    !!!       A(i3,i1:i2-1)=W(i4:i5-1)*rhm1 !!! dep à alpha
    !!!       A(i3,i4:i5-1)=W(i1:i2-1)*rhm1 !!! dep à pi

    !!!       ! dependance de la pression à vx
    !!!       A(i4:i5-1,i3)=W(i2:i3-1)*MF1D(ix)%F(1:Nl)%c0**2

    !!!       ! dependance de v à a,b,c
    !!!           !!! précalcul
    !!!       call init_Finger(MF1D(i)%F(iph))
    !!!       !!!------------------ vx -------------------
    !!!       !!! dep de vx à a
    !!!       A(i3,i5:i6-1)=-rhm1*W(1:Nl)*dS11dai(1,1:Nl)
    !!!       A(i3,i6:i7-1)=-rhm1*W(1:Nl)*dS11dai(2,1:Nl)
    !!!       A(i3,i7:i8-1)=-rhm1*W(1:Nl)*dS11dai(3,1:Nl)
    !!!       !!! dep de vx à b
    !!!       A(i3,i8 :i9-1) =-rhm1*W(1:Nl)*dS11dbi(1,1:Nl)
    !!!       A(i3,i9 :i10-1)=-rhm1*W(1:Nl)*dS11dbi(2,1:Nl)
    !!!       A(i3,i10:i11-1)=-rhm1*W(1:Nl)*dS11dbi(3,1:Nl)
    !!!       !!! dep de vx à c
    !!!       A(i3,i11:i12-1)=-rhm1*W(1:Nl)*dS11dci(1,1:Nl)
    !!!       A(i3,i12:i13-1)=-rhm1*W(1:Nl)*dS11dci(2,1:Nl)
    !!!       A(i3,i13:Nw)   =-rhm1*W(1:Nl)*dS11dci(3,1:Nl)

    !!!       !!!------------------ vy -------------------
    !!!       !!! dep de vy à a
    !!!       A(i3+1,i5:i6-1)=-rhm1*W(1:Nl)*dS12dai(1,1:Nl)
    !!!       A(i3+1,i6:i7-1)=-rhm1*W(1:Nl)*dS12dai(2,1:Nl)
    !!!       A(i3+1,i7:i8-1)=-rhm1*W(1:Nl)*dS12dai(3,1:Nl)
    !!!       !!! dep de vy à b
    !!!       A(i3+1,i8 :i9-1) =-rhm1*W(1:Nl)*dS12dbi(1,1:Nl)
    !!!       A(i3+1,i9 :i10-1)=-rhm1*W(1:Nl)*dS12dbi(2,1:Nl)
    !!!       A(i3+1,i10:i11-1)=-rhm1*W(1:Nl)*dS12dbi(3,1:Nl)
    !!!       !!! dep de vy à c
    !!!       A(i3+1,i11:i12-1)=-rhm1*W(1:Nl)*dS12dci(1,1:Nl)
    !!!       A(i3+1,i12:i13-1)=-rhm1*W(1:Nl)*dS12dci(2,1:Nl)
    !!!       A(i3+1,i13:Nw)   =-rhm1*W(1:Nl)*dS12dci(3,1:Nl)

    !!!       !!!------------------ vz -------------------
    !!!       !!! dep de vz à a
    !!!       A(i3+2,i5:i6-1)=-rhm1*W(1:Nl)*dS13dai(1,1:Nl)
    !!!       A(i3+2,i6:i7-1)=-rhm1*W(1:Nl)*dS13dai(2,1:Nl)
    !!!       A(i3+2,i7:i8-1)=-rhm1*W(1:Nl)*dS13dai(3,1:Nl)
    !!!       !!! dep de vz à b
    !!!       A(i3+2,i8 :i9-1) =-rhm1*W(1:Nl)*dS13dbi(1,1:Nl)
    !!!       A(i3+2,i9 :i10-1)=-rhm1*W(1:Nl)*dS13dbi(2,1:Nl)
    !!!       A(i3+2,i10:i11-1)=-rhm1*W(1:Nl)*dS13dbi(3,1:Nl)
    !!!       !!! dep de vz à c
    !!!       A(i3+2,i11:i12-1)=-rhm1*W(1:Nl)*dS13dci(1,1:Nl)
    !!!       A(i3+2,i12:i13-1)=-rhm1*W(1:Nl)*dS13dci(2,1:Nl)
    !!!       A(i3+2,i13:Nw)   =-rhm1*W(1:Nl)*dS13dci(3,1:Nl)

    !!!       ! dependance de a à vx
    !!!       A(i5:i8-1,i3)=W(i5:i8-1)
    !!!       ! dependance de a à vy
    !!!       A(i5:i8-1,i3)=W(i8:i11-1)
    !!!       ! dependance de a à vz
    !!!       A(i5:i8-1,i3)=W(i11:Nw)

    !!!       !!!--- dWdt
    !!!       !!!dWdt=-A*(Wip1-Wi)/dx

    !!!       WiL=WiL+0.5*dt/dx*A*(WiL-WiR)
    !!!       WiR=WiR+0.5*dt/dx*A*(WiL-WiR)

    !!!       !!!--- cellule suivante...
    !!!       Wim1=Wi
    !!!       Wi=Wip1
    !!!   enddo

    !!!---V1

    !!!--- PRIMITIVE TO CONSERVATIVE ---
    call primitive2conservative_1D

    !!!--- RECONSTRUCTION MUSCL MANDATORY---
    !call MUSCL

    !!---pour V2
    Wi(1:Nl)            = MF1D(1)%F(1:Nl)%f
    Wi(Nl+1:2*Nl)       = MF1D(1)%F(1:Nl)%rh
    Wi(2*Nl+1)          = MF1D(1)%vx
    Wi(2*Nl+2)          = MF1D(1)%vy
    Wi(2*Nl+3)          = MF1D(1)%vz
    Wi(2*Nl+4:3*Nl+3)   = MF1D(1)%F(1:Nl)%p
    Wi(3*Nl+4:4*Nl+3)   = MF1D(1)%F(1:Nl)%a(1)
    Wi(4*Nl+4:5*Nl+3)   = MF1D(1)%F(1:Nl)%a(2)
    Wi(5*Nl+4:6*Nl+3)   = MF1D(1)%F(1:Nl)%a(3)
    Wi(6*Nl+4:7*Nl+3)   = MF1D(1)%F(1:Nl)%b(1)
    Wi(7*Nl+4:8*Nl+3)   = MF1D(1)%F(1:Nl)%b(2)
    Wi(8*Nl+4:9*Nl+3)   = MF1D(1)%F(1:Nl)%b(3)
    Wi(9*Nl+4:10*Nl+3)  = MF1D(1)%F(1:Nl)%c(1)
    Wi(10*Nl+4:11*Nl+3) = MF1D(1)%F(1:Nl)%c(2)
    Wi(11*Nl+4:12*Nl+3) = MF1D(1)%F(1:Nl)%c(3)

    Wim1=Wi

    do i=1,N1D

        !!!---V2:

        !!!---Wip1
        ip1=min(i+1,N1D)
        !!!---Wip1
        Wip1(1:Nl)            = MF1D(ip1)%F(1:Nl)%f
        Wip1(Nl+1:2*Nl)       = MF1D(ip1)%F(1:Nl)%rh
        Wip1(2*Nl+1)          = MF1D(ip1)%vx
        Wip1(2*Nl+2)          = MF1D(ip1)%vy
        Wip1(2*Nl+3)          = MF1D(ip1)%vz
        Wip1(2*Nl+4:3*Nl+3)   = MF1D(ip1)%F(1:Nl)%p
        Wip1(3*Nl+4:4*Nl+3)   = MF1D(ip1)%F(1:Nl)%a(1)
        Wip1(4*Nl+4:5*Nl+3)   = MF1D(ip1)%F(1:Nl)%a(2)
        Wip1(5*Nl+4:6*Nl+3)   = MF1D(ip1)%F(1:Nl)%a(3)
        Wip1(6*Nl+4:7*Nl+3)   = MF1D(ip1)%F(1:Nl)%b(1)
        Wip1(7*Nl+4:8*Nl+3)   = MF1D(ip1)%F(1:Nl)%b(2)
        Wip1(8*Nl+4:9*Nl+3)   = MF1D(ip1)%F(1:Nl)%b(3)
        Wip1(9*Nl+4:10*Nl+3)  = MF1D(ip1)%F(1:Nl)%c(1)
        Wip1(10*Nl+4:11*Nl+3) = MF1D(ip1)%F(1:Nl)%c(2)
        Wip1(11*Nl+4:12*Nl+3) = MF1D(ip1)%F(1:Nl)%c(3)

        !!!--Deltai
        Deltai=(Wip1-Wi)
        do j=1,12*Nl+3
            if(abs(Deltai(j)).gt.tolmuscl)then
                r=(Wi(j)-Wim1(j))/(Wip1(j)-Wi(j))
                Deltai(j)=LIMITEUR(r)*Deltai(j)
            endif
        enddo

        !!!--- WL/R
        WL=Wi-0.5_PR*Deltai
        WR=Wi+0.5_PR*Deltai

        !!!---V1:

        ! !!!--- WL
        ! WL(1:Nl)            = MFrec_L(i)%F(1:Nl)%f
        ! WL(Nl+1:2*Nl)       = MFrec_L(i)%F(1:Nl)%rh
        ! WL(2*Nl+1)          = MFrec_L(i)%vx
        ! WL(2*Nl+2)          = MFrec_L(i)%vy
        ! WL(2*Nl+3)          = MFrec_L(i)%vz
        ! WL(2*Nl+4:3*Nl+3)   = MFrec_L(i)%F(1:Nl)%p
        ! WL(3*Nl+4:4*Nl+3)   = MFrec_L(i)%F(1:Nl)%a(1)
        ! WL(4*Nl+4:5*Nl+3)   = MFrec_L(i)%F(1:Nl)%a(2)
        ! WL(5*Nl+4:6*Nl+3)   = MFrec_L(i)%F(1:Nl)%a(3)
        ! WL(6*Nl+4:7*Nl+3)   = MFrec_L(i)%F(1:Nl)%b(1)
        ! WL(7*Nl+4:8*Nl+3)   = MFrec_L(i)%F(1:Nl)%b(2)
        ! WL(8*Nl+4:9*Nl+3)   = MFrec_L(i)%F(1:Nl)%b(3)
        ! WL(9*Nl+4:10*Nl+3)  = MFrec_L(i)%F(1:Nl)%c(1)
        ! WL(10*Nl+4:11*Nl+3) = MFrec_L(i)%F(1:Nl)%c(2)
        ! WL(11*Nl+4:12*Nl+3) = MFrec_L(i)%F(1:Nl)%c(3)

        ! !!!--- WR
        ! WR(1:Nl)            = MFrec_R(i)%F(1:Nl)%f
        ! WR(Nl+1:2*Nl)       = MFrec_R(i)%F(1:Nl)%rh
        ! WR(2*Nl+1)          = MFrec_R(i)%vx
        ! WR(2*Nl+2)          = MFrec_R(i)%vy
        ! WR(2*Nl+3)          = MFrec_R(i)%vz
        ! WR(2*Nl+4:3*Nl+3)   = MFrec_R(i)%F(1:Nl)%p
        ! WR(3*Nl+4:4*Nl+3)   = MFrec_R(i)%F(1:Nl)%a(1)
        ! WR(4*Nl+4:5*Nl+3)   = MFrec_R(i)%F(1:Nl)%a(2)
        ! WR(5*Nl+4:6*Nl+3)   = MFrec_R(i)%F(1:Nl)%a(3)
        ! WR(6*Nl+4:7*Nl+3)   = MFrec_R(i)%F(1:Nl)%b(1)
        ! WR(7*Nl+4:8*Nl+3)   = MFrec_R(i)%F(1:Nl)%b(2)
        ! WR(8*Nl+4:9*Nl+3)   = MFrec_R(i)%F(1:Nl)%b(3)
        ! WR(9*Nl+4:10*Nl+3)  = MFrec_R(i)%F(1:Nl)%c(1)
        ! WR(10*Nl+4:11*Nl+3) = MFrec_R(i)%F(1:Nl)%c(2)
        ! WR(11*Nl+4:12*Nl+3) = MFrec_R(i)%F(1:Nl)%c(3)

        !!!---alpha       
        do iph=1,Nl
           AdW(iph)=MF1D(i)%vx*(MFrec_L(i)%F(iph)%f-MFrec_R(i)%F(iph)%f)
        enddo

        !!!---rho     
        do iph=1,Nl
           AdW(iph+Nl)=MF1D(i)%vx*(MFrec_L(i)%F(iph)%rh-MFrec_R(i)%F(iph)%rh)+&
               MF1D(i)%F(iph)%rh*(MFrec_L(i)%vx-MFrec_R(i)%vx)    
           !MFrec_L(i)%F(iph)%rh=MFrec_L(i)%F(iph)%rh+0.5_PR*dt*SVL(i)*AdW
           !MFrec_R(i)%F(iph)%rh=MFrec_R(i)%F(iph)%rh+0.5_PR*dt*SVR(i)*AdW
        enddo

        !!!---vx  
        !!! dep de vx à la pression
        k=2*Nl+1 ; AdW(k)=0.0_PR 
        do iph=1,Nl
            AdW(k)=AdW(k)+MF1D(i)%F(iph)%p*(MFrec_L(i)%F(iph)%f-MFrec_R(i)%F(iph)%f)
            AdW(k)=AdW(k)+MF1D(i)%F(iph)%f*(MFrec_L(i)%F(iph)%p-MFrec_R(i)%F(iph)%p)   
        enddo

        !!! dependance de vx aux vecteurs de base
        !do iph=1,Nl
        !    !!! précalcul
        !    call init_Finger(MF1D(i)%F(iph))
        !    !!! dep de vx à al
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*dS11dai(k)*(MFrec_L(i)%F(iph)%a(k)-MFrec_R(i)%F(iph)%a(k))
        !    enddo
        !    !!! dep de vx à bl
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*dS11dbi(k)*(MFrec_L(i)%F(iph)%b(k)-MFrec_R(i)%F(iph)%b(k))
        !    enddo
        !    !!! dep de vx à cl
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*dS11dci(k)*(MFrec_L(i)%F(iph)%c(k)-MFrec_R(i)%F(iph)%c(k))
        !    enddo
        !    AdW(k)=AdW(k)-(MF1D(i)%F(iph)%sig(1,1)+MF1D(i)%p)*(MFrec_L(i)%F(iph)%f-MFrec_R(i)%F(iph)%f)
        !enddo

        AdW(k)=AdW(k)/MF1D(i)%rh
        AdW(k)=AdW(k)+MF1D(i)%vx*(MFrec_L(i)%vx-MFrec_R(i)%vx)

        !MFrec_L(i)%vx=MFrec_L(i)%vx+0.5_PR*dt*SVL(i)*AdW
        !MFrec_R(i)%vx=MFrec_R(i)%vx+0.5_PR*dt*SVR(i)*AdW

        !!!!---vy 
        k=2*Nl+2 ; AdW(k)=0.0_PR
        !!! dependance de vy aux vecteurs de base
        !do iph=1,Nl
        !    !!! précalcul
        !    call init_Finger(MF1D(i)%F(iph))
        !    !!! dep de vy à al
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*ds12dai(k)*(MFrec_L(i)%F(iph)%a(k)-MFrec_R(i)%F(iph)%a(k))
        !    enddo
        !    !!! dep de vy à bl
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*ds12dbi(k)*(MFrec_L(i)%F(iph)%b(k)-MFrec_R(i)%F(iph)%b(k))
        !    enddo
        !    !!! dep de vy à cl
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*ds12dci(k)*(MFrec_L(i)%F(iph)%c(k)-MFrec_R(i)%F(iph)%c(k))
        !    enddo
        !    AdW(k)=AdW(k)-(MF1D(i)%F(iph)%sig(1,2))*(MFrec_L(i)%F(iph)%f-MFrec_R(i)%F(iph)%f)
        !enddo

        AdW(k)=AdW(k)/MF1D(i)%rh
        AdW(k)=AdW(k)+MF1D(i)%vx*(MFrec_L(i)%vy-MFrec_R(i)%vy)

        !MFrec_L(i)%vy=MFrec_L(i)%vy+0.5_PR*dt*SVL(i)*AdW
        !MFrec_R(i)%vy=MFrec_R(i)%vy+0.5_PR*dt*SVR(i)*AdW

        !!!!---vz
        k=2*Nl+3 ; AdW(k)=0.0_PR
        !!! dependance de vz aux vecteurs de base
        !do iph=1,Nl
        !    !!! précalcul
        !    call init_Finger(MF1D(i)%F(iph))
        !    !!! dep de vz à al
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*ds13dai(k)*(MFrec_L(i)%F(iph)%a(k)-MFrec_R(i)%F(iph)%a(k))
        !    enddo
        !    !!! dep de vz à bl
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*ds13dbi(k)*(MFrec_L(i)%F(iph)%b(k)-MFrec_R(i)%F(iph)%b(k))
        !    enddo
        !    !!! dep de vz à cl
        !    do k=1,3
        !     AdW(k)=AdW(k)-MF1D(i)%F(iph)%f*ds13dci(k)*(MFrec_L(i)%F(iph)%c(k)-MFrec_R(i)%F(iph)%c(k))
        !    enddo
        !    AdW(k)=AdW(k)-(MF1D(i)%F(iph)%sig(1,3))*(MFrec_L(i)%F(iph)%f-MFrec_R(i)%F(iph)%f)
        !enddo

        AdW(k)=AdW(k)/MF1D(i)%rh
        AdW(k)=AdW(k)+MF1D(i)%vx*(MFrec_L(i)%vz-MFrec_R(i)%vz)

        !MFrec_L(i)%vz=MFrec_L(i)%vz+0.5_PR*dt*SVL(i)*AdW
        !MFrec_R(i)%vz=MFrec_R(i)%vz+0.5_PR*dt*SVR(i)*AdW

        !!!---p
        do iph=1,Nl
            if(MF1D(i)%F(iph)%rh.lt.0.0_PR)then
               print*, 'PROBLEM COMPACT: MF1D%rh<0'
               stop 
            endif
            c0=ph(iph)%soundspeed(MF1D(i)%F(iph),MF1D(i)%F(iph)%rh,MF1D(i)%F(iph)%p)
            AdW(k+iph)=MF1D(i)%vx*(MFrec_L(i)%F(iph)%p-MFrec_R(i)%F(iph)%p)+&
                MF1D(i)%F(iph)%rh*c0**2*(MFrec_L(i)%vx-MFrec_R(i)%vx)    
            !!!MFrec_L(i)%F(iph)%p=MFrec_L(i)%F(iph)%p+0.5_PR*dt*SVL(i)*AdW
            !!!MFrec_R(i)%F(iph)%p=MFrec_R(i)%F(iph)%p+0.5_PR*dt*SVR(i)*AdW
        enddo

        !!!---al
        j=3*Nl+3
        do iph=1,Nl
            do k=1,3
                j=j+1  
                AdW(j)=MF1D(i)%vx*(MFrec_L(i)%F(iph)%a(k)-MFrec_R(i)%F(iph)%a(k))&
                      +MF1D(i)%F(iph)%a(k)*(MFrec_L(i)%vx-MFrec_R(i)%vx)&
                      +MF1D(i)%F(iph)%b(k)*(MFrec_L(i)%vy-MFrec_R(i)%vy)&
                      +MF1D(i)%F(iph)%c(k)*(MFrec_L(i)%vz-MFrec_R(i)%vz) 
                !MFrec_L(i)%F(iph)%a(k)=MFrec_L(i)%F(iph)%a(k)+0.5_PR*dt*SVL(i)*AdW
                !MFrec_R(i)%F(iph)%a(k)=MFrec_R(i)%F(iph)%a(k)+0.5_PR*dt*SVR(i)*AdW
            enddo
        enddo
        !!!---bl
        do iph=1,Nl
            do k=1,3
                j=j+1
                AdW(j)=MF1D(i)%vx*(MFrec_L(i)%F(iph)%b(k)-MFrec_R(i)%F(iph)%b(k))
                !MFrec_L(i)%F(iph)%b(k)=MFrec_L(i)%F(iph)%b(k)+0.5_PR*dt*SVL(i)*AdW
                !MFrec_R(i)%F(iph)%b(k)=MFrec_R(i)%F(iph)%b(k)+0.5_PR*dt*SVR(i)*AdW
            enddo
        enddo
        !!!---cl
        do iph=1,Nl
            do k=1,3
                j=j+1
                AdW(j)=MF1D(i)%vx*(MFrec_L(i)%F(iph)%c(k)-MFrec_R(i)%F(iph)%c(k))
                !MFrec_L(i)%F(iph)%c(k)=MFrec_L(i)%F(iph)%c(k)+0.5_PR*dt*SVL(i)*AdW
                !MFrec_R(i)%F(iph)%c(k)=MFrec_R(i)%F(iph)%c(k)+0.5_PR*dt*SVR(i)*AdW
            enddo
        enddo

        if(j.ne.(12*Nl+3))then
            write(iout,*) '   PROBLEM COMPACT W : ', j, ' .ne. ', 12*Nl+3
            stop
        endif

        !!!--- MAJ W

        WL=WL+0.5_PR*dt*SVL(i)*AdW
        WR=WR+0.5_PR*dt*SVR(i)*AdW

        !!!---MAJ R --------------------------- 

        MFrec_R(i)%F(1:Nl)%f    =  WR(1:Nl)            
        MFrec_R(i)%F(1:Nl)%rh   =  WR(Nl+1:2*Nl)       
        MFrec_R(i)%vx           =  WR(2*Nl+1)          
        MFrec_R(i)%vy           =  WR(2*Nl+2)          
        MFrec_R(i)%vz           =  WR(2*Nl+3)          
        MFrec_R(i)%F(1:Nl)%p    =  WR(2*Nl+4:3*Nl+3)   
        MFrec_R(i)%F(1:Nl)%a(1) =  WR(3*Nl+4:4*Nl+3)   
        MFrec_R(i)%F(1:Nl)%a(2) =  WR(4*Nl+4:5*Nl+3)   
        MFrec_R(i)%F(1:Nl)%a(3) =  WR(5*Nl+4:6*Nl+3)   
        MFrec_R(i)%F(1:Nl)%b(1) =  WR(6*Nl+4:7*Nl+3)   
        MFrec_R(i)%F(1:Nl)%b(2) =  WR(7*Nl+4:8*Nl+3)   
        MFrec_R(i)%F(1:Nl)%b(3) =  WR(8*Nl+4:9*Nl+3)   
        MFrec_R(i)%F(1:Nl)%c(1) =  WR(9*Nl+4:10*Nl+3)  
        MFrec_R(i)%F(1:Nl)%c(2) =  WR(10*Nl+4:11*Nl+3) 
        MFrec_R(i)%F(1:Nl)%c(3) =  WR(11*Nl+4:12*Nl+3) 

        MFrec_R(i)%rh=sum(MFrec_R(i)%F(1:Nl)%f*MFrec_R(i)%F(1:Nl)%rh)
        do iph=1,Nl
            MFrec_R(i)%F(iph)%Y=MFrec_R(i)%F(iph)%f*MFrec_R(i)%F(iph)%rh/MFrec_R(i)%rh
        enddo
        do iph=1,Nl
            MFrec_R(i)%F(iph)%eh=ph(iph)%energie_hydro(MFrec_R(i)%F(iph),MFrec_R(i)%F(iph)%rh,MFrec_R(i)%F(iph)%rh)
        enddo

        MFrec_R(i)%u=sum(MFrec_R(i)%F(1:Nl)%Y*MFrec_R(i)%F(1:Nl)%eh)
 
        !!!---relaxation de p ?
        !do iph=1,Nl
        !   call relax_p(Frec_R(i))
        !enddo

        call MAJ_meca(MFrec_R(i))

        MFrec_R(i)%e=MFrec_R(i)%u+MFrec_R(i)%ee+0.5_PR*MFrec_R(i)%vx**2

        MFrec_R(i)%c=soundspeed_mixt(MFrec_R(i))

        !!!---MAJ L --------------------------- 

        MFrec_L(i)%F(1:Nl)%f    =  WL(1:Nl)           
        MFrec_L(i)%F(1:Nl)%rh   =  WL(Nl+1:2*Nl)      
        MFrec_L(i)%vx           =  WL(2*Nl+1)         
        MFrec_L(i)%vy           =  WL(2*Nl+2)         
        MFrec_L(i)%vz           =  WL(2*Nl+3)         
        MFrec_L(i)%F(1:Nl)%p    =  WL(2*Nl+4:3*Nl+3)  
        MFrec_L(i)%F(1:Nl)%a(1) =  WL(3*Nl+4:4*Nl+3)  
        MFrec_L(i)%F(1:Nl)%a(2) =  WL(4*Nl+4:5*Nl+3)  
        MFrec_L(i)%F(1:Nl)%a(3) =  WL(5*Nl+4:6*Nl+3)  
        MFrec_L(i)%F(1:Nl)%b(1) =  WL(6*Nl+4:7*Nl+3)  
        MFrec_L(i)%F(1:Nl)%b(2) =  WL(7*Nl+4:8*Nl+3)  
        MFrec_L(i)%F(1:Nl)%b(3) =  WL(8*Nl+4:9*Nl+3)  
        MFrec_L(i)%F(1:Nl)%c(1) =  WL(9*Nl+4:10*Nl+3) 
        MFrec_L(i)%F(1:Nl)%c(2) =  WL(10*Nl+4:11*Nl+3)
        MFrec_L(i)%F(1:Nl)%c(3) =  WL(11*Nl+4:12*Nl+3)

        MFrec_L(i)%rh=sum(MFrec_L(i)%F(1:Nl)%f*MFrec_L(i)%F(1:Nl)%rh)
        do iph=1,Nl
           MFrec_L(i)%F(iph)%Y=MFrec_L(i)%F(iph)%f*MFrec_L(i)%F(iph)%rh/MFrec_L(i)%rh
        enddo
        do iph=1,Nl
           MFrec_L(i)%F(iph)%eh=ph(iph)%energie_hydro(MFrec_L(i)%F(iph),MFrec_L(i)%F(iph)%rh,MFrec_L(i)%F(iph)%p)
        enddo

        MFrec_L(i)%u=sum(MFrec_L(i)%F(1:Nl)%Y*MFrec_L(i)%F(1:Nl)%eh)
 
        call MAJ_meca(MFrec_L(i))

        MFrec_L(i)%e=MFrec_L(i)%u+MFrec_L(i)%ee+0.5_PR*MFrec_L(i)%vx**2

        MFrec_L(i)%c=soundspeed_mixt(MFrec_L(i))

        !!!--- cellule suivante...
        Wim1=Wi
        Wi=Wip1

    enddo

    call SRC 

    !!!--- MAJ VARIABLES CONSERVATIVES ---
    U(1:Nl+4,1:N1D) = U(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt

    !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
    G(1:Nl,1:11,1:N1D) = G(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt

    !!!--- CONSERVATIVE TO PRIMITIVE
    call conservative2primitive_1D

    write(iout,*) 'coucouzzz ! 8=>'

    !!!--- PRESSURE RELAXATION ---
    call relaxation
    write(iout,*) 'relax finii ooo'

    !!!--- RESET INTERNAL ENERGY ---
    call reset_internal_energy_1D

end subroutine COMPACT


subroutine RK1(dt) !!!(Nl,N,SVL,SVR,ir,F,dt) !!!U,G,dt)

    implicit none
    !integer, intent(in) :: Nl, N
    !real(PR), intent(in) :: SVL(1:N),SVR(1:N),ir(1:N)
    !type(fluide), intent(inout) :: F(1:N)
    real(PR), intent(in) :: dt
    !real(PR) :: U(1:Nl+4,1:N), G(1:Nl,1:11,1:N)
    !real(PR) :: dUdt(1:Nl+4,1:N), dGdt(1:Nl,1:11,1:N)
    integer :: i, j, k
    logical :: check_det=.false.
    real(PR) :: t1,t2,t3,t4,t5,t6

    dUdt=0.0_PR ; dGdt=0.0_PR

    !!!--- PRIMITIVE TO CONSERVATIVE ---
    call primitive2conservative_1D
    
    !!!--- RECONSTRUCTION ? ---
    if(Reconstruct) call MUSCL

    !!!--- CALCUL DES SOURCES (src_euler1D)
    call SRC 

    call cpu_time(t1)  

    !!!--- MAJ VARIABLES CONSERVATIVES ---
    U(1:Nl+4,1:N1D) = U(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt

    !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
    G(1:Nl,1:11,1:N1D) = G(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt

    call cpu_time(t2) ; t_src=t_src+t2-t1

    !!!--- CONSERVATIVE TO PRIMITIVE
    call conservative2primitive_1D

    !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
    !!! call nonconservative2primitive_1D
    
    !!!--- PRESSURE RELAXATION ---
    call relaxation

    !!!--- RESET INTERNAL ENERGY ---
    call reset_internal_energy_1D

end subroutine RK1

subroutine check_rhodet(message)

   !!!  check_det
   implicit none

   character(len=*), intent(in) :: message
   real(PR) :: arho1, arho2
   integer :: i

   do i=1,N1D

   arho1=MF1D(i)%F(1)%f*MF1D(i)%F(1)%rh
   !rho2=ph(1)%rho0*MF1D(i)%al(1,1)

   arho2=MF1D(i)%F(1)%f0*MF1D(i)%F(1)%rh0*sqrt(MF1D(i)%F(1)%detG)

   if(2.0_PR*abs(arho1-arho2)/(arho1+arho2).gt.1.0e-3_PR)then
      print*, trim(adjustl(message))
      print*, 'alpha0*rho0*detG**1/2 .ne. alpha*rho:', i
      print*, arho1, 'ne',  arho2
      print*, 'rho1=',MF1D(i)%F(1)%rh,'rho0=',MF1D(i)%F(1)%rh0
      print*, 'alpha1=',MF1D(i)%F(1)%f, 'alpha0=', MF1D(i)%F(1)%f0
      print*, 'sqrt(detG)):', sqrt(MF1D(i)%F(1)%detG)
      stop 
   endif

   enddo

end subroutine check_rhodet

subroutine RK2(dt) !!!(Nl,N,SVL,SVR,ir,F,dt)

    implicit none
    
    real(PR), intent(in) :: dt
    integer :: i, j, k

    dUdt=0.0_PR ; dGdt=0.0_PR

    !!!---STEP 1 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!---Save current state
        U0=U ; G0=G

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = U0(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*0.5_PR*dt

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = G0(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*0.5_PR*dt

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

    !!!---STEP 2 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC
        
        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = U0(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = G0(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation 

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

end subroutine RK2

subroutine RK2_SSP(dt) !!!(Nl,N,SVL,SVR,ir,F,dt)

    implicit none
    
    real(PR), intent(in) :: dt
    integer :: i, j, k

    dUdt=0.0_PR ; dGdt=0.0_PR

    !!!---STEP 1 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!---Save current state
        U0=U ; G0=G

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = U0(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = G0(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

    !!!---STEP 2 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC
        
        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = 0.5_PR*( U(1:Nl+4,1:N1D) + U0(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt )

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = 0.5_PR*( G(1:Nl,1:11,1:N1D) + G0(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt )

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation 

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

end subroutine RK2_SSP

subroutine RK3(dt)   !!!(Nl,N,SVL,SVR,ir,F) !!!,U,G,dt)

    implicit none
    real(PR), intent(in) :: dt
    integer :: i, j, k

    !!!---STEP 1 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!---Save current state
        U0=U ; G0=G

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = U0(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = G0(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

    !!!---STEP 2 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!---Save current state
        U1=U ; G1=G

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = 0.75_PR*U0(1:Nl+4,1:N1D) + 0.25_PR*( U1(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt )

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = 0.75_PR*G0(1:Nl,1:11,1:N1D) + 0.25_PR*( G1(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt )

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

    !!!---STEP 3 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!---Save current state
        U1=U ; G1=G

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = r13*( U0(1:Nl+4,1:N1D) + 2.0_PR*( U1(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt ) )

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = r13*( G0(1:Nl,1:11,1:N1D) + 2.0_PR*( G1(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt ) )

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

end subroutine RK3

subroutine RK4(dt)   !!!(Nl,N,SVL,SVR,ir,F) !!!,U,G,dt)

    implicit none
    real(PR), intent(in) :: dt
    integer :: i, j, k

    !!!---STEP 1 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!---Save current state
        U0=U ; G0=G

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = U0(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*0.5_PR*dt

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = G0(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*0.5_PR*dt

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

    !!!---STEP 2 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!---Save current state
        U1=U ; G1=G

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = U0(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*0.5_PR*dt

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = G0(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*0.5_PR*dt

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

    !!!---STEP 3 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!---Save current state
        U2=U ; G2=G

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = U0(1:Nl+4,1:N1D) + dUdt(1:Nl+4,1:N1D)*dt

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = G0(1:Nl,1:11,1:N1D) + dGdt(1:Nl,1:11,1:N1D)*dt

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

    !!!---STEP 4 :

        !!!--- PRIMITIVE TO CONSERVATIVE ---
        call primitive2conservative_1D

        !!!--- RECONSTRUCTION ? ---
        if(Reconstruct) call MUSCL 

        !!!--- CALCUL DES SOURCES (src_euler1D)
        call SRC 

        !!!--- MAJ VARIABLES CONSERVATIVES ---
        U(1:Nl+4,1:N1D) = r13*( -U0(1:Nl+4,1:N1D) + U1(1:Nl+4,1:N1D) + 2.0_PR*U2(1:Nl+4,1:N1D)&
                                + U(1:Nl+4,1:N1D) + 0.5_PR*dUdt(1:Nl+4,1:N1D)*dt )

        !!!--- MAJ VARIABLES NON-CONSERVATIVES ---
        G(1:Nl,1:11,1:N1D) = r13*( -G0(1:Nl,1:11,1:N1D) + G1(1:Nl,1:11,1:N1D) + 2.0_PR*G2(1:Nl,1:11,1:N1D)&
                                   + G(1:Nl,1:11,1:N1D) + 0.5_PR*dGdt(1:Nl,1:11,1:N1D)*dt )

        !!!--- CONSERVATIVE TO PRIMITIVE
        call conservative2primitive_1D

        !!!--- NON-CONSERVATIVE TO PRIMITIVE ---
        !!! call nonconservative2primitive_1D
 
        !!!--- PRESSURE RELAXATION ---
        call relaxation

        !!!--- RESET INTERNAL ENERGY ---
        call reset_internal_energy_1D

end subroutine RK4


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
   !type(fluide), intent(in) :: flu
   real(PR) :: tol, fx(1:1), dfdx(1:1), dx
   integer :: it    

   call F(1,x,fx)
   tol=abs(fx(1))
   it=0

   do while(tol.gt.1.0e-6_PR.and.it.lt.10000)
 
      call F(1,x(1),fx(1:1))
      call DF(1,x(1),dfdx(1:1))

      dx=-fx(1)/dfdx(1)

      x(1)=x(1)+dx
      x(1)=max(x(1),1.0e-10_PR)
      !x(1)=max(x(1),xmin)
      !x(1)=min(x(1),xmax)

      it=it+1

      tol=abs(fx(1))

      if(verbnewton) print*, 'it', it, 'tol', tol, x(1)

      !if(it.ge.9900)then
      ! print*, 'it=', it, x(1), fx(1), dx, dfdx(1)
      !endif

   enddo

   Nit_Newton=it

   if(it.ge.10000)then
     print*, '********************************************'
     print*, 'PROBLEM NEWTON it>10000'
     print*, 'fl0:', fl0(:)
     !print*, 'flrhl0(:)', flrhl0(:)
     !print*, 'rhl0(:)', flrhl0(:)/fl0(:)
     !print*, 'el0(:)', el0(:)
     print*, 'pl0:', pl0(:)
     print*, 'xmin=',xmin
     print*, 'xmax=',xmax
     print*, '********************************************'
 
     err_relax=3

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

subroutine init_multifluide(MF)

    implicit none
    type(multifluide), intent(inout) :: MF
    integer :: i
    
    !!!---grandeurs de mélange 
    MF%rh=0.0_PR
    MF%e=0.0_PR
    MF%u=0.0_PR
    MF%ee=0.0_PR
    !!!---vitesses hydrodynamiques:
    MF%vx=0.0_PR
    MF%vy=0.0_PR
    MF%vz=0.0_PR
    
    !!!---grandeurs phasiques
    
    !!!---propriétés physiques des matériaux (SG EOS ):
    MF%Nl=Nl
    allocate(MF%F(1:Nl))
    
    !!!---matrice identité
    Id(:,:)=0.0_PR
    do i=1,3
       Id(i,i)=1.0_PR
    enddo 

end subroutine init_multifluide

subroutine init_phase(M)

    implicit none
    type(mesh), intent(inout) :: M
    integer :: i,j,k,iph,imat
    logical :: found
 
    !!!---présence de matériaux de type 2 ?
    found=.false.
    do iph=1,M%Nl
        if(ph(iph)%typ.eq.2)then
            found=.true.
        endif
    enddo
    if(found) allocate(mat(1:M%Nl))
 
    do iph=1,M%Nl
 
        if(ph(iph)%typ.eq.1)then
 
            if(trim(adjustl(ph(iph)%nom)).eq.'air')then
 
                ph(iph)%p_sge=p_sge_air
                ph(iph)%g_sge=g_sge_air
                ph(iph)%mu=mu_air
                ph(iph)%sigy=sigy_air
                ph(iph)%rho0=rho0_air 
                ph(iph)%cv=cv_air
                ph(iph)%sig=sig_air
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            elseif(trim(adjustl(ph(iph)%nom)).eq.'eau')then
 
                ph(iph)%p_sge=p_sge_eau
                ph(iph)%g_sge=g_sge_eau
                ph(iph)%mu=mu_eau
                ph(iph)%sigy=sigy_eau
                ph(iph)%rho0=rho0_eau 
                ph(iph)%cv=cv_eau
                ph(iph)%sig=sig_eau
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            elseif(trim(adjustl(ph(iph)%nom)).eq.'Fe')then
 
                ph(iph)%p_sge=p_sge_Fe
                ph(iph)%g_sge=g_sge_Fe
                ph(iph)%mu=mu_Fe
                ph(iph)%sigy=sigy_Fe
                ph(iph)%rho0=rho0_Fe 
                ph(iph)%cv=cv_Fe
                ph(iph)%sig=sig_Fe
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            elseif(trim(adjustl(ph(iph)%nom)).eq.'Al')then
 
                ph(iph)%p_sge=p_sge_Al
                ph(iph)%g_sge=g_sge_Al
                ph(iph)%mu=mu_Al
                ph(iph)%sigy=sigy_Al
                ph(iph)%rho0=rho0_Al 
                ph(iph)%cv=cv_Al
                ph(iph)%sig=sig_Al
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            elseif(trim(adjustl(ph(iph)%nom)).eq.'Ti')then
 
                ph(iph)%p_sge=p_sge_Ti
                ph(iph)%g_sge=g_sge_Ti
                ph(iph)%mu=mu_Ti
                ph(iph)%sigy=sigy_Ti
                ph(iph)%rho0=rho0_Ti 
                ph(iph)%cv=cv_Ti
                ph(iph)%sig=sig_Ti
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            elseif(trim(adjustl(ph(iph)%nom)).eq.'Cu')then
 
                ph(iph)%p_sge=p_sge_Cu
                ph(iph)%g_sge=g_sge_Cu
                ph(iph)%mu=mu_Cu
                ph(iph)%sigy=sigy_Cu
                ph(iph)%rho0=rho0_Cu 
                ph(iph)%cv=cv_Cu
                ph(iph)%sig=sig_Cu
                ph(iph)%Tvap=Tvap_Cu
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            elseif(trim(adjustl(ph(iph)%nom)).eq.'Cug')then
 
                ph(iph)%p_sge=p_sge_Cug
                ph(iph)%g_sge=g_sge_Cug
                ph(iph)%mu=mu_Cug
                ph(iph)%sigy=sigy_Cug
                ph(iph)%rho0=rho0_Cug 
                ph(iph)%cv=cv_Cug 
                ph(iph)%sig=sig_Cug
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            elseif(trim(adjustl(ph(iph)%nom)).eq.'epoxy')then
 
                ph(iph)%p_sge=p_sge_epoxy
                ph(iph)%g_sge=g_sge_epoxy
                ph(iph)%mu=mu_epoxy
                ph(iph)%sigy=sigy_epoxy
                ph(iph)%rho0=rho0_epoxy 
                ph(iph)%cv=cv_epoxy 
                ph(iph)%sig=sig_epoxy
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            elseif(trim(adjustl(ph(iph)%nom)).eq.'spinel')then
 
                ph(iph)%p_sge=p_sge_spinel
                ph(iph)%g_sge=g_sge_spinel
                ph(iph)%mu=mu_spinel
                ph(iph)%sigy=sigy_spinel
                ph(iph)%rho0=rho0_spinel 
                ph(iph)%cv=cv_spinel 
                ph(iph)%sig=sig_spinel
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            else !!! défaut
 
                ph(iph)%p_sge=p_sge_air
                ph(iph)%g_sge=g_sge_air
                ph(iph)%mu=mu_air
                ph(iph)%sigy=sigy_air
                ph(iph)%rho0=rho0_air 
                ph(iph)%cv=cv_air
                ph(iph)%sig=sig_air
                ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            endif
 
            ph(iph)%soundspeed => soundspeed_sge 
            write(iout,*) '     > soundspeed    => soundspeed_sge'
            ph(iph)%pressure => pressure_sge
            write(iout,*) '     > pressure      => pressure_sge'
            ph(iph)%energie_hydro => energie_hydro_sge
            write(iout,*) '     > energie_hydro => energie_hydro_sge'
            ph(iph)%pressure_star => pressure_star_sge 
            write(iout,*) '     > pressure_star => pressure_star_sge'
            ph(iph)%rh_from_PT => rh_from_PT_sgeT
            write(iout,*) '     > rh_from_PT     => rh_from_PT_sgeT'
   
        elseif(ph(iph)%typ.eq.2)then
 
            call read_table(iph,trim(adjustl(ph(iph)%filetab)))
  
            !ph(iph)%soundspeed => soundspeed_tab 
            !ph(iph)%pressure => pressure_tab 
            !ph(iph)%energie_hydro => energie_hydro_tab 
            !ph(iph)%pressure_star => pressure_star_sge 
 
            ph(iph)%soundspeed => soundspeed_sge 
            write(iout,*) '     > soundspeed    => soundspeed_sge'
            ph(iph)%pressure => pressure_sge
            write(iout,*) '     > pressure      => pressure_sge'
            ph(iph)%energie_hydro => energie_hydro_sge
            write(iout,*) '     > energie_hydro => energie_hydro_sge'
            ph(iph)%pressure_star => pressure_star_sge 
            write(iout,*) '     > pressure_star => pressure_star_sge'
            ph(iph)%rh_from_PT => rh_from_PT_sgeT
            write(iout,*) '     > rh_from_PT     => rh_from_PT_sgeT'

        elseif(ph(iph)%typ.eq.3)then
 
            imat=ph(iph)%imat
            ph(iph)%p_sge=Input%mat(imat)%pinf
            ph(iph)%g_sge=Input%mat(imat)%gam
            ph(iph)%mu   =Input%mat(imat)%mu
            ph(iph)%sigy =Input%mat(imat)%sigy
            ph(iph)%rho0 =Input%mat(imat)%rho0
            ph(iph)%cv   =Input%mat(imat)%cv
            ph(iph)%sig  =Input%mat(imat)%sigma
            ph(iph)%ener0=( 1.0e5_PR + ph(iph)%g_sge*ph(iph)%p_sge )/( (ph(iph)%g_sge - 1.0_PR)*ph(iph)%rho0 )
 
            ph(iph)%soundspeed => soundspeed_sge 
            write(iout,*) '     > soundspeed    => soundspeed_sge'
            ph(iph)%pressure => pressure_sge
            write(iout,*) '     > pressure      => pressure_sge'
            ph(iph)%energie_hydro => energie_hydro_sge
            write(iout,*) '     > energie_hydro => energie_hydro_sge'
            ph(iph)%pressure_star => pressure_star_sge 
            write(iout,*) '     > pressure_star => pressure_star_sge'
            ph(iph)%T_from_rP => T_from_rP_sge
            write(iout,*) '     > T_from_rP     => T_from_rP_sge'
            ph(iph)%rh_from_PT => rh_from_PT_sgeT
            write(iout,*) '     > rh_from_PT     => rh_from_PT_sgeT'

        elseif(ph(iph)%typ.eq.4)then
 
            imat=ph(iph)%imat
            ph(iph)%p_sge = Input%mat(imat)%pinf
            ph(iph)%g_sge = Input%mat(imat)%gam
            ph(iph)%q     = Input%mat(imat)%q
            ph(iph)%qp    = Input%mat(imat)%qp
            ph(iph)%mu    = Input%mat(imat)%mu
            ph(iph)%sigy  = Input%mat(imat)%sigy
            ph(iph)%rho0  = Input%mat(imat)%rho0
            ph(iph)%cv    = Input%mat(imat)%cv
            ph(iph)%sig   = Input%mat(imat)%sigma
 
            ph(iph)%soundspeed => soundspeed_sgeT 
            write(iout,*) '     > soundspeed    => soundspeed_sgeT'
            ph(iph)%pressure => pressure_sgeT
            write(iout,*) '     > pressure      => pressure_sgeT'
            ph(iph)%energie_hydro => energie_hydro_sgeT
            write(iout,*) '     > energie_hydro => energie_hydro_sgeT'
            ph(iph)%pressure_star => pressure_star_sgeT 
            write(iout,*) '     > pressure_star => pressure_star_sgeT'
            ph(iph)%T_from_rP => T_from_rP_sgeT
            write(iout,*) '     > T_from_rP     => T_from_rP_sgeT'
            ph(iph)%rh_from_PT => rh_from_PT_sgeT
            write(iout,*) '     > rh_from_PT     => rh_from_PT_sgeT'
 
        endif
 
    enddo 
 
    do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl
        M%MF(i,j,k)%F(iph)%p_sge = ph(iph)%p_sge
        M%MF(i,j,k)%F(iph)%g_sge = ph(iph)%g_sge
        M%MF(i,j,k)%F(iph)%mu    = ph(iph)%mu
        M%MF(i,j,k)%F(iph)%sigy  = ph(iph)%sigy
        M%MF(i,j,k)%F(iph)%rh0   = ph(iph)%rho0
        M%MF(i,j,k)%F(iph)%rh    = ph(iph)%rho0
        M%MF(i,j,k)%F(iph)%Cv    = ph(iph)%Cv
    enddo; enddo; enddo ; enddo
 
end subroutine init_phase

subroutine set_TP(M,T,P)

   implicit none
   type(mesh), intent(inout) :: M
   real(PR), intent(in) :: T, P
   real(PR) :: u, rho
   integer :: i,j,k,iph,ir,iT
   logical :: found

   M%MF(:,:,:)%T=T

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl    
   
   if(ph(iph)%typ.eq.1)then

      M%MF(i,j,k)%F(iph)%p_sge = ph(iph)%p_sge
      M%MF(i,j,k)%F(iph)%g_sge = ph(iph)%g_sge
      M%MF(i,j,k)%F(iph)%mu    = ph(iph)%mu
      M%MF(i,j,k)%F(iph)%sigy  = ph(iph)%sigy
      M%MF(i,j,k)%F(iph)%rh0   = ph(iph)%rho0
      M%MF(i,j,k)%F(iph)%cv    = ph(iph)%cv

      M%MF(i,j,k)%F(iph)%p=P
      M%MF(i,j,k)%F(iph)%T=T

      M%MF(i,j,k)%F(iph)%eh0=energie_hydro_sge(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh0,1.0e5_PR)
      u=u_from_TP_sge(M%MF(i,j,k)%F(iph),T,P)
      rho=r_from_uP_sge(M%MF(i,j,k)%F(iph),u,P)
      M%MF(i,j,k)%F(iph)%eh=u
      M%MF(i,j,k)%F(iph)%rh=rho

   elseif(ph(iph)%typ.eq.2)then

      M%MF(i,j,k)%F(iph)%p=P
      M%MF(i,j,k)%F(iph)%T=T

      !if(i.eq.1) print*, ' get_index_TP'      
      call get_index_TP(iph,T,P,ir,iT)
      !if(i.eq.1) print*, ' r_from_rT'      

      rho=r_from_TP(iph,T,P)
      if(i.eq.1) print*, ' coucou rho from table=', rho

      M%MF(i,j,k)%F(iph)%p_sge = pinf_from_rT(iph,rho,T)
      M%MF(i,j,k)%F(iph)%g_sge = g_from_rT(iph,rho,T)

      u=E_from_rT(iph,rho,T)
      if(i.eq.1) print*, ' coucou u from table=', u      
      u=ph(iph)%energie_hydro(M%MF(i,j,k)%F(iph),rho,P)
      if(i.eq.1) print*, ' coucou u sge=', u      


      M%MF(i,j,k)%F(iph)%rh=rho
      M%MF(i,j,k)%F(iph)%eh=u

      M%MF(i,j,k)%F(iph)%p=P_from_rT(iph,rho,T)
      if(i.eq.1) print*, ' coucou P from table=', M%MF(i,j,k)%F(iph)%p
      M%MF(i,j,k)%F(iph)%p=ph(iph)%pressure(M%MF(i,j,k)%F(iph),rho,u)
      if(i.eq.1) print*, ' coucou P sge=', M%MF(i,j,k)%F(iph)%p


      M%MF(i,j,k)%F(iph)%mu   = 0.0_PR
      M%MF(i,j,k)%F(iph)%sigy = 0.0_PR
      M%MF(i,j,k)%F(iph)%rh0  = rho
      M%MF(i,j,k)%F(iph)%eh0  =energie_hydro_sge(M%MF(i,j,k)%F(iph),M%MF(i,j,k)%F(iph)%rh0,1.0e5_PR)
      M%MF(i,j,k)%F(iph)%cv   = Cv_from_rT(iph,rho,T)
      M%MF(i,j,k)%F(iph)%ZTF  = Z_from_rT(iph,rho,T)

      if(i.eq.1.and.j.eq.1.and.k.eq.1)then
         print*, 'set_TP: iph=', iph, M%MF(i,j,k)%F(iph)%g_sge, M%MF(i,j,k)%F(iph)%p_sge, M%MF(i,j,k)%F(iph)%cv
      endif

   endif

   enddo; enddo; enddo ; enddo

   !!!---Grandeurs de mélange : 
   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
      M%MF(i,j,k)%rh=sum(M%MF(i,j,k)%F(1:Nl)%f*M%MF(i,j,k)%F(1:Nl)%rh)
   enddo ; enddo ; enddo

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx ; do iph=1,M%Nl 
      M%MF(i,j,k)%F(iph)%Y=M%MF(i,j,k)%F(iph)%f*M%MF(i,j,k)%F(iph)%rh/M%MF(i,j,k)%rh
   enddo ; enddo ; enddo ; enddo

   write(iout,*) ' MAJ meca'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
    call MAJ_meca(M%MF(i,j,k))
   enddo ; enddo ; enddo

   write(iout,*) ' MAJ mixt'

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_mixt(M%MF(i,j,k))
   enddo ; enddo ; enddo

   do k=1,M%Nz ; do j=1,M%Ny ; do i=1,M%Nx
   call MAJ_THERMO(M%MF(i,j,k))
   enddo ; enddo ; enddo

end subroutine set_TP

subroutine set_TP0(MF,T,P)

   implicit none
   type(multifluide), intent(inout) :: MF
   real(PR), intent(in) :: T, P
   real(PR) :: u, rho
   integer :: i,j,k,iph,ir,iT
   logical :: found

   MF%T=T

   do iph=1,Nl    
   
   if(ph(iph)%typ.eq.1)then

      MF%F(iph)%p_sge = ph(iph)%p_sge
      MF%F(iph)%g_sge = ph(iph)%g_sge
      MF%F(iph)%mu    = ph(iph)%mu
      MF%F(iph)%sigy  = ph(iph)%sigy
      MF%F(iph)%rh0   = ph(iph)%rho0
      MF%F(iph)%cv    = ph(iph)%cv

      MF%F(iph)%p=P
      MF%F(iph)%T=T

      u=u_from_TP_sge(MF%F(iph),T,P)
      rho=r_from_uP_sge(MF%F(iph),u,P)
      MF%F(iph)%eh=u
      MF%F(iph)%rh=rho
      MF%F(iph)%eh0=u

   elseif(ph(iph)%typ.eq.2)then

      MF%F(iph)%p=P
      MF%F(iph)%T=T

      !if(i.eq.1) print*, ' get_index_TP'      
      call get_index_TP(iph,T,P,ir,iT)
      !if(i.eq.1) print*, ' r_from_rT'      

      rho=r_from_TP(iph,T,P)
      if(i.eq.1) print*, ' coucou rho from table=', rho

      MF%F(iph)%p_sge = pinf_from_rT(iph,rho,T)
      MF%F(iph)%g_sge = g_from_rT(iph,rho,T)

      u=E_from_rT(iph,rho,T)
      if(i.eq.1) print*, ' coucou u from table=', u      
      u=ph(iph)%energie_hydro(MF%F(iph),rho,P)
      if(i.eq.1) print*, ' coucou u sge=', u      


      MF%F(iph)%rh=rho
      MF%F(iph)%eh=u

      MF%F(iph)%p=P_from_rT(iph,rho,T)
      if(i.eq.1) print*, ' coucou P from table=', MF%F(iph)%p
      MF%F(iph)%p=ph(iph)%pressure(MF%F(iph),rho,u)
      if(i.eq.1) print*, ' coucou P sge=', MF%F(iph)%p

      MF%F(iph)%mu   = 0.0_PR
      MF%F(iph)%sigy = 0.0_PR
      MF%F(iph)%rh0  = rho
      MF%F(iph)%eh0  = energie_hydro_sge(MF%F(iph),MF%F(iph)%rh0,1.0e5_PR)
      MF%F(iph)%cv   = Cv_from_rT(iph,rho,T)
      MF%F(iph)%ZTF  = Z_from_rT(iph,rho,T)

      if(i.eq.1.and.j.eq.1.and.k.eq.1)then
         print*, 'set_TP: iph=', iph, MF%F(iph)%g_sge, MF%F(iph)%p_sge, MF%F(iph)%cv
      endif

   endif

   enddo

   !!!---Grandeurs de mélange : 
   MF%rh=sum(MF%F(1:Nl)%f*MF%F(1:Nl)%rh)

   do iph=1,Nl 
      MF%F(iph)%Y=MF%F(iph)%f*MF%F(iph)%rh/MF%rh
   enddo 

   call MAJ_meca(MF)

   call MAJ_mixt(MF)

   call MAJ_THERMO(MF)

end subroutine set_TP0


subroutine set_rT(MF,T)

   implicit none
   type(multifluide), intent(inout) :: MF
   real(PR), intent(in) :: T
   integer :: i,j,k,iph,ir,iT
   logical :: found

   MF%T=T

   do iph=1,Nl    
   
   if(ph(iph)%typ.eq.1)then

      MF%F(iph)%p_sge = ph(iph)%p_sge
      MF%F(iph)%g_sge = ph(iph)%g_sge
      MF%F(iph)%mu    = ph(iph)%mu
      MF%F(iph)%sigy  = ph(iph)%sigy
      MF%F(iph)%rh0   = ph(iph)%rho0
      MF%F(iph)%cv    = ph(iph)%cv

      !F%rhl(iph)=rho
      MF%F(iph)%T=T
      
      MF%F(iph)%eh0=energie_hydro_sge(MF%F(iph),MF%F(iph)%rh0,1.0e5_PR)
   
      MF%F(iph)%eh=MF%F(iph)%eh0+MF%F(iph)%cv*(T-300.0_PR)

      MF%F(iph)%p=ph(iph)%pressure(MF%F(iph),MF%F(iph)%rh,MF%F(iph)%eh)

   endif

   enddo

   !!!---Grandeurs de mélange : 
   MF%rh=sum(MF%F(1:Nl)%f*MF%F(1:Nl)%rh)

   do iph=1,Nl 
      MF%F(iph)%Y=MF%F(iph)%f*MF%F(iph)%rh/MF%rh
   enddo 

   call MAJ_meca(MF)

   call MAJ_mixt(MF)

   call MAJ_THERMO(MF)

end subroutine set_rT



!!!-- MATHS RELATIVES AU TENSEUR DE FINGER

subroutine init_Finger(F)

   !!!  Pré-calcul certaines quantités 
   !!!  pour accélérer les calculs.
   !!!  A appeler avant d'accéder aux
   !!!  routines de calcul du tenseur de Finger 
   !!!  ou de ses dérivées !!! 

   implicit none
   type(fluide), intent(inout) :: F
   real(PR) :: coeff
   real(PR) :: t1,t2

   !call cpu_time(t1)

   !coeff=F%fl(iph)**(-r13)
   a(1:3)=F%a(1:3)   !*coeff
   b(1:3)=F%b(1:3)   !*coeff
   c(1:3)=F%c(1:3)   !*coeff

   !a(1:3)=F%al(iph,1:3)
   !b(1:3)=F%bl(iph,1:3)
   !c(1:3)=F%cl(iph,1:3)

   aa=a(1)**2+a(2)**2+a(3)**2
   bb=b(1)**2+b(2)**2+b(3)**2
   cc=c(1)**2+c(2)**2+c(3)**2

   ab=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
   ac=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)
   bc=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)

   ab2=ab**2
   ac2=ac**2
   bc2=bc**2

   call Finger_tensor(FT)
   call Square_Finger_tensor(FT2)

   detFT=det3x3(FT)

   !detFTm16=detFT**(-1.0_PR/6.0_PR)

   detFT12=sqrt(detFT)
   detFTm13=detFT**(-1.0_PR/3.0_PR)
   detFTm23=detFTm13**2
   !detFTm43=detFTm23**2
   !detFTm53=detFTm43*detFTm13
   murr0=F%mu*F%rh/F%rh0

   !!!---normalisation de FT et FT2
   FT(:,:)=FT(:,:)*detFTm13
   FT2(:,:)=FT2(:,:)*detFTm23


   J1=FT(1,1)+FT(2,2)+FT(3,3)
   J2=FT2(1,1)+FT2(2,2)+FT2(3,3)

   mu=F%mu
   
   !rho=F%rhl(iph)
   !alpham1=1.0_PR/F%fl(iph)
   !alpharhom1=1.0_PR/(F%fl(iph)*F%rhl(iph))

   !call cpu_time(t2)

   !t_init_Finger=t_init_Finger+t2-t1

end subroutine init_Finger

subroutine Finger_tensor(G)

   implicit none
   real(PR), intent(out) :: G(1:3,1:3)

   !G(1,1:3)=(/ aa   , a(2)   , a(3)   /)
   !G(2,1:3)=(/ a(2) , 1.0_PR , 0.0_PR /)
   !G(3,1:3)=(/ a(3) , 0.0_PR , 1.0_PR /)


   G(1,1)=aa
   G(1,2)=ab
   G(1,3)=ac

   G(2,1)=G(1,2)
   G(2,2)=bb
   G(2,3)=bc

   G(3,1)=G(1,3)
   G(3,2)=G(2,3)
   G(3,3)=cc

end subroutine Finger_tensor

subroutine Square_Finger_tensor(G2)

   implicit none
   real(PR), intent(out) :: G2(1:3,1:3)

   G2(1,1)=aa**2+ab2+ac2
   G2(2,2)=ab2+bb**2+bc2
   G2(3,3)=ac2+bc2+cc**2

   G2(1,2)=ab*(aa+bb)+ac*bc
   G2(1,3)=ac*(aa+cc)+ab*bc
   G2(2,3)=bc*(bb+cc)+ab*ac

   G2(2,1)=G2(1,2)
   G2(3,1)=G2(1,3)
   G2(3,2)=G2(2,3)

end subroutine Square_Finger_tensor

!!!---dérivées du tenseur de Finger

real(PR) function dG_11dai(i)
   implicit none
   integer, intent(in) :: i
   dG_11dai=2.0_PR*a(i)
end function dG_11dai

real(PR) function dG_11dbi(i)
   implicit none
   integer, intent(in) :: i
   dG_11dbi=0.0_PR
end function dG_11dbi

real(PR) function dG_11dci(i)
   implicit none
   integer, intent(in) :: i
   dG_11dci=0.0_PR
end function dG_11dci

real(PR) function dG_22dai(i)
   implicit none
   integer, intent(in) :: i
   dG_22dai=0.0_PR
end function dG_22dai

real(PR) function dG_22dbi(i)
   implicit none
   integer, intent(in) :: i
   dG_22dbi=2.0_PR*b(i)
end function dG_22dbi

real(PR) function dG_22dci(i)
   implicit none
   integer, intent(in) :: i
   dG_22dci=0.0_PR
end function dG_22dci

real(PR) function dG_33dai(i)
   implicit none
   integer, intent(in) :: i
   dG_33dai=0.0_PR
end function dG_33dai

real(PR) function dG_33dbi(i)
   implicit none
   integer, intent(in) :: i
   dG_33dbi=0.0_PR
end function dG_33dbi

real(PR) function dG_33dci(i)
   implicit none
   integer, intent(in) :: i
   dG_33dci=2.0_PR*c(i)
end function dG_33dci

real(PR) function dG_12dai(i)
   implicit none
   integer, intent(in) :: i
   dG_12dai=b(i)
end function dG_12dai

real(PR) function dG_12dbi(i)
   implicit none
   integer, intent(in) :: i
   dG_12dbi=a(i)
end function dG_12dbi

real(PR) function dG_12dci(i)
   implicit none
   integer, intent(in) :: i
   dG_12dci=0.0_PR
end function dG_12dci

real(PR) function dG_13dai(i)
   implicit none
   integer, intent(in) :: i
   dG_13dai=c(i)
end function dG_13dai

real(PR) function dG_13dbi(i)
   implicit none
   integer, intent(in) :: i
   dG_13dbi=0.0_PR
end function dG_13dbi

real(PR) function dG_13dci(i)
   implicit none
   integer, intent(in) :: i
   dG_13dci=a(i)
end function dG_13dci

real(PR) function dG_23dai(i)
   implicit none
   integer, intent(in) :: i
   dG_23dai=0.0_PR
end function dG_23dai

real(PR) function dG_23dbi(i)
   implicit none
   integer, intent(in) :: i
   dG_23dbi=c(i)
end function dG_23dbi

real(PR) function dG_23dci(i)
   implicit none
   integer, intent(in) :: i
   dG_23dci=b(i)
end function dG_23dci

!!!---dérivées du carré du tenseur Finger

real(PR) function dG2_11dai(i)
   implicit none
   integer, intent(in) :: i
   dG2_11dai=2.0_PR*(2.0_PR*aa*a(i)+ab*b(i)+ac*c(i))
end function dG2_11dai

real(PR) function dG2_11dbi(i)
   implicit none
   integer, intent(in) :: i
   dG2_11dbi=2.0_PR*ab*a(i)
end function dG2_11dbi

real(PR) function dG2_11dci(i)
   implicit none
   integer, intent(in) :: i
   dG2_11dci=2.0_PR*ac*a(i)
end function dG2_11dci

real(PR) function dG2_22dai(i)
   implicit none
   integer, intent(in) :: i
   dG2_22dai=2.0_PR*ab*b(i)
end function dG2_22dai

real(PR) function dG2_22dbi(i)
   implicit none
   integer, intent(in) :: i
   dG2_22dbi=2.0_PR*(2.0_PR*bb*b(i)+ab*a(i)+bc*c(i))
end function dG2_22dbi

real(PR) function dG2_22dci(i)
   implicit none
   integer, intent(in) :: i
   dG2_22dci=2.0_PR*bc*b(i)
end function dG2_22dci

real(PR) function dG2_33dai(i)
   implicit none
   integer, intent(in) :: i
   dG2_33dai=2.0_PR*ac*c(i)
end function dG2_33dai

real(PR) function dG2_33dbi(i)
   implicit none
   integer, intent(in) :: i
   dG2_33dbi=2.0_PR*bc*c(i)
end function dG2_33dbi

real(PR) function dG2_33dci(i)
   implicit none
   integer, intent(in) :: i
   dG2_33dci=2.0_PR*(2.0_PR*cc*c(i)+ac*a(i)+bc*b(i))
end function dG2_33dci

real(PR) function dG2_12dai(i)
   implicit none
   integer, intent(in) :: i
   dG2_12dai=2.0_PR*ab*a(i)+(aa+bb)*b(i)+bc*c(i)
end function dG2_12dai

real(PR) function dG2_12dbi(i)
   implicit none
   integer, intent(in) :: i
   dG2_12dbi=2.0_PR*ab*b(i)+(aa+bb)*a(i)+ac*c(i)
end function dG2_12dbi

real(PR) function dG2_12dci(i)
   implicit none
   integer, intent(in) :: i
   dG2_12dci=ac*b(i)+bc*a(i)
end function dG2_12dci

real(PR) function dG2_13dai(i)
   implicit none
   integer, intent(in) :: i
   dG2_13dai=2.0_PR*ac*a(i)+(aa+cc)*c(i)+bc*b(i)
end function dG2_13dai

real(PR) function dG2_13dbi(i)
   implicit none
   integer, intent(in) :: i
   dG2_13dbi=ab*c(i)+bc*a(i)
end function dG2_13dbi

real(PR) function dG2_13dci(i)
   implicit none
   integer, intent(in) :: i
   dG2_13dci=2.0_PR*ac*c(i)+(aa+cc)*a(i)+ab*b(i)
end function dG2_13dci

real(PR) function dG2_23dai(i)
   implicit none
   integer, intent(in) :: i
   dG2_23dai=ab*c(i)+ac*b(i)
end function dG2_23dai

real(PR) function dG2_23dbi(i)
   implicit none
   integer, intent(in) :: i
   dG2_23dbi=2.0_PR*bc*b(i)+(bb+cc)*c(i)+ac*a(i)
end function dG2_23dbi

real(PR) function dG2_23dci(i)
   implicit none
   integer, intent(in) :: i
   dG2_23dci=2.0_PR*bc*c(i)+(bb+cc)*b(i)+ab*a(i)
end function dG2_23dci

!!!---dérivées du determinant du tenseur Finger

real(PR) function ddetGdai(i)
   implicit none
   integer, intent(in) :: i
   ddetGdai=2.0_PR*(bb*cc*a(i)+bc*(ac*b(i)+ab*c(i))-bb*ac*c(i)-bc2*a(i)-cc*ab*b(i))
end function ddetGdai

real(PR) function ddetGdbi(i)
   implicit none
   integer, intent(in) :: i
   ddetGdbi=2.0_PR*(aa*cc*b(i)+ac*(ab*c(i)+bc*a(i))-ac2*b(i)-aa*bc*c(i)-cc*ab*a(i))
end function ddetGdbi

real(PR) function ddetGdci(i)
   implicit none
   integer, intent(in) :: i
   ddetGdci=2.0_PR*(aa*bb*c(i)+ab*(ac*b(i)+bc*a(i))-bb*ac*a(i)-aa*bc*b(i)-ab2*c(i))
end function ddetGdci

!!!---Dérivées des invariants mécaniques

real(PR) function dJ1dai(i)
   implicit none
   integer, intent(in) :: i
   dJ1dai=2.0_PR*a(i)
end function dJ1dai

real(PR) function dJ1dbi(i)
   implicit none
   integer, intent(in) :: i
   dJ1dbi=2.0_PR*b(i)
end function dJ1dbi

real(PR) function dJ1dci(i)
   implicit none
   integer, intent(in) :: i
   dJ1dci=2.0_PR*c(i)
end function dJ1dci

real(PR) function dJ2dai(i)
   implicit none
   integer, intent(in) :: i
   dJ2dai=4.0_PR*(aa*a(i)+ab*b(i)+ac*c(i))
end function dJ2dai

real(PR) function dJ2dbi(i)
   implicit none
   integer, intent(in) :: i
   dJ2dbi=4.0_PR*(bb*b(i)+ab*a(i)+bc*c(i))
end function dJ2dbi

real(PR) function dJ2dci(i)
   implicit none
   integer, intent(in) :: i
   dJ2dci=4.0_PR*(cc*c(i)+ac*a(i)+bc*b(i))
end function dJ2dci


!!!---Dérivées du tenseur de contrainte par rapport aux vecteurs de base

!!!---S11 :
real(PR) function ds11dai(i)
   implicit none
   integer, intent(in) :: i

   ds11dai=-murr0*( &
      (dG2_11dai(i)-r13*dJ2dai(i))*detFTm23&
     -r23*detFTm53*ddetGdai(i)*(FT2(1,1)-r13*J2*Id(1,1))&
     -(dG_11dai(i)-r13*dJ1dai(i))*detFTm13&
     +r13*detFTm43*ddetGdai(i)*(FT(1,1)-r13*J1) ) 

end function ds11dai

real(PR) function ds11dbi(i)
   implicit none
   integer, intent(in) :: i

   ds11dbi=-murr0*( &
      (dG2_11dbi(i)-r13*dJ2dbi(i))*detFTm23&
     -r23*detFTm53*ddetGdbi(i)*(FT2(1,1)-r13*J2*Id(1,1))&
     -(dG_11dbi(i)-r13*dJ1dbi(i))*detFTm13&
     +r13*detFTm43*ddetGdbi(i)*(FT(1,1)-r13*J1) ) 

end function ds11dbi

real(PR) function ds11dci(i)
   implicit none
   integer, intent(in) :: i

   ds11dci=-murr0*( &
      (dG2_11dci(i)-r13*dJ2dci(i))*detFTm23&
     -r23*detFTm53*ddetGdci(i)*(FT2(1,1)-r13*J2*Id(1,1))&
     -(dG_11dci(i)-r13*dJ1dci(i))*detFTm13&
     +r13*detFTm43*ddetGdci(i)*(FT(1,1)-r13*J1) ) 

end function ds11dci

!!!---S12:
real(PR) function ds12dai(i)
   implicit none
   integer, intent(in) :: i

   ds12dai=-murr0*( &
      (dG2_12dai(i)-r13*dJ2dai(i))*detFTm23&
     -r23*detFTm53*ddetGdai(i)*(FT2(1,2)-r13*J2*Id(1,2))&
     -(dG_12dai(i)-r13*dJ1dai(i))*detFTm13&
     +r13*detFTm43*ddetGdai(i)*(FT(1,2)-r13*J1) ) 

end function ds12dai

real(PR) function ds12dbi(i)
   implicit none
   integer, intent(in) :: i

   ds12dbi=-murr0*( &
      (dG2_12dbi(i)-r13*dJ2dbi(i))*detFTm23&
     -r23*detFTm53*ddetGdbi(i)*(FT2(1,2)-r13*J2*Id(1,2))&
     -(dG_12dbi(i)-r13*dJ1dbi(i))*detFTm13&
     +r13*detFTm43*ddetGdbi(i)*(FT(1,2)-r13*J1) ) 

end function ds12dbi

real(PR) function ds12dci(i)
   implicit none
   integer, intent(in) :: i

   ds12dci=-murr0*( &
      (dG2_12dci(i)-r13*dJ2dci(i))*detFTm23&
     -r23*detFTm53*ddetGdci(i)*(FT2(1,2)-r13*J2*Id(1,2))&
     -(dG_12dci(i)-r13*dJ1dci(i))*detFTm13&
     +r13*detFTm43*ddetGdci(i)*(FT(1,2)-r13*J1) ) 

end function ds12dci

!!!---S13:
real(PR) function ds13dai(i)
   implicit none
   integer, intent(in) :: i

   ds13dai=-murr0*( &
      (dG2_13dai(i)-r13*dJ2dai(i))*detFTm23&
     -r23*detFTm53*ddetGdai(i)*(FT2(1,3)-r13*J2*Id(1,3))&
     -(dG_13dai(i)-r13*dJ1dai(i))*detFTm13&
     +r13*detFTm43*ddetGdai(i)*(FT(1,3)-r13*J1) ) 

end function ds13dai

real(PR) function ds13dbi(i)
   implicit none
   integer, intent(in) :: i

   ds13dbi=-murr0*( &
      (dG2_13dbi(i)-r13*dJ2dbi(i))*detFTm23&
     -r23*detFTm53*ddetGdbi(i)*(FT2(1,3)-r13*J2*Id(1,3))&
     -(dG_13dbi(i)-r13*dJ1dbi(i))*detFTm13&
     +r13*detFTm43*ddetGdbi(i)*(FT(1,3)-r13*J1) ) 

end function ds13dbi

real(PR) function ds13dci(i)
   implicit none
   integer, intent(in) :: i

   ds13dci=-murr0*( &
      (dG2_13dci(i)-r13*dJ2dci(i))*detFTm23&
     -r23*detFTm53*ddetGdci(i)*(FT2(1,3)-r13*J2*Id(1,3))&
     -(dG_13dci(i)-r13*dJ1dci(i))*detFTm13&
     +r13*detFTm43*ddetGdci(i)*(FT(1,3)-r13*J1) ) 

end function ds13dci




real(PR) function det3x3(G)

    implicit none
    real(PR), intent(in) :: G(1:3,1:3)

    det3x3=G(1,1)*G(2,2)*G(3,3)+G(2,1)*G(3,2)*G(1,3)+G(3,1)*G(1,2)*G(2,3)&
          -G(3,1)*G(2,2)*G(1,3)-G(1,1)*G(3,2)*G(2,3)-G(2,1)*G(1,2)*G(3,3)

end function det3x3



subroutine hpsort(n, ra, ind)

  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps" in hslt.
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  implicit none  
  !-input/output variables
  integer, intent(in)   :: n  
  integer :: ind (n)  
  real(PR) :: ra (n)
  !-local variables
  integer :: i, ir, j, l, iind  
  real(PR) :: rra  

  ! initialize index array
  IF (ind (1) .eq.0) then  
     DO i = 1, n  
        ind (i) = i  
     ENDDO
  ENDIF
  ! nothing to order
  IF (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    

end subroutine hpsort

logical function hslt( a, b )
    !  compare two real number and return the result
    REAL(PR) :: a, b
    IF( abs(a-b) <  1.0e-12_PR ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
end function hslt


end module mod_euler
