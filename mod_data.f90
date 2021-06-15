module mod_data


implicit none

integer, parameter, public :: PR=selected_real_kind(8)
real(PR), parameter, public :: Pi=2.0_PR*Asin(1.0_PR)
!!!!========== CONSTANTES PHYSIQUES
real(PR), parameter :: Avogadro=6.0221409e23_PR
real(PR), parameter :: kb=1.38064852e-23_PR !!! J K-1
real(PR), parameter :: qe=1.602176e-19_PR !!! C
real(PR), parameter :: eV=qe/kb !!!! K
real(PR), parameter :: uma=1.660538921e-24_PR !!!g
real(PR), parameter :: mp=1.672621898e-24_PR  !!! g
real(PR), parameter :: mn=1.674927471e-24_PR  !!! g
real(PR), parameter :: me=9.10938356e-28_PR !!!g
real(PR), parameter :: meSI=9.10938356e-31_PR !!!g
real(PR), parameter :: Rg=8.31447_PR !!! J mol-1 K-1
real(PR), parameter :: h=6.626e-34_PR !!!! J s
real(PR), parameter :: hb=h/(2.0_PR*Pi) !!!! J s
real(PR), parameter :: clight=299792458.0_PR
real(PR), parameter :: eps0=8.854187817e-12_PR
real(PR), parameter :: RBohr=h**2*eps0/(pi*meSI*qe**2) !!! Bohr radius

!!!!============ MATERIAUX
!real(PR), dimension(92) :: A !!! g/mol ou nucléon/atome 
!real(PR), dimension(92) :: ge !!! mJ/mol/K
!real(PR), dimension(92) :: n,m !!! Constants de Lennard Jones
!real(PR), dimension(92) :: Tms !!! Temperature de fusion a la densité de reference K
!character(len=2), dimension(92) :: symbol


!!!!============ MATERIAUX
!character (len=2) :: element_name(92) !< the element abbreviated name
!real(PR) :: cohesive_energy(92) !< The cohesive energy for each element, eV
!real(PR) :: melting_point(92)   !< The melting point of each element, eV
!real(PR) :: boiling_point(92)   !< The boiling point of each element, eV
!real(PR) :: rhos(92)   !< Mass density of the element in the solid phase, g cm^-3
!real(PR) :: debye_temp(92)      !< Debeye Temperature at solid density, eV
!real(PR) :: gruneisen(92)       !< Grüneisen parameter at solid density
!real(PR) :: bulk_modulus(92)    !< Bulk modulus, GPa
!real(PR) :: bm_pderiv(92)       !< Pressure derivative of the Bulk Modulus
!!!!--critical parameter : just data for post-processing
!real(PR) :: Tcrit(92)       !< critical temperature 
!real(PR) :: Pcrit(92)       !< critical pressure
!real(PR) :: rcrit(92)     !< critical density 
!!!!---Fermi energies
!real(PR) :: E_Fermi(92)
!real(PR) :: v_Fermi(92)
!!!!---for Desjarlais Z model
!real(PR) :: Eio1(92) ! First ionisation energy (eV)
!real(PR) :: pola(92) ! atomic polarizability RBohr**3
!real(PR) :: sigSL(92) ! ratio of conductivities sig_solid/sig_liquid at Tm
!real(PR) :: sig0(92) ! reference conductivities for solids at T=300 K
!real(PR) :: sigmax(92) ! max conductivity before scaling at T=300K andrho=rhosolide 
!real(PR) :: lamb0(92) ! reference thermal conductivities for solids at T=300 K
!real(PR) :: lambmax(92) ! max thermal conductivity before scaling 

integer :: iout=6

integer :: it

logical :: solid=.true.

!logical :: SGE=.true.

logical :: MUSCL_REC=.true.

logical :: radial=.false.

logical :: initT=.false.

!integer :: Nl !!! nombre de phases

!!!---maillage
!integer :: Nx,Ny,Nz


type fluide !!! conteneur avec toutes les grandeurs primitives

   !!!---variables phasiques
   real(PR) :: f=1.0e-12_PR   !!! volume fractions (alphal)
   real(PR) :: rh=1.0_PR      !!! densité des phases
   real(PR) :: Y=0.0_PR       !!! fraction massique
   real(PR) :: p=0.0_PR       !!! pression de phase
   real(PR) :: eh=0.0_PR      !!! énergie hydro de phase
   real(PR) :: ee=0.0_PR      !!! énergie elastique de phase
   real(PR) :: c0=0.0_PR      !!! vitesse du son

   !!!--- EOS:
   integer :: iph
   !!!---Stiffened gas Equation constants (Nl)
   real(PR) :: p_sge, g_sge
   real(PR) :: q=0.0_PR, qp=0.0_PR

   !!!---Mecanical properties (Nl)
   real(PR) :: sig(1:3,1:3)=0.0_PR !!! tenseur de contrainte
   real(PR) :: mu=0.0_PR   !!! coefficient de 
   real(PR) :: sigy=1.e20_PR !!! limite d'elasticité
 
   !!!---local basis vectors
   real(PR), dimension(3) :: a=(/1.0_PR,0.0_PR,0.0_PR/) !!! composantes x des vecteurs de base
   real(PR), dimension(3) :: b=(/0.0_PR,1.0_PR,0.0_PR/) !!! composantes y des vecteurs de base
   real(PR), dimension(3) :: c=(/0.0_PR,0.0_PR,1.0_PR/) !!! composantes z des vecteurs de base
   real(PR) :: detG=0.0_PR
 
   !!!---termes sources
   real(PR) :: Qe=0.0_PR
   real(PR) :: Qpx=0.0_PR
   real(PR) :: Qpy=0.0_PR
   real(PR) :: Qpz=0.0_PR

   !!!---etat précédent
   real(PR) :: rh0=0.0_PR
   real(PR) :: eh0=0.0_PR
   real(PR) :: f0=0.0_PR
   real(PR) :: frh0=0.0_PR

   !!!---conductivité
   real(PR) :: ZTF=0.0_PR
   real(PR) :: sigma=1.0e-6_PR
   real(PR) :: J=0.0_PR !!! A m^-2

   !!!---propriétés thermiques 
   real(PR) :: T=300.0_PR
   real(PR) :: cv=1.0e3_PR !!! moyenne volumique
 
end type fluide

type multifluide

   integer :: Nl !!! nombre de fluides
   type(fluide), allocatable :: F(:) !!! Nl   

   !!!---Grandeurs physiques globales:
   real(PR) :: p,vx,vy,vz,rh,e,c,u,ee
   real(PR) :: sig(1:3,1:3) !!! Pa
   real(PR) :: sigma=1.0e-6_PR !!! S m^-1
   real(PR) :: J=0.0_PR !!! A m^-2
   real(PR) :: B=0.0_PR !!! T
   real(PR) :: T=300.0_PR
   real(PR) :: cv=1.0e3_PR !!! moyenne volumique

   !!!---termes sources
   real(PR) :: Qe=0.0_PR
   real(PR) :: Qpx=0.0_PR
   real(PR) :: Qpy=0.0_PR
   real(PR) :: Qpz=0.0_PR

   !!!---variable pour sorties
   real(PR) :: sonde(1:2)

   !!!---tag
   logical :: tag=.false.

endtype multifluide


type mesh
   integer :: Nx, Ny, Nz
   real(PR) :: Lx,Ly,Lz
   real(PR), allocatable :: x(:), y(:), z(:)
   real(PR), allocatable :: xm(:), ym(:), zm(:)
   real(PR), allocatable :: dx(:), dy(:), dz(:)
   real(PR), allocatable :: dxm(:), dym(:), dzm(:)
   !!!--- problem 1D x :
   real(PR), allocatable :: SVLx(:),SVRx(:),Sgeox(:)
   !real(PR), allocatable :: Ux(:,:)
   !real(PR), allocatable :: Gx(:,:,:)

   !!!---Fluide:
   integer :: Nl
   type(multifluide), allocatable :: MF(:,:,:) !!! Nx, Ny, Nz

   !!!---temps
   integer :: Nt
   real(PR) :: dt
   real(PR) :: dtmin=1.0_PR

end type mesh

!!!----- Objet matériau

type materiaux

    character(len=20) :: nom='--'
    integer :: Z,A,id
    integer :: typ=1 
    character(len=100) :: filemat=''
    integer :: matfiletype=1
    !---rho,T,P
    real(PR) :: rho=1.0_PR, T=300.0_PR, P=101325.0_PR, U=3.0e5_PR
    !---etat de reference
    real(PR) :: rho0=1.0_PR, T0=300.0_PR, P0=101325.0_PR, U0=3.0e5_PR
    !---propriétés thermiques
    real(PR) :: lambda=0.1_PR !!! (W/m/K)
    real(PR) :: Cv=1000.0_PR
    !---propriétés electromag
    real(PR) :: sigma=1.0e-6_PR  !!! (S/m)
    real(PR) :: epsr=1.0_PR !!! permittivitÃ© relative
    real(PR) :: mur=1.0_PR
    !---propriétés radiatives
    real(PR) :: kappa=1.0e-20_PR
    real(PR) :: eps=1.0_PR  !!! emissivité
    !---EOS SGE
    real(PR) :: pinf=0.0_PR
    real(PR) :: gam=1.4_PR
    real(PR) :: q=0.0_PR
    real(PR) :: qp=0.0_PR
    !---propriétés méca
    real(PR) :: mu=0.0_PR   !!! module de cisaillement
    real(PR) :: sigy=0.0_PR !!! limite d'elasticité
    !---propriétés geometriques
    real(PR) :: R, ep, sec
    !--- numero de table (typ=2)
    integer :: itab=0
end type materiaux

type boite
   integer :: id=0 !!! pour verifier si elle existe
   integer :: forme=0 !!! 1: rectangle, 2: sphere
   real(PR) :: xmin, xmax, ymin, ymax, zmin, zmax
   real(PR) :: xc, yc, zc, R, Ri, D, Di, L
   real(PR) :: teta0, teta1, teta2
   character(len=1) :: dir='x' !!! x, y, z
end type boite

type init_field
    type(boite) :: box
    integer :: Nvar=0
    character(len=10) :: nom(1:30)=''
    real(PR)          :: val(1:30)=0.0_PR
endtype init_field

type input_data

    integer  :: Ndim=1

    character(len=200) :: nom='no name'

    !!!---Maillage:
    integer  :: Nx=1
    integer  :: Ny=1
    integer  :: Nz=1
    real(PR) :: Lx=1.0_PR
    real(PR) :: Ly=1.0_PR
    real(PR) :: Lz=1.0_PR

    !!!---Temps:
    integer  :: Noutput=1
    real(PR) :: dtoutput=1.0e-6_PR
    real(PR) :: CFL=0.5_PR

    !!!---fluides:
    integer :: Nl=1

    !!!---EOS:
    integer :: EOS=1 ! 1: SGE, 2: SGET

    !!!---time-scheme:
    character(len=10) :: tscheme='RK1'

    !!!---soundspeed_mixt
    integer :: soundspeed_mixt=1 ! 1: frozen 2: Wood

    !!!---reconstruction
    character(len=10) :: reco='NOREC'

    !!!---Limiteur
    character(len=10) :: Limiteur='MINMOD'

    !!!---Just thermo
    integer :: just_thermo=0

    !!!---Just thermo
    integer :: isoT=0 !!! active le modèle de phases isothermes de Dr. Tholin.:
                      !!! toutes les phases dans une cellule donnée ont la même température 

    !!!---matériaux:
    integer :: Nmat
    type(materiaux), allocatable :: mat(:)

    !!!---affectations fluides -> materiaux
    integer, allocatable :: f2m(:) !!! Nl

    !!!---initialization
    integer :: Nfield=0
    type(init_field), allocatable :: field(:)

end type input_data

type(input_data) :: Input

!contains
!
!
!!!!========================================   INITIALIZATION DATAS  =============================================
!
!
!
!
!subroutine INIT_EOS
!implicit none
!
!!!!!!============ noms
!
!      element_name( 1:15) = (/ ' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne', 'Na', 'Mg', 'Al', 'Si', ' P' /)
!      element_name(16:30) = (/ ' S', ' C', 'Ar', ' K', 'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn' /)
!      element_name(31:45) = (/ 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh' /)
!      element_name(46:60) = (/ 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd' /)
!      element_name(61:75) = (/ 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', ' W', 'Re' /)
!      element_name(76:86) = (/ 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn' /)
!
!
!!!!!!============ Gruneisen coefficients
!      gruneisen( 1: 5) = (/ 1.79146_PR, 2.12563_PR, 1.67765_PR, 1.20888_PR, 1.16175_PR /)
!      gruneisen( 6:10) = (/ 1.22269_PR, 1.62934_PR, 1.37132_PR, 1.56789_PR, 1.61263_PR /)
!      gruneisen(11:15) = (/ 1.86850_PR, 1.59939_PR, 2.13600_PR, 1.51886_PR, 1.67780_PR /)
!      gruneisen(16:20) = (/ 1.65461_PR, 1.68219_PR, 1.83847_PR, 2.14301_PR, 1.86535_PR /)
!      gruneisen(21:25) = (/ 1.59661_PR, 1.43857_PR, 1.34148_PR, 1.29223_PR, 1.29483_PR /)
!      gruneisen(26:30) = (/ 1.28196_PR, 1.25924_PR, 1.25737_PR, 1.28197_PR, 1.36904_PR /)
!      gruneisen(31:35) = (/ 1.46800_PR, 1.52829_PR, 1.50856_PR, 1.60737_PR, 1.82110_PR /)
!      gruneisen(36:40) = (/ 1.89290_PR, 2.19694_PR, 1.94694_PR, 1.69033_PR, 1.52959_PR /)
!      gruneisen(41:45) = (/ 1.42424_PR, 1.36844_PR, 1.33680_PR, 1.32271_PR, 1.32631_PR /)
!      gruneisen(46:50) = (/ 1.34886_PR, 1.40126_PR, 1.49116_PR, 1.57043_PR, 1.58419_PR /)
!      gruneisen(51:55) = (/ 1.63400_PR, 1.68750_PR, 1.79625_PR, 1.97781_PR, 2.28434_PR /)
!      gruneisen(56:60) = (/ 2.00220_PR, 1.72999_PR, 1.69308_PR, 1.69829_PR, 1.68300_PR /)
!      gruneisen(61:65) = (/ 1.66813_PR, 1.67849_PR, 1.84531_PR, 1.66448_PR, 1.65060_PR /)
!      gruneisen(66:70) = (/ 1.64207_PR, 1.63544_PR, 1.62752_PR, 1.61903_PR, 1.79243_PR /)
!      gruneisen(71:75) = (/ 1.60946_PR, 1.49241_PR, 1.41479_PR, 1.37144_PR, 1.34784_PR /)
!      gruneisen(76:80) = (/ 1.33266_PR, 1.33626_PR, 1.36141_PR, 1.39305_PR, 1.52937_PR /)
!      gruneisen(81:86) = (/ 1.59162_PR, 1.61599_PR, 1.68526_PR, 1.71311_PR, 1.84537_PR, 2.10319_PR /)
!
!!!!!!============ Debye Temperatures eV
!      debye_temp( 1: 5) = (/ 0.00395978_PR, 0.000868131_PR, 0.0233086_PR, 0.0592689_PR, 0.0741264_PR /)
!      debye_temp( 6:10) = (/ 0.0871836_PR, 0.00801577_PR, 0.00875263_PR, 0.00729161_PR, 0.00466196_PR /)
!      debye_temp(11:15) = (/ 0.0148618_PR, 0.0273893_PR, 0.0301668_PR, 0.0373163_PR, 0.0141415_PR /)
!      debye_temp(16:20) = (/ 0.0155253_PR, 0.00990004_PR, 0.00609248_PR, 0.00967749_PR, 0.0208953_PR /)
!      debye_temp(21:25) = (/ 0.031323_PR, 0.0357402_PR, 0.0403533_PR, 0.0414447_PR, 0.033887_PR /)
!      debye_temp(26:30) = (/ 0.0368149_PR, 0.0365957_PR, 0.0356441_PR, 0.0302851_PR, 0.0195756_PR /)
!      debye_temp(31:35) = (/ 0.0117083_PR, 0.0219696_PR, 0.0208158_PR, 0.0128074_PR, 0.00797565_PR /)
!      debye_temp(36:40) = (/ 0.00493837_PR, 0.00647994_PR, 0.0139088_PR, 0.0213249_PR, 0.0257068_PR /)
!      debye_temp(41:45) = (/ 0.0314038_PR, 0.0334277_PR, 0.0311502_PR, 0.0322989_PR, 0.0294256_PR /)
!      debye_temp(46:50) = (/ 0.0256731_PR, 0.0198242_PR, 0.0125638_PR, 0.00990966_PR, 0.010498_PR /)
!      debye_temp(51:55) = (/ 0.0133653_PR, 0.0113588_PR, 0.00761392_PR, 0.00430148_PR, 0.00468924_PR /)
!      debye_temp(56:60) = (/ 0.0102879_PR, 0.013345_PR, 0.0128289_PR, 0.013399_PR, 0.0138865_PR /)
!      debye_temp(61:65) = (/ 0.0142989_PR, 0.0138945_PR, 0.0110479_PR, 0.014918_PR, 0.0151129_PR /)
!      debye_temp(66:70) = (/ 0.0153073_PR, 0.015503_PR, 0.0155385_PR, 0.0156902_PR, 0.0106347_PR /)
!      debye_temp(71:75) = (/ 0.0159904_PR, 0.0198009_PR, 0.0241115_PR, 0.0264354_PR, 0.0259907_PR /)
!      debye_temp(76:80) = (/ 0.0256124_PR, 0.0230175_PR, 0.0191656_PR, 0.0148823_PR, 0.00545389_PR /)
!      debye_temp(81:86) = (/ 0.00806864_PR, 0.00800644_PR, 0.00717512_PR, 0.00685952_PR, 0.00647832_PR, 0.00320344_PR /)
!
!!!!!!============ Solid density g / cc !!! 25°C 298.15 K
!      rhos( 1: 5) = (/ 0.0760_PR, 0.1248_PR, 0.535_PR, 1.848_PR, 2.460_PR /)
!      rhos( 6:10) = (/ 2.260_PR, 1.026_PR, 2.000_PR, 1.516_PR, 1.444_PR /)
!      rhos(11:15) = (/ 0.968_PR, 1.738_PR, 2.698_PR, 2.330_PR, 1.823_PR /)
!      rhos(16:20) = (/ 1.960_PR, 2.030_PR, 1.656_PR, 0.856_PR, 1.550_PR /)
!      rhos(21:25) = (/ 2.985_PR, 4.507_PR, 6.110_PR, 7.140_PR, 7.470_PR /)
!      rhos(26:30) = (/ 7.874_PR, 8.9_PR, 8.908_PR, 8.920_PR, 7.140_PR /)
!      rhos(31:35) = (/ 5.904_PR, 5.323_PR, 5.727_PR, 4.819_PR, 3.120_PR /)
!      rhos(36:40) = (/ 2.82_PR, 1.532_PR, 2.630_PR, 4.472_PR, 6.511_PR /)
!      rhos(41:45) = (/ 8.570_PR, 10.280_PR, 11.5_PR, 12.370_PR, 12.450_PR /)
!      rhos(46:50) = (/ 12.023_PR, 10.490_PR, 8.650_PR, 7.310_PR, 7.310_PR /)
!      rhos(51:55) = (/ 6.697_PR, 6.240_PR, 4.940_PR, 3.540_PR, 1.879_PR /)
!      rhos(56:60) = (/ 3.510_PR, 6.146_PR, 6.689_PR, 6.640_PR, 7.010_PR /)
!      rhos(61:65) = (/ 7.264_PR, 7.353_PR, 5.244_PR, 7.901_PR, 8.219_PR /)
!      rhos(66:70) = (/ 8.551_PR, 8.795_PR, 9.066_PR, 9.321_PR, 6.570_PR /)
!      rhos(71:75) = (/ 9.841_PR, 13.310_PR, 16.650_PR, 19.250_PR, 21.020_PR /)
!      rhos(76:80) = (/ 22.59_PR, 22.56_PR, 21.090_PR, 19.3_PR, 13.534_PR /)
!      rhos(81:86) = (/ 11.850_PR, 11.340_PR, 9.780_PR, 9.196_PR, 7.000_PR, 4.40_PR /)
!
!      
!!!!!!============ Cohesive energy eV / atom
!      cohesive_energy( 1: 5) = (/ 9.3146e-3_PR, 8.60243e-4_PR, 1.52356_PR, 3.07822_PR, 5.25474_PR /)
!      cohesive_energy( 6:10) = (/ 7.41053_PR, 0.0289166_PR, 0.0353425_PR, 0.0338915_PR, 0.0181377_PR /)
!  !!! cohesive_energy(11:15) = (/ 1.0126_PR, 1.32664_PR, 3.03676_PR, 3.72081_PR, 0.128518_PR /) !!!
!      cohesive_energy(11:15) = (/ 1.0126_PR, 1.32664_PR, 3.39000_PR, 3.72081_PR, 0.128518_PR /) !!!
!      cohesive_energy(16:20) = (/ 0.101571_PR, 0.105717_PR, 0.0673684_PR, 0.79702_PR, 1.60648_PR /)
!      cohesive_energy(21:25) = (/ 3.29587_PR, 4.40486_PR, 4.69506_PR, 3.51352_PR, 2.28016_PR /)
!      cohesive_energy(26:30) = (/ 3.59644_PR, 3.88664_PR, 3.91773_PR, 3.10931_PR, 1.23336_PR /)
!      cohesive_energy(31:35) = (/ 2.65328_PR, 3.4617_PR, 0.335806_PR, 0.269474_PR, 0.153393_PR /)
!      cohesive_energy(36:40) = (/ 0.0934867_PR, 0.746235_PR, 1.41992_PR, 3.93846_PR, 6.01134_PR /)
!      cohesive_energy(41:45) = (/ 7.15142_PR, 6.21862_PR, 5.70041_PR, 6.01134_PR, 5.13036_PR /)
!      cohesive_energy(46:50) = (/ 3.93846_PR, 2.64292_PR, 1.03644_PR, 2.38381_PR, 3.00567_PR /)
!      cohesive_energy(51:55) = (/ 0.704777_PR, 0.49749_PR, 0.216615_PR, 0.131006_PR, 0.673684_PR /)
!      cohesive_energy(56:60) = (/ 1.45101_PR, 4.14575_PR, 3.62753_PR, 3.42024_PR, 2.95385_PR /)
!      cohesive_energy(61:65) = (/ 3.00567_PR, 1.81377_PR, 1.81377_PR, 3.16113_PR, 3.05749_PR /)
!      cohesive_energy(66:70) = (/ 2.90202_PR, 2.74656_PR, 2.95385_PR, 2.59109_PR, 1.6583_PR /)
!  !!! cohesive_energy(71:75) = (/ 4.30122_PR, 6.52956_PR, 7.61781_PR, 8.2915_PR, 7.30688_PR /)
!      cohesive_energy(71:75) = (/ 4.30122_PR, 6.52956_PR, 8.10000_PR, 8.2915_PR, 7.30688_PR /)
!      cohesive_energy(76:80) = (/ 6.52956_PR, 5.80405_PR, 5.07854_PR, 3.42024_PR, 0.613571_PR /)
!      cohesive_energy(81:86) = (/ 1.71012_PR, 1.84486_PR, 1.6583_PR, 1.03644_PR, 0.414575_PR, 0.176194_PR /)
!
!!!!!!============ melting point eV
!      melting_point( 1: 5) = (/ 0.00121763_PR, 0.0000825656_PR, 0.0394307_PR, 0.135594_PR, 0.20408_PR /)
!      melting_point( 6:10) = (/ 0.332274_PR, 0.00547975_PR, 0.00476708_PR, 0.00465409_PR, 0.00213454_PR /)
!      melting_point(11:15) = (/ 0.0322327_PR, 0.0802321_PR, 0.081129_PR, 0.146632_PR, 0.0275813_PR /)
!      melting_point(16:20) = (/ 0.0337528_PR, 0.0149183_PR, 0.0072875_PR, 0.0292482_PR, 0.096919_PR /)
!      melting_point(21:25) = (/ 0.15767_PR, 0.168708_PR, 0.18974_PR, 0.189479_PR, 0.132031_PR /)
!      melting_point(26:30) = (/ 0.157409_PR, 0.153672_PR, 0.150196_PR, 0.118005_PR, 0.0602016_PR /)
!      melting_point(31:35) = (/ 0.0263263_PR, 0.105289_PR, 0.0947462_PR, 0.0429472_PR, 0.0231053_PR /)
!      melting_point(36:40) = (/ 0.0100634_PR, 0.0271563_PR, 0.0912698_PR, 0.156366_PR, 0.18496_PR /)
!      melting_point(41:45) = (/ 0.239019_PR, 0.251708_PR, 0.211207_PR, 0.22659_PR, 0.194433_PR /)
!      melting_point(46:50) = (/ 0.158878_PR, 0.107329_PR, 0.0516444_PR, 0.0373501_PR, 0.0438971_PR /)
!      melting_point(51:55) = (/ 0.0785486_PR, 0.0628072_PR, 0.0336216_PR, 0.0140231_PR, 0.0262115_PR /)
!      melting_point(56:60) = (/ 0.0869242_PR, 0.103698_PR, 0.0930949_PR, 0.104654_PR, 0.112476_PR /)
!      melting_point(61:65) = (/ 0.119342_PR, 0.116909_PR, 0.0951808_PR, 0.137854_PR, 0.141591_PR /)
!      melting_point(66:70) = (/ 0.146458_PR, 0.151847_PR, 0.153846_PR, 0.158018_PR, 0.09492_PR /)
!      melting_point(71:75) = (/ 0.168273_PR, 0.217812_PR, 0.285951_PR, 0.32115_PR, 0.300639_PR /)
!      melting_point(76:80) = (/ 0.287341_PR, 0.238063_PR, 0.177425_PR, 0.116229_PR, 0.020365_PR /)
!      melting_point(81:86) = (/ 0.0501608_PR, 0.0521997_PR, 0.0473188_PR, 0.0458152_PR, 0.049987_PR, 0.0175691_PR /)
! 
!
!
!!!!!!============ boiling point eV
!      boiling_point( 1: 5) = (/ 0.00176256_PR, 0.000366765_PR, 0.140375_PR, 0.23841_PR, 0.371384_PR /)
!      boiling_point( 6:10) = (/ 0.373731_PR, 0.00672345_PR, 0.00784373_PR, 0.00739006_PR, 0.00235269_PR /)
!      boiling_point(11:15) = (/ 0.100482_PR, 0.118473_PR, 0.242669_PR, 0.275782_PR, 0.0481184_PR /)
!      boiling_point(16:20) = (/ 0.0623909_PR, 0.0207813_PR, 0.00759169_PR, 0.0897054_PR, 0.152716_PR /)
!      boiling_point(21:25) = (/ 0.269698_PR, 0.309417_PR, 0.319846_PR, 0.25588_PR, 0.202864_PR /)
!      boiling_point(26:30) = (/ 0.272393_PR, 0.278129_PR, 0.276912_PR, 0.278129_PR, 0.102568_PR /)
!      boiling_point(31:35) = (/ 0.215292_PR, 0.268829_PR, 0.0771033_PR, 0.0832739_PR, 0.0288675_PR /)
!      boiling_point(36:40) = (/ 0.0104233_PR, 0.0835347_PR, 0.143851_PR, 0.314458_PR, 0.406931_PR /)
!      boiling_point(41:45) = (/ 0.436046_PR, 0.426921_PR, 0.394416_PR, 0.384421_PR, 0.344877_PR /)
!      boiling_point(46:50) = (/ 0.281258_PR, 0.211642_PR, 0.0904007_PR, 0.20382_PR, 0.249883_PR /)
!      boiling_point(51:55) = (/ 0.161668_PR, 0.109608_PR, 0.0397575_PR, 0.0143534_PR, 0.0820572_PR /)
!      boiling_point(56:60) = (/ 0.186264_PR, 0.3248_PR, 0.315761_PR, 0.309678_PR, 0.293164_PR /)
!      boiling_point(61:65) = (/ 0.284473_PR, 0.180441_PR, 0.156453_PR, 0.306201_PR, 0.304463_PR /)
!      boiling_point(66:70) = (/ 0.246841_PR, 0.2584_PR, 0.273001_PR, 0.193217_PR, 0.127686_PR /)
!      boiling_point(71:75) = (/ 0.319412_PR, 0.423792_PR, 0.498101_PR, 0.506531_PR, 0.510095_PR /)
!      boiling_point(76:80) = (/ 0.459339_PR, 0.408582_PR, 0.356175_PR, 0.271958_PR, 0.0547436_PR /)
!      boiling_point(81:86) = (/ 0.15176_PR, 0.175747_PR, 0.159669_PR, 0.107348_PR, 0.0530289_PR, 0.0183774_PR /)
!
!!!!!!============ Bulk Modulus GPa
!      !bulk_modulus( 1: 5) = (/ 0.166_PR, 0.47_PR, 13.10_PR, 116.8_PR, 210._PR /)
!      bulk_modulus( 1: 5) = (/ 0.166_PR, 0.0047_PR, 13.10_PR, 116.8_PR, 210._PR /)
!      !bulk_modulus( 6:10) = (/ 444.3_PR, 301._PR, 2.97_PR, 15.0_PR, 1.07_PR /)
!      bulk_modulus( 6:10) = (/ 444.3_PR, 301._PR, 2.97_PR, 1.50_PR, 1.07_PR /) ! Fix for flourine
!      bulk_modulus(11:15) = (/ 6.43_PR, 33._PR, 75.8_PR, 100.1_PR, 36.0_PR /)
!      bulk_modulus(16:20) = (/ 14.5_PR, 15.18_PR, 2.86_PR, 4.25_PR, 17.4_PR /)
!      bulk_modulus(21:25) = (/ 60.0_PR, 110.02_PR, 154._PR, 193._PR, 158._PR /)
!      bulk_modulus(26:30) = (/ 163.4_PR, 203._PR, 161._PR, 137.4_PR, 65.0_PR /)
!      bulk_modulus(31:35) = (/ 55.0_PR, 87.4_PR, 78.107_PR, 48.1_PR, 14.1_PR /)
!      bulk_modulus(36:40) = (/ 3.34_PR, 5.83_PR, 11.7_PR, 34.0_PR, 92.0_PR /)
!      bulk_modulus(41:45) = (/ 173._PR, 265._PR, 344._PR, 267._PR, 318._PR /)
!      bulk_modulus(46:50) = (/ 195._PR, 109._PR, 42.0_PR, 36.8_PR, 55.8_PR /)
!      bulk_modulus(51:55) = (/ 46.7_PR, 24.0_PR, 30.38_PR, 3.63_PR, 1.723_PR /)
!      bulk_modulus(56:60) = (/ 7.65_PR, 24.8_PR, 23.08_PR, 27.16_PR, 25.38_PR /)
!      bulk_modulus(61:65) = (/ 38.0_PR, 30.7_PR, 11.7_PR, 34.0_PR, 37.0_PR /)
!      bulk_modulus(66:70) = (/ 36.3_PR, 38.9_PR, 44.0_PR, 44.5_PR, 14.6_PR /)
!      bulk_modulus(71:75) = (/ 47.0_PR, 122._PR, 194._PR, 307._PR, 372._PR /)
!      bulk_modulus(76:80) = (/ 412._PR, 383._PR, 288.4_PR, 180._PR, 292._PR /)
!      bulk_modulus(81:86) = (/ 35.3_PR, 43.2_PR, 54.7_PR, 39.45_PR, 30.00_PR , 3.50_PR /)
!
!!!!!!============ Bulk Modulus derivativ with p
!      bm_pderiv( 1: 5) = (/ 7.102_PR, 4.01_PR, 2.80_PR, 4.6_PR, 2.23_PR /)
!      bm_pderiv( 6:10) = (/ 4.04_PR, 4.02_PR, 7.78_PR, 3.00_PR, 8.40_PR /)
!      bm_pderiv(11:15) = (/ 3.84_PR, 4.32_PR, 4.109_PR, 4.16_PR, 4.50_PR /)
!      bm_pderiv(16:20) = (/ 7.00_PR, 1.59_PR, 7.20_PR, 3.63_PR, 3.22_PR /)
!      bm_pderiv(21:25) = (/ 2.80_PR, 3.59_PR, 4.20_PR, 4.80_PR, 4.60_PR /)
!      bm_pderiv(26:30) = (/ 5.38_PR, 3.60_PR, 7.55_PR, 5.52_PR, 4.60_PR /)
!      bm_pderiv(31:35) = (/ 4.70_PR, 4.31_PR, 4.317_PR, 4.33_PR, 4.60_PR /)
!      bm_pderiv(36:40) = (/ 7.20_PR, 2.08_PR, 2.41_PR, 5.00_PR, 4.00_PR /)
!      bm_pderiv(41:45) = (/ 4.10_PR, 4.70_PR, 4.54_PR, 4.50_PR, 4.86_PR /)
!      bm_pderiv(46:50) = (/ 5.42_PR, 5.87_PR, 6.50_PR, 5.25_PR, 4.08_PR /)
!      bm_pderiv(51:55) = (/ 4.90_PR, 2.30_PR, 6.13_PR, 8.08_PR, 3.87_PR /)
!      bm_pderiv(56:60) = (/ 3.45_PR, 2.80_PR, 4.945_PR, 0.096_PR, 3.119_PR /)
!      bm_pderiv(61:65) = (/ 1.50_PR, 2.50_PR, 3.00_PR, 4.20_PR, 2.70_PR /)
!      bm_pderiv(66:70) = (/ 3.71_PR, 2.83_PR, 0.8332_PR, 2.71_PR, 1.205_PR /)
!      bm_pderiv(71:75) = (/ 3.16_PR, 3.38_PR, 3.55_PR, 4.30_PR, 4.05_PR /)
!      bm_pderiv(76:80) = (/ 4.40_PR, 3.10_PR, 5.05_PR, 5.61_PR, 5.50_PR /)
!      bm_pderiv(81:86) = (/ 4.11_PR, 4.90_PR, 4.90_PR, 4.89_PR, 6.00_PR, 7.00_PR /)      
!!!!!!============ First ionization energy (eV) NIST
!   Eio1(1)=13.59843449_PR ; Eio1(2)=24.58738880_PR ;  Eio1(3)=5.39171495_PR
!   Eio1(4)=9.322699_PR    ; Eio1(5)=8.298019_PR    ;  Eio1(6)=11.2602880_PR
!   Eio1(7)=14.53413_PR    ; Eio1(8)=13.618055_PR   ;  Eio1(9)=17.42282_PR
!   Eio1(10)=21.564540_PR  ; Eio1(11)=5.1390769_PR  ;  Eio1(12)=7.646236_PR
!   Eio1(13)=5.985769_PR   ; Eio1(14)=8.15168_PR    ;  Eio1(15)=10.486686_PR
!   Eio1(16)=10.36001_PR   ; Eio1(17)=12.967632_PR  ;  Eio1(18)=15.7596117_PR
!   Eio1(19)=4.34066369_PR ; Eio1(20)=6.1131554_PR  ;  Eio1(21)=6.56149_PR
!   Eio1(22)=6.828120_PR   ; Eio1(23)=6.746187_PR   ;  Eio1(24)=6.76651_PR
!   Eio1(25)=7.4340379_PR  ; Eio1(26)=7.9024681_PR  ;  Eio1(27)=7.88101_PR
!   Eio1(28)=7.639878_PR   ; Eio1(29)=7.726380_PR   ;  Eio1(30)=9.394197_PR
!   Eio1(31)=5.9993020_PR  ; Eio1(32)=7.899435_PR   ;  Eio1(33)=9.78855_PR
!   Eio1(34)=9.752392_PR   ; Eio1(35)=11.81381_PR   ;  Eio1(36)=13.9996053_PR
!   Eio1(37)=4.1771280_PR  ; Eio1(38)=5.69486740_PR ;  Eio1(39)=6.21726_PR
!   Eio1(40)=6.63412_PR    ; Eio1(41)=6.75885_PR    ;  Eio1(42)=7.09243_PR
!   Eio1(43)=7.11938_PR    ; Eio1(44)=7.36050_PR    ;  Eio1(45)=7.45890_PR
!   Eio1(46)=8.336839_PR   ; Eio1(47)=7.576234_PR   ;  Eio1(48)=8.993820_PR
!   Eio1(49)=5.7863556_PR  ; Eio1(50)=7.343918_PR   ;  Eio1(51)=8.608389_PR
!   Eio1(52)=9.00966_PR    ; Eio1(53)=10.451260_PR  ;  Eio1(54)=12.1298436_PR
!   Eio1(55)=3.893905695_PR; Eio1(56)=5.2116646_PR  ;  Eio1(57)=5.5769_PR
!   Eio1(58)=5.5386_PR     ; Eio1(59)=5.4702_PR     ;  Eio1(60)=5.5250_PR
!   Eio1(61)=5.577_PR      ; Eio1(62)=5.64371_PR    ;  Eio1(63)=5.670385_PR
!   Eio1(64)=6.14980_PR    ; Eio1(65)=5.8638_PR     ;  Eio1(66)=5.93905_PR
!   Eio1(67)=6.0215_PR     ; Eio1(68)=6.1077_PR     ;  Eio1(69)=6.18431_PR
!   Eio1(70)=6.254160_PR   ; Eio1(71)=5.425871_PR   ;  Eio1(72)=6.825070_PR
!   Eio1(73)=7.549571_PR   ; Eio1(74)=7.86403_PR    ;  Eio1(75)=7.83352_PR
!   Eio1(76)=8.43823_PR    ; Eio1(77)=8.96702_PR    ;  Eio1(78)=8.95883_PR
!   Eio1(79)=9.225554_PR   ; Eio1(80)=10.437504_PR  ;  Eio1(81)=6.1082873_PR
!   Eio1(82)=7.4166799_PR  ; Eio1(83)=7.285516_PR   ;  Eio1(84)=8.414_PR
!   Eio1(85)=9.31751_PR    ; Eio1(86)=10.74850_PR   ;  Eio1(87)=4.0727410_PR
!   Eio1(88)=5.2784239_PR  ; Eio1(89)=5.380226_PR   ;  Eio1(90)=6.30670_PR
!   Eio1(91)=5.89_PR       ; Eio1(92)=6.19405_PR
!
!!!!!!============ Polarizability Schwerdtfeger 2014 (p=alpa*E) 
!!!!    "Table of experimental and calculated static dipole polarizabilities
!!!!     for the electronic ground states of the neutral elements"
!!!!     in atomic unit: 1 a.u. = 0.14818474 A^3 = 1.6487773e-41 C m^2 / V 
!
!   pola(1 )=4.5_PR    ; pola(2 )=1.383191_PR ; pola(3 )=164.0_PR
!   pola(4 )=37.755_PR ; pola(5 )=20.5_PR     ; pola(6 )=11.0_PR
!   pola(7 )=7.6_PR    ; pola(8 )=6.04_PR     ; pola(9 )=3.76_PR
!   pola(10)=2.670_PR  ; pola(11)=162.6_PR    ; pola(12)=71.7_PR
!   pola(13)=46_PR     ; pola(14)=36.7_PR     ; pola(15)=24.7_PR
!   pola(16)=19.6_PR   ; pola(17)=14.7_PR     ; pola(18)=11.10_PR
!   pola(19)=290.6_PR  ; pola(20)=169_PR      ; pola(21)=120_PR
!   pola(22)=99_PR     ; pola(23)=84_PR       ; pola(24)=78_PR
!   pola(25)=63_PR     ; pola(26)=57_PR       ; pola(27)=51_PR
!   pola(28)=46_PR     ; pola(29)=53.44_PR    ; pola(30)=38.8_PR
!   pola(31)=54.9_PR   ; pola(32)=39.43_PR    ; pola(33)=29.1_PR
!   pola(34)=26.24_PR  ; pola(35)=21.9_PR     ; pola(36)=17.075_PR
!   pola(37)=316_PR    ; pola(38)=186_PR      ; pola(39)=153_PR
!   pola(40)=121_PR    ; pola(41)=106_PR      ; pola(42)=86_PR
!   pola(43)=77_PR     ; pola(44)=65_PR       ; pola(45)=58_PR
!   pola(46)=32_PR     ; pola(47)=52.2_PR     ; pola(48)=49.65_PR
!   pola(49)=68.7_PR   ; pola(50)=42.4_PR     ; pola(51)=45_PR
!   pola(52)=37_PR     ; pola(53)=35.1_PR     ; pola(54)=27.815_PR
!   pola(55)=401.0_PR  ; pola(56)=268_PR      ; pola(57)=210_PR
!   pola(58)=200_PR    ; pola(59)=190_PR      ; pola(60)=212_PR
!   pola(61)=203_PR    ; pola(62)=194_PR      ; pola(63)=187_PR
!   pola(64)=159_PR    ; pola(65)=172_PR      ; pola(66)=165_PR
!   pola(67)=159_PR    ; pola(68)=153_PR      ; pola(69)=147_PR
!   pola(70)=142_PR    ; pola(71)=148_PR      ; pola(72)=109_PR
!   pola(73)=88_PR     ; pola(74)=75_PR       ; pola(75)=65_PR
!   pola(76)=57_PR     ; pola(77)=51_PR       ; pola(78)=44_PR
!   pola(79)=35.1_PR   ; pola(80)=33.91_PR    ; pola(81)=51_PR
!   pola(82)=47.1_PR   ; pola(83)=50_PR       ; pola(84)=46_PR
!   pola(85)=45.6_PR   ; pola(86)=33.18_PR    ; pola(87)=317.8_PR
!   pola(88)=246.2_PR  ; pola(89)=217_PR      ; pola(90)=217_PR
!   pola(91)=171_PR    ; pola(92)=137_PR
!
!   !!!! pola=pola*1.6487773e-41_PR !!! Cm/(V/m)
!   pola=pola*0.1481847e-30_PR !!! m^3
!!!!---numéro atomique A
!!A(13)=26.9815386_PR
!!A(26)=55.845_PR
!!A(73)=180.94788_PR
!!A(80)=200.59_PR
!!A(82)=207.2_PR
!
!!!!======== Mass number
!
!A(1:92)= (/ 1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,30,40,&
!            45,48,51,52,55,56,58,58,64,65,70,73,75,79,80,84,85,88,89,91,&
!            93,96,98,101,103,106,108,112,115,119,122,128,127,131,133,137,139,140,141,144,&
!            145,150,152,157,159,163,165,167,169,173,175,178,181,184,186,190,192,195,197,201,&
!            204,207,209,209,210,222,223,226,227,232,231,238/)
!            !237,244,243,247,247,251,252,257,&
!            !258,259,262,261,268,263,264,269,268,272,273,277,286,289,288,292,292,293 /)
!
!
!symbol(1:92)=(/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',&
!               'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',&
!               'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',&
!               'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',&
!               'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U '/)
!                 !'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','&
!                 ! Md ','No ','Lr ','Rf ','Db ','Sg ','Bh ','Hs ','Mt ','Ds ','Rg','Uub','Uut','Uuq','Uup','Uuh','Uus','Uuo'/)
!
!!!!---densité,de,référence,du,solide,g,/,cc
!!rhos(14)=2.6989_PR !!!! 2.773_PR 
!!rhos(80)=14.38_PR !!!! 13.546,14.38,16.39 
!
!!!!---Températures de fusion du solide K
!Tms(:)=melting_point(:)*eV
!
!!Tms(14)=934.0_PR !!!!K
!!Tms(80)=234.3_PR !!!!K
!
!!!!---énergie de cohesion
!!Ecoh(14)=3.39_PR !!! eV/atom
!!Ecoh(80)=0.67_PR !!! eV/atom
!
!!!!---Electronic heat capacity constant : C = ge*T + A*T^3 Kittel p157(djvu:173) ge en mJ/mol/K**2
!ge(13)=1.35_PR
!ge(26)=4.98_PR
!ge(29)=0.695_PR
!ge(73)=5.9_PR
!ge(80)=1.79_PR
!ge(82)=2.98_PR
!
!!!!---critical parameters : Ray 2017 EOS for metals in the liquid-vapor region
!!!!                         MPa, K, g cm^-3
!!!!---Na: 
!Tcrit(11)=2448.0_PR ; rcrit(11)=0.20_PR ; Pcrit(13)=0.48e2_PR 
!!!!---Al: 
!Tcrit(13)=5700.0_PR ; rcrit(13)=0.32_PR ; Pcrit(13)=1.87e2_PR 
!!!!---Fe:
!Tcrit(26)=6900.0_PR ; rcrit(26)=1.77_PR ; Pcrit(26)=8.82e2_PR
!!!!---Cu:
!Tcrit(29)=7800.0_PR ; rcrit(29)=2.31_PR ; Pcrit(29)=8.94e2_PR  
!!!!---Mo:
!Tcrit(42)=9500.0_PR ; rcrit(42)=2.82_PR ; Pcrit(42)=9.64e2_PR  
!!!!---Hg:
!Tcrit(80)=1755.0_PR ; rcrit(80)=5.99_PR ; Pcrit(80)=1.65e2_PR  
!!!!---Pb:
!Tcrit(82)=5400.0_PR ; rcrit(82)=2.80_PR ; Pcrit(82)=4.01e2_PR  
!!!!---Ta:
!Tcrit(73)=9900.0_PR ; rcrit(73)=4.19_PR ; Pcrit(73)=10.6e2_PR  
!!!!---W:
!Tcrit(74)=12000.0_PR; rcrit(74)=4.04_PR ; Pcrit(74)=10.21e2_PR  
!!!!---U:
!Tcrit(92)=7000.0_PR ; rcrit(92)=8.11_PR ; Pcrit(92)=7.11e2_PR  
!
!!!!!---Likalter 2002 (bar, K, g/cm^3) !!!!---Likalter 2002: autre ref
!!Pcrit(13)=3120.0_PR                     !Pcrit(13)=4470.0_PR 
!!Tcrit(13)=8860.0_PR                     !Tcrit(13)=8000.0_PR 
!!rcrit(13)=0.28_PR                     !rcrit(13)=0.64_PR
!
!!!!---Cu:
!!!!---Likalter 2002 (bar, K, g/cm^3) 
!!Pcrit(29)  =6510.0_PR ; Tcrit(29)  =8440.0_PR ; rcrit(29)=1.940_PR                   
!
!!!!!---Likalter 2002: autre ref
!   !Pcrit(29)=7460.0_PR ; Tcrit(29)=8390.0_PR ; !rcrit(29)=2.39.0_PR
!
!!!!---Lennard Jones parameters for metals RAY 2017
!m(11)=5.4823_PR ; n(11)=0.5898_PR !!! Na
!m(13)=1.9508_PR ; n(13)=0.4531_PR !!! Al 
!m(26)=1.8988_PR ; n(26)=0.6950_PR !!! Fe
!m(29)=2.3857_PR ; n(29)=0.6727_PR !!! Cu
!m(42)=3.2131_PR ; n(42)=0.4887_PR !!! Mo
!m(73)=2.6315_PR ; n(73)=0.5568_PR !!! Ta
!m(74)=2.6278_PR ; n(74)=0.4635_PR !!! W
!m(80)=4.6765_PR ; n(80)=0.747_PR  !!! Hg
!m(82)=1.4414_PR ; n(82)=1.3008_PR !!! Pb
!m(92)=6.5917_PR ; n(92)=0.3989_PR !!! U
!
!!m(13)=2.0_PR    ; n(13)=0.5_PR !!! Al ---RAY
!!m(42)=2.0_PR    ; n(42)=0.65_PR !!! Mo
!!m(73)=3.0_PR ; n(73)=0.5_PR !!! Ta
!!m(80)=3.4_PR    ; n(80)=1.05_PR   !!! Hg
!!m(82)=2.3_PR   ; n(82)=1.0_PR !!! Pb
!
!!!!---Default Lenard Jones
!!n(:)=m(:)=6.0_PR
!
!!!!---Fermi energy of metals (eV) : 
!E_Fermi=0.0_PR
!!!! Ashcroft, N. W. and Mermin, N. D., Solid State Physics, Saunders, 1976.
!E_Fermi(3 ) = 4.74_PR !!! Li 
!E_Fermi(4 ) = 14.3_PR !!! Be
!E_Fermi(11) = 3.24_PR !!! Na
!E_Fermi(12) = 7.08_PR !!! Mg
!E_Fermi(13) = 11.7_PR !!! Al
!E_Fermi(19) = 2.12_PR !!! K 
!E_Fermi(20) = 4.69_PR !!! Ca
!E_Fermi(25) = 10.9_PR !!! Mn
!E_Fermi(26) = 11.1_PR !!! Fe
!E_Fermi(29) = 7.00_PR !!! Cu
!E_Fermi(30) = 9.47_PR !!! Zn
!E_Fermi(31) = 10.4_PR !!! Ga
!E_Fermi(37) = 1.85_PR !!! Rb
!E_Fermi(38) = 3.93_PR !!! Sr
!E_Fermi(41) = 5.32_PR !!! Nb
!E_Fermi(47) = 5.49_PR !!! Ag
!E_Fermi(48) = 7.47_PR !!! Cd
!E_Fermi(49) = 8.63_PR !!! In
!E_Fermi(50) = 10.2_PR !!! Sn
!E_Fermi(51) = 10.9_PR !!! Sb
!E_Fermi(55) = 1.59_PR !!! Cs
!E_Fermi(56) = 3.64_PR !!! Ba
!E_Fermi(79) = 5.53_PR !!! Au
!E_Fermi(80) = 7.13_PR !!! Hg
!E_Fermi(81) = 8.15_PR !!! Tl
!E_Fermi(82) = 9.47_PR !!! Pb
!E_Fermi(83) = 9.90_PR !!! Bi
!!!! Panda 2012 : Electronic structure and equilibrium properties
!!!! of hcp titanium and zirconium
!E_Fermi(40) = 6.9389_PR !!! Zr
!E_Fermi(22) = 8.8437_PR !!! Ti
!
!
!v_Fermi(:)=sqrt(2.0_PR*E_Fermi(:)*qe/(me*1.0e-3_PR))
!
!
!!!!======== Ratio of conductivities sig_solid/sig_liquid from
!!!! "the electronic properties of liguid metals" N. E. CUSACK
!sigSL(:)=1.0_PR
!sigSL(3 )=1.64_PR
!sigSL(11)=1.451_PR
!sigSL(12)=1.78_PR
!sigSL(13)=2.20_PR
!sigSL(14)=0.034_PR
!sigSL(19)=1.56_PR
!sigSL(25)=0.61_PR
!sigSL(26)=1.09_PR
!sigSL(27)=1.05_PR
!sigSL(28)=1.3_PR
!sigSL(29)=2.04_PR
!sigSL(30)=2.24_PR
!sigSL(31)=3.12_PR
!sigSL(32)=0.053_PR
!sigSL(34)=1.0e3_PR
!sigSL(37)=1.60
!sigSL(42)=1.23_PR
!sigSL(47)=2.09_PR
!sigSL(48)=1.97_PR
!sigSL(49)=2.18_PR
!sigSL(50)=2.10_PR
!sigSL(51)=0.61_PR
!sigSL(52)=0.091_PR
!sigSL(55)=1.66_PR
!sigSL(56)=1.62_PR
!sigSL(74)=1.08_PR
!sigSL(78)=1.40_PR
!sigSL(79)=2.28_PR
!sigSL(80)=4.94_PR
!sigSL(81)= 2.06_PR
!sigSL(82)=1.94_PR
!sigSL(83)=0.35_PR
!
!!!!======== Ref conductivity (S/m) for solids at 300 K
!sig0(:)=1.0e7_PR
!sig0(13)=37.7e6_PR 
!sig0(26)=9.93e6_PR
!sig0(29)=59.6e6_PR
!
!!!!======== Max conductivity for solids at 300 K obtained before scaling from
!!!!         Lee more algorithm
!sigmax(:)=sig0(:)
!sigmax(13)=86191303.0_PR
!sigmax(26)=268630338.0_PR
!sigmax(29)=208077178.0_PR
!
!!!!======== Ref thermal conductivity (S/m) for solids at 300 K
!lamb0(:)=1.0e2_PR
!lamb0(13)=237.0_PR 
!lamb0(26)=80.2_PR
!lamb0(29)=401.0_PR
!
!!!!======== Max thermal conductivity for solids at 300 K obtained before scaling from
!!!!         Lee more algorithm
!
!lambmax(:)=lamb0(:)
!lambmax(13)=632.289_PR
!lambmax(26)=1970.64_PR
!lambmax(29)=1526.43_PR
!
!
!
!end subroutine INIT_EOS
!



end module mod_data

