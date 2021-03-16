!!
!> Equilibrium Air Chemistry (EACh)
!!
!! This model implements the curve-fits for equilibrium air plasma properties
!! given by D'Angola et. al. in Ref. 1.  All properties provided are valid for
!! pressures from 0.01 to 100.0 atm and temperatures from 50 to 60,000 K.  The
!! species considered are the following:
!!
!!    N2    N2+    N      N+     N++    N+++  N++++
!!    O2    O2+    O2-    O      O-     O+    O++
!!    O+++  O++++  NO     NO+    e-
!!
!! A number of properties are provided as function of T and P.  The following
!! are provided via curve-fits in Ref. 1:
!!
!!    1) species mole fractions
!!    2) mean molar mass
!!    3) specific enthalpy
!!    4) specific heat
!!    5) specific entropy
!!    6) electric conductivity
!!    7) thermal conductivity
!!    8) viscosity
!!
!! In addition to the above, the library is extended to provide other properties
!! that can be derived from the above such as density and species number 
!! densities.  See the public declaration below for a list of all available 
!! variables and procedures.
!!
!! References:
!!    [1]   A. D'Angola, G. Colonna, C. Gorse, and M. Capitelli, "Thermodynamic
!!          and transport properties in equilibrium air plasmas in a wide
!!          pressure and temperature range.", Eur. Phys. J. D, vol 46, 2008
!!
!! Author:
!!    J.B. Scoggins (jbscoggi@gmail.com)
!!
!!                                                                     6.15.2011
module air_ETL
    implicit none    
    private
    
    public :: NSPECIES,              &
              RU,                    &
              ONEATM,                &
              species_names,         &
              mole_fraction,         &
              mole_fractions,        &
              molar_mass,            &  ! [kg/mol]
              mass_density,          &  ! [kg/m^3]
              molar_density,         &  ! [mol/m^3]
              number_density,        &  ! [molecule/m^3]
              number_densities,      &  ! [molecule/m^3]
              specific_enthalpy,     &  ! [J/kg]
              specific_heat,         &  ! [J/kg-K]
              specific_entropy,      &  ! [J/kg-K]
              electric_conductivity, &  ! [S/m]
              thermal_conductivity,  &  ! [W/m-K]
              viscosity                 ! [kg/m-s]
              
    double precision, parameter :: RU     = 8.314472D0    ! Universal gas constant in J/mol-K
    double precision, parameter :: ONEATM = 101325.0D0    ! One atmosphere in Pa
    double precision, parameter :: NA     = 6.0221415D23  ! Avagodro's number in molecule/mol
    
    integer, parameter :: NSPECIES = 19
    character(len=5), dimension(NSPECIES), parameter :: species_names = [ &
        'N2   ', 'N2+  ', 'N    ', 'N+   ', 'N++  ',                      &
        'N+++ ', 'N++++', 'O2   ', 'O2+  ', 'O2-  ',                      & 
        'O    ', 'O-   ', 'O+   ', 'O++  ', 'O+++ ',                      &
        'O++++', 'NO   ', 'NO+  ', 'e-   '                                &
    ]
    
    integer, parameter :: MF_TYPE_1 = 1
    integer, parameter :: MF_TYPE_2 = 2
    integer, parameter :: MF_TYPE_3 = 3
    
    integer, parameter :: MAX_SIGMOIDS = 11
    
    double precision, dimension(0:MAX_SIGMOIDS) :: c
    double precision, dimension(0:MAX_SIGMOIDS) :: d
    double precision, dimension(0:MAX_SIGMOIDS) :: a
    
    double precision :: cm
    double precision :: dm
    double precision :: sm
    double precision :: Patm
    double precision :: w0
    
    interface mole_fraction
        module procedure mole_fraction_id
        module procedure mole_fraction_name
    end interface
    
    contains
    
    !!
    !> Gas density in kg/m^3.
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function mass_density(T, P) result(rho)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
    
        rho = molar_density(T, P) * molar_mass(T, P)
    end function mass_density
    
    !!
    !> Molar density in mol/m^3.
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function molar_density(T, P) result(cm)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        cm = P / (RU * T)
    end function molar_density
    
    !!
    !> Total number density in molecules/m^3.
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function number_density(T, P) result(n)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        n = molar_density(T, P) * NA
    end function number_density    
    
    !!
    !> Computes individual number densities for each species.
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!  @param n holds individual number densities on return
    !!
    subroutine number_densities(T, P, n)
        implicit none
        double precision,               intent(in) :: T
        double precision,               intent(in) :: P
        double precision, dimension(:), intent(out) :: n
        
        double precision :: n_total
        
        ! Compute total number density
        n_total = number_density(T, P)
        
        ! Multiply by mole fractions to get individual number densities
        call mole_fractions(T, P, n)
        n(1:NSPECIES) = n(1:NSPECIES) * n_total
    end subroutine
    
    !!
    !> Mean molar mass in kg/mol. (Eq. 53, Table 21)
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function molar_mass(T, P) result(Mw)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        integer :: i
        
        Patm = P / ONEATM
        
        c(0) = 0.028811D0
        c(1) = log_polynomial([8.170734D0, 5.708244D-2, 1.293374D-3], Patm)
        c(2) = log_polynomial([8.805680D0, 5.468057D-2, 1.121881D-3], Patm)
        c(3) = log_polynomial([9.525862D0, 6.639994D-2, 7.836529D-4,-2.447910D-4,-2.415297D-5], Patm)
        c(4) = log_polynomial([1.055726D1, 8.397717D-3, 9.849480D-4, 3.539965D-4,-4.236150D-5], Patm)
        c(5) = log_polynomial([1.020784D1, 2.553473D-2,-3.549988D-3], Patm)
        c(6) = log_polynomial([1.096136D1, 2.887564D-2,-3.621097D-4], Patm)
        
        d(1) = log_polynomial([6.380594D0, 1.046470D-1, 8.553860D-4,-1.572857D-4], Patm)
        d(2) = log_polynomial([7.080690D0, 1.142540D-1, 6.869247D-4,-2.257365D-4], Patm)
        d(3) = log_polynomial([7.888211D0, 9.954169D-2,-1.327510D-4,-2.926560D-4,-4.717532D-5], Patm)
        d(4) = log_polynomial([8.707609D0, 3.713173D-2,-1.670186D-2,-5.094908D-4, 4.248200D-4], Patm)
        d(5) = log_polynomial([8.422438D0, 1.125955D-1,-3.204629D-3,-1.655103D-3,-2.051312D-4], Patm)
        d(6) = log_polynomial([9.253817D0, 1.341329D-2,-6.004835D-3, 1.860800D-3,-1.229602D-4], Patm)
        
        a(1) = log_polynomial([-5.452539D0,-2.762076D-2,-3.327630D-3,-2.453118D-4,-6.332107D-6], Patm)
        a(2) = log_polynomial([-4.595514D0, 1.328152D-2, 9.294096D-4,-8.243998D-5,-9.490079D-6], Patm)
        a(3) = log_polynomial([-4.971377D0,-1.668833D-2,-2.409638D-3,-2.840529D-4,-2.934495D-5], Patm)
        a(4) = log_polynomial([-6.720756D0, 7.203127D-2, 6.766486D-3,-1.019894D-3, 9.196578D-5], Patm)
        a(5) = log_polynomial([-6.218117D0,-7.145834D-2, 6.529894D-4, 1.599394D-3, 1.981881D-5], Patm)
        a(6) = log_polynomial([-6.611171D0, 8.990124D-2,-5.418532D-3], Patm)
        
        Mw = c(0)
        do i = 1,6
            Mw = Mw - a(i)*sigmoid(T, c(i), d(i))
        end do        
    end function molar_mass
    
    !!
    !> Specific enthalpy in J/kg. (Eq. 55, Table 22)
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function specific_enthalpy(T, P) result(H)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        integer :: i
        
        Patm = P / ONEATM
        
        ! First summation
        c(1) = polynomial([2.350912D-1,-1.120236D-3,-2.508755D-5], log(Patm))
        c(2) = polynomial([1.542966D-5, 6.556647D-7], log(Patm))
        H = T * (c(1) + c(2) * T)
        
        ! Second summation
        c(1) = log_polynomial([8.164839D0, 5.283021D-2, 4.741812D-4,-1.276598D-4,-9.877950D-6], Patm)
        c(2) = log_polynomial([8.856133D0, 5.964702D-2, 1.745638D-3, 2.343688D-5,-3.102821D-6], Patm)
        c(3) = log_polynomial([9.593196D0, 7.089945D-2, 1.640521D-3,-1.055407D-4,-1.510653D-5], Patm)
        c(4) = log_polynomial([1.030572D1, 6.607308D-2, 1.512694D-3,-5.009486D-5,-5.192881D-6, 1.116840D-6], Patm)
        c(5) = log_polynomial([1.076031D1, 6.404003D-2, 9.621465D-4,-1.883920D-5], Patm)
        c(6) = log_polynomial([1.109244D1, 6.026294D-2, 1.125935D-3,-2.170126D-5,-3.141895D-6], Patm)
        c(7) = log_polynomial([1.314544D1, 2.079129D0, 9.992304D-1, 2.223931D-1, 1.963954D-2,-1.622592D-4,-1.094608D-5, 2.304744D-5, 1.817656D-6], Patm)
        
        d(1) = log_polynomial([6.513247D0, 1.040239D-1,-8.104042D-4,-2.991537D-4, 4.348437D-5, 6.258153D-6], Patm)
        d(2) = log_polynomial([6.981907D0, 1.119408D-1, 4.185626D-3,-2.499247D-4,-5.209456D-5], Patm)
        d(3) = log_polynomial([7.910995D0, 1.006930D-1,-1.608832D-3,-2.526731D-4], Patm)
        d(4) = log_polynomial([8.320951D0, 7.474585D-2, 1.789257D-3, 5.273341D-4, 3.755570D-5, 3.425485D-6], Patm)
        d(5) = log_polynomial([8.846750D0, 1.307197D-1,-2.943134D-4,-6.425060D-4], Patm)
        d(6) = log_polynomial([8.942747D0, 8.687938D-2, 1.554323D-2, 3.584506D-5,-2.447405D-4], Patm)
        d(7) = log_polynomial([-1.743314D0,-1.807206D1,-1.393980D1,-5.232064D0,-7.607736D-1, 8.529592D-2, 4.967101D-2, 7.733746D-3, 5.507513D-4, 1.527569D-5], Patm)
        
        a(1) = log_polynomial([6.587335D0,-6.112145D-2,-9.108114D-3,-9.569561D-4,-1.128838D-4,-8.757988D-6], Patm)
        a(2) = log_polynomial([8.740885D0, 3.050736D-3, 1.599171D-3,-2.859059D-4,-5.371695D-5], Patm)
        a(3) = log_polynomial([1.014496D1,-1.833015D-2,-4.265166D-3,-8.321612D-4,-6.481810D-5], Patm)
        a(4) = log_polynomial([1.082665D1,-4.777223D-2,-4.682547D-3], Patm)
        a(5) = log_polynomial([1.145937D1, 5.122940D-4,-8.805300D-3,-1.193042D-3], Patm)
        a(6) = log_polynomial([1.172458D1,-5.461477D-2, 3.413385D-3, 7.407737D-4,-1.644220D-4,-1.878043D-5], Patm)
        a(7) = log_polynomial([-1.011841D1,-2.295953D1,-1.220667D1,-3.504472D0,-4.373233D-1, 1.127311D-2, 6.598926D-3,-2.119755D-4,-1.369506D-4,-8.311253D-6], Patm)
        
        do i = 1,7
            H = H + a(i) * sigmoid(T, c(i), d(i))
        end do
        
        ! Convert from cal/g to J/kg
        H = H * 4.184D3
    end function specific_enthalpy
    
    !!
    !> Specific heat in J/kg-K. (Eq. 56, Table 23)
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function specific_heat(T, P) result(Cp)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        integer :: i
        
        Patm = P / ONEATM
        
        ! First summation
        c(0) = polynomial([ 4.303513D-3, 6.487052D-2,-6.616517D-3, 2.391262D-4], log(Patm))
        c(1) = polynomial([-2.472201D-5, 1.865503D-5,-1.298963D-6], log(Patm))
        Cp = c(0) + c(1) * T
        
        ! Second summation
        c(1) = log_polynomial([1.008338D1, 8.730410D-2, 8.590102D-3, 5.892083D-4], Patm)
        c(2) = log_polynomial([7.803107D0, 6.576559D-2, 1.214098D-4,-2.773380D-4], Patm)
        c(3) = log_polynomial([1.174382D1, 1.351866D-1,-5.421755D-2,-1.623227D-2, 7.438041D-4, 4.977237D-4, 3.456374D-5], Patm)
        c(4) = log_polynomial([1.052223D1, 8.741211D-3,-2.198211D-4], Patm)
        c(5) = log_polynomial([8.854075D0, 5.131262D-2, 1.507223D-3, 3.158892D-4], Patm)
        
        d(1) = log_polynomial([8.043134D0, 3.294414D-1, 4.681080D-2,-1.509745D-3,-4.534410D-4], Patm)
        d(2) = log_polynomial([6.212416D0, 1.085758D-1,-1.459860D-2, 4.297049D-4], Patm)
        d(3) = log_polynomial([6.596434D0,-5.025506D0,-3.238969D0,-6.901133D-1,-1.855573D-2, 8.925939D-3, 7.161553D-4], Patm)
        d(4) = log_polynomial([9.241400D0,-6.373646D-2,-7.339952D-3, 5.024652D-4], Patm)
        d(5) = log_polynomial([8.973554D0, 1.209818D-2, 2.753489D-3, 2.401117D-4], Patm)
        
        a(1) = log_polynomial([-3.497916D-1,-2.900058D-1,-3.839544D-2,-6.284876D-3,-4.130292D-4], Patm)
        a(2) = log_polynomial([-2.305816D0, 9.286290D-2,-1.095463D-2,-1.929857D-3,-2.358095D-4], Patm)
        a(3) = log_polynomial([1.885717D1, 1.564732D1, 3.946648D0, 2.094257D-1,-4.423704D-2,-5.854376D-3,-1.750164D-4], Patm)
        a(4) = log_polynomial([9.844680D-1,-2.553591D-1,-8.898889D-3, 1.493946D-3, 6.005988D-5], Patm)
        a(5) = log_polynomial([4.664598D-1,-2.233574D-1,-1.441672D-2,-1.177062D-3,-6.026800D-5], Patm)
        
        do i = 1,5
            Cp = Cp + a(i) * sigmoid(T, c(i), d(i))
        end do
        
        ! Third summation
        c( 6) = log_polynomial([8.164620D0, 5.272624D-2, 5.356645D-4,-4.303413D-5], Patm)
        c( 7) = log_polynomial([8.857074D0, 5.974192D-2, 1.621499D-3, 2.811880D-5], Patm)
        c( 8) = log_polynomial([9.598608D0, 7.030841D-2, 9.720862D-4,-8.467979D-5], Patm)
        c( 9) = log_polynomial([1.030698D1, 6.396773D-2, 1.387554D-3, 5.277379D-5], Patm)
        c(10) = log_polynomial([1.076627D1, 6.602355D-2, 1.098331D-3,-2.395208D-5], Patm)
        c(11) = log_polynomial([1.109308D1, 5.876202D-2, 1.243864D-3, 3.958414D-5], Patm)
        
        d( 6) = log_polynomial([6.414342D0, 7.141268D-2,-3.184188D-3,-3.896052D-4], Patm)
        d( 7) = log_polynomial([7.031326D0, 9.966653D-2, 2.637695D-3,-3.740228D-5], Patm)
        d( 8) = log_polynomial([7.915774D0, 9.011657D-2,-7.629395D-4,-1.579088D-4], Patm)
        d( 9) = log_polynomial([8.306511D0, 8.714253D-2, 6.812322D-4], Patm)
        d(10) = log_polynomial([8.764870D0, 1.276501D-1,-5.083689D-4,-7.452322D-4, 2.885332D-5], Patm)
        d(11) = log_polynomial([8.959570D0, 1.014329D-1, 1.073510D-2,-1.155143D-3,-2.731432D-4], Patm)
        
        a( 6) = log_polynomial([-8.639771D-1,-2.135237D-1,-1.735545D-2,-1.885139D-3,-1.226041D-4], Patm)
        a( 7) = log_polynomial([9.596448D-1,-1.130686D-1,-2.461674D-3, 4.743607D-5], Patm)
        a( 8) = log_polynomial([1.431534D0,-1.255579D-1,-5.407784D-3,-4.894608D-4], Patm)
        a( 9) = log_polynomial([1.658985D0,-1.098660D-1,-7.382403D-3,-1.597338D-3,-1.259823D-4], Patm)
        a(10) = log_polynomial([1.638978D0,-1.238859D-1,-3.036868D-3,-1.130285D-3,-1.070291D-4], Patm)
        a(11) = log_polynomial([1.933029D0,-1.248750D-1,-1.646256D-2, 5.253210D-4, 3.143929D-4], Patm)
        
        do i = 6,11
            Cp = Cp + a(i) * gaussian(T, c(i), d(i))
        end do
        
        ! Convert from cal/g-K to J/kg-K
        Cp = Cp * 4.184D3         
    end function specific_heat
    
    !!
    !> Specific entropy in J/kg-K. (Eq. 57, Table 24)
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function specific_entropy(T, P) result(S)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        integer :: i
        double precision    :: w0
        
        Patm = P / ONEATM
        
        ! Psi term
        a(0) = polynomial([7.247773D-1,-5.579293D-2, 1.246960D-3], log(Patm))
        w0   = log_polynomial([-1.949369D0, 4.114017D-2,-1.494867D-4], Patm)
        S = phi(T, a(0), w0) 
        
        ! Sum over sigmoids
        c(1) = log_polynomial([8.148014D0, 5.310698D-2, 1.031964D-3], Patm)
        c(2) = log_polynomial([8.839993D0, 5.776858D-2, 1.370456D-3], Patm)
        c(3) = log_polynomial([9.574258D0, 6.744128D-2, 1.024908D-3,-4.207616D-5], Patm)
        c(4) = log_polynomial([1.029506D1, 6.492224D-2, 1.056643D-3,-1.643501D-5], Patm)
        c(5) = log_polynomial([1.073845D1, 6.189030D-2, 1.130363D-3], Patm)
        c(6) = log_polynomial([1.108196D1, 6.043375D-2, 1.256963D-3], Patm)
        c(7) = log_polynomial([1.188475D1, 8.409849D-2, 1.321048D-3], Patm)
        
        d(1) = log_polynomial([6.487571D0, 1.023051D-1, 2.174761D-3], Patm)
        d(2) = log_polynomial([7.037019D0, 1.068246D-1, 2.382441D-3], Patm)
        d(3) = log_polynomial([8.030117D0, 1.278105D-1, 2.875819D-3, 5.995551D-5], Patm)
        d(4) = log_polynomial([8.459989D0, 1.137276D-1, 3.515754D-3], Patm)
        d(5) = log_polynomial([8.822908D0, 1.163591D-1, 3.457547D-3], Patm)
        a(6) = log_polynomial([9.251052D0, 1.085228D-1, 1.802991D-3], Patm)
        a(7) = log_polynomial([1.082135D1, 1.048287D-1,-5.310563D-3], Patm)
        
        a(1) = polynomial([2.181324D-1,-2.219875D-2,-3.107110D-4], log(Patm))
        a(2) = polynomial([9.599015D-1,-5.505086D-2, 1.666018D-5], log(Patm))
        a(3) = log_polynomial([6.970847D-1,-6.594736D-2,-2.375941D-3,-6.719048D-5], Patm)
        a(4) = log_polynomial([7.418974D-1,-6.340965D-2,-1.805794D-3,-5.053043D-5], Patm)
        a(5) = log_polynomial([7.657208D-1,-5.775822D-2,-7.902876D-4], Patm)
        a(6) = polynomial([2.767445D0,-1.126949D-1, 2.483520D-3], log(Patm))
        a(7) = log_polynomial([1.503593D0, 6.710278D-2, 1.417762D-3], Patm)
        
        do i = 1,7
            S = S + a(i) * sigmoid(T, c(i), d(i))
        end do
        
        ! Convert from cal/g-K to J/kg-K
        S = S * 4.184D3 
    end function specific_entropy
    
    !!
    !> Electric conductivity in S/m. (Eq. 58, Table 25)
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function electric_conductivity(T, P) result(sigma)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        integer :: i
        
        Patm = P / ONEATM
        
        ! Xi term
        a(0) = log_polynomial([1.635045D0, 4.450390D-2,-5.928863D-4], Patm)
        c(0) = log_polynomial([5.748398D0, 6.411299D-2], Patm)
        d(0) = log_polynomial([1.786355D0,-1.212690D-2,-2.375673D-4], Patm)
        w0   = log_polynomial([1.419925D0,-3.875497D-2], Patm)
        sigma = xi(log(T), a(0), c(0), d(0), w0)
        
        ! Sum the sigmoids
        c(1) = log_polynomial([8.930803D0, 5.718843D-2, 1.093759D-3], Patm)
        c(2) = log_polynomial([8.576847D0, 1.004174D-1, 7.406762D-3,-1.095186D-3], Patm)
        c(3) = log_polynomial([1.023493D1, 6.651575D-2, 1.090308D-3,-6.576415D-5, 4.715318D-7], Patm)
        c(4) = log_polynomial([1.072380D1, 5.671452D-2, 1.468210D-4, 2.608196D-5, 6.511719D-6], Patm)
        c(5) = log_polynomial([1.106431D1, 5.578774D-2, 6.639655D-4], Patm)
        c(6) = log_polynomial([1.023203D1, 8.703300D-2, 5.007602D-3], Patm)
        c(7) = log_polynomial([2.946755D1,-4.289010D0,-3.224136D-1, 9.371814D-2], Patm)
        
        d(1) = log_polynomial([7.014976D0, 7.625175D-2, 3.011941D-4], Patm)
        d(2) = log_polynomial([9.113182D0,-8.202725D-2, 6.299430D-3, 9.099828D-4], Patm)
        d(3) = log_polynomial([8.039563D0, 1.435966D-1, 8.862611D-3,-3.478227D-4,-3.614711D-5], Patm)
        d(4) = log_polynomial([8.556977D0, 2.227207D-1,-2.773160D-3,-1.684219D-3, 1.878188D-4], Patm)
        d(5) = log_polynomial([9.309043D0, 1.366510D-1,-2.599317D-3], Patm)
        d(6) = log_polynomial([1.130562D1,-2.184155D-2,-1.865543D-4], Patm)
        d(7) = log_polynomial([2.430324D1,-2.653523D0,-3.309222D-1, 4.769061D-2], Patm)
        
        a(1) = log_polynomial([4.493934D-2,-9.063708D-3,-2.489500D-3], Patm)
        a(2) = polynomial([1.593153D0, 4.137850D-2, 1.430491D-2,-4.403957D-7], log(Patm))
        a(3) = -polynomial([2.627897D-1, 2.917558D-3, 3.036205D-3,-1.926651D-4,-2.917018D-5], log(Patm))
        a(4) = -polynomial([1.707216D-1, 2.035164D-2, 1.809127D-3,-9.630175D-5, 1.781249D-5], log(Patm))
        a(5) = -polynomial([2.480007D-1, 2.217818D-2, 9.094614D-4], log(Patm))
        a(6) = polynomial([3.636707D0,-1.436268D-1,-2.934926D-3], log(Patm))
        a(7) = a(3) + a(4) + a(5) - a(1) - a(2) - a(6)
        
        do i = 1,7
            sigma = sigma + a(i) * sigmoid(T, c(i), d(i))
        end do
        
        ! Exponential
        sigma = exp(sigma)
    end function electric_conductivity
    
    !!
    !> Thermal conductivity in W/K-m. (Eq. 59, Table 26)
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function thermal_conductivity(T, P) result(lambda)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        integer :: i
        
        Patm = P / ONEATM
        
        ! First term
        a(0) = -1.283401D1
        lambda = a(0)
        
        ! Second term
        c(1) = 6.622230D0
        c(2) = log_polynomial([2.080280D0, 7.242844D-3, 1.959358D-4], Patm)
        c(3) = log_polynomial([2.114417D0, 6.588084D-3, 1.041527D-4], Patm)
        c(4) = log_polynomial([2.275814D0, 2.789634D-3, 1.613876D-4], Patm)
        c(5) = log_polynomial([2.361178D0, 6.072448D-3,-1.995121D-5], Patm)
        c(6) = log_polynomial([2.467469D0, 5.822255D-3], Patm)
        
        d(1) = 1.184624D1
        d(2) = log_polynomial([-1.421111D0, 7.326017D-2, 1.685275D-3], Patm)
        d(3) = log_polynomial([-1.873926D0, 7.669056D-2, 4.311158D-3], Patm)
        d(4) = log_polynomial([-1.078759D0, 1.962265D-2,-8.795026D-3, 5.277830D-4], Patm)
        d(5) = log_polynomial([-1.820338D0, 1.075866D-1], Patm)
        d(6) = log_polynomial([-2.830928D-2,-2.218935D-2,-2.718760D-3], Patm)
        
        a(1) = 1.991839D1
        a(2) = log_polynomial([7.133981D-1,-2.282818D-2, 5.491632D-4], Patm)
        a(3) = -log_polynomial([8.309337D-1,-2.699607D-3, 2.836175D-3], Patm)
        a(4) = log_polynomial([5.566144D-1, 1.402546D-1,-4.355200D-3, 1.422302D-4], Patm)
        a(5) = log_polynomial([-1.893687D0, 3.628971D-2, 9.796743D-3], Patm)
        a(6) = log_polynomial([1.153927D0,-1.647523D-2,-2.502041D-3], Patm)
        
        do i = 1,6
            lambda = lambda + a(i) * sigmoid(log(T), c(i), d(i))
        end do
        
        ! Third term
        c( 7) = log_polynomial([2.061427D0, 1.117607D-3,-3.916231D-4], Patm)
        c( 8) = log_polynomial([2.205458D0, 6.659429D-3, 1.324918D-4], Patm)
        c( 9) = log_polynomial([2.183883D0, 5.938113D-3, 1.877191D-4, 4.341127D-6], Patm)
        c(10) = log_polynomial([2.255570D0, 5.826924D-3,-5.486194D-5,-1.143664D-5], Patm)
        
        d( 7) = log_polynomial([-3.290998D1,-3.353576D0,-5.634466D-1], Patm)
        d( 8) = log_polynomial([-1.819779D0, 3.825355D-3, 1.202891D-3], Patm)
        d( 9) = log_polynomial([-9.494270D-1, 3.609984D-2, 1.528015D-3,-9.686251D-5], Patm)
        d(10) = log_polynomial([-1.374699D0, 2.577156D-2,-1.763376D-3], Patm)
        
        a( 7) = log_polynomial([-1.700917D2,-2.131620D1,-3.099200D0], Patm)
        a( 8) = -log_polynomial([-1.456072D-1,-1.437036D-1,-1.480764D-3], Patm)
        a( 9) = log_polynomial([1.055279D0,-2.677491D-2, 2.446759D-3], Patm)
        a(10) = log_polynomial([2.885339D-1,-7.133722D-2, 2.612269D-4, 2.585150D-4], Patm)
        
        do i = 7,10
            lambda = lambda + a(i) * gaussian(log(T), c(i), d(i))
        end do
        
        ! Exponential
        lambda = exp(lambda)
    end function thermal_conductivity
    
    !!
    !> Dynamic viscosity in kg/m-s. (Eq. 60, Table 27)
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function viscosity(T, P) result(eta)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: P
        
        integer :: i
        
        Patm = P / ONEATM
        
        ! Xi term
        a(0) = -2.490318D-3
        c(0) =  1.138022D-5
        d(0) =  2.525477D2
        w0   =  1.0D0
        eta = xi(T, a(0), c(0), d(0), w0)
        
        ! First sigmoid sum
        c(1) = -1.158238D4
        c(2) = log_polynomial([8.758933D0, 5.609203D-2, 7.113878D-4], Patm)
        c(3) = log_polynomial([8.183393D0, 5.531418D-2, 1.161696D-3], Patm)
        c(4) = log_polynomial([9.196899D0, 6.227176D-2, 1.047858D-3,-1.062417D-4,-1.844923D-6], Patm)
        c(5) = log_polynomial([1.054992D1, 6.447025D-2,-3.834145D-4,-3.294639D-5,-3.605812D-6], Patm)
        
        d(1) = 8.370196D3
        d(2) = log_polynomial([7.621521D0, 6.802267D-2,-4.173943D-3], Patm)
        d(3) = log_polynomial([6.544301D0, 1.395158D-1, 1.269937D-3], Patm)
        d(4) = log_polynomial([7.345945D0, 9.087033D-2,-2.859605D-3,-1.787083D-4, 1.598906D-4], Patm)
        d(5) = log_polynomial([8.500778D0, 7.811525D-2, 4.703012D-3,-1.262204D-4,-1.791684D-5], Patm)
        
        a(1) = 2.658346D-3
        a(2) = log_polynomial([-9.146259D0, 9.214388D-2,-7.532526D-3], Patm)
        a(3) = log_polynomial([-1.077843D1, 8.010183D-2, 5.530383D-3], Patm)
        a(4) = log_polynomial([-9.136467D0, 4.321416D-2,-1.415683D-2,-9.589284D-4, 2.418933D-4, 5.834458D-6], Patm)
        a(5) = log_polynomial([-1.924773D-2,-1.929031D-1,-7.597965D-2, 1.232504D-3, 2.797944D-4], Patm)
        
        do i = 1,5
            eta = eta + a(i) * sigmoid(T, c(i), d(i))
        end do
        eta = log(eta)
        
        ! Second sigmoid sum
        c( 6) = log_polynomial([9.648995D0, 6.284331D-2, 8.307533D-4,-5.453268D-6], Patm)
        c( 7) = log_polynomial([1.029898D1, 6.646081D-2, 9.291080D-4,-2.151764D-5], Patm)
        c( 8) = log_polynomial([1.077964D1, 6.865954D-2, 1.085963D-3,-3.640453D-5], Patm)
        c( 9) = log_polynomial([1.108799D1, 5.677599D-2, 4.945738D-4,-2.418338D-5], Patm)
        c(10) = log_polynomial([1.144639D1, 6.234266D-2, 3.785377D-3], Patm)
        
        d( 6) = log_polynomial([8.298063D0, 8.885346D-2,-2.901675D-3,-4.450595D-4], Patm)
        d( 7) = log_polynomial([9.012095D0, 9.149373D-2,-3.140624D-3,-3.285520D-6], Patm)
        d( 8) = log_polynomial([8.301383D0, 3.547869D-2, 3.053608D-3, 1.705129D-3, 4.357310D-5], Patm)
        d( 9) = log_polynomial([9.032754D0, 1.718233D-1,-1.352010D-3,-2.482520D-4, 1.256822D-4], Patm)
        d(10) = polynomial([2.879605D4, 2.066908D3, 1.929331D2,-6.651374D1,-7.803606D0], log(Patm))
        
        a( 6) = -polynomial([3.551938D0,-3.852851D-1,-1.698205D-2, 7.712558D-4, 1.558067D-4], log(Patm))
        a( 7) = -log_polynomial([2.202713D0,-5.805578D-3,-8.393797D-3,-1.542107D-4,-2.149336D-5,-3.876960D-7], Patm)
        a( 8) = -log_polynomial([-9.551600D-1,-1.743228D-1,-2.627017D-3, 2.020135D-3, 1.148529D-4], Patm)
        a( 9) = -log_polynomial([-4.892131D-1, 3.979950D-2, 2.397782D-3,-2.138908D-4, 1.140375D-5], Patm)
        a(10) = polynomial([1.134500D0,-5.153304D-2, 6.888543D-3], log(Patm))
        
        do i = 6,10
            eta = eta + a(i) * sigmoid(T, c(i), d(i))
        end do
        eta = exp(eta)
    end function viscosity
    
    !!
    !> Fills X with the equilibrium mole fractions of each species.
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    subroutine mole_fractions(T, P, X)
        implicit none
        double precision,               intent(in)  :: T
        double precision,               intent(in)  :: P
        double precision, dimension(:), intent(out) :: X
        
        integer :: i
        
        do i = 1,NSPECIES
            X(i) = mole_fraction_id(T, P, i)
        end do
    end subroutine mole_fractions
    
    !!
    !> Returns the equilibrium mole fraction of species sp_id as a function of
    !! T and P.
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function mole_fraction_id(T, P, sp_id) result(X)
        implicit none
        double precision,    intent(in) :: T
        double precision,    intent(in) :: P
        integer, intent(in) :: sp_id
        
        if (sp_id < 1 .or. sp_id > NSPECIES) then
            write(*,*) 'EACh:mole_fraction: Bad species index!'
            X = 0.0D0
        end if
        
        X = mole_fraction_name(T, P, species_names(sp_id))
    end function mole_fraction_id
    
    !!
    !> Returns the equilibrium mole fraction of species sp_name as a function of
    !! T and P.
    !!
    !!  @param T temperature in K
    !!  @param P pressure in Pa
    !!
    double precision function mole_fraction_name(T, P, sp_name) result(X)
        implicit none
        double precision,             intent(in) :: T
        double precision,             intent(in) :: P
        character(len=*), intent(in) :: sp_name
        
        Patm = P / ONEATM
        
        select case (trim(sp_name))
        case ('N2')
            call get_mole_fraction_N2(T, Patm, X)
        case ('N2+')
            call get_mole_fraction_N2P1(T, Patm, X)
        case ('N')
            call get_mole_fraction_N(T, Patm, X)
        case ('N+')
            call get_mole_fraction_NP1(T, Patm, X)
        case ('N++')
            call get_mole_fraction_NP2(T, Patm, X)
        case ('N+++')
            call get_mole_fraction_NP3(T, Patm, X)
        case ('N++++')
            call get_mole_fraction_NP4(T, Patm, X)
        case ('O2')
            call get_mole_fraction_O2(T, Patm, X)
        case ('O2+')
            call get_mole_fraction_O2P1(T, Patm, X)
        case ('O2-')
            call get_mole_fraction_O2M1(T, Patm, X)
        case ('O')
            call get_mole_fraction_O(T, Patm, X)
        case ('O-')
            call get_mole_fraction_OM1(T, Patm, X)
        case ('O+')
            call get_mole_fraction_OP1(T, Patm, X)
        case ('O++')
            call get_mole_fraction_OP2(T, Patm, X)
        case ('O+++')
            call get_mole_fraction_OP3(T, Patm, X)
        case ('O++++')
            call get_mole_fraction_OP4(T, Patm, X)
        case ('NO')
            call get_mole_fraction_NO(T, Patm, X)
        case ('NO+')
            call get_mole_fraction_NOP1(T, Patm, X)
        case ('e-')
            call get_mole_fraction_e(T, Patm, X)
        case default
            write(*,*) 'EACh:mole_fraction: Bad species name!'
            X = 0.0D0
        end select
    end function mole_fraction_name
    
    !!
    !! X_N2 (Table 2)
    !!
    subroutine get_mole_fraction_N2(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(0) = 0.8D0
        c(1) = log_polynomial([8.110148D0,4.553237D-2,-8.193725D-4,-2.156896D-4], P)
        c(2) = log_polynomial([8.812799D0,5.665883D-2, 1.293767D-3], P)
        
        d(1) = log_polynomial([6.561427D0,1.422222D-1, 7.476314D-4,-8.715905D-4], P)
        d(2) = log_polynomial([7.016774D0,1.058804D-1, 3.292541D-3, 2.267238D-4], P)
        
        a(2) = log_polynomial([-4.037422D-1,-7.147343D-4,4.492235D-4,9.648313D-5,-1.284083D-8], P)
        a(1) = c(0) - a(2)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 2, MF_TYPE_1)
    end subroutine get_mole_fraction_N2
    
    !!
    !! X_N2+ (Table 3)
    !!
    subroutine get_mole_fraction_N2P1(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([8.884490D0,5.573065D-2,1.616982D-3,6.738352D-5], P)
        c(2) = log_polynomial([9.203463D0,7.494796D-2,2.541069D-3,7.257196D-5,6.051419D-6], P)
        c(3) = log_polynomial([9.449201D0,6.238298D-2,1.564727D-3,5.575073D-5], P)
        
        d(1) = log_polynomial([6.552069D0, 1.058201D-1,3.989777D-3, 1.801416D-4], P)
        d(2) = log_polynomial([7.294752D0,-1.099569D-3,4.040325D-3, 2.717526D-3,-5.081078D-5,-3.474609D-5], P)
        d(3) = log_polynomial([7.762006D0, 1.260807D-1,2.223845D-3,-1.231135D-4], P)
        
        a(1) =  log_polynomial([-9.746298D0,4.199007D-1,3.143417D-3, 3.882378D-4], P)
        a(3) = -log_polynomial([-9.992503D0,4.689265D-1,1.182887D-3,-1.176687D-4], P)
        a(2) = -a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_2)       
    end subroutine get_mole_fraction_N2P1
    
    !!
    !! X_N (Table 4)
    !!
    subroutine get_mole_fraction_N(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([8.812279D0,5.474146D-2,1.019131D-3], P)
        c(2) = log_polynomial([9.516473D0,6.520807D-2,1.270979D-3,-3.857140D-5,-5.540006D-6], P)
        c(3) = log_polynomial([9.864059D0,3.659617D-2,7.907898D-3], P)
        cm   = log_polynomial([8.405373D0,4.371184D-2,1.893389D-3, 1.927737D-4], P)
        
        d(1) = log_polynomial([7.051737D0,1.128378D-1,2.407727D-3,-1.247502D-5], P)
        d(2) = log_polynomial([7.949412D0,1.206502D-1,1.785666D-3,-3.344976D-5], P)
        d(3) = log_polynomial([8.814892D0,5.421480D-2,1.537056D-3], P)
        dm   = log_polynomial([6.923056D0,1.987863D-1,3.645361D-3,-5.777817D-4], P)
        
        a(1) = polynomial([8.188731D-1,2.581889D-3,1.395103D-4], log(P))
        a(3) = -log_polynomial([-3.552214D0,4.085111D-1,-2.961084D-2], P)
        a(2) = -a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_3)
    end subroutine get_mole_fraction_N
    
    !!
    !! X_N+ (Table 5)
    !!
    subroutine get_mole_fraction_NP1(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([9.494309D0,5.588021D-2, 2.479295D-3, 5.228102D-4, 5.047984D-5,-1.606423D-6,-8.671283D-7,-5.919943D-8], P)
        c(2) = log_polynomial([9.511003D0,8.651691D-2,-5.130145D-5,-2.847046D-4], P)
        c(3) = log_polynomial([1.037880D1,6.497169D-2, 3.027193D-3, 1.559114D-4,-2.230902D-7, 3.440378D-6], P)
        c(4) = log_polynomial([1.025494D1,6.494578D-2, 1.277401D-3], P)
        cm   = log_polynomial([9.276675D0,8.451791D-2,-7.509912D-3, 1.762683D-3,-2.856325D-4, 3.392616D-5,-5.010435D-6, 3.875277D-7], P)
        
        d(1) = log_polynomial([8.228344D0, 2.288911D-1,-7.989931D-4,-1.145501D-3], P)
        d(2) = log_polynomial([7.645166D0, 8.574186D-2, 2.708947D-4, 6.210369D-4], P)
        d(3) = log_polynomial([8.810742D0, 1.305064D-1,-1.083168D-3, 4.025862D-5, 1.348428D-4,-2.273123D-5], P)
        d(4) = log_polynomial([8.187912D0, 1.182600D-1, 6.307194D-3, 2.948945D-4, 1.136590D-6], P)
        dm   = log_polynomial([7.931270D0,-4.388112D-2, 2.643605D-2,-1.501361D-3,-2.178943D-4, 2.476492D-5], P)
        
        a(1) =  log_polynomial([-1.211184D0, 2.634222D-4, 2.560470D-3], P)
        a(2) =  log_polynomial([-2.230927D0, 2.047906D-2,-2.220684D-3], P)
        a(4) = -log_polynomial([-1.200592D0,-3.074481D-2, 4.780888D-3,8.341989D-4,6.160353D-6,-2.708386D-6], P)
        a(3) = -a(4) - a(2) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 4, MF_TYPE_3)
    end subroutine get_mole_fraction_NP1
    
    !!
    !! X_N++ (Table 6)
    !!
    subroutine get_mole_fraction_NP2(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([1.018105D1,6.182886D-2,4.542717D-4, 1.665348D-4,-1.688929D-5], P)
        c(2) = log_polynomial([1.071778D1,6.267958D-2,1.384143D-3, 1.803319D-5], P)
        c(3) = log_polynomial([1.081164D1,6.929471D-2,3.005312D-3, 5.422861D-5], P)
        cm   = log_polynomial([1.020635D1,6.787015D-2,2.930559D-3,-2.387278D-5,-1.580874D-5], P)
        
        d(1) = log_polynomial([8.328213D0,7.134338D-2, 8.440573D-3,-1.913632D-4], P)
        d(2) = log_polynomial([8.558877D0,1.280075D-1, 7.408166D-3,-6.068102D-5,-3.499092D-5], P)
        d(3) = log_polynomial([9.008121D0,1.059058D-1, 3.835047D-3,-5.778232D-4], P)
        dm   = log_polynomial([8.637046D0,1.730873D-1,-2.312739D-3,-1.253255D-4, 6.870714D-5], P)
        
        a(1) =  log_polynomial([-1.320561D0,4.613513D-3, 1.563146D-3,9.805924D-5], P)
        a(3) = -log_polynomial([-2.441955D0,1.600937D-2,-1.796504D-2,4.445771D-5], P)
        a(2) = -a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_3)
    end subroutine get_mole_fraction_NP2
    
    !!
    !! X_N+++ (Table 7)
    !!
    subroutine get_mole_fraction_NP3(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([1.070665D1,6.722548D-2,6.769799D-5, 4.111595D-5], P)
        c(2) = log_polynomial([1.105085D1,5.890335D-2,1.918852D-3, 9.521033D-5], P)
        cm   = log_polynomial([1.066404D1,5.711793D-2,1.063676D-3,-1.137507D-6], P)
        
        d(1) = log_polynomial([9.340050D0,5.929963D-2, 1.505109D-3, 2.034159D-4], P)
        d(2) = log_polynomial([9.258763D0,1.273121D-1,-6.021997D-4,-2.540618D-4], P)
        dm   = log_polynomial([8.726521D0,1.521811D-1, 2.430293D-3,-4.716643D-4], P)
        
        a(1) = log_polynomial([-1.339800D0,1.954622D-2,-3.939015D-3,-4.170049D-4], P)
        a(2) = -a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 2, MF_TYPE_3)        
    end subroutine get_mole_fraction_NP3
    
    !!
    !! X_N++++ (Table 8)
    !!
    subroutine get_mole_fraction_NP4(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([1.100960D1, 7.368509D-2, 1.075589D-3], P)
        c(2) = log_polynomial([1.206372D1,-1.734608D-3,-1.447988D-2, 1.590266D-3], P)
        c(3) = log_polynomial([1.280436D1,-1.896326D-1, 2.801196D-2], P)
        cm   = log_polynomial([1.100986D1, 4.882927D-2, 3.853047D-4,-1.475967D-6], P)
        
        d(1) = log_polynomial([9.329126D0, 7.704932D-2, 2.666225D-3], P)
        d(2) = log_polynomial([1.019997D1,-1.423777D-1,-4.095877D-2, 2.180861D-3, 2.368183D-4], P)
        d(3) = log_polynomial([1.103058D1,-2.553162D-1, 2.330651D-2], P)
        dm   = log_polynomial([9.006971D0, 1.074664D-1,-1.472426D-3,-2.722012D-4], P)
        
        a(1) =  log_polynomial([-1.849635D0,-4.491118D-3,-3.702617D-4], P)
        a(3) = -log_polynomial([-6.074622D-1,6.073274D-1, 9.963043D-2, 5.415504D-3], P)
        a(2) = -a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_3)    
    end subroutine get_mole_fraction_NP4
    
    !!
    !! X_O2 (Table 9)
    !!
    subroutine get_mole_fraction_O2(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(0)   = 0.2D0
        c(1) = log_polynomial([7.851965D0,-4.971670D-2,-1.438515D-2,-8.639710D-4], P)
        c(2) = log_polynomial([8.148167D0, 4.575379D-2, 1.841872D-4], P)
        
        d(1) = log_polynomial([6.500431D0, 7.318423D-2,-2.704126D-3,-2.824658D-4], P)
        d(2) = log_polynomial([6.459154D0, 1.486515D-1, 5.919587D-3,-3.159509D-5,-4.048213D-5], P)
        
        a(2) = log_polynomial([-1.685730D0, 3.728595D-2,-5.172240D-3, 2.021941D-4, 6.195083D-5,-5.999263D-6], P)
        a(1) = c(0) - a(2)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 2, MF_TYPE_1)    
    end subroutine get_mole_fraction_O2
    
    !!
    !! X_O2+ (Table 10)
    !!
    subroutine get_mole_fraction_O2P1(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([8.794853D0, 4.659480D-2, 5.610403D-4, 1.044006D-4,-1.835079D-5], P)
        c(2) = log_polynomial([8.991604D0, 5.142449D-2, 1.298498D-3, 4.051458D-4, 1.170299D-5], P)
        c(3) = log_polynomial([9.563817D0, 7.340431D-2, 7.915772D-4,-1.592330D-4,-1.027704D-5], P)
        c(4) = log_polynomial([8.900254D0, 3.563862D-2, 1.399785D-3, 1.003372D-4,-3.618984D-5], P)
        
        d(1) = log_polynomial([7.268996D0, 9.440745D-2,-2.146537D-3, 4.167152D-5, 3.077941D-5], P)
        d(2) = log_polynomial([7.456563D0, 1.277214D-1, 8.479742D-3, 8.341173D-4,-1.597360D-4], P)
        d(3) = log_polynomial([7.834428D0, 1.245447D-1, 4.949361D-3, 3.875066D-5,-2.966365D-5], P)
        d(4) = log_polynomial([7.450971D0, 9.288765D-2,-1.491663D-3, 7.510663D-4,-9.458429D-5], P)
        
        a(1) =  log_polynomial([-1.373444D1, 6.627381D-1,-1.950471D-2, 7.469315D-4, 1.358278D-4], P)
        a(2) =  log_polynomial([-1.419853D1, 4.889623D-1,-6.123742D-3, 5.940913D-4, 9.783232D-5], P)
        a(4) = -log_polynomial([-1.342851D1, 6.025406D-1,-1.482459D-2, 1.461126D-4, 1.408990D-4], P)
        a(3) = -a(4) - a(2) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 4, MF_TYPE_2)    
    end subroutine get_mole_fraction_O2P1
    
    !!
    !! X_O2- (Table 11)
    !!
    subroutine get_mole_fraction_O2M1(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([8.151744D0, 5.269669D-2, 1.328087D-3, 9.918314D-5, 6.931618D-6], P)
        c(2) = log_polynomial([8.327753D0, 6.884887D-2, 2.843931D-3, 1.083879D-4], P)
        c(3) = log_polynomial([8.601320D0, 7.342289D-2,-9.411900D-4,-1.339663D-4, 3.126379D-5], P)
        c(4) = log_polynomial([9.428115D0, 5.014640D-2,-3.340382D-4, 7.998702D-6], P)
        
        d(1) = log_polynomial([6.140093D0, 1.051897D-1, 2.939827D-3, 1.422812D-4], P)
        d(2) = log_polynomial([6.644117D0, 1.374513D-1, 4.095263D-3, 8.402722D-5,-1.242256D-5,-7.990825D-6, 8.075101D-7, 2.001120D-7], P)
        d(3) = log_polynomial([6.985630D0, 8.256947D-2, 9.999196D-3,-2.953396D-5,-1.526330D-4], P)
        d(4) = log_polynomial([7.530896D0, 1.558330D-1, 4.905502D-3,-8.411242D-4], P)
        
        a(1) =  log_polynomial([-2.009128D1, 1.218472D0,-1.023713D-2,-3.693487D-4], P)
        a(3) = -log_polynomial([-2.169571D1, 1.231117D0, 1.792651D-3, 2.558252D-4,-1.732401D-4, 8.498995D-6, 1.264359D-6], P)
        a(4) = -log_polynomial([-2.472870D1, 1.526884D0, 1.203852D-3, 1.430794D-4], P)
        a(2) = -a(4) - a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 4, MF_TYPE_2)    
    end subroutine get_mole_fraction_O2M1
    
    !!
    !! X_O (Table 12)
    !!
    subroutine get_mole_fraction_O(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([8.145639D0, 5.431612D-2, 2.023998D-3, 1.003745D-4], P)
        c(2) = log_polynomial([9.811744D0, 7.436247D-2,-1.239267D-4, 6.132060D-4], P)
        c(3) = log_polynomial([8.819734D0, 5.805213D-2, 1.501067D-3, 2.511693D-5], P)
        c(4) = log_polynomial([9.554277D0, 6.746571D-2, 8.910292D-4,-4.496226D-5], P)
        cm   = log_polynomial([7.940783D0, 6.741169D-2,-2.087042D-3, 3.972481D-4,-3.481686D-5, 1.485858D-6], P)
        
        d(1) = log_polynomial([6.576786D0, 1.491330D-1, 3.724868D-3,-1.382563D-4, 1.947915D-6, 1.082756D-6], P)
        d(2) = log_polynomial([9.044853D0, 5.997097D-3, 4.532508D-4, 6.756744D-4], P)
        d(3) = log_polynomial([6.918048D0, 9.326905D-2, 2.506390D-3, 1.395474D-4], P)
        d(4) = log_polynomial([8.033301D0, 1.233674D-1, 1.651217D-3,-3.811131D-5], P)
        dm   = log_polynomial([6.664764D0, 4.575484D-2, 4.557480D-3], P)
        
        a(1) =  log_polynomial([-1.139281D0,-1.050647D-2,-1.022007D-3,-4.830320D-5,-3.531305D-6,-2.296630D-7], P)
        a(2) = -log_polynomial([-4.979500D0, 2.665257D-1, 1.458327D-2,-2.533456D-3,-3.704428D-4, 2.339924D-5], P)
        a(4) = -log_polynomial([-1.615594D0,-5.157778D-3,-1.550658D-3,-1.264223D-4, 2.343404D-5, 3.184705D-6], P)
        a(3) = -a(4) - a(2) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 4, MF_TYPE_3)    
    end subroutine get_mole_fraction_O
    
    !!
    !! X_O- (Table 13)
    !!
    subroutine get_mole_fraction_OM1(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([8.415326D0, 5.157258D-2, 2.024706D-3, 1.312425D-4, 1.315036D-5], P)
        c(2) = log_polynomial([9.270258D0, 5.316281D-2, 1.482070D-3,-5.476676D-5,-9.733849D-6], P)
        c(3) = log_polynomial([9.598507D0, 6.569448D-2, 5.303147D-4,-9.613381D-5,-7.330576D-6], P)
        
        d(1) = log_polynomial([6.462668D0, 6.272626D-2,-6.193918D-3, 6.376014D-4, 2.245471D-4], P)
        d(2) = log_polynomial([7.724023D0, 9.838486D-2, 4.215920D-3,-6.990084D-5,-3.230965D-5], P)
        d(3) = log_polynomial([7.809041D0, 1.423877D-1, 4.366188D-3,-8.184536D-5,-1.524608D-5], P)
        
        a(1) = log_polynomial([-1.492297D1, 9.064321D-1,-8.724265D-3,-2.165125D-4, 1.166368D-4], P)
        a(2) = log_polynomial([-1.175041D1, 7.618857D-1, 1.501595D-3, 1.781504D-4, 9.215991D-6], P)
        a(3) = -a(1) - a(2)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_2)    
    end subroutine get_mole_fraction_OM1
    
    !!
    !! X_O+ (Table 14)
    !!
    subroutine get_mole_fraction_OP1(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([9.588569D0, 6.997026D-2, 9.769379D-4,-6.246775D-5,-4.877947D-6], P)
        c(2) = log_polynomial([1.020364D1, 6.299762D-2,-1.091887D-3, 3.702998D-5], P)
        c(3) = log_polynomial([1.027215D1, 4.672465D-2, 1.597850D-4, 9.311678D-6], P)
        cm   = log_polynomial([9.115221D0, 6.168847D-2, 2.270858D-3, 1.412631D-4], P)
        
        d(1) = log_polynomial([8.044970D0, 1.175891D-1, 1.645336D-3,-9.489377D-5,-9.694619D-6], P)
        d(2) = log_polynomial([8.680331D0, 1.325526D-1, 2.754338D-3,-7.964755D-5], P)
        d(3) = log_polynomial([8.696369D0, 1.339624D-1, 1.995427D-3,-3.323281D-5], P)
        dm   = log_polynomial([7.651684D0, 1.477558D-1,-1.967294D-3,-9.075769D-4], P)
        
        a(1) =  log_polynomial([-2.319093D0,-7.610174D-3,-1.953269D-3,-3.002482D-4,-1.751192D-5], P)
        a(3) = -log_polynomial([-1.436906D0, 2.872384D-1, 2.978317D-2, 1.769679D-4,-9.414001D-5], P)
        a(2) = -a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_3)    
    end subroutine get_mole_fraction_OP1
    
    !!
    !! X_O++ (Table 15)
    !!
    subroutine get_mole_fraction_OP2(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([1.029003D1, 4.517420D-2,-1.618224D-5, 2.245678D-4, 3.130833D-5,-2.423868D-6,-3.903368D-7], P)
        c(2) = log_polynomial([1.082680D1, 7.388982D-2, 9.267668D-4], P)
        c(3) = log_polynomial([1.078471D1, 5.999115D-2, 1.044468D-3], P)
        cm   = log_polynomial([1.029386D1, 8.048612D-2,-4.497818D-4, 3.852087D-5], P)
        
        d(1) = log_polynomial([8.449025D0, 1.233942D-1,-3.128794D-3,-5.456369D-4, 5.445584D-5, 6.520078D-6], P)
        d(2) = log_polynomial([9.267200D0, 5.532633D-2, 2.362320D-3, 6.299569D-4, 1.122230D-5,-2.869166D-6,-4.451869D-7], P)
        d(3) = log_polynomial([8.785646D0, 9.165132D-2, 9.925663D-4], P)
        dm   = log_polynomial([8.843594D0, 4.195145D-2, 1.187095D-2, 1.964457D-4,-4.989937D-5,-9.711143D-7], P)
        
        a(1) =  polynomial([7.063013D-2,-5.187789D-4,-9.288238D-6], log(P))
        a(3) = -log_polynomial([-2.991458D0,-5.757422D-2,-3.835760D-3], P)
        a(2) = -a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_3)    
    end subroutine get_mole_fraction_OP2
    
    !!
    !! X_O+++ (Table 16)
    !!
    subroutine get_mole_fraction_OP3(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([1.074207D1, 5.260080D-2, 4.936255D-4,-4.405321D-5,-3.025027D-6,-5.425422D-7], P)
        c(2) = log_polynomial([1.111558D1, 5.973321D-2, 2.038965D-3, 9.054082D-5], P)
        cm   = log_polynomial([1.077506D1, 6.587529D-2, 2.491665D-4, 1.077355D-4], P)
        
        d(1) = log_polynomial([8.835975D0, 1.411710D-1, 2.773994D-3,-6.211959D-4, 3.813517D-6, 1.323357D-5,-1.119305D-6,-3.062376D-7], P)
        d(2) = log_polynomial([9.317898D0, 1.146590D-1,-4.219919D-4,-1.986513D-4,-9.501572D-6], P)
        dm   = log_polynomial([9.367809D0, 3.868631D-2,-7.976461D-4, 6.108727D-4], P)
        
        a(1) = log_polynomial([-2.760009D0, 3.495500D-3,-5.357864D-3,-2.144466D-4, 9.251860D-6,-9.005345D-7], P)
        a(2) = -a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 2, MF_TYPE_3)    
    end subroutine get_mole_fraction_OP3
    
    !!
    !! X_O++++ (Table 17)
    !!
    subroutine get_mole_fraction_OP4(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([1.114079D1, 6.128099D-2, 1.305781D-3,-4.745385D-5,-1.294845D-5,-6.416314D-7], P)
        c(2) = log_polynomial([1.133963D1, 5.445065D-2,-3.976441D-4, 1.251159D-4], P)
        c(3) = log_polynomial([1.473199D1,-3.158041D-1,-3.070674D-2, 6.776443D-3], P)
        cm   = log_polynomial([1.097387D1, 5.385207D-2, 3.454294D-5,-9.334055D-5], P)
        
        d(1) = log_polynomial([9.124558D0, 1.015232D-1,-1.452067D-3,-4.363441D-4,-9.737843D-6, 1.643326D-6], P)
        d(2) = log_polynomial([9.165912D0, 3.362575D-2, 1.118630D-3,-3.084012D-4,-7.665827D-5], P)
        d(3) = log_polynomial([1.306288D1,-3.228563D-1,-3.275522D-2, 6.750116D-3], P)
        dm   = log_polynomial([9.008289D0, 5.266326D-2,-2.558320D-4, 3.532844D-5], P)
        
        a(1) =  log_polynomial([-3.273424D0,-9.222532D-3,-2.546540D-3,-6.142466D-4,-6.803461D-5,-1.480622D-6], P)
        a(3) = -log_polynomial([-3.227410D0,-4.108171D-3,-6.841752D-4,-3.928651D-5], P)
        a(2) = -a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_3)    
    end subroutine get_mole_fraction_OP4
    
    !!
    !! X_NO (Table 18)
    !!
    subroutine get_mole_fraction_NO(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([7.942600D0, 2.917164D-2, 6.775381D-4, 2.209082D-5], P)
        c(2) = log_polynomial([8.274503D0, 6.655621D-2, 2.214534D-3, 3.856329D-5], P)
        c(3) = log_polynomial([8.364477D0, 7.365241D-2, 2.771836D-3,-5.013391D-6,-5.293616D-6], P)
        
        d(1) = log_polynomial([6.780323D0, 6.029139D-2, 4.276063D-4], P)
        d(2) = log_polynomial([6.495225D0, 7.930874D-2,-1.952605D-3,-7.384374D-4,-5.231985D-5], P)
        d(3) = log_polynomial([7.549495D0, 9.399569D-2], P)
        
        a(1) =  log_polynomial([-2.397641D0, 9.644207D-2], P)
        a(3) = -log_polynomial([-2.923272D0, 1.671984D-1], P)
        a(2) = -a(3) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 3, MF_TYPE_2)    
    end subroutine get_mole_fraction_NO
    
    !!
    !! X_NO+ (Table 19)
    !!
    subroutine get_mole_fraction_NOP1(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([8.740893D0, 4.144123D-2, 3.456197D-4], P)
        c(2) = log_polynomial([8.817743D0, 4.865084D-2, 7.358462D-4], P)
        c(3) = log_polynomial([8.899564D0, 6.228872D-2, 1.910295D-3, 5.292903D-5], P)
        c(4) = log_polynomial([9.221935D0, 8.005371D-2, 3.728793D-3,-1.235847D-4,-6.058282D-6], P)
        
        d(1) = log_polynomial([6.996599D0, 6.789593D-2, 1.320085D-3, 2.143434D-5, 6.597691D-6, 1.625852D-7], P)
        d(2) = log_polynomial([6.260938D0, 9.417073D-2, 7.841151D-3], P)
        d(3) = log_polynomial([7.246371D0, 1.012940D-1, 4.389279D-3,-2.344414D-5,-1.533963D-5], P)
        d(4) = log_polynomial([7.940040D0, 1.021609D-1, 5.411563D-3, 1.592304D-4,-7.583651D-5], P)
        
        a(1) =  log_polynomial([-7.135266D0, 4.617651D-2,-7.097386D-4], P)
        a(2) =  log_polynomial([-8.753925D0, 1.392942D-1, 1.317873D-2], P)
        a(4) = -log_polynomial([-8.772596D0, 1.382471D-1,-1.513718D-3, 1.822779D-3, 5.774867D-5], P)
        a(3) = -a(4) - a(2) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 4, MF_TYPE_2)    
    end subroutine get_mole_fraction_NOP1
    
    !!
    !! X_e- (Table 20)
    !!
    subroutine get_mole_fraction_e(T, P, X)
        implicit none
        double precision, intent(in)  :: T
        double precision, intent(in)  :: P
        double precision, intent(out) :: X
        
        ! Evaluate logarithmic polynomials for a, c, and d
        c(1) = log_polynomial([9.514823D0, 6.426626D-2, 3.538392D-4,-1.093881D-4], P)
        c(2) = log_polynomial([1.025313D1, 6.613035D-2, 2.106960D-3, 1.249059D-4,-3.254728D-6,-1.073094D-6,-4.149968D-7,-4.918145D-8], P)
        c(3) = log_polynomial([1.074247D1, 6.026184D-2, 6.834881D-4, 1.412968D-6], P)
        c(4) = log_polynomial([1.106632D1, 5.734452D-2, 1.326880D-3, 4.870977D-5], P)
        c(5) = log_polynomial([1.009244D1, 5.691765D-2, 2.642057D-3, 3.297719D-5], P)
        c(6) = log_polynomial([1.219498D2,-3.565001D0, 7.046916D-1, 3.062083D-1,-2.940975D-2], P)
        cm   = log_polynomial([6.343867D0, 1.473478D0,-2.628976D-1, 2.653667D-2,-1.170989D-3], P)
        
        d(1) = log_polynomial([7.931006D0, 1.174176D-1,-5.369256D-4,-1.640676D-4, 2.876393D-5], P)
        d(2) = log_polynomial([8.461864D0, 1.033435D-1,-6.800325D-3,-2.171111D-3, 8.042855D-5, 3.126866D-5, 3.548083D-6, 1.732832D-7], P)
        d(3) = log_polynomial([8.457103D0, 1.570495D-1, 2.577271D-2,-4.699755D-4,-7.340190D-4,-1.521958D-6, 7.337098D-6,-1.937258D-8], P)
        d(4) = log_polynomial([9.134358D0, 1.817063D-1, 8.463508D-3], P)
        d(5) = log_polynomial([9.041428D0, 9.809302D-2, 1.899235D-3,-1.329754D-4,-2.357106D-5], P)
        d(6) = log_polynomial([1.163952D2,-3.232407D0, 6.981116D-1, 2.997466D-1,-2.749064D-2], P)
        dm   = log_polynomial([1.029159D1, 3.502366D-2,-1.043994D-2,-7.498040D-4, 1.464646D-4, 1.031691D-5,-3.878009D-7], P)
        
        a(1) =  log_polynomial([-3.932487D-1, 7.116035D-4, 4.083493D-4, 3.307562D-4, 2.215248D-5,-4.020145D-6], P)
        a(2) =  log_polynomial([-1.599518D0,-3.681243D-2,-1.499672D-2,-4.875898D-3,-9.278204D-5, 8.792226D-5, 1.273088D-5], P)
        a(3) =  log_polynomial([-3.031217D0,-1.236964D-2, 4.999807D-3, 4.130827D-4,-5.879976D-5,-5.643378D-6,-2.118479D-7,-8.835667D-8], P)
        a(4) =  log_polynomial([-3.096101D0, 5.690833D-2, 1.063005D-2, 8.066239D-4], P)
        a(6) = -log_polynomial([-3.465436D-1,-2.831472D-3,-1.021467D-3,-7.753035D-5], P)
        a(5) = -a(6) - a(2) - a(1)
        
        ! Mole fraction
        X = evaluate_mole_fraction(T, 6, MF_TYPE_3)    
    end subroutine get_mole_fraction_e
    
    !!
    !> @brief Computes the mole fraction of a species at a given temperature T 
    !!        using precomputed sigmoid coefficients as functions of pressure.  
    !!
    !! Three mole fraction functions are defined as follows:
    !!
    !!    MF_TYPE_1:  
    !!
    !!        X(T) = c0 - sum_{j=1}^n a_j*s_j(T)
    !!
    !!    MF_TYPE_2:  
    !!
    !!        X(T) = max[ 0.0, sum_{j=1}^n a_j*s_j(T) ]
    !!
    !!    MF_TYPE_3:
    !!
    !!        X(T) = max[ 0.0, a_1*s_m(T)*s_1(T) + sum_{j=2}^n a_j*s_j(T) ] 
    !!
    !! where s(T) is the sigmoid function that takes c and d as parameters.
    !! Note that all of the necessary a, c, and d's must be computed prior to 
    !! calling this function.
    !!
    !! @param T    temperature in K
    !! @param n    summation limit as defined in the mole fraction functions
    !! @param func function type (MF_TYPE_1, MF_TYPE_2, or MF_TYPE_3)
    !!
    double precision function evaluate_mole_fraction(T, n, func) result(X)
        implicit none
        double precision,    intent(in) :: T
        integer, intent(in) :: n
        integer, intent(in) :: func
        
        integer :: i
        
        select case (func)
        case (MF_TYPE_1)
            X = c(0)
            do i = 1,n
                X = X - a(i) * sigmoid(T, c(i), d(i))
            end do
        case (MF_TYPE_2)
            X = 0.0D0
            do i = 1,n
                X = X + a(i) * sigmoid(T, c(i), d(i))
            end do
            X = max(0.0D0, X)
        case (MF_TYPE_3)
            X = a(1) * sigmoid(T, cm, dm) * sigmoid(T, c(1), d(1))
            do i = 2,n
                X = X + a(i) * sigmoid(T, c(i), d(i))
            end do
            X = max(0.0D0, X)
        end select
    end function evaluate_mole_fraction
    
    !!
    !> Evaluates the sigmoid function s(T, c, d) where
    !!
    !!    $s = exp(q) / (exp(q) + exp(-q)), and q = (T - c) / d.
    !!
    double precision function sigmoid(T, c, d)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: c
        double precision, intent(in) :: d
        
        double precision :: q
        
        q = (T - c) / d
        sigmoid = exp(q)/(exp(q) + exp(-q))
    end function sigmoid
    
    double precision function gaussian(T, c, d)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: c
        double precision, intent(in) :: d
        
        double precision :: q
        
        q = (T - c) / d
        gaussian = exp(-q*q)
    end function gaussian
    
    double precision function xi(T, a, c, d, w)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: a
        double precision, intent(in) :: c
        double precision, intent(in) :: d
        double precision, intent(in) :: w
        
        xi = a - c * exp(-(T / d)**w)
    end function xi
    
    double precision function phi(T, a, w)
        implicit none
        double precision, intent(in) :: T
        double precision, intent(in) :: a
        double precision, intent(in) :: w
        
        phi = a * T**w
    end function phi
    
    !!
    !! Evaluates the polynomial
    !!
    !!    p(x) = a(1) + a(2)*x + a(3)*x**2 + ...
    !!
    double precision function polynomial(a, x)
        implicit none
        double precision, dimension(:), intent(in) :: a
        double precision,               intent(in) :: x
        
        integer :: i, n
        n = size(a)
        
        polynomial = a(n)
        do i = n-1,1,-1
            polynomial = a(i) + polynomial*x
        end do 
    end function polynomial
    
    !!
    !! Evaluates the logarithmic polynomial
    !!
    !!    lp(x) = exp( a(1) + a(2)*log(x) + a(3)*log(x)**2 + ... )
    !!
    double precision function log_polynomial(a, x)
        implicit none
        double precision, dimension(:), intent(in) :: a
        double precision,               intent(in) :: x
        
        log_polynomial = exp(polynomial(a, log(x)))
    end function log_polynomial
end module air_ETL
