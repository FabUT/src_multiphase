module mod_elec

use mod_data !!!, only : mesh, fluide, PR, pi, radial

use air_ETL, only : electric_conductivity 

use mod_euler, only : phase, ph

implicit none


real(PR) :: mu0=4.0e-7_PR*Pi
real(PR) :: current=0.0_PR
real(PR) :: U=0.0_PR
real(PR) :: R=0.0_PR



contains 


subroutine SOLVE_elec(M,time)

   implicit none
   
   type(mesh), intent(inout) :: M
   real(PR), intent(in) :: time
   real(PR) :: sigl,sig,G,P,T,S,curr,Te
   !real(PR) :: Jl(1:M%Nl)
   integer :: i,j,k,iph
   
   if(time.lt.-20.0e-9_PR)then 

      do i=1,M%Nx
         M%MF(i,1,1)%Qe=1.5e15_PR*exp(-((M%x(i)-0.0e-2_PR)/2.0e-3_PR)**4)
         do iph=1,M%Nl
           !!!---Terme source initial
           M%MF(i,1,1)%F(iph)%Qe=M%MF(i,1,1)%Qe
         enddo
      enddo 

   else

      !!!---Effet Joule

      current=biexp(time)

      !!!---calcul conductance

      G=0.0_PR

      do i=1,M%Nx
 
        sig=0.0_PR
        T=min(M%MF(i,1,1)%T,70000.0_PR)

        do iph=1,M%Nl
           P=min(M%MF(i,1,1)%F(iph)%p,1.0e8_PR)
           !sigl=electric_conductivity(T,P)

           if(ph(iph)%typ.eq.1)then
             sigl=ph(iph)%sig
             if(M%MF(i,1,1)%F(iph)%T.gt.ph(iph)%Tvap)then
             sigl=1.0e4_PR
             endif

           elseif(ph(iph)%typ.eq.2)then
             Te=M%MF(i,1,1)%F(iph)%T
             sigl=1.0e7_PR
             !call transport_parameters(im,M%F(i,1,1)%rhl(iph),M%F(i,1,1)%Tl(iph),Te,sigl,lambda,Sth)                 
           endif


           M%MF(i,1,1)%F(iph)%sigma=sigl

           !!!---modèle avec limite de percolation
           !if(M%F(i,1,1)%fl(1).ge.0.5_PR)then
           !   !!!---distribution parallèle
           !   sig=sig+M%F(i,1,1)%fl(iph)*sigl
           !else
           !   !!!---distribution série
           !   if(iph.eq.1)then
           !     sig=sigl/M%F(i,1,1)%fl(iph)
           !   else
           !     sig=1.0_PR/(1.0_PR/sig+M%F(i,1,1)%fl(iph)/sigl)
           !   endif
           !endif

           !!!---modèle parallèle
           sig=sig+M%MF(i,1,1)%F(iph)%f*sigl

        enddo



        M%MF(i,1,1)%sigma=sig

        if(radial)then
           S=Pi*M%xm(i)**2        
           if(i.gt.1) S=S-Pi*M%xm(i-1)**2
        else
           S=M%dx(i)*M%dy(1)
        endif
 

        G=G+sig*S

      enddo

      R=1.0_PR/G

      !!!---calcul tension
      U=R*current


      !!!---Calcul densité de l'effet Joule
      do i=1,M%Nx

         M%MF(i,1,1)%J=0.0_PR

         !!!---Jl différent dans chaque phase

         !!!---distribution parallèle
            do iph=1,M%Nl
              M%MF(i,1,1)%F(iph)%J=M%MF(i,1,1)%F(iph)%sigma*U/M%dz(1)
            enddo


         !!!---Modèle avec percolation    
         !if(M%F(i,1,1)%fl(1).ge.0.5_PR)then
         !   do iph=1,M%Nl
         !     M%F(i,1,1)%Jl(iph)=M%F(i,1,1)%sigmal(iph)*U/M%dz(1)
         !   enddo
         !else
         !   do iph=1,M%Nl
         !     M%F(i,1,1)%Jl(iph)=0.0_PR !!!!   M%F(i,1,1)%sigma*U/M%dz(1)
         !   enddo
         !endif

         M%MF(i,1,1)%cv=0.0_PR
         do iph=1,M%Nl
            M%MF(i,1,1)%cv=M%MF(i,1,1)%cv+M%MF(i,1,1)%F(iph)%Y*M%MF(i,1,1)%F(iph)%cv
         enddo
 
         !!!---V1 : effet Joule différent dans chaque phase
         !do iph=1,M%Nl
         !  M%F(i,1,1)%Qel(iph)=M%F(i,1,1)%Jl(iph)**2/M%F(i,1,1)%sigmal(iph)
         !enddo
         !M%F(i,1,1)%J= sum(M%F(i,1,1)%fl(1:M%Nl)*M%F(i,1,1)%Jl(1:M%Nl))
         !M%F(i,1,1)%Qe=sum(M%F(i,1,1)%fl(1:M%Nl)*M%F(i,1,1)%Qel(1:M%Nl))

         !!!---V2 : Chaufage uniforme
         !M%F(i,1,1)%J=sum(M%F(i,1,1)%fl(1:M%Nl)*M%F(i,1,1)%Jl(1:M%Nl))
         M%MF(i,1,1)%J=M%MF(i,1,1)%sigma*U/M%dz(1)
         M%MF(i,1,1)%Qe=M%MF(i,1,1)%J**2/M%MF(i,1,1)%sigma
        do iph=1,M%Nl
           M%MF(i,1,1)%F(iph)%Qe=M%MF(i,1,1)%Qe*M%MF(i,1,1)%F(iph)%rh*M%MF(i,1,1)%F(iph)%cv/(M%MF(i,1,1)%rh*M%MF(i,1,1)%cv)
        enddo

      enddo 

      !!!---Calcul du champ magnetique et des forces de Laplace

      curr=0.0_PR

      do i=1,M%Nx

         if(radial)then

            S=Pi*M%xm(i)**2        
            if(i.gt.1) S=S-Pi*M%xm(i-1)**2

            curr=curr+M%MF(i,1,1)%J*S

            M%MF(i,1,1)%B=0.5_PR*mu0*curr/(Pi*M%x(i))

         else

            S=M%dx(i)*M%dy(1)
           
            curr=curr+M%MF(i,1,1)%J*S

            M%MF(i,1,1)%B=0.5_PR*mu0*curr/M%dy(1)          
 
         endif



         !M%F(i,1,1)%Qpx=-M%F(i,1,1)%J*M%F(i,1,1)%B

      enddo 

      !if(.not.radial)then

      !     Bmax=maxval(M%F(:,1,1)%B)
      !     M%F(:,1,1)%B=M%F(:,1,1)%B-0.5_PR*Bmax

      !endif


   endif

end subroutine solve_elec



real(PR) function biexp(t)
   implicit none
   real(PR), intent(in) :: t
   real(PR) :: A, alpha, beta, coef


   coef=1.0_PR
   
   ! D:  
   A=109405.0_PR ; alpha=22708.0_PR ; beta=1294530.0_PR 
   
   ! A:
   !A=218810.0_PR ; alpha=11354.0_PR ; beta=647265.0_PR

   ! fit cas test plaque_trou
   !beta=20000.000000000000000      
   !alpha=215233.43577297761624      
   !A=-998.43814979270929472   

   !!!---Grifon 100 kA 
   !A=-1044746.2664139411842
   !alpha=81169.502440601603134
   !beta=1.d0/(0.08d0*200.d-6)
   

   biexp = coef*A*( exp(-alpha*t) - exp(-beta*t) ) 
   
end function biexp


!!!!=============================  CONDUCTIVITY MODEL  ====================================
!!!!
!!!! Lee More 1984 + Desjarlais 2001
!!!!
!!!!=======================================================================================
!
!
!subroutine transport_properties(Z,rho,T,Zbar,sig,K,S,&
!                                sig1,sig2,K1,K2,na,ne,nn,ni,&
!                                lnl,tau,tau1,tau2,tau_en,cross_en,mu)
!
!    use mod_Fermi_Dirac, only : fd1h, fdm1h
!
!   implicit none
!   integer, intent(in) :: Z
!   real(PR), intent(in) :: rho, T
!   real(PR), intent(in) :: Zbar
!   real(PR), intent(out) :: sig,K,S,sig1,sig2,K1,K2,lnl,tau
!   real(PR), intent(out) :: na, nn, ni, ne !!! density of atoms, neutrals, ions, electrons
!   real(PR) :: F12, x, ln_Lambda, Te, Ti, Zi
!   real(PR) :: Tmelt, R0, bmin, bmax, lambda_dense
!   real(PR), intent(out) :: tau1, tau2, tau_en, cross_en
!   real(PR) :: VF, vth, TF, EF, mu
!   real(PR) :: lb, ldbc, ldbg, lcca
!   !!!---for Redmer e-n model:
!   real(PR) :: alphaD,lambda_e,Fm12,kappa,kwave_e,kar,kr0,Ak,Bk,Ck,Dk,Ek
!
!   Te=T ; Ti=T
!
!   na=rho/A(Z)*Avogadro*1.0e6_PR !!! m^-3
!
!   ne=na*Zbar
!
!   if(Zbar.ge.1.0_PR)then
!     ni=na
!     nn=0.0_PR
!     Zi=Zbar
!   else
!     ni=ne
!     nn=na-ni
!     Zi=1.0_PR
!   endif
!
!   EF=hb**2/(2.0_PR*meSI)*(3.0_PR*Pi**2*ne)**(2.0_PR/3.0_PR) !!!J
!
!   !mu=EF
!
!   !print*, 'ne=', ne, 'Zbar=', Zbar, 'na=', na, 'Te=', Te, 'EF=', EF
!
!   mu=chemical_potential(EF,Te)
!
!   !print*, 'mu=', mu, 'Zbar=', Zbar
!
!   TF=EF/kb
!
!   x=mu/(kb*Te)
! 
!   F12=fd1h(x)
!
!   !bmax=sqrt(1.0_PR/(&
!   !         4.0_PR*Pi*ne*qe**2/(kb*sqrt(Te**2+TF**2))&
!   !        +4.0_PR*Pi*ni*Zbar**2*qe**2/(kb*Ti)&
!   !        ))
!
!   !bmax=sqrt(1.0_PR/(&
!   !         ne*qe**2/(eps0*kb*sqrt(Te**2+TF**2))&
!   !        +ni*Zbar**2*qe**2/(eps0*kb*Ti)&
!   !        ))
!
!   bmax=sqrt(1.0_PR/(&
!            ne*qe**2/(eps0*kb*sqrt(Te**2+TF**2))&
!           +ni*Zi**2*qe**2/(eps0*kb*Ti)&
!           ))
!
!
!   !bmax=sqrt((eps0*kb*Te)/(ne*qe**2))
!
!
!   !!!---electron velocity 
!   vF=sqrt(2.0_PR*EF/meSI)
!   !vF=v_Fermi(Z)
!
!   vth=max(sqrt(3.0_PR*kB*Te/meSI),vF)
!
!   !vth=max(sqrt(3.0_PR*kB*Te/meSI),0.0_PR)
!
!
!   !bmin=max( Zbar*qe**2/( 4.0_PR*Pi*eps0*meSI*vth**2 ), h/(2.0_PR*meSI*vth))
!
!   bmin=max( Zi*qe**2/( 4.0_PR*Pi*eps0*meSI*vth**2 ), h/(2.0_PR*meSI*vth))
!
!   !bmin=sqrt( (Zbar*qe**2/( 4.0_PR*Pi*eps0*3.0_PR*kb*Te))**2 + (hb/(2.0_PR*sqrt(3.0_PR*me*kb*Te)))**2 )
!
!
!   !!!---papier de Lee More 1984
!   ln_lambda=0.5_PR*log(1.0_PR+bmax**2/bmin**2)
!
!   !!!---Hayes 2016 (sur l'implementation de 84)
!   !ln_lambda=0.5_PR*log(exp(1.0_PR)+bmax**2/bmin**2)
!
!   !!!---Hayes 2016 (sur l'implementation de 84)
!   !lb=exp(1.0_PR)
!   !ldbc=sqrt((eps0*kb*Te)/(ne*qe**2))
!   !ldbg=hb/(2.0_PR*sqrt(3.0_PR*meSI*kb*Te))
!   !Lcca=1.0_PR/(4.0_PR*Pi*eps0)*Zbar*qe**2/(3.0_PR*kb*Te)
!
!   !ln_lambda=0.5_PR*log(lb**2+Ldbc**2/(Ldbg**2+Lcca**2))
!
!   !ln_lambda=max(0.5_PR*log(lb**2+Ldbc**2/(Ldbg**2+Lcca**2)),2.0_PR)
!
!
!   lnl=ln_lambda
!
!   !!!====Coulomb cross section for e-i with Debye-Hückel cut-off
!   !!!    and degeneracy correction
!
!  ! tau1=( (3.0_PR*sqrt(meSI)*(kb*Te)**(1.5_PR))/&
!  !        (2.0_PR*sqrt(2.0_PR)*Pi*Zbar**2*ni*qe**4*ln_Lambda) )*&
!  !        ( 1.0_PR+exp(-x) )*F12
!
!   !!!--version 0 de Lee et More
!   !tau1=( (3.0_PR*(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR))/&
!   !       (2.0_PR*sqrt(2.0_PR)*Pi*Zi**2*ni*qe**4*max(ln_Lambda,2.0_PR)) )*&
!   !       ( 1.0_PR + exp(-x) )*F12
!
!   !!!--version de Nanagan 2015 avec correction pour e-e scattering
!   tau1=( (3.0_PR*(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR))/&
!          (2.0_PR*sqrt(2.0_PR)*Pi*Zi**2*ni*qe**4*max(ln_Lambda,2.0_PR)) )*&
!          Fc_alpha(Zi)*( 1.0_PR + exp(-x) )*F12
!
!  
!
!
!   !tau1=( (3.0_PR*(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR))/&
!   !       (2.0_PR*sqrt(2.0_PR)*Pi*Zbar**2*na*qe**4*max(ln_Lambda,2.0_PR)) )*&
!   !       ( 1.0_PR + exp(-x) )*F12
!
!
!   !tau=(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR)/&
!   !    (Pi*Zbar**2*ni*qe**4*ln_lambda)
!
!   !!!===Sphere dure ?
!
!   if(nn.gt.0.0_PR)then
!
!      !!!---Lee More (2e-15 cm^2 from desjarlais 2001)--------------------
!      !cross_en=2.0e-19_PR
!      !!!---Max value from Desjarlais (30e-15 cm^2 from desjarlais 2001)--------------------
!      !cross_en=30.0e-19_PR
!      !!! Desjarlais 2001: model from Redmer 1999 -------------------
!      alphaD=pola(Z)
!      lambda_e=sqrt(2.0_PR*Pi*hb**2/(kb*Te*meSI)) !!! thermal wavelength (cf eq 4)
!      Fm12=fdm1h(x)     
!      !!! eq 14 in Redmer1999
!      kappa=sqrt(qe**2/(kb*Te*eps0)*2.0_PR/lambda_e**3*Fm12)     
!      kwave_e=meSI*vth/hb !!! electron wave number Desjarlais 2001
!      r0=(alphaD*RBohr/(2.0_PR*real(Z,PR)**(1.0_PR/3.0_PR)))**0.25_PR
!      kar=kappa*r0
!      kr0=kwave_e*r0
!
!      Ak=1.0_PR+2.0_PR*kar+7.0_PR/pi**2*kar**2+pi/7.0_PR*kar**3
!      Bk=exp(-18.0_PR*kar)
!      Ck=(1.0_PR+22.0_PR*kar-11.3_PR*kar**2+33.0_PR*kar**4)/&
!         (1.0_PR+6.0_PR*kar+4.7_PR*kar**2+2.0_PR*kar**4)
!      Dk=(1.0_PR+28.0_PR*kar+13.8_PR*kar**2+3.2_PR*kar**3)/&
!         (1.0_PR+8.0_PR*kar+10.0_PR*kar**2+kar**3)
!      Ek=1.0_PR+0.1_PR*kar+0.3665_PR*kar**2
!
!      cross_en=pi**3*(alphaD/(2.0_PR*r0*RBohr))**2/&
!      (Ak**2+3.0_PR*Bk*kr0+7.5_PR*Ck*kr0**2-3.4_PR*Dk*kr0**3+10.6668_PR*Ek*kr0**4)
!
!      !print*,  'cross_en=', T, cross_en, 'r0**2', r0**2
!      !!!--------------------------------------------------------------
!         
!      tau_en=1.0_PR/(nn*vth*cross_en)
!      tau=1.0_PR/(1.0_PR/tau1+1.0_PR/tau_en)
!      !tau =tau1    
! 
!   else
!
!      tau=tau1
!
!   endif
!
!
!   !!!====Bloch-Grüneisen - Zimann formula for solid-liquid conductivities
!  
!   Tmelt=Tm(Z,rho)
!   R0=(1.0_PR/na)**(1.0/3.0_PR)
!
!   if(T.le.Tmelt)then
!     lambda_dense=50.0_PR*R0*(Tmelt/T)
!   else
!     lambda_dense=50.0_PR*R0*(Tmelt/T)/sigSL(Z)
!   endif
!
!   !tau2=max(lambda_dense/vF,R0/vF)
!   tau2=lambda_dense/vF
!
!   !taumin=min(p2a+p2b/
!   !tau=max(tau,taumin)
!
!
!   !!!===Final tau
!
!   !tau=tau1
!   !write(6,"(4(ES14.7,' '))"), ne, tau,A_alpha(x), x
! 
!   sig1=ne*qe**2*tau/meSI*A_alpha(x)
!   sig2=ne*qe**2*tau2/meSI
!   !!!---scaling of sig2 conductivity for solids
!   sig2=sig2*sig0(Z)/sigmax(Z)
!
!
!   sig=sig1+sig2   !!! max(sig1,sig2)
!   sig=max(sig,0.001_PR)
!
!   K1=ne*kb*(kb*Te)*tau/meSI*A_beta(x)
!   K2=ne*kb*(kb*Te)*tau2/meSI*A_beta(x)
!   !!!---scaling of K2 with reference value
!   K2=K2*lamb0(Z)/lambmax(Z)
!
!   K=K1+K2
!   K=max(K,0.001_PR)
!
!   S=(kb/qe)*A_gamma(x)
!
!   !sigS=sig1
!   !sigM=sig2
!
!   !sigS=sig1
!   !sigM=sig2
!
!
!   !tau=( 3.0_PR/4.0_PR*sqrt(meSI/(2.0_PR*Pi))*(kb*Te)**(1.5_PR) )/&
!   !    (Zbar**2*ni*qe**4*ln_lambda)
!
!
!   !!!---Spitzer
!   !bmax=sqrt(eps0*kb*Te/(ne*qe**2*(1.0_PR+Zbar)))
!   !vth=sqrt(3.0_PR*kB*Te/meSI)
!   !bmin=1.0_PR/(4.0_PR*Pi*eps0)*Zbar*qe**2/(meSI*vth**2)
!   !ln_lambda=log(bmax/bmin)
!
!   !bmin=Zbar*qe**2/(4.0_PR*Pi*eps0*meSI*vth**2)
!!   bmin=max(Zbar*qe**2/(4.0_PR*Pi*eps0*meSI*vth**2), h/(2.0_PR*meSI*vth))
!   !bmin=max(Zbar*qe**2/(4.0_PR*Pi*eps0*meSI*vth**2), h/(2.0_PR*meSI*vth))
!
!  ! ln_lambda=0.5_PR*log(1.0_PR+bmax**2/bmin**2)
!
!   !tau=(4.0_PR*Pi*eps0)**2*sqrt(meSI)*(kb*Te)**(1.5_PR)/&
!   !    (Pi*Zbar**2*ni*qe**4*ln_lambda)
!
!   !tau=3.0_PR/4.0_PR*sqrt(meSI)/sqrt(2.0_PR*pi)*(4.0_PR*Pi*eps0)**2*(kb*Te)**(1.5_PR)/(Zbar**2*ni*qe**4*ln_lambda)
!
!   !sigS=ne*qe**2*tau/meSI*32.0_PR/(3.0_PR*Pi)
!
!
!   !tau2=3.0_PR*Pi*hb**3/(4.0_PR*meSI*Zbar*qe**4*ln_lambda)
!
!
!   !sigM=ne*qe**2*tau2/meSI
!
!   !!!!!print*, 'ne=', ne, 'Zb=', Zbar, 'tau=', tau, 'lnL=', ln_lambda, 'tau2=', tau2
!   !print*, 'ne=', ne, 'Zb=', Zbar, 'lnL=', ln_lambda, 'bmin/max=', bmin, bmax
!
!end subroutine transport_properties
!
!
!subroutine compute_ZFS(Z,rho,T,ZTF,ZS,ZFS)
!
!   !!!---modification de Z par la formule (3) de Desjarlais 2001
!
!   implicit none
!   integer, intent(in) :: Z
!   real(PR), intent(in) :: rho, T
!   real(PR), intent(in) :: ZTF
!   real(PR), intent(out) :: ZFS,ZS
!   real(PR) :: na, Io, g0, g1, Ra, K, fe
!
!   na=rho/A(Z)*Avogadro*1.0e6_PR !!! m^-3
! 
!   Io=Eio1(Z)*qe
!   g0=2.0_PR
!   g1=1.0_PR
!   Ra=(3.0_PR/(4.0_PR*Pi*na))**(1.0_PR/3.0_PR)
!
!   !!!-----Version Originale
!
!   !K=2.0_PR*g1/g0*1.0_PR/na*(2.0_PR*Pi*me*kb*T/h**2)**1.5_PR*&
!   !  exp(-Io/(kb*T)*(1.0_PR-min(0.0_PR,(1.5_PR*qe**2/(Io*4.0_PR*pi*eps0*Ra))**1.5_PR)) )
!
!   !fe=0.5_PR*(sqrt(K**2+4.0_PR*K)-K)
!
!   !!!!!!fe=max(fe,1.0e-10_PR)
!
!
!   !ZFS=fe**(2.0_PR/ZTF**2)*ZTF+(1.0_PR-fe**(2.0_PR/ZTF**2))*fe 
!
!   !ZS=fe
!
!   !!!-----Version Fabien
! 
!   K=2.0_PR*g1/g0*(2.0_PR*Pi*me*kb*T/h**2)**1.5_PR*&
!     exp(-Io/(kb*T)*(1.0_PR-min(1.0_PR,(1.5_PR*qe**2/(Io*4.0_PR*pi*eps0*Ra))**1.5_PR)) )
!
!   fe=0.5_PR*(sqrt(K**2+4.0_PR*K*na)-K)/na
!  
!   ZS=fe
! 
!   ZFS=fe**(2.0_PR/ZTF**2)*ZTF+(1.0_PR-fe**(2.0_PR/ZTF**2))*fe
!
!   !print*, 'fe=', fe, 'K=', K, 'ZTF=', ZTF, 'ZFS=', ZFS
!
!end subroutine compute_ZFS
!
!real(PR) function chemical_potential(EF,T)
!
!   implicit none
!   real(PR), intent(in) :: EF, T
!   real(PR) :: Rm, denom, a(1:3), b(1:4), xi, chi
!   integer :: i
!
!
!   !!!---fit Rm3 de Managan 2015
!   !a(1)=0.19972_PR
!   !a(2)=0.17258_PR
!   !a(3)=0.145_PR
!
!   !b(1)=0.25829_PR
!   !b(2)=0.28756_PR
!   !b(3)=0.16842_PR
!   !b(4)=0.145_PR
!
!   !xi=sqrt(EF/(kb*T))
!
!   !Rm=4.0_PR/sqrt(3.0_PR) 
!   !do i=1,3
!   !   Rm=Rm+a(i)*xi**i
!   !enddo
!
!   !denom=1.0_PR
!   !do i=1,4   
!   !   denom=denom+b(i)*xi**i
!   !enddo
!
!   !chi=Rm/denom*xi**3
!
!   !chemical_potential=kb*T*log( exp(chi)-1.0_PR )
!
!   !!!---Zimmerman's form (see Managan 2015) :
!
!   xi=sqrt(EF/(kb*T))
!
!   chi=(0.7531_PR+0.1679_PR*xi+0.3108*xi**2)/&
!       (1.0_PR+0.2676_PR*xi+0.2280*xi**2+0.3099*xi**3)*xi**3
!
!   !write(6,*) 'xi=', xi, 'chi=', chi
!   !write(6,*) 'exp(chi)', exp(chi)
!
!
!   chemical_potential=kb*T*log( max(exp(min(chi,100.0_PR))-1.0_PR,1.0e-100_PR) )
!
!   if(chemical_potential.ne.chemical_potential)then
!     print*, 'PROBLEM CHEMICAL POTENTIAL'
!     stop
!   endif
!
!
!end function chemical_potential
!
!real(PR) function Fc_alpha(Z)
!
!    implicit none
!    real(PR), intent(in) :: Z
!    real(PR) :: x
!
!    x=1.0_PR/Z
!
!    Fc_alpha=0.295_PR/(1.0_PR-&
!    (0.0678_PR+0.4924_PR*x+0.9760_PR*x**2+0.3008_PR*x**3)/&
!    (0.0961_PR+0.7778_PR*x+1.5956_PR*x**2+1.3008_PR*x**3) )
!
!end function Fc_alpha
!
!
!
!real(PR) function A_alpha(x)
!
!    use mod_Fermi_Dirac, only : fd1h, fd6h, fd4h
!
!    implicit none
!    real(PR), intent(in) :: x
!    real(PR) :: F12, F3, F2
!
!    F12=fd1h(x)
!    F2=fd4h(x)
!    !F3=fd6h(x)
!
!    !A_alpha=4.0_PR/3.0_PR*F3/( (1.0_PR+exp(-1.0_PR*x))*F12**2 )
!
!    A_alpha=4.0_PR/3.0_PR*F2/( (1.0_PR+exp(-x))*F12**2 )
!
!end function A_alpha
!
!
!real(PR) function A_alpha2(x)
!
!    implicit none
!    real(PR), intent(in) :: x
!    real(PR) :: F12, F3, a(1:3), b(1:3), denom, y
!    integer :: i
!
!     y=log(1.0_PR+exp(x))
!     a(1)=3.39_PR
!     a(2)=0.347_PR
!     a(3)=0.129_PR
!     b(1)=0.0_PR
!     b(2)=0.511_PR
!     b(3)=0.124_PR
!
!     A_alpha2=0.0_PR 
!     do i=1,3
!        A_alpha2=A_alpha2+a(i)*y**(i-1)
!     enddo
!
!     denom=1.0_PR
!     do i=2,3
!        denom=denom+b(i)*y**(i-1)
!     enddo
!
!     A_alpha2=A_alpha2/denom
!
!end function A_alpha2
!
!real(PR) function A_beta(x)
!
!    use mod_Fermi_Dirac, only : fd1h, fd4h, fd6h, fd8h
!
!    implicit none
!    real(PR), intent(in) :: x
!    real(PR) :: F12, F2, F3, F4
!
!    F12=fd1h(x) 
!    F2=fd4h(x) 
!    F3=fd6h(x) 
!    F4=fd8h(x)
!
!    A_beta=20.0_PR/9.0_PR*F4*(1.0_PR-16.0_PR*F3**2/(15.0_PR*F4*F2))/&
!          ( (1.0_PR+exp(-x))*F12**2 )
!
!end function A_beta
!
!real(PR) function A_gamma(x)
!
!    use mod_Fermi_Dirac, only : fd1h, fd3h, fd4h, fd6h
!
!    implicit none
!    real(PR), intent(in) :: x
!    real(PR) :: F12, F32, F2, F3
!
!    F12=fd1h(x) 
!    F32=fd3h(x)
!    F2=fd4h(x) 
!    F3=fd6h(x) 
!
!    A_gamma=5.0_PR/3.0_PR*F32/F12-4.0_PR/3.0_PR*F3/F2
!
!end function A_gamma
!




end module mod_elec
