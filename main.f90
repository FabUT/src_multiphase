program multi

                     
use mod_data

use mod_euler, only : Nl, Nx,init_euler, solve_euler, dtcfl,&
                      ph, soundspeed_mixt,&
                      test_convection, test_SOD_air, test_eau_air,&
                      test_eau_air_strong, test_rar_rar, test_cavitation,&
                      test_cavitation2, test_relax,test_fil,test_fil2,&
                      Teuler, Trelax, test_relax2, test_protection,&
                      test_relax3, test_Ndanou2015, test_meca,&
                      test_solid_solid_shock_tube, test_cavitation_Cu,&
                      test_advection_Cu,test_solid_gas_shock_tube,&
                      test_epoxy_spinel, test_epoxy_spinel_strong


use mod_elec, only : U,R,current,solve_elec

use mod_mat_tab, only : mat, mat_tab, read_table, test_lecture_tab

use air_ETL

use mod_output, only : output

implicit none

type(mesh) :: M
integer :: iph, i,j,k,i2,j2,Nitloc
real(PR) :: t, dtloc, cp, cv
real(PR) :: CFL=0.5_PR

!call test_lecture_tab

!stop

!cp=specific_heat(9100.0_PR,0.15e6_PR)
!cv=cp-RU/molar_mass(9100.0_PR,0.15e6_PR)
!
!
!
!print*, cp, cv, cp/cv
!
!cp=specific_heat(9200.0_PR,0.15e6_PR)
!cv=cp-RU/molar_mass(9200.0_PR,0.15e6_PR)
!
!
!print*, cp, cv, cp/cv
!
!
!stop

!!!---initialization des variables primitives
   write(iout,*) ""
   write(iout,fmt='(A)') "INITIALISATION"
   write(iout,*) ""

      !call init_euler(100,1,1,2)

      !call test_convection(M)
      
      !call test_SOD_air(M)

      !call test_snes

      !call test_relax(M)
      !call test_relax2(M)

      !call test_relax3(M)

      !call test_eau_air(M)

      !call test_eau_air_strong(M) !!?????

      !call test_epoxy_spinel(M)

      !call test_epoxy_spinel_strong(M)

      call test_rar_rar(M)

      !call test_cavitation(M)

      !call test_cavitation2

      !call test_fil(M)

      !call test_protection(M)

     !call test_fil2(M)

     !!!---meca :
   
     !call test_meca(M)

     !call test_advection_Cu(M)

     !call test_solid_solid_shock_tube(M) 

     !call test_solid_gas_shock_tube(M) 

     !call test_cavitation_Cu(M)

     !call test_Ndanou2015(M)

      write(iout,*) 'sound speed:'
      do iph=1,M%Nl
        write(iout,*) 'phase:', iph, 'c=: ',ph(iph)%soundspeed(M%MF(1,1,1)%F(iph),M%MF(1,1,1)%F(iph)%rh,M%MF(1,1,1)%F(iph)%p)
      enddo
      write(iout,*) 'cf=: ',M%MF(1,1,1)%c

   !!!---boucle en temps
   
   write(iout,*) ""
   write(iout,fmt='(A)') "DEBUT DE LA BOUCLE EN TEMPS"
   write(iout,*) ""

   !!!---sortie it=0 : 
   call output(t,M)
  
   t=0.0_PR
 
   DO it=1,M%Nt
   
      !!!-----------------------  Sorties
   
      dtloc=dtcfl(M,CFL)
      Nitloc=floor(M%dt/dtloc)+1
   
      write(iout,fmt='(A13,I6,A3,I6,A3,ES14.7,A5,ES14.7,A5,I6)')&
      "     ---- it=", it,' / ', M%Nt, ' t=', t, ' dtl=', dtloc, ' Nit~',Nitloc  

      !call SOLVE_elec(M,t)

      !print*, 'current=', current
 
      !if(t.lt.20.0e-9_PR)then 
      !   do i=1,Nx
      !      M%F(i,1,1)%Qe=1.5e15_PR*exp(-((M%x(i)-0.0e-2_PR)/2.0e-3_PR)**2)
      !      do iph=1,Nl
      !      !!!---Effet Joule
      !         !M%F(i,1,1)%Qel(iph)=M%F(i,1,1)%Qe*M%F(i,1,1)%sigmal(iph)/M%F(i,1,1)%sigma
      !         M%F(i,1,1)%Qel(iph)=M%F(i,1,1)%Qe
      !      enddo
      !   enddo 
      !else
      !   do i=1,Nx
      !      M%F(i,1,1)%Qe=0.0_PR
      !      M%F(i,1,1)%Qel(:)=0.0_PR
      !   enddo
      !endif


      call SOLVE_euler(M,M%dt,CFL)

 
      t=t+M%dt
   
      call output(t,M)

      !write(iout,fmt="(A3,3(A4,ES14.7,' '))") '   ', 'dt1=', dt1, 'dt2=', dt2, 'dt3=', dt3

   ENDDO   !!!---it


write(iout,fmt="(A3,2(A7,ES14.7,' '))") '   ', 'Teuler=', Teuler, 'Trelax=', Trelax

print*, 'maxvy=', maxval(M%MF(:,:,:)%vy), 'minvy=', minval(M%MF(:,:,:)%vy)

close(2)

write(iout,'(A)') ""
write(iout,'(A)') ' CALCUL MULTI FINI :) !'

end program multi
