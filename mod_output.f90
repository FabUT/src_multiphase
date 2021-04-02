module mod_output

use mod_data

use mod_euler, only : Nl, Nx, t_init_Finger,t_euler, t_EM, t_src, t_relaxp, t_reset,&
                      t_p2c, t_c2p, t_nc2p, t_tot, t_reco

implicit none

character(len=50) :: FOLDER_OUT='output/'


contains

subroutine output(t,M)

implicit none

   !!!---grandeurs de mélange 
   real(PR), intent(in) :: t
   type(mesh), intent(in) :: M
   integer :: i, id, iph
   character(len=5) :: phasenum=''
   integer, save :: itwrite=0
   logical :: sonde=.true.
   logical :: outtime=.true.

   if(itwrite.eq.0)then

      open(unit=2,file=trim(adjustl(FOLDER_OUT))//'out.dat', status='replace')
    
      do iph=1,Nl

         write(phasenum,fmt='(I2)') iph
         id=10+iph
         open(unit=id,file=trim(adjustl(FOLDER_OUT))//'out_'//trim(adjustl(phasenum))//'.dat', status='replace')
        
      enddo
   
      !!!---sondes
      if(sonde)then
      open(unit=3,file=trim(adjustl(FOLDER_OUT))//'sonde.dat', status='replace')
      endif

      itwrite=1

      !!!---sondes
      if(outtime)then
      open(unit=4,file=trim(adjustl(FOLDER_OUT))//'time.dat', status='replace')
      write(4,*) '# t_tot t_euler t_src t_p2c t_nc2p t_c2p t_relaxp t_reset t_reco'
      endif

      itwrite=1

   endif


   !!!----sorties sur les temps de calcul
   if(outtime)then
   write(4,*) t_tot, t_euler, t_src, t_p2c, t_nc2p, t_c2p, t_relaxp, t_reset, t_reco
   endif


   !!!----sorties par phase------------------


   do iph=1,Nl

      id=10+iph

      do i=1,M%Nx

      write(id,fmt="(13(ES14.7,' '))") M%x(i), M%MF(i,1,1)%F(iph)%f, M%MF(i,1,1)%F(iph)%Y ,&
                            M%MF(i,1,1)%F(iph)%p,&   
                            M%MF(i,1,1)%F(iph)%rh,&  
                            M%MF(i,1,1)%F(iph)%eh,&
                            M%MF(i,1,1)%F(iph)%J,&
                            M%MF(i,1,1)%F(iph)%T,&
                            M%MF(i,1,1)%F(iph)%Cv,&
                            M%MF(i,1,1)%F(iph)%a(1),&
                            M%MF(i,1,1)%F(iph)%ee,&
                            M%MF(i,1,1)%F(iph)%detG,&
                            M%MF(i,1,1)%F(iph)%rh0

      enddo

      write(id,*) ""
      write(id,*) ""

   enddo

   !!!----sorties globales------------------

   do i=1,M%Nx

        write(2,fmt="(14(ES14.7,' '))") M%x(i), M%MF(i,1,1)%rh, M%MF(i,1,1)%p, M%MF(i,1,1)%vx, M%MF(i,1,1)%e, M%MF(i,1,1)%c,& 
                         M%MF(i,1,1)%sig(1,1), M%MF(i,1,1)%T,&
                         M%MF(i,1,1)%sigma, M%MF(i,1,1)%J, M%MF(i,1,1)%B, M%MF(i,1,1)%Cv,&
                         M%MF(i,1,1)%u, M%MF(i,1,1)%ee

   enddo

   write(2,*) ""
   write(2,*) ""


   !!!----sorties sondes------------------
   if(sonde)then

   do i=1,M%Nx

        write(3,fmt="(3(ES14.7,' '))") M%x(i), M%MF(i,1,1)%sonde(1), M%MF(i,1,1)%sonde(2)

   enddo

   write(3,*) ""
   write(3,*) ""

   endif


   !!!----------------------------------------

   !do i=1,M%Nx
   !  !write(2,fmt="(18(ES14.7,' '))") x(i), rh(i,1,1), p(i,1,1), vx(i,1,1), e(i,1,1), c(i,1,1), fl(1:Nl,i,1,1),&
   !  !                              sig(1:3,1,i,1,1), sig(1:3,2,i,1,1), sig(1:3,3,i,1,1)

   !  


   !  elseif(M%Nl.eq.2)then

   !  write(2,fmt="(26(ES14.7,' '))") M%x(i), M%F(i,1,1)%rh, M%F(i,1,1)%p, M%F(i,1,1)%vx, M%F(i,1,1)%e, M%F(i,1,1)%c,& 
   !                         M%F(i,1,1)%fl(1) ,  M%F(i,1,1)%fl(2) ,&
   !                         M%F(i,1,1)%Yl(1) ,  M%F(i,1,1)%Yl(2) ,&
   !                         M%F(i,1,1)%pl(1) ,  M%F(i,1,1)%pl(2) ,&
   !                         M%F(i,1,1)%rhl(1),  M%F(i,1,1)%rhl(2),&
   !                         M%F(i,1,1)%elh(1),  M%F(i,1,1)%elh(2),&
   !                         M%F(i,1,1)%sig(1,1), M%F(i,1,1)%T,&
   !                         M%F(i,1,1)%sigma, M%F(i,1,1)%J,M%F(i,1,1)%B,&
   !                         M%F(i,1,1)%Jl(1),M%F(i,1,1)%Jl(2),&
   !                         M%F(i,1,1)%Tl(1),M%F(i,1,1)%Tl(2),&
   !                         M%F(i,1,1)%Cv

   !  endif

   !  if(M%Nl.eq.1)then
   !  write(2,fmt="(26(ES14.7,' '))") M%x(i), M%F(i,1,1)%rh, M%F(i,1,1)%p, M%F(i,1,1)%vx, M%F(i,1,1)%e, M%F(i,1,1)%c,& 
   !                         M%F(i,1,1)%fl(1),0.0_PR,&
   !                         M%F(i,1,1)%Yl(1),0.0_PR,&
   !                         M%F(i,1,1)%pl(1),0.0_PR,&
   !                         M%F(i,1,1)%rhl(1),0.0_PR,&
   !                         M%F(i,1,1)%elh(1),0.0_PR,&
   !                         M%F(i,1,1)%sig(1,1), M%F(i,1,1)%T,&
   !                         M%F(i,1,1)%sigma, M%F(i,1,1)%J,M%F(i,1,1)%B,&
   !                         M%F(i,1,1)%Jl(1),0.0_PR,&
   !                         M%F(i,1,1)%Tl(1),0.0_PR,&
   !                         M%F(i,1,1)%Cv

   !  endif

   !enddo

   !write(2,*) ""
   !write(2,*) ""

   !i=1
   !write(iout,*) '----------------------------------'
   !write(iout,fmt="(3(ES14.7,' '))")  sig(1:3,1,i,1,1)
   !write(iout,fmt="(3(ES14.7,' '))")  sig(1:3,2,i,1,1)
   !write(iout,fmt="(3(ES14.7,' '))")  sig(1:3,3,i,1,1)
   !write(iout,*) '----------------------------------'


end subroutine output

end module mod_output
