module mod_input

use mod_data, only : iout, Input, materiaux, PR, Pi

implicit none

contains

subroutine LOAD_INPUT

implicit none
integer :: i,j,k,ierr,id
logical :: fin
character(len=200) :: dump

id=1
open(unit=id,file='input.dat',status='old',action='read',iostat=ierr) 

write(iout,*) ''
write(iout,'(A)') '   ------------- LECTURE input.dat: -----------------'
write(iout,*) ''

fin=.false.

do while(.not.fin)
   
    read(id,'(A)') dump
    
    if(scan(dump,'!#').eq.0)THEN
      
        !!!___nom
        if(index(dump,'nom=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%nom
            write(iout,'(A,A)') '   nom=', trim(adjustl(Input%nom))
        endif
 
        !!!___dim
        if(index(dump,'Ndim=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Ndim
            write(iout,'(A,I2)') '   Ndim=', Input%Ndim
        endif
    
        !!!___Nx
        if(index(dump,'Nx=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Nx
            write(iout,'(A,I6)') '   Nx=', Input%Nx
        endif
 
        !!!___Ny
        if(index(dump,'Ny=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Ny
            write(iout,'(A,I6)') '   Ny=', Input%Ny
        endif
 
        !!!___Nz
        if(index(dump,'Nz=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Nz
            write(iout,'(A,I6)') '   Nz=', Input%Nz
        endif
 
        !!!___Lx
        if(index(dump,'Lx=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Lx
            write(iout,'(A,ES12.5)') '   Lx=', Input%Lx
        endif
 
        !!!___Ly
        if(index(dump,'Ly=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Ly
            write(iout,'(A,ES12.5)') '   Ly=', Input%Ly
        endif
 
        !!!___Lz
        if(index(dump,'Lz=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Lz
            write(iout,'(A,ES12.5)') '   Lz=', Input%Lz
        endif

        !!!___Noutput
        if(index(dump,'Noutput=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Noutput
            write(iout,'(A,I6)') '   Noutput=', Input%Noutput
        endif
 
        !!!___dtoutput
        if(index(dump,'dtoutput=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%dtoutput
            write(iout,'(A,ES12.5)') '   dtoutput=', Input%dtoutput
        endif

        !!!___Nl : nombre de fluides
        if(index(dump,'Nl=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Nl
            write(iout,'(A,I3)') '   Nl=', Input%Nl
            allocate(Input%f2m(1:Input%Nl))
            !!!---lecture des affectations fluide -> materiaux
            do k=1,Input%Nl
                i=0
                do while(i.eq.0)
                   read(id,'(A)') dump
                   if(scan(dump,'!#').eq.0) i=index(dump,'->')
                enddo
                read(dump(i+2:),*) Input%f2m(k)
                write(iout,'(A,I3,A,I3)') '   fluid ', k, ' -> mat ', Input%f2m(k)
            enddo
        endif
 
        !!!___MATERIAUX
        if(index(dump,'MATERIAUX').ne.0)then

            write(iout,*) ''
            write(iout,'(A)') '   --- MATERIAUX :'
            !lecture N
            i=0 ; do while(i.eq.0) ; read(id,'(A)') dump ;  i=index(dump,'Nmat=') ; enddo
            i=index(dump,'=') ; read(dump(i+1:),*) Input%Nmat 
            write(iout,'(A9,I2)') '   Nmat: ', Input%Nmat
            allocate(Input%mat(1:Input%Nmat))
    
            if(Input%Nmat.gt.0)then
     
                do k=1,Input%Nmat
         
                    write(iout,*) ""
                    write(iout,'(A12,I2,A41)') '   Materiau ', k , ' ========================================'
                    write(iout,*) ""

                    !lecture typ
                    i=0
                    do while(i.eq.0)
                       read(id,'(A)') dump
                       i=index(dump,'typ=')
                    enddo
                    i=index(dump,'=') ; read(dump(i+1:),*) Input%mat(k)%typ
                    write(iout,'(A12,I2)') '      > typ=', Input%mat(k)%typ
     
                    if(Input%mat(k)%typ.eq.1)then
                        read(id,'(A)') dump ;  i=index(dump,'filemat=')
                        if(i.eq.0)then
                           write(iout,*) '   NO FILE FOR MAT 1 :( -> stop'
                           stop
                        else
                           j=index(dump,'=')
                           read(dump(j+1:),'(A)') Input%mat(k)%filemat
                           write(iout,'(A16,A)') '      > filemat=', trim(adjustl(Input%mat(k)%filemat))
                           call read_mat(trim(adjustl(Input%mat(k)%filemat)),Input%mat(k))
                        endif
                        !!!--modif du materiau?
                        i=0
                        do while(i.eq.0)
                            read(id,'(A)') dump
                            IF(scan(dump,'!#').eq.0)THEN
                            i=index(dump,'endmat')
                            if(i.eq.0)then
                                !!!---modif de rho ?
                                j=index(dump,'rho=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%rho 
                                   write(iout,*) '     > rho=', Input%mat(k)%rho
                                endif
                                !!!---modif de T ?
                                j=index(dump,'T=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%T
                                   write(iout,*) '     > T=', Input%mat(k)%T
                                endif
                                !!!---modif de P ?
                                j=index(dump,'P=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%P
                                   write(iout,*) '     > P=', Input%mat(k)%P
                                endif
                                !!!---modif de lambda ?
                                j=index(dump,'lambda=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%lambda
                                   write(iout,*) '     > lambda=', Input%mat(k)%lambda
                                endif
                                !!!---modif de Cv ?
                                j=index(dump,'Cv=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%Cv
                                   write(iout,*) '     > Cv=', Input%mat(k)%Cv
                                endif
                                !!!---modif de sigma ?
                                j=index(dump,'sigma=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%sigma
                                   write(iout,*) '     > sigma=', Input%mat(k)%sigma
                                endif
                                !!!---modif de epsr ?
                                j=index(dump,'epsr=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%epsr
                                   write(iout,*) '     > epsr=', Input%mat(k)%epsr
                                endif
                                !!!---modif de mur ?
                                j=index(dump,'mur=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%mur
                                   write(iout,*) '     > mur=', Input%mat(k)%mur
                                endif
                                !!!---modif de kappa ?
                                j=index(dump,'kappa=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%kappa
                                   write(iout,*) '     > kappa=', Input%mat(k)%kappa
                                endif
                                !!!---modif de eps ?
                                j=index(dump,'eps=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%eps
                                   write(iout,*) '     > eps=', Input%mat(k)%eps
                                endif
                                !!!---modif de rho0 ?
                                j=index(dump,'rho0=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%rho0
                                   write(iout,*) '     > rho=', Input%mat(k)%rho0
                                endif
                                !!!---modif de T0 ?
                                j=index(dump,'T0=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%T0
                                   write(iout,*) '     > T0=', Input%mat(k)%T0
                                endif
                                !!!---modif de P0 ?
                                j=index(dump,'P0=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%P0
                                   write(iout,*) '     > P0=', Input%mat(k)%P0
                                endif
                                !!!---modif de gam ?
                                j=index(dump,'gam=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%gam
                                   write(iout,*) '     > gam=', Input%mat(k)%gam
                                endif
                                !!!---modif de P0 ?
                                j=index(dump,'pinf=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%pinf
                                   write(iout,*) '     > pinf=', Input%mat(k)%pinf
                                endif
                                !!!---modif de R ?
                                j=index(dump,'R=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%R
                                   write(iout,*) '     > R=', Input%mat(k)%R
                                   Input%mat(k)%sec=Pi*Input%mat(k)%R**2
                                endif
                                !!!---modif de ep ?
                                j=index(dump,'ep=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%ep
                                   write(iout,*) '     > ep=', Input%mat(k)%ep
                                endif
                                !!!---modif de sec ?
                                j=index(dump,'sec=')
                                if(j.gt.0)then
                                   j=index(dump,'=')
                                   read(dump(j+1:),*) Input%mat(k)%sec
                                   write(iout,*) '     > sec=', Input%mat(k)%sec
                                   Input%mat(k)%R=sqrt(Input%mat(k)%sec/Pi)
                                endif
                                ENDIF !!! scan(!#)
                            endif
                        enddo !!! while(endmat)
                        !!!---calcul de U0
                        Input%mat(k)%U0=Input%mat(k)%Cv*Input%mat(k)%T0
                        Input%mat(k)%U=Input%mat(k)%Cv*Input%mat(k)%T

                    endif
                enddo
            endif !!! Nmat>0
        endif !!! MATERIAUX

        IF(index(dump,'FIN').ne.0)then
            FIN=.true.
        ENDIF

    endif !!! scan #!
    
enddo !!! do while(.not.fini)

write(iout,*) ''
write(iout,'(A)') '   --------------------------------------------------'
write(iout,*) ''

end subroutine LOAD_INPUT

subroutine read_mat(file_in,mat)

   implicit none
   character(len=*), intent(in) :: file_in
   type(materiaux), intent(inout) :: mat
   integer :: id, i, j, ierr
   character(len=200) :: dump    
   logical :: FIN
   logical :: verb_mat=.true.

   id = 2

   open( unit=id, file=file_in, status='old',action='read',&
             form='formatted', IOSTAT=ierr)
   if( ierr .ne. 0 )then
   write(iout,*) "   > There was an error while opening  ", file_in
   else
   write(iout,*) "     -------------------------------------------------"
   write(iout,*) "     read_mat: open ", file_in
   endif

   FIN=.false.

   do while(.not.FIN)
   
      read(id,'(A)') dump
   
      IF(scan(dump,'!#').eq.0)THEN
 
         IF(index(dump,'FIN').ne.0)then
            FIN=.true.
         ENDIF 
 
         !!!--- nom ---
         j=index(dump,'nom=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%nom 
             if(verb_mat) write(iout,*) '       > nom=', mat%nom
         endif
       
         !!!--- rho ---
         j=index(dump,'rho=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%rho 
             if(verb_mat) write(iout,*) '       > rho=', mat%rho
         endif
         !!!--- T ---
         j=index(dump,'T=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%T
             if(verb_mat) write(iout,*) '       > T=', mat%T
         endif
         !!!--- P ---
         j=index(dump,'P=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%P
             if(verb_mat) write(iout,*) '       > P=', mat%P
         endif
         !!!--- lambda ---
         j=index(dump,'lambda=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%lambda
             if(verb_mat) write(iout,*) '       > lambda=', mat%lambda
         endif
         !!!--- Cv ---
         j=index(dump,'Cv=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%Cv
             if(verb_mat) write(iout,*) '       > Cv=', mat%Cv
             mat%U=mat%T*mat%Cv
             if(verb_mat) write(iout,*) '       > U=', mat%U
         endif
         !!!--- sigma ---
         j=index(dump,'sigma=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%sigma
             if(verb_mat) write(iout,*) '       > sigma=', mat%sigma
         endif
         !!!--- epsr ---
         j=index(dump,'epsr=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%epsr
             if(verb_mat) write(iout,*) '       > epsr=', mat%epsr
         endif
         !!!--- mur ---
         j=index(dump,'mur=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%mur
             if(verb_mat) write(iout,*) '       > mur=', mat%mur
         endif
         !!!--- kappa ---
         j=index(dump,'kappa=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%kappa
             if(verb_mat) write(iout,*) '       > kappa=', mat%kappa
         endif
         !!!--- eps ---
         j=index(dump,'eps=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%eps
             if(verb_mat) write(iout,*) '       > eps=', mat%eps
         endif
         !!!--- rho0 ---
         j=index(dump,'rho0=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%rho0 
             if(verb_mat) write(iout,*) '       > rho0=', mat%rho0
         endif
         !!!--- T0 ---
         j=index(dump,'T0=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%T0
             if(verb_mat) write(iout,*) '       > T0=', mat%T0
             mat%U0=mat%T0*mat%Cv
             if(verb_mat) write(iout,*) '       > U0=', mat%U0
         endif
         !!!--- P0 ---
         j=index(dump,'P0=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%P0
             if(verb_mat) write(iout,*) '       > P0=', mat%P0
         endif
         !!!--- eps ---
         j=index(dump,'gam=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%gam
             if(verb_mat) write(iout,*) '       > gam=', mat%gam
         endif
         !!!--- eps ---
         j=index(dump,'pinf=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%pinf
             if(verb_mat) write(iout,*) '       > gam=', mat%pinf
         endif
         !!!--- R ---
         j=index(dump,'R=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%R
             if(verb_mat) write(iout,*) '       > R=', mat%R
             mat%sec=Pi*mat%R**2
         endif
         !!!--- ep ---
         j=index(dump,'ep=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%ep
             if(verb_mat) write(iout,*) '       > ep=', mat%ep
         endif
         !!!--- sec ---
         j=index(dump,'sec=')
         if(j.gt.0)then
             j=index(dump,'=')
             read(dump(j+1:),*) mat%sec
             if(verb_mat) write(iout,*) '       > sec=', mat%sec
             mat%R=sqrt(mat%sec/Pi)
         endif

      ENDIF !!! scan(!,#)

   
   ENDDO !!! DO WHILE(.not.FIN)
   
   close(id)
   write(iout,'(A)') "      read_mat fini"
   write(iout,'(A)') "      -------------------------------------------------"

end subroutine read_mat



end module mod_input
