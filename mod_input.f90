module mod_input

use mod_data, only : iout, Input, materiaux, PR, Pi, boite

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


        !!!___boites
        if(index(dump,'Nboite=').ne.0)then
            i=index(dump,'=')
            read(dump(i+1:),*) Input%Nboite
            write(iout,'(A,I2)') '   Nboite=', Input%Nboite
            if(Input%Nboite.gt.0)then
              allocate(Input%box(1:Input%Nboite))
              k=0
              do while(k.lt.Input%Nboite)
                  read(id,'(A)') dump
                  if(scan(dump,'!#').eq.0)THEN
                    i=index(dump,'box=')
                    if(i.gt.0)then
                       k=k+1
                       call readbox(dump(i:),Input%box(k))
                    endif
                  endif
              enddo
            endif
        endif

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

subroutine readbox(dump,box)

   implicit none
   character(len=*), intent(in) :: dump
   type(boite), intent(out) :: box
   integer :: j1,j2,j3

   j1=index(dump,'box=')
   j2=index(dump,'[')
   j3=index(dump,']')

   if(j1.eq.0)then
      write(iout,'(A)') 'error reading box: j1'
      write(iout,'(A)') dump
      call flush(iout)
   elseif(j2.eq.0)then
      write(iout,'(A)') 'error reading box: j2'
      write(iout,'(A)') dump
      call flush(iout)
   elseif(j3.eq.0)then
      write(iout,'(A)') 'error reading box: j3'
      write(iout,'(A)') dump
      call flush(iout)
   else

      box%id=1

      read(dump(j1+4:j2-1),*) box%forme
      write(iout,'(A,I2)') '        > forme=', box%forme

      if(box%forme.eq.1)then !!! boite rectangulaire
         read(dump(j2+1:j3-1),*)  box%xmin, box%xmax, box%ymin, box%ymax, box%zmin, box%zmax
         write(iout,"(A,6(ES12.5,' '),A2)") '        > [',box%xmin, box%xmax, box%ymin, box%ymax, box%zmin, box%zmax,' ]' 
      elseif(box%forme.eq.2)then !!! sphere de centre x,y,z de rayon R
         read(dump(j2+1:j3-1),*)  box%xc, box%yc, box%zc, box%R
         write(iout,"(A,4(ES12.5,' '),A2)") '        > [',box%xc, box%yc, box%zc, box%R,' ]' 
      elseif(box%forme.eq.3)then !!! cylindre de diametre D, de longueur L, centre du cercle de base en x,y,z
         read(dump(j2+1:j3-1),*)  box%dir, box%xc, box%yc, box%zc, box%D, box%L
         write(iout,"(A,A1,5(ES12.5,' '),A2)") '        > [', box%dir, box%xc, box%yc, box%zc, box%D, box%L, ' ]' 
         box%R=0.5_PR*box%D 
      elseif(box%forme.eq.4)then !!! cylindre creux entre Di et D de longueur L, centre du cercle de base en x,y,z
         read(dump(j2+1:j3-1),*)  box%dir, box%xc, box%yc, box%zc, box%Di, box%D, box%L
         write(iout,"(A,A1,6(ES12.5,' '),A2)") '        > [', box%dir, box%xc, box%yc, box%zc, box%Di, box%D, box%L, ' ]' 
         box%R=0.5_PR*box%D
         box%Ri=0.5_PR*box%Di
      elseif(box%forme.eq.5)then !!! stripes: cylindre creux Di-D de longueur L et stripes entre teta1 et teta2:
         read(dump(j2+1:j3-1),*)  box%dir, box%xc, box%yc, box%zc, box%Di, box%D, box%L, box%teta0, box%teta1, box%teta2
         write(iout,"(A,A1,9(ES12.5,' '),A2)") '        > [', box%dir, box%xc, box%yc, box%zc, box%Di, box%D, box%L, &
                                                           box%teta0, box%teta1, box%teta2, ' ]' 
         box%R=0.5_PR*box%D
         box%Ri=0.5_PR*box%Di
         box%teta0=box%teta0*Pi/180.0_PR
         box%teta1=box%teta1*Pi/180.0_PR
         box%teta2=box%teta2*Pi/180.0_PR
      endif
 

   endif 

end subroutine readbox

logical function insidebox(x,y,z,box)

   implicit none
   real(PR), intent(in) :: x,y,z
   type(boite), intent(in) :: box
   real(PR) :: R, dx, dy, dz, teta

   insidebox=.false.

   if(box%id.eq.0)then

      insidebox=.true.

   elseif(box%forme.eq.1)then !!! paralellepipede rectange
   
      if( (x.ge.box%xmin.and.x.le.box%xmax) .and. &
          (y.ge.box%ymin.and.y.le.box%ymax) .and. &
          (z.ge.box%zmin.and.z.le.box%zmax) )then

         insidebox=.true.

      endif

   elseif(box%forme.eq.2)then !!! sphere

      R=sqrt((x-box%xc)**2+(y-box%yc)**2+(z-box%zc)**2)
 
      if(R.lt.box%R)then

         insidebox=.true.

      endif

   elseif(box%forme.ge.3.and.box%forme.le.6)then 
   !!! cylindre plein(3)/creux(4)/creux avec stripes(5) 

      R=1.0e30_PR ; teta=0.0_PR

      if(box%dir.eq.'x')then
         if(box%L.gt.0.0_PR)then 
            if(x.gt.box%xc.and.x.lt.(box%xc+box%L))then
               dy=(y-box%yc) ; dz=(z-box%zc)
               R=sqrt(dy**2+dz**2)
            endif
         else
            if(x.lt.box%xc.and.x.gt.(box%xc+box%L))then
               dy=(y-box%yc) ; dz=(z-box%zc)
               R=sqrt(dy**2+dz**2)
            endif
         endif
      elseif(box%dir.eq.'y')then
         if(box%L.gt.0.0_PR)then 
            if(y.gt.box%yc.and.y.lt.(box%yc+box%L))then
               dx=(x-box%xc) ; dz=(z-box%zc)
               R=sqrt(dx**2+dz**2)
            endif
         else
            if(y.lt.box%yc.and.y.gt.(box%yc+box%L))then
               dx=(x-box%xc) ; dz=(z-box%zc)
               R=sqrt(dx**2+dz**2)
            endif
         endif
      elseif(box%dir.eq.'z')then
         if(box%L.gt.0.0_PR)then 
            if(z.gt.box%zc.and.z.lt.(box%zc+box%L))then
               dx=(x-box%xc) ; dy=(y-box%yc) 
               R=sqrt(dx**2+dy**2)
            endif
         else
            if(z.lt.box%zc.and.z.gt.(box%zc+box%L))then
               dx=(x-box%xc) ; dy=(y-box%yc) 
               R=sqrt(dx**2+dy**2)
            endif
         endif
      endif
      if(R.lt.box%R)then
         if(box%forme.eq.3)then !!! cylindre plein
            insidebox=.true.
         elseif(box%forme.eq.4)then !!! cylindre creux
            if(R.gt.box%Ri)then 
               insidebox=.true.            
            endif
         elseif(box%forme.eq.5)then !!! stripes
            if(R.gt.box%Ri)then
               if(box%dir.eq.'x')then 
                  if(dy.ne.0.0_PR.or.dz.ne.0.0_PR) teta=atan2(dz,dy)
               elseif(box%dir.eq.'y')then
                  if(dx.ne.0.0_PR.or.dz.ne.0.0_PR) teta=atan2(dz,dx)
               elseif(box%dir.eq.'z')then
                  if(dx.ne.0.0_PR.or.dy.ne.0.0_PR) teta=atan2(dy,dx)
               endif
               if(teta.lt.0.0_PR) teta=teta+2.0_PR*Pi
               teta=teta-box%teta0
               insidebox=inside_dteta(teta,box%teta1,box%teta2)
            endif 
         endif
      endif

   endif !!! forme

end function insidebox 

recursive function inside_dteta(teta0, teta1, teta2) result(inside)

   !!!---
   !!! Check if there is k such that teta0 is between [k*teta1,k*teta2] 
   !!!---
   implicit none
   logical :: inside
   real(PR), intent(in) :: teta0
   real(PR), intent(in) :: teta1, teta2
   real(PR) :: teta

   teta=teta0
   inside=.false.

   if(teta.gt.teta2)then
      teta=teta-teta2
      inside=inside_dteta(teta,teta1,teta2)
   else
      if(teta.gt.teta1)then
         inside=.true.
      else
         inside=.false.
      endif
   endif

end function inside_dteta

end module mod_input
