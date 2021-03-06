!***********************************************************************
!*                                                                     *
!*                   NONLINEAR EULER EQUATION-AEROACUSTIC              *
!*                                                                     *
!* UPDATE            : 29-08/2013                                      *
!* MODIFICATION      : IMPLEMENTAÇÃO DA PML NAS  EQUAÇÕES NAO LINEARES *
!*                     DE EULER EM VARIAVEIS PRIMITIVAS.               *
!                                                                      *
!* DATE              : 16-08/2013                                      *
!* MODIFY  BY        : JHONATAN                                        *
!* BASED ON          : PML ABC FOR NONLINEAR EULER EQUATIONS IN        *
!*                     PRIMITIVE VARIABLES,AIAA 2009-6 - HU, LIN, LI *
!* GRID              : GENERAL GRID(M,N)                               *
!* ALGORITHMS        : EXPLICIT METHOD- COMPLETE EULER EQUATION        *
!*                     2D-STRECHING NA PML EN X E Y, FILTRO 10 ORDEM   *
!* RESOLT METHOD     : RUNGE-KUTA 4 ORDER                              *
!* BOUNDARY CONDITION:                                                 * 
!* FORTRAN 90        : MODULE , THOMA'S ALGORITHM                      *
!*                                                                     *
!***********************************************************************

program EULER

use global
use globalq
use initialize
use outt
use filtering
use rk

!........................................................................
      implicit none

      integer:: iter,icount,cf,i,j
      
      real ::timei,timef
      real(kind=ip)::t
      real(kind=ip),dimension(4,im,jm):: U
      character(40) :: temporalfile1,temporalfile2,temporalfile3
      
      write (temporalfile1,'("temp1_D60.dat")') 
      open(unit=1000,file=temporalfile1,status='unknown')

      write (temporalfile2,'("temp2_D60.dat")') 
      open(unit=1001,file=temporalfile2,status='unknown')
        
      write (temporalfile3,'("temp3_D60.dat")') 
      open(unit=1002,file=temporalfile3,status='unknown')
      
      call CPU_TIME(timei)
      
      pi=acos(-1.d0)  

!................Inicialização das variaveis..............................
      call init(U)
!.........................................................................
!Contador Graficas
      icount = 0
!Contador filtro 
      cf=0
!.........................................................................
      write(*,*)m(411),n(221)

      do iter = 0,maxit

         t=iter*dt
       
         call rkpml(U,t)

         if(cf==50)then

            write(*,*)iter 

!        do j =1,jmax
!           write(200,*)n(j),u(1,111,j)
!        end do
         
           !call filtercompact1(U(1,imaxpml+1:imax-D,:))
!         call filtercompactx(U)
 !        call filtercompacty(U)

        call output(U,icount)
!        do j =1,jmax
!           write(201,*)n(j),u(1,111,j)
!        end do

             icount = icount+1
             cf=0
         end if

         write(1000,*)dt*iter,U(4,206,26),u(3,206,26)
         write(1001,*)dt*iter,U(4,411,221),u(3,411,221)
         write(1002,*)dt*iter,U(4,206,56),u(3,206,56)

         cf=cf+1

      end do 

         close(unit=1000)
         close(unit=1001)
         close(unit=1002)

      call CPU_TIME(timef) 
      
      open(unit=10,file="tempototal_D60.dat",status='unknown')
      write(*,*) 'Tempo :',timef-timei,'segundos'
      write(10,*)'Tempo :',timef-timei,'segundos',',dx=dy=',dm
      close(unit=10)

end program EULER

