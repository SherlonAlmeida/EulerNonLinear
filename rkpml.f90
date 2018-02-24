module rk

use global
use globalq
use rhspml
use rhsq
use cc

contains

subroutine rkpml(U,t)
     
      implicit none         
      integer:: k, i, j

      real(kind=ip)::t,dt2, dt6
      real(kind=ip),dimension(4,im,jm):: U
      real(kind=ip),dimension(4,im,jm):: Up1,Up2,Up3,Up4
      real(kind=ip),dimension(4,im,jm):: dUdt,dUdtp1,dUdtp2,dUdtp3
      real(kind=ip),dimension(4,im,jm):: q1p1,q1p2,q1p3
      real(kind=ip),dimension(4,im,jm):: q2p1,q2p2,q2p3
   
!......................................................................
      dt2 = dt * 0.5d0
      dt6 = dt / 6.d0

!...................Parametros Fonte...................................
      w  = 0.5d0*pi
      r0 = 0.03d0  

!....................fonte...................................................
     do i=1,imax
     do j=1,jmax
     !s(i,j)= sin(w*(t+dt2))*exp(-log(2.d0)*((m(i)+0.5d0)**2+(n(j))**2)/(r0*r0))
     s(i,j)= sin(w*(t))*exp(-log(2.d0)*((m(i)+0.5d0)**2+(n(j))**2)/(r0*r0))
     end do
     end do
!.......................................................................

      call rhs(dUdt,U,q1,q2)
      call rhsq1(u,q1,q1dtp0)
      call rhsq2(u,q2,q2dtp0)

!.................. First Step .....................................

      do k=1,4 
        do j = 2,jmax-1
            do i = 2,imax-1
                Up1(k,i,j) = U(k,i,j) + dt2*dUdt(k,i,j)
            end do
        end do
      end do 

      call ccpml(Up1)

      call q1fun(q1p1,q1dtp0,1,imaxpml,1,jmax,dt2)                     
      call q2fun(q2p1,q2dtp0,1,imaxpml,1,jmax,dt2)                     

      call q1fun(q1p1,q1dtp0,imax-D+1,imax,1,jmax,dt2)                   
      call q2fun(q2p1,q2dtp0,imax-D+1,imax,1,jmax,dt2)                   

      call q1fun(q1p1,q1dtp0,imaxpml+1,imax-D,1,jmaxpml,dt2)             
      call q2fun(q2p1,q2dtp0,imaxpml+1,imax-D,1,jmaxpml,dt2)             

      call q1fun(q1p1,q1dtp0,imaxpml+1,imax-D,jmax-D+1,jmax,dt2)           
      call q2fun(q2p1,q2dtp0,imaxpml+1,imax-D,jmax-D+1,jmax,dt2)           

!....................Fonte...................................................
      do i=1,imax
         do j=1,jmax
            !s(i,j)= sin(w*(t+dt2))*exp(-log(2.d0)*((m(i)+0.5d0)**2+(n(j))**2)/(r0*r0))
            s(i,j)= sin(w*(t))*exp(-log(2.d0)*((m(i)+0.5d0)**2+(n(j))**2)/(r0*r0))
         end do
      end do
!.......................................................................

      call rhs(dUdtp1,Up1,q1p1,q2p1)
      call rhsq1(up1,q1p1,q1dtp1)
      call rhsq2(up1,q2p1,q2dtp1)

!..................  Second Step.....................................

      do k=1,4 
        do j = 2,jmax-1
            do i = 2,imax-1
                Up2(k,i,j) = U(k,i,j) + dt2*dUdtp1(k,i,j)
            end do
        end do
      end do 

      call ccpml(Up2)

      call q1fun(q1p2,q1dtp1,1,imaxpml,1,jmax,dt2)                     
      call q1fun(q1p2,q1dtp1,imax-D+1,imax,1,jmax,dt2)                   
      call q1fun(q1p2,q1dtp1,imaxpml+1,imax-D,1,jmaxpml,dt2)             
      call q1fun(q1p2,q1dtp1,imaxpml+1,imax-D,jmax-D+1,jmax,dt2)           

      call q2fun(q2p2,q2dtp1,1,imaxpml,1,jmax,dt2)                     
      call q2fun(q2p2,q2dtp1,imax-D+1,imax,1,jmax,dt2)                   
      call q2fun(q2p2,q2dtp1,imaxpml+1,imax-D,1,jmaxpml,dt2)             
      call q2fun(q2p2,q2dtp1,imaxpml+1,imax-D,jmax-D+1,jmax,dt2)           

!....................fonte...................................................
        do i=1,imax
           do j=1,jmax
              !s(i,j)= sin(w*(t+dt))*exp(-log(2.d0)*((m(i)+0.5d0)**2+(n(j))**2)/(r0*r0))
              s(i,j)= sin(w*(t))*exp(-log(2.d0)*((m(i)+0.5d0)**2+(n(j))**2)/(r0*r0))
           end do
        end do
!.......................................................................

      call rhs(dUdtp2,Up2,q1p2,q2p2)
      call rhsq1(up2,q1p2,q1dtp2)
      call rhsq2(up2,q2p2,q2dtp2)

!..................  Third Step.....................................

      do k=1,4 
        do j = 2,jmax-1
            do i = 2,imax-1
                Up3(k,i,j) = U(k,i,j) + dt*dUdtp2(k,i,j)
            end do
        end do
      end do 

      call ccpml(Up3)

      call q1fun(q1p3,q1dtp2,1,imaxpml,1,jmax,dt)                     
      call q1fun(q1p3,q1dtp2,imax-D+1,imax,1,jmax,dt)                   
      call q1fun(q1p3,q1dtp2,imaxpml+1,imax-D,1,jmaxpml,dt)             
      call q1fun(q1p3,q1dtp2,imaxpml+1,imax-D,jmax-D+1,jmax,dt)           

      call q2fun(q2p3,q2dtp2,1,imaxpml,1,jmax,dt)                     
      call q2fun(q2p3,q2dtp2,imax-D+1,imax,1,jmax,dt)                   
      call q2fun(q2p3,q2dtp2,imaxpml+1,imax-D,1,jmaxpml,dt)             
      call q2fun(q2p3,q2dtp2,imaxpml+1,imax-D,jmax-D+1,jmax,dt)           

!....................fonte...................................................
        do i=1,imax
           do j=1,jmax
              !s(i,j)= sin(w*(t+dt6))*exp(-log(2.d0)*((m(i)+0.5d0)**2+(n(j))**2)/(r0*r0))
              s(i,j)= sin(w*(t))*exp(-log(2.d0)*((m(i)+0.5d0)**2+(n(j))**2)/(r0*r0))
           end do
        end do
!.......................................................................

      call rhs(dUdtp3,Up3,q1p3,q2p3)
      call rhsq1(up3,q1p3,q1dtp3)
      call rhsq2(up3,q2p3,q2dtp3)

!.......................Fourt Step ..................................
!
      call q1f(1,imaxpml,1,jmax,dt6)                    
      call q1f(imax-D+1,imax,1,jmax,dt6)                
      call q1f(imaxpml+1,imax-D,1,jmaxpml,dt6)          
      call q1f(imaxpml+1,imax-D,jmax-D+1,jmax,dt6)      
                                                 
      call q2f(1,imaxpml,1,jmax,dt6)                    
      call q2f(imax-D+1,imax,1,jmax,dt6)                
      call q2f(imaxpml+1,imax-D,1,jmaxpml,dt6)          
      call q2f(imaxpml+1,imax-D,jmax-D+1,jmax,dt6)      


      do k=1,4 
        do j = 2,jmax-1
          do i = 2,imax-1
               U(k,i,j) = U(k,i,j) + dt6 * (dUdt(k,i,j) +   &     
           &       2.d0*dUdtp1(k,i,j) + 2.d0*dUdtp2(k,i,j) + dUdtp3(k,i,j))
     
          end do
        end do
       end do

       call ccpml(U)
    
end subroutine

subroutine q1fun(qp,qt,ii,fi,ij,fj,dt)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: qp
      real(kind=ip),dimension(4,im,jm):: qt
      real(kind=ip):: dt

      do k=1,4
        do j = ij,fj
          do i = ii,fi

              qp(k,i,j) = q1(k,i,j) + dt*qt(k,i,j)

          end do
        end do
     end do

end subroutine

subroutine q2fun(qp,qt,ii,fi,ij,fj,dt)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: qp
      real(kind=ip),dimension(4,im,jm):: qt
      real(kind=ip):: dt

      do k=1,4
        do j = ij,fj
          do i = ii,fi

              qp(k,i,j) = q2(k,i,j) + dt*qt(k,i,j)

          end do
        end do
     end do

end subroutine

subroutine q1f(ii,fi,ij,fj,dt)

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip):: dt


      do k=1,4
        do i = ii,fi
           do j=ij,fj

               q1(k,i,j) = q1(k,i,j) + dt * (q1dtp0(k,i,j) +          &
     &           2.d0*q1dtp1(k,i,j) + 2.d0*q1dtp2(k,i,j) + q1dtp3(k,i,j))
     
 
          end do
        end do
      end do

end subroutine

subroutine q2f(ii,fi,ij,fj,dt)

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip):: dt


      do k=1,4
        do i = ii,fi
           do j=ij,fj

               q2(k,i,j) = q2(k,i,j) + dt * (q2dtp0(k,i,j) +          &
     &           2.d0*q2dtp1(k,i,j) + 2.d0*q2dtp2(k,i,j) + q2dtp3(k,i,j))
     
 
          end do
        end do
      end do

end subroutine


end module

