module  eqeuler 

use global
use global2 

contains

subroutine nonlinearmatriz(U,dE,A,B)

implicit none                     
integer::k,l,i,j 
real(kind=ip),dimension(4,im,jm)::U,dE
real(kind=ip),dimension(4,4,im,jm)::A,B
real(kind=ip),dimension(im,jm)::temp,temp2


temp2(imaxpml+1:imax-D,jmaxpml+1:jmax-D) = 0.d0

do k=1,4
 do l=1,4

  do i=imaxpml+1,imax-D
    do j=jmaxpml+1,jmax-D
        temp(i,j) = A(k,l,i,j)*dUdm(l,i,j) + B(k,l,i,j)*dUdn(l,i,j)
    end do
  end do
  temp2(imaxpml+1:imax-D,jmaxpml+1:jmax-D) = temp(imaxpml+1:imax-D,jmaxpml+1:jmax-D)+temp2(imaxpml+1:imax-D,jmaxpml+1:jmax-D)

end do
       dE(k,imaxpml+1:imax-D,jmaxpml+1:jmax-D)=temp2(imaxpml+1:imax-D,jmaxpml+1:jmax-D)
       temp2(imaxpml+1:imax-D,jmaxpml+1:jmax-D) = 0.d0
end do

end subroutine 

subroutine nonlinear(U,dE,A,B)

implicit none                     
integer::k,l,i,j 
real(kind=ip),dimension(4,im,jm)::U,dE
real(kind=ip),dimension(4,4,im,jm)::A,B
real(kind=ip),dimension(im,jm)::temp,temp2

do i=imaxpml+1,imax-D
     do j=jmaxpml+1,jmax-D
       dE(1,i,j) = A(1,1,i,j)*dUdm(1,i,j) + A(1,2,i,j)*dUdm(2,i,j) + &
  &                B(1,1,i,j)*dUdn(1,i,j) + B(1,3,i,j)*dUdn(3,i,j)
       dE(2,i,j) = A(2,2,i,j)*dUdm(2,i,j) + A(2,4,i,j)*dUdm(4,i,j) + &
  &                B(2,2,i,j)*dUdn(2,i,j)
       dE(3,i,j) = A(3,3,i,j)*dUdm(3,i,j) + &
  &                B(3,3,i,j)*dUdn(3,i,j) + B(3,4,i,j)*dUdn(4,i,j) 
       dE(4,i,j) = A(4,2,i,j)*dUdm(2,i,j) + A(4,4,i,j)*dUdm(4,i,j) + &
  &                B(4,3,i,j)*dUdn(3,i,j) + B(4,4,i,j)*dUdn(4,i,j)
     end do
end do


end subroutine 
end module
