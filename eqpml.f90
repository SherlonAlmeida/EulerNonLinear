module eqpml

use global
use global2

contains

subroutine pml(U,q1,q2,dpml,A,B)

implicit none 
integer:: i, j
real(kind=ip),dimension(4,im,jm):: U,q1,q2,dpml
real(kind=ip),dimension(4,4,im,jm):: A,B


call corner(U,A,B,q1,q2,dpml,1,imaxpml,1,jmaxpml) 
call ylayer(U,A,B,q1,q2,dpml,imaxpml+1,imax-D,1,jmaxpml) 
call corner(U,A,B,q1,q2,dpml,imax-D+1,imax,1,jmaxpml) 

call xlayer(U,A,B,q1,q2,dpml,1,imaxpml,jmaxpml+1,jmax-D) 
call xlayer(U,A,B,q1,q2,dpml,imax-D+1,imax,jmaxpml+1,jmax-D) 

call corner(U,A,B,q1,q2,dpml,1,imaxpml,jmax-D+1,jmax) 
call ylayer(U,A,B,q1,q2,dpml,imaxpml+1,imax-D,jmax-D+1,jmax) 
call corner(U,A,B,q1,q2,dpml,imax-D+1,imax,jmax-D+1,jmax) 

end subroutine              

subroutine xlayer(U,A,B,q1,q2,dpml,ii,fi,ij,fj)
implicit none                     
integer::k,l,i,j,ii,fi,ij,fj
real(kind=ip),dimension(4,im,jm)::U,q1,q2,dpml
real(kind=ip),dimension(4,4,im,jm)::A,B
real(kind=ip),dimension(im,jm)::temp,temp2

do i=ii,fi
   do j=ij,fj

        dpml(1,i,j) = A(1,1,i,j)*(dUdm(1,i,j)/meshx(i) + Beta*sigmax(i)*(U(1,i,j)-Ub(1,i,j)) - sigmax(i)*q1(1,i,j)) + &
                      A(1,2,i,j)*(dUdm(2,i,j)/meshx(i) + Beta*sigmax(i)*(U(2,i,j)-Ub(2,i,j)) - sigmax(i)*q1(2,i,j)) + &
                      B(1,1,i,j)*dUdn(1,i,j)                                                                        + &
                      B(1,3,i,j)*dUdn(3,i,j)                                                                        - & 
                      Bb(1,1,i,j)*(dUbdn(1,i,j))                                                                    - &
                      Bb(1,3,i,j)*(dUbdn(3,i,j))

        dpml(2,i,j) = A(2,2,i,j)*(dUdm(2,i,j)/meshx(i) + Beta*sigmax(i)*(U(2,i,j)-Ub(2,i,j)) - sigmax(i)*q1(2,i,j)) + &
                      A(2,4,i,j)*(dUdm(4,i,j)/meshx(i) + Beta*sigmax(i)*(U(4,i,j)-Ub(4,i,j)) - sigmax(i)*q1(4,i,j)) + &
                      B(2,2,i,j)*dUdn(2,i,j)                                                                        - & 
                      Bb(2,2,i,j)*(dUbdn(2,i,j))                                                                     

        dpml(3,i,j) = A(3,3,i,j)*(dUdm(3,i,j)/meshx(i) + Beta*sigmax(i)*(U(3,i,j)-Ub(3,i,j)) - sigmax(i)*q1(3,i,j)) + &
                      B(3,3,i,j)*dUdn(3,i,j)                                                                        + &
                      B(3,4,i,j)*dUdn(4,i,j)                                                                        - & 
                      Bb(3,3,i,j)*(dUbdn(3,i,j))                                                                    - &
                      Bb(3,4,i,j)*(dUbdn(4,i,j))
 
        dpml(4,i,j) = A(4,2,i,j)*(dUdm(2,i,j)/meshx(i) + Beta*sigmax(i)*(U(2,i,j)-Ub(2,i,j)) - sigmax(i)*q1(2,i,j)) + &
                      A(4,4,i,j)*(dUdm(4,i,j)/meshx(i) + Beta*sigmax(i)*(U(4,i,j)-Ub(4,i,j)) - sigmax(i)*q1(4,i,j)) + &
                      B(4,3,i,j)*dUdn(3,i,j)                                                                        + &
                      B(4,4,i,j)*dUdn(4,i,j)                                                                        - & 
                      Bb(4,3,i,j)*(dUbdn(3,i,j))                                                                    - &
                      Bb(4,4,i,j)*(dUbdn(4,i,j))
  
   end do
end do

end subroutine 

subroutine ylayer(U,A,B,q1,q2,dpml,ii,fi,ij,fj)
implicit none                     
integer::k,l,i,j,ii,fi,ij,fj
real(kind=ip),dimension(4,im,jm)::U,q1,q2,dpml
real(kind=ip),dimension(4,4,im,jm)::A,B
real(kind=ip),dimension(im,jm)::temp,temp2

  do i=ii,fi
    do j=ij,fj

       dpml(1,i,j) = A(1,1,i,j)*(dUdm(1,i,j))                                  + &
                     A(1,2,i,j)*(dUdm(2,i,j))                                  + &
                     B(1,1,i,j)*(dUdn(1,i,j)/meshy(j) -sigmay(j)*q2(1,i,j))    + &
                     B(1,3,i,j)*(dUdn(3,i,j)/meshy(j) -sigmay(j)*q2(3,i,j))    - & 
                     Bb(1,1,i,j)*(dUbdn(1,i,j)/meshy(j))                       - &
                     Bb(1,3,i,j)*(dUbdn(3,i,j)/meshy(j))

       dpml(2,i,j) = A(2,2,i,j)*(dUdm(2,i,j))                                  + &
                     A(2,4,i,j)*(dUdm(4,i,j))                                  + &
                     B(2,2,i,j)*(dUdn(2,i,j)/meshy(j) -sigmay(j)*q2(2,i,j))    - & 
                     Bb(2,2,i,j)*(dUbdn(2,i,j)/meshy(j))    

       dpml(3,i,j) = A(3,3,i,j)*(dUdm(3,i,j))                                  + &
                     B(3,3,i,j)*(dUdn(3,i,j)/meshy(j) -sigmay(j)*q2(3,i,j))    + &
                     B(3,4,i,j)*(dUdn(4,i,j)/meshy(j) -sigmay(j)*q2(4,i,j))    - & 
                     Bb(3,3,i,j)*(dUbdn(3,i,j)/meshy(j))                       - &
                     Bb(3,4,i,j)*(dUbdn(4,i,j)/meshy(j))

       dpml(4,i,j) = A(4,2,i,j)*(dUdm(2,i,j))                                  + &
                     A(4,4,i,j)*(dUdm(4,i,j))                                  + &
                     B(4,3,i,j)*(dUdn(3,i,j)/meshy(j) -sigmay(j)*q2(3,i,j))    + &
                     B(4,4,i,j)*(dUdn(4,i,j)/meshy(j) -sigmay(j)*q2(4,i,j))    - & 
                     Bb(4,3,i,j)*(dUbdn(3,i,j)/meshy(j))                       - &
                     Bb(4,4,i,j)*(dUbdn(4,i,j)/meshy(j))
    end do
  end do

end subroutine 

subroutine corner(U,A,B,q1,q2,dpml,ii,fi,ij,fj)
implicit none                     
integer::k,l,i,j,ii,fi,ij,fj
real(kind=ip),dimension(4,im,jm)::U,q1,q2,dpml
real(kind=ip),dimension(4,4,im,jm)::A,B
real(kind=ip),dimension(im,jm)::temp,temp2

do i=ii,fi
   do j=ij,fj

       dpml(1,i,j) = A(1,1,i,j)*(dUdm(1,i,j)/meshx(i) + Beta*sigmax(i)*(U(1,i,j)-Ub(1,i,j)) - sigmax(i)*q1(1,i,j)) + &
                     A(1,2,i,j)*(dUdm(2,i,j)/meshx(i) + Beta*sigmax(i)*(U(2,i,j)-Ub(2,i,j)) - sigmax(i)*q1(2,i,j)) + &
                     B(1,1,i,j)*(dUdn(1,i,j)/meshy(j) -sigmay(j)*q2(1,i,j))                                        + &
                     B(1,3,i,j)*(dUdn(3,i,j)/meshy(j) -sigmay(j)*q2(3,i,j))                                        - & 
                     Bb(1,1,i,j)*(dUbdn(1,i,j)/meshy(j))                                                           - &
                     Bb(1,3,i,j)*(dUbdn(3,i,j)/meshy(j))

       dpml(2,i,j) = A(2,2,i,j)*(dUdm(2,i,j)/meshx(i) + Beta*sigmax(i)*(U(2,i,j)-Ub(2,i,j)) - sigmax(i)*q1(2,i,j)) + &
                     A(2,4,i,j)*(dUdm(4,i,j)/meshx(i) + Beta*sigmax(i)*(U(4,i,j)-Ub(4,i,j)) - sigmax(i)*q1(4,i,j)) + &
                     B(2,2,i,j)*(dUdn(2,i,j)/meshy(j) -sigmay(j)*q2(2,i,j))                                        - & 
                     Bb(2,2,i,j)*(dUbdn(2,i,j)/meshy(j))                                                          

       dpml(3,i,j) = A(3,3,i,j)*(dUdm(3,i,j)/meshx(i) + Beta*sigmax(i)*(U(3,i,j)-Ub(3,i,j)) - sigmax(i)*q1(3,i,j)) + &
                     B(3,3,i,j)*(dUdn(3,i,j)/meshy(j) -sigmay(j)*q2(3,i,j))                                        + &
                     B(3,4,i,j)*(dUdn(4,i,j)/meshy(j) -sigmay(j)*q2(4,i,j))                                        - & 
                     Bb(3,3,i,j)*(dUbdn(3,i,j)/meshy(j))                                                           - &
                     Bb(3,4,i,j)*(dUbdn(4,i,j)/meshy(j))

       dpml(4,i,j) = A(4,2,i,j)*(dUdm(2,i,j)/meshx(i) + Beta*sigmax(i)*(U(2,i,j)-Ub(2,i,j)) - sigmax(i)*q1(2,i,j)) + &
                     A(4,4,i,j)*(dUdm(4,i,j)/meshx(i) + Beta*sigmax(i)*(U(4,i,j)-Ub(4,i,j)) - sigmax(i)*q1(4,i,j)) + &
                     B(4,3,i,j)*(dUdn(3,i,j)/meshy(j) -sigmay(j)*q2(3,i,j))                                        + &
                     B(4,4,i,j)*(dUdn(4,i,j)/meshy(j) -sigmay(j)*q2(4,i,j))                                        - & 
                     Bb(4,3,i,j)*(dUbdn(3,i,j)/meshy(j))                                                           - &
                     Bb(4,4,i,j)*(dUbdn(4,i,j)/meshy(j))

   end do
end do

end subroutine 
!*******************************MATRIZ**************************************
subroutine xlayermatriz(U,A,B,q1,q2,dpml,ii,fi,ij,fj)
implicit none                     
integer::k,l,i,j,ii,fi,ij,fj
real(kind=ip),dimension(4,im,jm)::U,q1,q2,dpml
real(kind=ip),dimension(4,4,im,jm)::A,B
real(kind=ip),dimension(im,jm)::temp,temp2

temp2(ii:fi,ij:fj)= 0.d0

do k=1,4
 do l=1,4
    do i=ii,fi
      do j=ij,fj
  
         temp(i,j) = A(k,l,i,j)*(dUdm(l,i,j)/meshx(i) & 
        &          + Beta*sigmax(i)*(U(l,i,j)-Ub(l,i,j))-sigmax(i)*q1(l,i,j)) &
        &          + B(k,l,i,j)*(dUdn(l,i,j))-Bb(k,l,i,j)*(dUbdn(l,i,j))
  
      end do
    end do
    temp2(ii:fi,ij:fj) = temp(ii:fi,ij:fj)+temp2(ii:fi,ij:fj)
  end do
  dpml(k,ii:fi,ij:fj) = temp2(ii:fi,ij:fj)
  temp2(ii:fi,ij:fj)  = 0.d0
end do


end subroutine 

subroutine ylayermatriz(U,A,B,q1,q2,dpml,ii,fi,ij,fj)
implicit none                     
integer::k,l,i,j,ii,fi,ij,fj
real(kind=ip),dimension(4,im,jm)::U,q1,q2,dpml
real(kind=ip),dimension(4,4,im,jm)::A,B
real(kind=ip),dimension(im,jm)::temp,temp2

temp2(ii:fi,ij:fj)= 0.d0

do k=1,4
 do l=1,4
  do i=ii,fi
    do j=ij,fj

       temp(i,j) = A(k,l,i,j)*(dUdm(l,i,j))&
      &          + B(k,l,i,j)*(dUdn(l,i,j)/meshy(j)-sigmay(j)*q2(l,i,j)) &
      &          - Bb(k,l,i,j)*(dUbdn(l,i,j)/meshy(j))

    end do
  end do
  temp2(ii:fi,ij:fj) = temp(ii:fi,ij:fj)+temp2(ii:fi,ij:fj)
 end do
 dpml(k,ii:fi,ij:fj) = temp2(ii:fi,ij:fj)
 temp2(ii:fi,ij:fj)  = 0.d0
end do

end subroutine 

subroutine cornermatriz(U,A,B,q1,q2,dpml,ii,fi,ij,fj)
implicit none                     
integer::k,l,i,j,ii,fi,ij,fj
real(kind=ip),dimension(4,im,jm)::U,q1,q2,dpml
real(kind=ip),dimension(4,4,im,jm)::A,B
real(kind=ip),dimension(im,jm)::temp,temp2

temp2(ii:fi,ij:fj)= 0.d0

do k=1,4
 do l=1,4
  do i=ii,fi
    do j=ij,fj

       temp(i,j) = A(k,l,i,j)*(dUdm(l,i,j)/meshx(i) & 
      &          + Beta*sigmax(i)*(U(l,i,j)-Ub(l,i,j))-sigmax(i)*q1(l,i,j))&
      &          + B(k,l,i,j)*(dUdn(l,i,j)/meshy(j)-sigmay(j)*q2(l,i,j)) & 
      &          - Bb(k,l,i,j)*(dUbdn(l,i,j)/meshy(j))

    end do
  end do
  temp2(ii:fi,ij:fj) = temp(ii:fi,ij:fj)+temp2(ii:fi,ij:fj)
 end do
 dpml(k,ii:fi,ij:fj) = temp2(ii:fi,ij:fj)
 temp2(ii:fi,ij:fj)  = 0.d0
end do

end subroutine 
end module
