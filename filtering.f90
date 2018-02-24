module filtering

use global, only : ip,im,jm,imax,jmax,D
use globalq

contains


subroutine filtercompactx(U)

implicit none
integer::k,i,j,n
real(kind=ip),dimension(4,im,jm)::U
real(kind=ip),dimension(im)::a,b,c,g
real(kind=ip)::a0,a1,a2,a3,a4,a5,alphaf

alphaf = 0.250d0

a5 = (1.d0-2.d0*alphaf)/512.d0
a4 = -1.d0/128.d0+alphaf/64.d0-6.d0*a5
a3 = 1.d0/32.d0 - alphaf/16.d0-4.d0*a4-11.d0*a5
a2 = -1.d0/8.d0 + alphaf/4.d0-2.d0*a3-4.d0*a4-6.d0*a5
a1 = 0.5d0+alphaf-a3-a5
a0 = 0.5d0+alphaf-a2-a4



a=alphaf
b=1.d0
c=alphaf

a(1)    = 0.d0
b(1)    = 1.d0
c(1)    = 0.d0

a(2)    = 0.d0
b(2)    = 1.d0
c(2)    = 0.d0

a(3)    = 0.d0
b(3)    = 1.d0
c(3)    = 0.d0

a(4)    = 0.d0
b(4)    = 1.d0
c(4)    = 0.d0

a(imax)    = 0.d0
b(imax)    = 1.d0
c(imax)    = 0.d0

a(imax-1)    = 0.d0
b(imax-1)    = 1.d0
c(imax-1)    = 0.d0

a(imax-2)    = 0.d0
b(imax-2)    = 1.d0
c(imax-2)    = 0.d0

a(imax-3)    = 0.d0
b(imax-3)    = 1.d0
c(imax-3)    = 0.d0

do k=1,4
do j=2,jmax-1

g(1)    = U(k,1,j)
g(2)    = U(k,2,j)!a0*U(k,2,j)+a1/2.d0*(U(k,3,j)+U(k,1,j))
g(3)    = U(k,3,j)!a0*u(1,3,j)+a1/2.d0*(U(k,4,j)+U(k,2,j)) + a2/2.d0*(U(k,5,j)+U(k,1,j))
g(4)    = U(k,4,j)!a0*u(1,4,j)+a1/2.d0*(U(k,5,j)+U(k,3,j)) + a2/2.d0*(U(k,6,j)+U(k,2,j)) + a3/2.d0*(U(k,7,j)+U(k,1,j))

g(imax)      = U(k,imax,j)
g(imax-1)    = U(k,imax-1,j)!a0*U(k,imax-1,j)+a1/2.d0*(U(k,imax,j)  +U(k,imax-2,j))
g(imax-2)    = U(k,imax-2,j)!a0*U(k,imax-2,j)+a1/2.d0*(U(k,imax-1,j)+U(k,imax-3,j)) + a2/2.d0*(U(k,imax,j)  +U(k,imax-4,j))
g(imax-3)    = U(k,imax-3,j)!a0*U(k,imax-3,j)+a1/2.d0*(U(k,imax-2,j)+U(k,imax-4,j)) + a2/2.d0*(U(k,imax-1,j)+U(k,imax-5,j))&
               !&  + a3/2.d0*(U(k,imax,j)+U(k,imax-6,j))

do i=5,imax-4

  g(i)   =  a0*U(k,i,j)+a1/2.d0*(U(k,i+1,j)+U(k,i-1,j))+a2/2.d0*(U(k,i+2,j)+U(k,i-2,j))&
        &  +a3/2.d0*(U(k,i+3,j)+U(k,i-3,j))+a4/2.d0*(U(k,i+4,j)+U(k,i-4,j))

end do

call thomas(a,b,c,g,U(k,:,j),imax)

end do
end do

end subroutine

subroutine filtercompacty(U)

implicit none
integer::k,i,j,n
real(kind=ip),dimension(4,im,jm)::U
real(kind=ip),dimension(jm)::a,b,c,g
real(kind=ip)::a0,a1,a2,a3,a4,a5,alphaf

alphaf = 0.250d0

a5 = (1.d0-2.d0*alphaf)/512.d0
a4 = -1.d0/128.d0+alphaf/64.d0-6.d0*a5
a3 = 1.d0/32.d0 - alphaf/16.d0-4.d0*a4-11.d0*a5
a2 = -1.d0/8.d0 + alphaf/4.d0-2.d0*a3-4.d0*a4-6.d0*a5
a1 = 0.5d0+alphaf-a3-a5
a0 = 0.5d0+alphaf-a2-a4



a=alphaf
b=1.d0
c=alphaf

a(1)    = 0.d0
b(1)    = 1.d0
c(1)    = 0.d0

a(2)    = 0.d0
b(2)    = 1.d0
c(2)    = 0.d0

a(3)    = 0.d0
b(3)    = 1.d0
c(3)    = 0.d0

a(4)    = 0.d0
b(4)    = 1.d0
c(4)    = 0.d0

a(jmax)    = 0.d0
b(jmax)    = 1.d0
c(jmax)    = 0.d0

a(jmax-1)    = 0.d0
b(jmax-1)    = 1.d0
c(jmax-1)    = 0.d0

a(jmax-2)    = 0.d0
b(jmax-2)    = 1.d0
c(jmax-2)    = 0.d0

a(jmax-3)    = 0.d0
b(jmax-3)    = 1.d0
c(jmax-3)    = 0.d0

do k=1,4
do i=2,imax-1

g(1)    = U(k,i,1)
g(2)    = U(k,i,2)!a0*U(k,2,j)+a1/2.d0*(U(k,3,j)+U(k,1,j))
g(3)    = U(k,i,3)!a0*u(1,3,j)+a1/2.d0*(U(k,4,j)+U(k,2,j)) + a2/2.d0*(U(k,5,j)+U(k,1,j))
g(4)    = U(k,i,4)!a0*u(1,4,j)+a1/2.d0*(U(k,5,j)+U(k,3,j)) + a2/2.d0*(U(k,6,j)+U(k,2,j)) + a3/2.d0*(U(k,7,j)+U(k,1,j))

g(jmax)      = U(k,i,jmax)
g(jmax-1)    = U(k,i,jmax-1)!a0*U(k,jmax-1,j)+a1/2.d0*(U(k,jmax,j)  +U(k,jmax-2,j))
g(jmax-2)    = U(k,i,jmax-2)!a0*U(k,jmax-2,j)+a1/2.d0*(U(k,jmax-1,j)+U(k,jmax-3,j)) + a2/2.d0*(U(k,jmax,j)  +U(k,jmax-4,j))
g(jmax-3)    = U(k,i,jmax-3)!a0*U(k,jmax-3,j)+a1/2.d0*(U(k,jmax-2,j)+U(k,jmax-4,j)) + a2/2.d0*(U(k,jmax-1,j)+U(k,jmax-5,j))&
               !&  + a3/2.d0*(U(k,jmax,j)+U(k,jmax-6,j))

do j=5,jmax-4

  g(j)   =  a0*U(k,i,j)+a1/2.d0*(U(k,i,j+1)+U(k,i,j-1))+a2/2.d0*(U(k,i,j+2)+U(k,i,j-2))&
        &  +a3/2.d0*(U(k,i,j+3)+U(k,i,j-3))+a4/2.d0*(U(k,i,j+4)+U(k,i,j-4))

end do

call thomas(a,b,c,g,U(k,i,:),jmax)

end do
end do



end subroutine

subroutine filtercompact2(U)

implicit none
integer::k,i,j,n
real(kind=ip),dimension(4,im,jm)::U
real(kind=ip),dimension(im)::a,b,c,g
real(kind=ip)::a0,a1,a2,a3,a4,a5,alphaf

alphaf = 0.25d0

a5 = (1.d0-2.d0*alphaf)/512.d0
a4 = -1.d0/128.d0+alphaf/64.d0-6.d0*a5
a3 = 1.d0/32.d0 - alphaf/16.d0-4*a4-11.d0*a5
a2 = -1.d0/8.d0 + alphaf/4.d0-2.d0*a3-4.d0*a4-6.d0*a5
a1 = 0.5d0+alphaf-a3-a5
a0 = 0.5d0+alphaf-a2-a4



a=alphaf
b=1.d0
c=alphaf

a(1)    = 0.d0
b(1)    = 1.d0
c(1)    = 0.d0

a(imax)    = 0.d0
b(imax)    = 1.d0
c(imax)    = 0.d0


!do k=1,4
do j=2,jmax-1

g(1)    = U(1,1,j)
g(2)    = a0*U(1,2,j)+a1/2.d0*(U(1,3,j)+U(1,1,j))
g(3)    = a0*u(1,3,j)+a1/2.d0*(U(1,4,j)+U(1,2,j)) + a2/2.d0*(U(1,5,j)+U(1,1,j))
g(4)    = a0*u(1,4,j)+a1/2.d0*(U(1,5,j)+U(1,3,j)) + a2/2.d0*(U(1,6,j)+U(1,2,j)) + a3/2.d0*(U(1,7,j)+U(1,1,j))

g(imax)      = U(1,imax,j)
g(imax-1)    = a0*U(1,imax-1,j)+a1/2.d0*(U(1,imax,j)  +U(1,imax-2,j))
g(imax-2)    = a0*U(1,imax-2,j)+a1/2.d0*(U(1,imax-1,j)+U(1,imax-3,j)) + a2/2.d0*(U(1,imax,j)  +U(1,imax-4,j))
g(imax-3)    = a0*U(1,imax-3,j)+a1/2.d0*(U(1,imax-2,j)+U(1,imax-4,j)) + a2/2.d0*(U(1,imax-1,j)+U(1,imax-5,j))&
            &  + a3/2.d0*(U(1,imax,j)+U(1,imax-6,j))

do i=5,imax-4

  g(i)   =  a0*U(1,i,j)+a1/2.d0*(U(1,i+1,j)+U(1,i-1,j))+a2/2.d0*(U(1,i+2,j)+U(1,i-2,j))&
        &  +a3/2.d0*(U(1,i+3,j)+U(1,i-3,j))+a4/2.d0*(U(1,i+4,j)+U(1,i-4,j))

end do

call thomas(a,b,c,g,U(1,:,j),imax)

end do
!end do


end subroutine

subroutine filtercompact(U)

implicit none
integer::k,i,j,n
real(kind=ip),dimension(4,im,jm)::U
real(kind=ip),dimension(im)::a,b,c,d
real(kind=ip)::a1,b1,c1,d1,alpha,beta

alpha = 4.75d0
beta  = 0.d0

!a1 = 1.d0/4.d0*(2.d0+3.d0*alpha) 
!b1 = 1.d0/16.d0*(9.d0+16.d0*alpha+10.d0*beta)
!c1 = 1.d0/4.d0*(alpha+4.d0*beta)
!d1 = 0.d0!1.d0/16.d0*(6*beta-1)

!a1 = 1.d0/8.d0*(5.d0+6.d0*alpha-6.d0*beta +16.d0*d)
!b1 = 1.d0/2.d0*(1.d0+2.d0*alpha+2.d0*beta-2.d0*d)
!c1 = -1.d0/8.d0*(1.d0-2.d0*alpha-14.d0*beta+16.d0*d)
!d1 = 0.d0!1.d0/16.d0*(6*beta-1)

a1 = 1.d0/16.d0*(11.d0+10.d0*alpha-10.d0*beta)
b1 = 1.d0/32.d0*(15.d0+34.d0*alpha+30.d0*beta)
c1 = 1.d0/16.d0*(-3.d0+ 6.d0*alpha+26.d0*beta)
d1 = 0.d0!1.d0/16.d0*(6*beta-1)

do j=1,jmax
do i=4,imax-3

  a(i)   = alpha
  b(i)   = 1.d0
  c(i)   = alpha
    
  d(i)   = a1*U(1,i,j) + c1/2.d0*(U(1,i+2,j)+U(1,i-2,j)) + b1/2.d0*(U(1,i+1,j) + U(1,i-1,j))
write(*,*)a(i),b(i),c(i),d(i)
end do

call thomas(a(4:imax-3),b(4:imax-3),c(4:imax-3),d(4:imax-3),U(1,4:imax-3,j),imax-6)
write(*,*)U(1,4:imax-3,j)

end do


U(1,1,:) = 15.d0/16.d0*U(1,1,:)+ 1.d0/16.d0*(4.d0*U(1,2,:)-6.d0*U(1,3,:)+4.d0*U(1,4,:)-U(1,5,:))
U(1,2,:) = 3.d0/4.d0*U(1,2,:)  + 1.d0/16.d0*(U(1,1,:)+6.d0*U(1,3,:)-4.d0*U(1,4,:)+U(1,5,:))
U(1,3,:) = 5.d0/8.d0*U(1,3,:)  + 1.d0/16.d0*(-U(1,1,:)+4.d0*U(1,2,:)+4.d0*U(1,4,:)-U(1,5,:))

U(1,imax,:)   = 15.d0/16.d0*U(1,imax,:)  + 1.d0/16.d0*(4.d0*U(1,imax-1,:)-6.d0*U(1,imax-2,:)+4.d0*U(1,imax-3,:)-U(1,imax-4,:))
U(1,imax-1,:) = 3.d0/4.d0*U(1,imax-1,:)  + 1.d0/16.d0*(U(1,imax,:)+6.d0*U(1,imax-2,:)-4.d0*U(1,imax-3,:)+U(1,imax-4,:))
U(1,imax-2,:) = 5.d0/8.d0*U(1,imax-2,:)  + 1.d0/16.d0*(-U(1,imax,:)+4.d0*U(1,imax-1,:)+4.d0*U(1,imax-3,:)-U(1,imax-4,:))

end subroutine

subroutine thomas(a,b,c,d,x,n)

implicit none
integer::i,n
real(kind=ip),dimension(n)::a,b,c,d
real(kind=ip),dimension(n)::cp,dp
real(kind=ip),dimension(n)::x
real(kind=ip)::m

cp(1)=c(1)/b(1)
dp(1)=d(1)/b(1)

do i=2,n

        m      = b(i)-cp(i-1)*a(i)        
        cp(i)  = c(i)/m
        dp(i)  = (d(i)-dp(i-1)*a(i))/m

end do
   
x(n)=dp(n)

do i=n-1,1,-1

     x(i)=dp(i)-cp(i)*x(i+1)
  
end do

end subroutine

subroutine filter(U)

integer::i,j,k
integer,dimension(11)::l
real(kind=ip),dimension(4,im,jm)::U
real(kind=ip),dimension(11)::F
real(kind=ip)::suma1,suma2,suma3,suma4
real(kind=ip)::sumaq11,sumaq12,sumaq13,sumaq14
real(kind=ip)::sumaq21,sumaq22,sumaq23,sumaq24

!D(1)  = -1.d0/1024.d0
!D(2)  =  5.d0/512.d0
!D(3)  = -45.d0/1024.d0
!D(4)  =  15.d0/128.d0
!D(5)  = -105.d0/512.d0
!D(6)  =  63.d0/512.d0
!D(7)  = -105.d0/512.d0
!D(8)  =  15.d0/128.d0
!D(9)  = -45.d0/1024.d0
!D(10) =  5.d0/512.d0
!D(11) = -1.d0/1024.d0

F(1)  = -1.d0/64.d0
F(2)  =  3.d0/32.d0
F(3)  = -15.d0/64.d0
F(4)  =  5.d0/14.d0
F(5)  = -15.d0/64.d0
F(6)  =  3.d0/32.d0
F(7)  = -1.d0/64.d0

do k=1, 7
  l(k)=-3+(k-1)
end do 

suma1=0.d0
suma2=0.d0
suma3=0.d0
suma4=0.d0

sumaq11=0.d0
sumaq12=0.d0
sumaq13=0.d0
sumaq14=0.d0

sumaq21=0.d0
sumaq22=0.d0
sumaq23=0.d0
sumaq24=0.d0


do j=jmaxpml+2,jmax-D-1
      do i=imaxpml+10,imax-D-10
      do k=1,7 

           suma1   = suma1+F(k)*U(1,i+l(k),j)
           suma2   = suma2+F(k)*U(2,i+l(k),j)
           suma3   = suma3+F(k)*U(3,i+l(k),j)
           suma4   = suma4+F(k)*U(4,i+l(k),j)

!           sumaq11 = sumaq11+F(k)*q1(1,i+l(k),j)
!           sumaq12 = sumaq12+F(k)*q1(2,i+l(k),j)
!           sumaq13 = sumaq13+F(k)*q1(3,i+l(k),j)
!           sumaq14 = sumaq14+F(k)*q1(4,i+l(k),j)
!
!           sumaq21 = sumaq21+F(k)*q2(1,i+l(k),j)
!           sumaq22 = sumaq22+F(k)*q2(2,i+l(k),j)
!           sumaq23 = sumaq23+F(k)*q2(3,i+l(k),j)
!           sumaq24 = sumaq24+F(k)*q2(4,i+l(k),j)
 
       end do

     U(1,i,j)     =   U(1,i,j)-suma1  
     U(2,i,j)     =   U(2,i,j)-suma2  
     U(3,i,j)     =   U(3,i,j)-suma3
     U(4,i,j)     =   U(4,i,j)-suma4  

!!   q1(1,i,j)   =   q1(1,i,j)-sumaq11
!!   q1(2,i,j)   =   q1(2,i,j)-sumaq12
!!   q1(3,i,j)   =   q1(3,i,j)-sumaq13
!!   q1(4,i,j)   =   q1(4,i,j)-sumaq14
!!            
!!   q2(1,i,j)   =   q2(1,i,j)-sumaq21
!!   q2(2,i,j)   =   q2(2,i,j)-sumaq22
!!   q2(3,i,j)   =   q2(3,i,j)-sumaq23
!!   q2(4,i,j)   =   q2(4,i,j)-sumaq24

   suma1=0.d0
   suma2=0.d0
   suma3=0.d0
   suma4=0.d0

!   sumaq11=0.d0
!   sumaq12=0.d0
!   sumaq13=0.d0
!   sumaq14=0.d0
!
!   sumaq21=0.d0
!   sumaq22=0.d0
!   sumaq23=0.d0
!   sumaq24=0.d0

   end do
end do

do i=imaxpml+2,imax-D-1
   do j=jmaxpml+10,jmax-D-10
      do k=1,7 

           suma1   = suma1+F(k)*U(1,i,j+l(k))
           suma2   = suma2+F(k)*U(2,i,j+l(k))
           suma3   = suma3+F(k)*U(3,i,j+l(k))
           suma4   = suma4+F(k)*U(4,i,j+l(k))

!           sumaq11 = sumaq11+F(k)*q1(1,i+l(k),j)
!           sumaq12 = sumaq12+F(k)*q1(2,i+l(k),j)
!           sumaq13 = sumaq13+F(k)*q1(3,i+l(k),j)
!           sumaq14 = sumaq14+F(k)*q1(4,i+l(k),j)
!
!           sumaq21 = sumaq21+F(k)*q2(1,i+l(k),j)
!           sumaq22 = sumaq22+F(k)*q2(2,i+l(k),j)
!           sumaq23 = sumaq23+F(k)*q2(3,i+l(k),j)
!           sumaq24 = sumaq24+F(k)*q2(4,i+l(k),j)
 
       end do

   U(1,i,j)     =   U(1,i,j)-suma1  
   U(2,i,j)     =   U(2,i,j)-suma2  
   U(3,i,j)     =   U(3,i,j)-suma3
   U(4,i,j)     =   U(4,i,j)-suma4  

!!   q1(1,i,j)   =   q1(1,i,j)-sumaq11
!!   q1(2,i,j)   =   q1(2,i,j)-sumaq12
!!   q1(3,i,j)   =   q1(3,i,j)-sumaq13
!!   q1(4,i,j)   =   q1(4,i,j)-sumaq14
!!            
!!   q2(1,i,j)   =   q2(1,i,j)-sumaq21
!!   q2(2,i,j)   =   q2(2,i,j)-sumaq22
!!   q2(3,i,j)   =   q2(3,i,j)-sumaq23
!!   q2(4,i,j)   =   q2(4,i,j)-sumaq24

   suma1=0.d0
   suma2=0.d0
   suma3=0.d0
   suma4=0.d0

!   sumaq11=0.d0
!   sumaq12=0.d0
!   sumaq13=0.d0
!   sumaq14=0.d0
!
!   sumaq21=0.d0
!   sumaq22=0.d0
!   sumaq23=0.d0
!   sumaq24=0.d0

   end do
end do
end subroutine 


end module
