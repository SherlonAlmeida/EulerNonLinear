module initialize
use global
use globalq
use diff
contains

subroutine init(U)
implicit none

integer::i,j

	real(kind=ip),dimension(4,im,jm)::U
	real(kind=ip):: u1b,u2b,T1,T2,gama,alpha,lamda 
	real(kind=ip):: Dx,Dy 
	real(kind=ip):: A,alpham 
        real(kind=ip),dimension(4,4,im,jm)::Ab,Bb
        real(kind=ip),dimension(jm)::drhobdn,dubadn

!................................ input data ..........................
namelist /grid/ imax, jmax, maxit,dt
namelist /freeadm/  mmin, mmax, nmin, nmax
namelist /pmlgrid/ imaxpml,jmaxpml,D,sigmamx,sigmamy,sigmamxc,sigmamyc,alpha,c0,delta

!................................ read input data ......................
open(unit=1,file='dados.in',status='old')

!................................ grid .................................
read (1,nml=grid)
write(*,nml=grid) 

!................................ free stream conditions ...............
read (1,nml=freeadm)
write(*,nml=freeadm)

!.................................PML GRIDe..............................
read (1,nml=pmlgrid)
write(*,nml=pmlgrid)
close(unit=1)

!.......................................................................
 dm       =    (mmax-mmin)/real(imax-1)
 dn       =    (nmax-nmin)/real(jmax-1)
 write(*,*)'dx=',dm
 write(*,*)'dy=',dn

! Definição do sigma 
 Dx       =    dm*D
 Dy       =    dn*D
! Começo da Pml ezquerda,direita,inferior,superior
 xli= mmin+Dx
 xld= mmax-Dx
 yli= nmin+Dy
 yls= nmax-Dy
!***********************************************

do i = 1,imax
   m(i)  = mmin + real(i-1)*dm
end do

do j = 1,jmax
   n(j) = nmin + real(j-1)*dn
end do

!Grid compress
!***********x
A      = 2.d0
alpham = 2.d0
meshx  = 1.d0
do i = imaxpml+1,1,-1
   meshx(i) =1.d0+A*(abs((m(i)-(xli))/(Dx)))**alpham
end do
do i = imax-D,imax
   meshx(i) =1.d0+A*(abs((m(i)-xld)/(Dx)))**alpham
end do
!*********y
meshy   = 1.d0
do j =jmaxpml+1,1,-1
   meshy(j) =1.d0+A*(abs((n(j)-(yli))/(Dy)))**alpham
end do
do j = jmax-D,jmax
   meshy(j) =1.d0+A*(abs((n(j)-yls)/(Dy)))**alpham
end do
!*******************

!SIGMAS
do i = imaxpml,1,-1
   sigmax(i) =sigmamx*(abs((m(i)-(xli))/(Dx)))**alpha
end do

do i = imax-D+1,imax
   sigmax(i) =sigmamx*(abs((m(i)-xld)/(Dx)))**alpha
end do

do i =imaxpml,1,-1
   sigmaxc(i) =sigmamxc*(abs((m(i)-(xli))/Dx))**alpha
end do

do i = imax-D+1,imax
   sigmaxc(i) =sigmamxc*(abs((m(i)-xld)/(Dx)))**alpha
end do


do j =jmaxpml,1,-1
   sigmay(j) =sigmamy*(abs((n(j)-(yli))/(Dy)))**alpha
end do

do j = jmax-D+1,jmax
   sigmay(j) =sigmamy*(abs((n(j)-yls)/(Dy)))**alpha
end do

do j = jmaxpml,1,-1
   sigmayc(j) =sigmamyc*(abs((n(j)-(yli))/(Dy)))**alpha
end do

do j = jmax-D+1,jmax
   sigmayc(j) =sigmamyc*(abs((n(j)-yls)/(Dy)))**alpha
end do

beta=-1.d0/c0

u1b  = 0.8d0
u2b  = 0.2d0
T1   = 1.0d0
T2   = 0.8d0
gama = 0.4d0
lamda= 1.4d0

do j=1,jmax
   uba(j)   = 0.5d0*(u1b+u2b+(u1b-u2b)*tanh(2.d0*n(j)/gama))
!  uba(j)   = 0.5d0*(u1b+u2b+(u1b-u2b)*tanh(2.d0*n(j)/gama))-0.5d0*(1.d0-(tanh((2.d0*n(j)/gama)))**2)
   Tb(j)   = T1*(uba(j)-u2b)/(u1b-u2b)+T2*(u1b-uba(j))/(u1b-u2b)+(lamda-1.d0)/2.d0*(u1b-Uba(j))*(uba(j)-u2b)
   rhob(j) = 1.d0/(Tb(j))
   !write(99999,*)n(j),uba(j) !Sherlon: Comentei isso pois nao tinha no com MPI
end do
 !  write(*,*)n(1),n(jmax),uba(1),uba(jmax)
!............................zerando as variaveis......................

!U  = 0.d0
q1 = 0.d0
q2 = 0.d0

!........................... inicialização ............................
do i=1,imax
  do j=1,jmax

   U(1,i,j)   = rhob(j)   
   U(2,i,j)   = uba(j)    
   U(3,i,j)   = 0.d0      
   U(4,i,j)   = 1.d0/lamda

  end do 
end do

do i=1,imax
  do j=1,jmax

   Ub(1,i,j)   = rhob(j) 
   Ub(2,i,j)   = uba(j)
   Ub(3,i,j)   = 0.d0
   Ub(4,i,j)   = 1.d0/lamda 

  end do 
end do

call dernb(dubadn,uba)
call dernb(drhobdn,rhob)

do i=1,imax
  do j=1,jmax

   dUbdn(1,i,j)   = drhobdn(j) 
   dUbdn(2,i,j)   = dubadn(j)
   dUbdn(3,i,j)   = 0.d0
   dUbdn(4,i,j)   = 0.d0 

  end do 
end do

Bb(1,1,:,:)  = Ub(3,:,:)
Bb(1,2,:,:)  = 0.d0
Bb(1,3,:,:)  = Ub(1,:,:)
Bb(1,4,:,:)  = 0.d0

Bb(2,1,:,:)  = 0.d0 
Bb(2,2,:,:)  = Ub(3,:,:)
Bb(2,3,:,:)  = 0.d0
Bb(2,4,:,:)  = 0.d0

Bb(3,1,:,:)  = 0.d0 
Bb(3,2,:,:)  = 0.d0
Bb(3,3,:,:)  = Ub(3,:,:)
Bb(3,4,:,:)  = 1.d0/Ub(1,:,:)

Bb(4,1,:,:)  = 0.d0 
Bb(4,2,:,:)  = 0.d0
Bb(4,3,:,:)  = lamda*Ub(4,:,:)
Bb(4,4,:,:)  = Ub(3,:,:)

end subroutine

end module
