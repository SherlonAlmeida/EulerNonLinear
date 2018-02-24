module global

      implicit none
      integer::im,jm,ip
      parameter(im=441,jm=441,ip=10)      
      integer:: maxit
      real(kind=ip):: mmin,mmax,nmin,nmax
      integer::imax,jmax
      real(kind=ip):: dm, dn, dt,pi,MA

      real(kind=ip),dimension(im):: meshx
      real(kind=ip),dimension(jm):: meshy
      real(kind=ip),dimension(jm):: uba,rhob,Tb
    
      integer::impml,jmpml
      parameter(impml=20,jmpml=20)
      integer::imaxpml,jmaxpml,D

      real(kind=ip):: xli,xld,yli,yls,sigmamx,delta,sigmamy,beta,sigmamxc,sigmamyc,c0

      real(kind=ip),dimension(im):: sigmax
      real(kind=ip),dimension(im):: sigmaxc
      real(kind=ip),dimension(jm):: sigmay
      real(kind=ip),dimension(jm):: sigmayc

      real(kind=ip),dimension(im,jm):: s
      real(kind=ip)::w,r0
   
      real(kind=ip),dimension(im):: m
      real(kind=ip),dimension(jm):: n

      real(kind=ip):: lamda,gama
      parameter(lamda=1.4d0)

      real(kind=ip),dimension(4,im,jm)::Ub,dUbdn

      real(kind=ip),dimension(4,4,im,jm)::Ab,Bb

end module



