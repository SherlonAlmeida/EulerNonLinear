module outt
use global
use diff
use global2
contains
 
subroutine output(U,icount)

      implicit none
      integer::i,j,icount

      real(kind=ip),dimension(4,im,jm):: U

      real(kind=ip):: mmin, mmax

      call derm(dudm,u,1,imax,1,jmax) 
      call dern(dudn,u,1,imax,1,jmax) 

      do j=1,jmax
         do i=1,imax
  
              write(09+icount,*) m(i),n(j),U(1,i,j),U(2,i,j),U(3,i,j),U(4,i,j),dUdm(3,i,j)-dudn(2,i,j)
              !write(09+icount,*) m(i),n(j),U(1,i,j),U(2,i,j),U(3,i,j),U(4,i,j)
          
          enddo
         write(09+icount,*)
      enddo
end subroutine

end module
