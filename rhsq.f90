Module rhsq
use global
use global2
 contains

subroutine  rhsq1(u,q1,q1dt)

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,q1,q1dt

     
      call q1corner(U,q1,q1dt,1,imaxpml,1,jmaxpml)      
      call q1ylayer(U,q1,q1dt,imaxpml+1,imax-D,1,jmaxpml)     
      call q1corner(U,q1,q1dt,imax-D+1,imax,1,jmaxpml)    

      call q1xlayer(U,q1,q1dt,1,imaxpml,jmaxpml+1,jmax-D) 
      call q1xlayer(U,q1,q1dt,imax-D+1,imax,jmaxpml+1,jmax-D)

      call q1corner(U,q1,q1dt,1,imaxpml,jmax-D+1,jmax)  
      call q1ylayer(U,q1,q1dt,imaxpml+1,imax-D,jmax-D+1,jmax) 
      call q1corner(U,q1,q1dt,imax-D+1,imax,jmax-D+1,jmax)    

 end subroutine

subroutine  rhsq2(u,q2,q2dt)

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,q2,q2dt

     
      call q2corner(U,q2,q2dt,1,imaxpml,1,jmaxpml)      
      call q2ylayer(U,q2,q2dt,imaxpml+1,imax-D,1,jmaxpml)     
      call q2corner(U,q2,q2dt,imax-D+1,imax,1,jmaxpml)    

      call q2xlayer(U,q2,q2dt,1,imaxpml,jmaxpml+1,jmax-D) 
      call q2xlayer(U,q2,q2dt,imax-D+1,imax,jmaxpml+1,jmax-D)

      call q2corner(U,q2,q2dt,1,imaxpml,jmax-D+1,jmax)  
      call q2ylayer(U,q2,q2dt,imaxpml+1,imax-D,jmax-D+1,jmax) 
      call q2corner(U,q2,q2dt,imax-D+1,imax,jmax-D+1,jmax)    

end subroutine


subroutine q1xlayer(U,q1,q1dt,ii,fi,ij,fj)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,q1,q1dt

       do k=1,4
         do j = ij,fj
            do i = ii,fi

                q1dt(k,i,j) = dUdm(k,i,j)/meshx(i) &
                 &          + beta*sigmax(i)*(U(k,i,j)-Ub(k,i,j)) &
                 &          - q1(k,i,j)*sigmax(i)

            end do
         end do
       end do 

end subroutine


subroutine q1ylayer(U,q1,q1dt,ii,fi,ij,fj)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,q1,q1dt


       do k=1,4
          do j = ij,fj
             do i = ii,fi

                q1dt(k,i,j) = dUdm(k,i,j)   

             end do
          end do
       end do

end subroutine

subroutine q1corner(U,q1,q1dt,ii,fi,ij,fj)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,q1,q1dt
      real(kind=ip):: dt2


       do k=1,4
          do j = ij,fj
             do i = ii,fi

                q1dt(k,i,j) = dUdm(k,i,j)/meshx(i) &
                 &          + beta*sigmax(i)*(U(k,i,j)-Ub(k,i,j)) &
                 &          - q1(k,i,j)*sigmax(i)  

             end do
          end do
       end do
   
end subroutine

subroutine q2xlayer(U,q2,q2dt,ii,fi,ij,fj)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,q2,q2dt
      real(kind=ip):: dt2


       do k=1,4
          do j = ij,fj
             do i = ii,fi

                q2dt(k,i,j) = dUdn(k,i,j) - dUbdn(k,i,j)   

             end do
          end do
       end do
   
end subroutine

subroutine q2ylayer(U,q2,q2dt,ii,fi,ij,fj)  

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,q2,q2dt
      real(kind=ip):: dt2


       do k=1,4
          do j = ij,fj
             do i = ii,fi

                q2dt(k,i,j) = dUdn(k,i,j)/meshy(j) &
                 &          - dUbdn(k,i,j)/meshy(j)   &
                 &          - q2(k,i,j)*sigmay(j)  

             end do
          end do
       end do
   
end subroutine

subroutine q2corner(U,q2,q2dt,ii,fi,ij,fj)    

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,q2,q2dt
      real(kind=ip):: dt2


       do k=1,4
          do j = ij,fj
             do i = ii,fi

                q2dt(k,i,j) = dUdn(k,i,j)/meshy(j) &
                 &          - dUbdn(k,i,j)/meshy(j)   &
                 &          - q2(k,i,j)*sigmay(j)  

             end do
          end do
       end do
   
end subroutine

end module
