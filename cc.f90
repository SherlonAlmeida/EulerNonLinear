module cc
use global
contains
!************Boundary conditions***********

subroutine ccpml(U)

implicit none
real(kind=ip),dimension(4,im,jm)::U 

!  rho=G(1,:,:)
!  u=G(2,:,:)
!  v=G(3,:,:)
!  P=G(4,:,:)

!Inflow Boundary!!!!!!!

    U(1,1,:)     = 0.d0
    U(1,1,:)     = rhob(:)

    U(2,1,:)     = 2.d0*u(2,2,:)-u(2,3,:)
    U(3,1,:)     = 0.d0

    U(4,1,:)     = 0.d0 
    U(4,1,:)     = 1.d0/lamda


! outflow boundary

   U(1,imax,:)   = 2.d0*U(1,imax-1,:)-U(1,imax-2,:)
   U(2,imax,:)   = 2.d0*U(2,imax-1,:)-U(2,imax-2,:)
   U(3,imax,:)   = 2.d0*U(3,imax-1,:)-U(3,imax-2,:)

   U(4,imax,:)   = 0.d0 
   U(4,imax,:)   = 1.d0/lamda

   
! upper boundary 
 
   U(1,:,jmax)   = 2.d0*U(1,:,jmax-1)-U(1,:,jmax-2)
   U(2,:,jmax)   = 2.d0*U(2,:,jmax-1)-U(2,:,jmax-2) 
   U(3,:,jmax)   = 2.d0*U(3,:,jmax-1)-U(3,:,jmax-2)

   U(4,:,jmax)   = 0.d0 
   U(4,:,jmax)   =  1.d0/lamda
 

! lower boundary

    U(1,:,1)     = 2.d0*U(1,:,2)-U(1,:,3)
    U(2,:,1)     = 2.d0*U(2,:,2)-U(2,:,3)
    U(3,:,1)     = 2.d0*U(3,:,2)-U(3,:,3)

    U(4,:,1)     = 0.d0 
    U(4,:,1)     = 1.d0/lamda
    
   
end subroutine     

!subroutine ccq(q1,q2,q3,q5)
!
!implicit none
!
!real(kind=ip),dimension(impml,jmpml)::q1,q2,q3,q5
!
!!Inflow Boundary!!!!!!!
!
!  q1(1,:)   = 2.d0*q1(2,:)-q1(3,:)
!  q2(1,:)   = 2.d0*q2(2,:)-q2(3,:)
!  q3(1,:)   = 2.d0*q3(2,:)-q3(3,:)
!  q5(1,:)   = 2.d0*q5(2,:)-q5(3,:)
!
!
!! oq2tflow boq2ndary
!
!   q1(imaxpml,:)      = 2.d0*q1(imaxpml-1,:)-q1(imaxpml-2,:) 
!   q2(imaxpml,:)      = 2.d0*q2(imaxpml-1,:)-q2(imaxpml-2,:)
!   q3(imaxpml,:)      = 2.d0*q3(imaxpml-1,:)-q3(imaxpml-2,:)
!   q5(imaxpml,:)      = 2.d0*q5(imaxpml-1,:)-q5(imaxpml-2,:)
!
!   
!! q2pper boq2ndary 
!
!  q1(:,jmaxpml)      = 2.d0*q1(:,jmaxpml-1)-q1(:,jmaxpml-2)
!  q2(:,jmaxpml)      = 2.d0*q2(:,jmaxpml-1)-q2(:,jmaxpml-2)
!  q3(:,jmaxpml)      = 2.d0*q3(:,jmaxpml-1)-q3(:,jmaxpml-2)
!  q5(:,jmaxpml)      = 2.d0*q5(:,jmaxpml-1)-q5(:,jmaxpml-2)
!
!
!! lower boq2ndary
!
!  q1(:,1)       = 2.d0*q1(:,2)-q1(:,3)
!  q2(:,1)       = 2.d0*q2(:,2)-q2(:,3)
!  q3(:,1)       = 2.d0*q3(:,2)-q3(:,3)
!  q5(:,1)       = 2.d0*q5(:,2)-q2(:,3)
!
!    
!   
!end subroutine     
end module
