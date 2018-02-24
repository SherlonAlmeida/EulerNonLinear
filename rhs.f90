module rhspml

use global
use derivs
use eqeuler
use eqpml

contains   

subroutine rhs(dUdt,U,q1,q2)

implicit none
     
integer::k,i,j
real(kind=ip),dimension(4,im,jm)::dUdt,U,q1,q2
real(kind=ip),dimension(4,im,jm)::dE,dpml
real(kind=ip),dimension(4,4,im,jm)::A,B

A(1,1,:,:)  = U(2,:,:)
A(1,2,:,:)  = U(1,:,:)
A(1,3,:,:)  = 0.0d0
A(1,4,:,:)  = 0.0d0

A(2,1,:,:)  = 0.d0 
A(2,2,:,:)  = U(2,:,:)
A(2,3,:,:)  = 0.d0
A(2,4,:,:)  = 1.d0/U(1,:,:)

A(3,1,:,:)  = 0.d0 
A(3,2,:,:)  = 0.d0
A(3,3,:,:)  = U(2,:,:)
A(3,4,:,:)  = 0.d0

A(4,1,:,:)  = 0.d0 
A(4,2,:,:)  = lamda*U(4,:,:)
A(4,3,:,:)  = 0.d0
A(4,4,:,:)  = U(2,:,:)


B(1,1,:,:)  = U(3,:,:)
B(1,2,:,:)  = 0.d0
B(1,3,:,:)  = U(1,:,:)
B(1,4,:,:)  = 0.d0

B(2,1,:,:)  = 0.d0 
B(2,2,:,:)  = U(3,:,:)
B(2,3,:,:)  = 0.d0
B(2,4,:,:)  = 0.d0

B(3,1,:,:)  = 0.d0 
B(3,2,:,:)  = 0.d0
B(3,3,:,:)  = U(3,:,:)
B(3,4,:,:)  = 1.d0/U(1,:,:)

B(4,1,:,:)  = 0.d0 
B(4,2,:,:)  = 0.d0
B(4,3,:,:)  = lamda*U(4,:,:)
B(4,4,:,:)  = U(3,:,:)

!........................ derivadas das variaveis dependentes ..........
   call deriv(U,q1,q2)

!........................ RHS da eq de Euler ...........................
   call nonlinear(U,dE,A,B) 

!........................ RHS da eq com PML ............................
   call pml(U,q1,q2,dpml,A,B)

!.......................................................................
do k=1,4
  do i=1,imax
    do j=1,jmax

        if(((i>imaxpml).and.(i<imax-D+1)).and.((j>jmaxpml).and.(j<jmax-D+1)))then

            if(k<4)then
               dUdt(k,i,j)  = - dE(k,i,j)  
            else 
               dUdt(k,i,j)  = - dE(k,i,j) + s(i,j) 
            end if 
        else 

            if(k<4)then
               dUdt(k,i,j)  = - dpml(k,i,j)
            else 
               dUdt(k,i,j)  = - dpml(k,i,j) + s(i,j)
            end if 

       end if    

    end do
  end do 
end do
end subroutine
      
end module
