module derivs

use global
use global2
use diff

contains

subroutine deriv(U,q1,q2)

  real(kind=ip),dimension(4,im,jm):: U   
  real(kind=ip),dimension(4,im,jm):: q1,q2


  call derm(dUdm,U,1,imax,1,jmax) 
  call dern(dUdn,U,1,imax,1,jmax) 

  call derm(dq1dm,q1,1,imaxpml,1,jmax)               
  call dern(dq1dn,q1,1,imaxpml,1,jmax)               
  call derm(dq2dm,q2,1,imaxpml,1,jmax)               
  call dern(dq2dn,q2,1,imaxpml,1,jmax)               

  call derm(dq1dm,q1,imax-D+1,imax,1,jmax)           
  call dern(dq1dn,q1,imax-D+1,imax,1,jmax)           
  call derm(dq2dm,q2,imax-D+1,imax,1,jmax)           
  call dern(dq2dn,q2,imax-D+1,imax,1,jmax)           

  call dern(dq1dn,q1,imaxpml+1,imax-D,1,jmaxpml)     
  call derm(dq1dm,q1,imaxpml+1,imax-D,1,jmaxpml)     
  call derm(dq2dm,q2,imaxpml+1,imax-D,1,jmaxpml)     
  call dern(dq2dn,q2,imaxpml+1,imax-D,1,jmaxpml)     

  call derm(dq1dm,q1,imaxpml+1,imax-D,jmax-D+1,jmax) 
  call dern(dq1dn,q1,imaxpml+1,imax-D,jmax-D+1,jmax) 
  call derm(dq2dm,q2,imaxpml+1,imax-D,jmax-D+1,jmax) 
  call dern(dq2dn,q2,imaxpml+1,imax-D,jmax-D+1,jmax) 

end subroutine

end module 
