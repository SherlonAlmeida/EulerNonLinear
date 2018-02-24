module diff
use global
contains
!.................................................................
subroutine dernb(dudn,u)

!.......................................................................
!     subrotina para calculo da derivada primeira de uma funcao u
!     utilizando diferencas centradas de quarta ordem nos pontos
!     internos do dominio, diferencas centradas de segunda ordem
!     nos pontos vizinhos 'a fronteira e diferencas unilateral de
!     segunda ordem nos pontos de fronteira.
!.......................................................................
      implicit none
      
      integer:: i, j
      
      real(kind=ip):: dv4th, dv2nd, onesp, onesm
      real(kind=ip):: ap2, ap1, am1, am2, del, a 
      real(kind=ip),dimension(jm):: u, dudn
     
    
      
!......................................................function statement
      dv4th(ap2,ap1,am1,am2,del) =                                      &
     &          (-ap2 + 8.d0*ap1 - 8.d0*am1 + am2)     / 12.d0 / del
      dv2nd(ap1,am1,del)   = ( ap1 - am1 )             * 0.5d0 / del
      onesp(a,ap1,ap2,del) = (-3.d0*a + 4.d0*ap1 - ap2) * 0.5d0 / del
      onesm(a,am1,am2,del) = ( 3.d0*a - 4.d0*am1 + am2) * 0.5d0 / del


!.................. Dentro do Dominio ..................................
      do j=3,jmax-2
         dudn(j) = dv4th(u(j+2),u(j+1),u(j-1),u(j-2),dn)
      enddo

!..........................Linha Especiais .............................
      j=2
         dudn(j) = dv2nd(u(j+1),u(j-1),dn)

!.......................................................................
      j=jmax-1
         dudn(j) = dv2nd(u(j+1),u(j-1),dn)

!.......................................................................
      j=jmax
         dudn(j) = onesm(u(j),u(j-1),u(j-2),dn)

!.......................................................................
      j=1
         dudn(j) = onesp(u(j),u(j+1),u(j+2),dn)
      
end subroutine


subroutine derm(dU,U,ii,fi,ij,fj)

!.......................................................................
!     subrotina para calculo da derivada primeira de uma funcao u
!     utilizando diferencas centradas de quarta ordem nos pontos
!     internos do dominio, diferencas centradas de segunda ordem
!     nos pontos vizinhos 'a fronteira e diferencas unilateral de
!     segunda ordem nos pontos de fronteira.
!.......................................................................
      implicit none
      
      integer:: k, i, j, ii,fi,ij,fj
      
      real(kind=ip):: dv4th, dv2nd, onesp, onesm
      real(kind=ip):: ap2, ap1, am1, am2, del, a 
      real(kind=ip),dimension(4,im,jm):: U,dU
    
     
      
!......................................................function statement
      dv4th(ap2,ap1,am1,am2,del) =                                      &
     &          (-ap2 + 8.d0*ap1 - 8.d0*am1 + am2)      / 12.d0 / del
      dv2nd(ap1,am1,del)   = ( ap1 - am1 )              * 0.5d0 / del
      onesp(a,ap1,ap2,del) = (-3.d0*a + 4.d0*ap1 - ap2) * 0.5d0 / del
      onesm(a,am1,am2,del) = ( 3.d0*a - 4.d0*am1 + am2) * 0.5d0 / del

do k=1,4
!.................. Dentro do Dominio ..................................
      do j=ij,fj
      do i=ii+2,fi-2      
         du(k,i,j)  = dv4th(u(k,i+2,j),u(k,i+1,j),u(k,i-1,j),u(k,i-2,j),dm)
      enddo
      enddo
!.......................................................................
      i=ii
      do j=ij,fj
         du(k,i,j) = onesp(u(k,i,j),u(k,i+1,j),u(k,i+2,j),dm)
      enddo
!......................................................................      
      
      i=ii+1
      do j=ij,fj
         du(k,i,j) = dv2nd(u(k,i+1,j),u(k,i-1,j),dm)
      enddo
!.......................................................................     
      
      i=fi
      do j=ij,fj
         du(k,i,j) = onesm(u(k,i,j),u(k,i-1,j),u(k,i-2,j),dm)
      enddo

!.......................................................................
      i=fi-1
      do j=ij,fj
         du(k,i,j) = dv2nd(u(k,i+1,j),u(k,i-1,j),dm)
      enddo

!.......................................................................
end do      
end subroutine

subroutine dern(dU,U,ii,fi,ij,fj)

!.......................................................................
!     subrotina para calculo da derivada primeira de uma funcao u
!     utilizando diferencas centradas de quarta ordem nos pontos
!     internos do dominio, diferencas centradas de segunda ordem
!     nos pontos vizinhos 'a fronteira e diferencas unilateral de
!     segunda ordem nos pontos de fronteira.
!.......................................................................
      implicit none
      
      integer:: ii,fi,ij,fj,k, i,j
      
      real(kind=ip):: dv4th, dv2nd, onesp, onesm
      real(kind=ip):: ap2, ap1, am1, am2, del, a 
      real(kind=ip),dimension(4,im,jm):: U,dU
     
     
      
!......................................................function statement
      dv4th(ap2,ap1,am1,am2,del) =                                      &
     &          (-ap2 + 8.d0*ap1 - 8.d0*am1 + am2)     / 12.d0 / del
      dv2nd(ap1,am1,del)   = ( ap1 - am1 )             * 0.5d0 / del
      onesp(a,ap1,ap2,del) = (-3.d0*a + 4.d0*ap1 - ap2) * 0.5d0 / del
      onesm(a,am1,am2,del) = ( 3.d0*a - 4.d0*am1 + am2) * 0.5d0 / del

do  k=1,4
!.................. Dentro do Dominio ..................................
      do j=ij+2,fj-2
      do i=ii,fi
         du(k,i,j) = dv4th(u(k,i,j+2),u(k,i,j+1),u(k,i,j-1),u(k,i,j-2),dn)
      enddo
      enddo

!..........................Linha Especiais .............................
      j=ij+1
      do i=ii,fi
         du(k,i,j) = dv2nd(u(k,i,j+1),u(k,i,j-1),dn)
      enddo

!.......................................................................
      j=fj-1
      do i=ii,fi
         du(k,i,j) = dv2nd(u(k,i,j+1),u(k,i,j-1),dn)
      enddo

!.......................................................................
      j=fj
      do i=ii,fi
         du(k,i,j) = onesm(u(k,i,j),u(k,i,j-1),u(k,i,j-2),dn)
      enddo

!.......................................................................
      j=ij
      do i=ii,fi
         du(k,i,j) = onesp(u(k,i,j),u(k,i,j+1),u(k,i,j+2),dn)
      enddo
!.......................................................................
end do      
end subroutine

end module
