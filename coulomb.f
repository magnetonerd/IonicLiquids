      implicit none
      
      double precision F(1000),r(1000),x(1000),y(1000),z(1000)
      double precision e0,q,pi,delta,r0,x0,y0,z0,kb,sqlaw(1000)
      integer N1,N2,i,j,k

      delta=1d-3
      pi=acos(-1d0)
      e0=pi*4d-7
      q=1.602d-19
      N1=1000
      N2=1000

      x0=10000d0
      y0=3d0
      z0=1d0

      do i=1,1000
         x(i)=x0+i*delta
         y(i)=y0+i*delta
         z(i)=z0+i*delta
      end do

      r0=0d0
      
      do i=1,1000
         r(i)=sqrt(x(i)**2+y(i)**2+z(i)**2)
      end do

      kb=1/(4*pi*e0)
      write(*,*)kb

      open(1,file='Coulomb.txt')
      do i=1,1000
         F(i)=kb*(q**2)*(N1*N2)/(r(i)**2)
         write(1,*)r(i),F(i)
      end do
      close(1)

      do i=2,999
         F(i)=log(F(i))
         F(i)=F(i)-log(kb*q**2*N1*N2)
         sqlaw(i)=F(i)/log(r(i))
         write(*,*)sqlaw(i)
      end do

      end
