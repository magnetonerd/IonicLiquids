      implicit none

      double precision E(1000,1000),Grad(1000),a(1000),v(1000)
      double precision q,u0,e0,pi,w,k,Efield,x(1000),t(1000)
      double precision delta,x0,t0,Lap(1000),Tder2(1000)
      double precision left(1000),right(1000),lamda
      integer n(1000),nx,ny,nz,total
      !parameter(nx=1000d0,ny=1000d0,nz=1000d0)
      integer i,j,l

      pi=acos(-1d0)
      q=1.602d-19 !Elementary charge unit (C)
      e0=8.854d-12 !Electric permitivity (F/m)
      u0=pi*4d-7 !Permiability (H/m)
      Efield=1d4 !Initial Electric Field (V/m)
      lamda=1d3 !Wavelength (m)
      w=1d5 !Angular Freq (1/s)
      k=sqrt((w**2/(u0*e0))-(q*1000/Efield)) !Wave number (1/m)
      delta=1d-3 !Step size

      write(*,*)'Charge',q
      write(*,*)'Permitivity',e0
      write(*,*)'Permiability',u0
      write(*,*)'Speed of Light',1/sqrt(u0*e0)
      write(*,*)'Initial Efield Strength',Efield
      write(*,*)'Wave Number',k
      write(*,*)'Angular Frequency',w
      write(*,*)'Step Size',delta
      
      !Initial conditions of position and time
      x0=0d0
      t0=0d0

      !Values position runs over
      do i=1,1000
         x(i)=x0+i*delta
      end do

      !values time runs over
      do j=1,1000
         t(j)=t0+j*delta
      end do

      !The Electric Field
      do i=1,1000
         do j=1,1000
            E(j,i)=Efield*cos(k*x(i)-w*t(j))
         end do
      end do

      !Laplacean of The Electric Field
      do i=2,999
         do j=1,1000
            Lap(i)=(E(j,i+1)+E(j,i-1)-(2*E(j,i)))/delta**2
         end do
      end do

      !number of charge carriers at site i
     
      total=0
      
      do i=1,1000
         n(i)=x(i)*1000
         total=total+n(i)
      end do

      write(*,*)'Number of Charges',total

      !Gradiant of Charge Carriers
      do i=2,999
         Grad(i)=(q/e0)*(n(i+1)-n(i-1))/(2*delta)
      end do
      
      open(1,file='Leftside.txt')
      do i=2,999
         left(i)=Lap(i)-Grad(i)
         write(1,*)x(i),left(i)
      end do
      close(1)

      do i=2,999
         v(i)=(x(i+1)-x(i-1))/(2*delta)
      end do
      
      do i=2,999
         a(i)=(v(i+1)-v(i-1))/(2*delta)
      end do

      do i=1,1000
         do j=2,999
            Tder2(j)=e0*u0*(E(j+1,i)+E(j-1,i)-(2*E(j,i)))/delta**2
         end do 
      end do

      open(2,file='Rightside.txt')
      do i=2,999
         right(i)=Tder2(i)+(u0*q*n(i)*a(i))
         write(2,*)x(i),right(i)
      end do
      close(2)

      open(3,file='Master.txt')
      do i=2,999
         write(3,*)x(i),left(i),right(i)
      end do
      close(3)

      end
