program DFT
implicit none
integer*8::i,j,k,N1,lmax1,step
integer*8,parameter::N=300,lmax=32
real*8::x1,x2,dx
real*8::summ,summ1,summ2,Egs
complex*16::zsumm
real*8,parameter::Lmap=20.0d0,Rmax=60
real*8,parameter::pi=3.141592653589790d0,small=10D-12
complex*16,parameter::IM=(0.D0,1.D0)

real*8,allocatable::r(:),r_p(:),xLGL(:),wLGL(:),xLG(:),wLG(:),d_r(:)
real*8,allocatable::D2(:,:),H(:,:),H1(:,:),H2(:,:),WORK(:),W(:)
real*8,allocatable::Vext(:),Vh(:),Vx1(:),Vx2(:),den1(:),den2(:),den(:)
real*8,allocatable::Vsic1(:),Vsic2(:),Vs1(:),Vs2(:)
complex*16,allocatable::Yml(:,:,:,:)

integer*8,parameter::LWORK=(4*N+2)*N,LWMAX=INT(34,kind=8)*N
real*8,parameter::ALPHA=2.0d0*Lmap/Rmax


INTEGER*8::INFO=0

allocate(r(N),r_p(N),xLGL(N),wLGL(N),xLG(lmax),wLG(lmax),D2(N,N),H1(N,N),H2(N,N),W(N),WORK(LWMAX))
allocate(Vext(N),den1(N),den2(N),den(N),d_r(N),Vh(N),Vx1(N),Vx2(N))
allocate(Vs1(N),Vs2(N),Yml(lmax,lmax,2,3))
N1=N
lmax1=lmax
open(2,file='vh.dat')
call LGL_root
call D_square
call Sph_Harmonics

!   External potential Vext
DO j=1,N
   Vext(j)=-2.D0/r(j)
END DO   


!      KS LOOP
DO step=1,60
  
  DO j=1,N
  !summ1=0
  !summ2=0
  !DO i=1,Z
  !  summ1=summ1+(d_r(j)*H1(j,i))**2
  !  summ2=summ2+(d_r(j)*H2(j,i))**2
  !  den1(j)=summ1
  !  den2(j)=summ2
    den1(j)=(d_r(j)*H1(j,1))**2*DCONJG(Ylm
    den2(j)=(d_r(j)*H2(j,1))**2
  !END DO  
  den(j)=den1(j)+den2(j)
  END DO
    
  !    Hartree potential , exchange potential
  DO j=1,N
    summ1=0
    summ2=0
    DO i=1,N
      if (i<=j) then
      summ1=summ1+4.D0*pi*1.0d0/r(j)*den(i)*r(i)**2*r_p(i)*wLGL(i)
      else
      summ2=summ2+4.D0*pi*den(i)*r(i)*r_p(i)*wLGL(i)
      end if
    END DO
    Vh(j)=0.7D0*Vh(j)+0.3D0*(summ1+summ2)
    Vx1(j)=0.7D0*Vx1(j)+0.3D0*(6.D0/pi*den1(j))**(1.D0/3.D0)
    Vx2(j)=0.7D0*Vx2(j)+0.3D0*(6.D0/pi*den2(j))**(1.D0/3.D0)
    Vs1(j)=Vext(j)+Vh(j)-Vx1(j)
    Vs2(j)=Vext(j)+Vh(j)-Vx2(j)
  END DO  
    
    
  !    Set up Hamiltonian
  DO j=1,N
    DO i=1,N   
       if (i==j) then
         H1(i,j)=D2(i,j)+Vs1(i)
         H2(i,j)=D2(i,j)+Vs2(i) 
       else
         H1(i,j)=D2(i,j)
         H2(i,j)=D2(i,j)
       end if 
    END DO
  END DO

  !     Solve TISE for spin-up electron
  CALL DSYEV( 'Vectors', 'Upper', N, H1 ,N, W, WORK, LWORK, INFO )
  !     Check for convergence.
  IF( INFO.GT.0 ) THEN
     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
     STOP  
  END IF
  

  write(*,*)'up',W(1)
  
  !     Solve TISE for spin-down electron
  CALL DSYEV( 'Vectors', 'Upper', N, H2 ,N, W, WORK, LWORK, INFO )
  !     Check for convergence.
  IF( INFO.GT.0 ) THEN
     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
     STOP  
  END IF


  write(*,*)'down',W(1)
  
END DO 




!      Testing real wavefunctions
!summ=0
!DO j=1,N
!  summ=summ+4*pi*(d_r(j)*H1(j,1))**2*r(j)**2*r_p(j)*wLGL(j)
!  write(1,*)r(j),d_r(j)*H1(j,1)
!ENd DO  
!write(*,*)'prob=',summ


contains
subroutine LGL_root()           
dx=(1.0d0/(N1+150))**2
!write(*,*)'dx=',dx
x1=-1.D0
DO i=1,N  !ith root finding
   x2=x1+dx
  DO while(Legendre_p(n1+1,x1)*Legendre_p(n1+1,x2)>0)
    x1=x2
    x2=x1+dx
    !write(*,*)'a=',a,'b=',b
  end DO
  x1=(x1+x2)/2.0d0
  x2=x2
  DO while(abs(x1-x2)>=small)
    x2=x1
    x1=x1-Legendre_p(n1+1,x1)/Legendre_pp(n1+1,x1)
  end DO
  xLGL(i)=x1
  x1=x1+dx
end DO
do i=1,n
  wLGL(i)=2.0d0/((n1+2.0d0)*(n1+1.0d0)*(legendre(n1+1,xLGL(i)))**2)
  r(i)=Lmap*(1.0d0+xLGL(i))/(1.0d0-xLGL(i)+ALPHA)
  r_p(i)=Lmap*(2.0d0+ALPHA)/((1.0d0-xLGL(i)+ALPHA)**2)
end do

!   factor brings eigenvector to real wavefunction ( Cn --->  Rn(r)*Y00(theta,phi) )
DO j=1,N
  d_r(j)=DSQRT((N1+1)*(N1+2)/2.0d0/r_p(j))*Legendre(N1+1,xLGL(j))/r(j)/DSQRT(4.D0*pi) 
END DO  
end subroutine

subroutine LG_root()
dx=(1.0d0/(N1+50.0d0))**2
!write(*,*)'h=',h,'order=',np+1
x1=-1.0d0
do i=1,lmax
   x2=x1+dx
  do while(Legendre(lmax1,x1)*Legendre(lmax1,x2)>0)
    x1=x2
    x2=x1+dx
  end do
  x1=(x1+x2)/2.0d0
  x2=x2
  do while(abs(x1-x2)>=small)
    x2=x1
    x1=x1-Legendre(lmax1,x1)/Legendre_p(lmax1,x1)
  end do
  xLG(i)=x1
  x1=x1+dx
end do

do i=1,lmax
  wLG(i)=2.0d0/((1.0d0-(xLG(i))**2)*(Legendre_p(lmax1,xLG(i)))**2)
end do
end subroutine

subroutine D_square()
!   Prepare D2 and external potential outside the main LOOP
DO j=1,N
   DO i=1,N
    if(i==j) then
      D2(i,j)=real((N1+1)*(N1+2),kind=8)/(r_p(i)**2*6.0d0*(1.0d0-xLGL(i)**2))
    else 
      D2(i,j)=1.0d0/(r_p(i)*r_p(j)*(xLGL(i)-xLGL(j))**2)
    end if
   end do
END DO
end subroutine

subroutine Sph_harmonics()
DO j=1,lmax
  DO k=1,lmax
    Ylm(j,k,1,1)=1.D0/DSQRT(4.D0*pi)
    Ylm(j,k,2,1)=DSQRT(3.D0/8.D0/pi)*DSIN(xLG(j))*DEXP(-IM*xLG(k))
    Ylm(j,k,2,2)=DSQRT(3.D0/4.D0/pi)*DCOS(xLG(j))
    Ylm(j,k,2,3)=-DSQRT(3.D0/8.D0/pi)*DSIN(xLG(j))*DEXP(IM*xLG(k))
  END DO
END DO
end subroutine

function Legendre(l,x)
implicit none
integer*8::k,l
real*8::x,Legendre
real*8,allocatable::s(:)
allocate(s(801))
s(1)=1.0d0
s(2)=x
do k=1,799
   s(k+2)=(real(2*k+1,kind=8)*x*s(k+1)-real(k,kind=8)*s(k))/real(k+1,kind=8)
end do
Legendre=s(l+1)
end function


function Legendre_p(l,x)
implicit none
integer*8::l,j,k
real*8::x,Legendre_p
real*8,allocatable::s(:,:)
allocate(s(801,2))
s(1,1)=1.0d0

do k=1,799
  do j=1,2
    if (j==1) then           !!zero derivative
      s(2,j)=x
      s(k+2,j)=(real(2*k+1,kind=8)*x*s(k+1,j)-real(k,kind=8)*s(k,j))/real(k+1,kind=8)
    else
      s(1,j)=0.0d0
      if (j==2) then          !!1st order derivative
        s(2,j)=1.0d0
        s(k+2,j)=real(2*k+1,kind=8)*s(k+1,j-1)+s(k,j)
      else
        s(2,j)=0.0d0
      end if
    end if
  end do
end do
Legendre_p=s(l+1,2)
end function

function Legendre_pp(l,x)
implicit none
integer*8::j,k,l
real*8::x,Legendre_pp
real*8::s(801,3)

s(1,1)=1.0d0

do k=1,799
  do j=1,3
    if (j==1) then           !!zero derivative
      s(2,j)=x
      s(k+2,j)=(real(2*k+1,kind=8)*x*s(k+1,j)-real(k,kind=8)*s(k,j))/real(k+1,kind=8)
    else
      s(1,j)=0
      if (j==2) then          !!1st order derivative
        s(2,j)=1.0d0
        s(k+2,j)=real(2*k+1,kind=8)*s(k+1,j-1)+s(k,j)
      else
        s(2,j)=0.0d0
        if(j==3) then         !!2nd order derivative
          s(k+2,j)=real(2*k+1,kind=8)*s(k+1,j-1)+s(k,j)
        else
          s(3,j)=0.0d0
        end if
      end if
    end if
  end do
end do
Legendre_pp=s(l+1,3)
end function

end program
