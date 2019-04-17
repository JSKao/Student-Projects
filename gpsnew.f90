!Legendre(n) n=1,2,...,700
!0.1688032D-01 ! 1*10**13 W/cm**2 !            !CHECK Transformation matrix 
!0.2387238D-01 ! 2*10**13 W/cm**2 !
!0.2923758D-01 ! 3*10**13 W/cm**2 !
!0.3774555D-01 ! 5*10**13 W/cm**2 !
!0.05879129D0 = 775 nm omega
!0.056954D0 = 800 nm  omega
!0.0428226D0 = 1064 nm omega
program gps
implicit none
integer*8::i,j,k,l,t,N1,lmax1,T_STEP1,Ncyc1,W1,W2
!  Independent parameter
integer*8,parameter::N=300,lmax=20,Ncyc=40,NHM=100
real*8,parameter::Lmap=40.0d0,Rmax=150,dt=0.10d0,w0=0.05879129D0,INTENS=0.029237580D0
real*8,parameter::pi=3.141592653589790d0,small=10D-12
INTEGER*8::INFO=0

!  Dependent parameter
real*8,parameter::ALPHA=2.0d0*Lmap/Rmax,OPTCL=2*pi/w0,Ro=100
real*8,parameter::Tlength=Ncyc*OPTCL
real*8,parameter::dw=2*pi/Tlength
INTEGER*8,parameter::LWMAX=INT(34,kind=8)*N,LWORK=(4*N+2)*N,T_STEP=FLOOR(Tlength/dt)
INTEGER*8,parameter::WP=FLOOR(NHM*w0/dw)


REAL*8,ALLOCATABLE::W(:),WORK(:)
real*8::x1,x2,dx
real*8::summ,summ1,prob,time05
real*8,allocatable::xLGL(:),xLG(:),wLGL(:),wLG(:),r(:),r_p(:),D2(:,:),V(:,:),H(:,:)
REAL*8,ALLOCATABLE::E(:,:)

REAL*8,ALLOCATABLE::TIME(:),F(:),powr(:),powra(:)
COMPLEX*16::ZSUMM,ZSUMM1
COMPLEX*16,PARAMETER::IM=(0.0d0,1.0d0),one=(1.0d0,0.0d0),zero=0.0d0
COMPLEX*16,ALLOCATABLE::gl(:,:),gl2(:,:),grth(:,:),grth1(:,:),grth2(:,:),Cnl(:,:,:),Rnl(:,:,:)
COMPLEX*16,ALLOCATABLE::absorb(:),S(:,:,:),B(:,:,:),T_l_th(:,:),T_th_l(:,:),dipl(:),dipa(:)

REAL*8,ALLOCATABLE::powr_peak(:),powra_peak(:)
real*8::Tlengthr,dwr


allocate(xLGL(N),xLG(lmax),wLGL(N),wLG(lmax),r(N),r_p(N),D2(N,N),V(N,lmax),H(N,N))
allocate(W(N),WORK(LWMAX),E(N,lmax),Cnl(N,N,lmax),Rnl(N,N,lmax))
ALLOCATE(gl(N,lmax),gl2(N,lmax),S(N,N,lmax),B(N,lmax,T_STEP),T_l_th(lmax,lmax),T_th_l(lmax,lmax))
ALLOCATE(grth(N,lmax),grth1(N,lmax),grth2(N,lmax),TIME(T_step),F(T_STEP),absorb(N))
ALLOCATE(dipl(T_STEP),dipa(T_STEP),powr_peak(WP),powra_peak(WP),powr(WP),powra(WP))

N1=N
lmax1=lmax
T_STEP1=T_STEP
!open(1,file='r.dat')
!open(2,file='g20.dat')
open(3,file='laser.dat')
open(4,file='absorb.dat')
open(5,file='dip.dat')
open(6,file='sp.dat')
!open(7,file='probwave.dat')
!open(8,file='prob.dat')
!open(9,file='profile.dat')
!open(10,file='dip2.dat')
!open(11,file='T_test.dat')
open(12,file='peak.dat')
!open(13,file='energy.dat')
write(9,*)'N=',N,'Rmax=',Rmax,'lmax=',lmax,'cycles=',Ncyc,'dt=',dt,'omega=',w0,'intensity',INTENS,'Absorber=',Ro

call LGL_root
call LG_root


DO t=1,T_STEP
  TIME(t)=t*dt
  time05=dt*(t-0.5D0)
  if(time(t)<=10*OPTCL) then
  F(t)=INTENS*(Dsin(pi*time(t)/20.0d0/OPTCL))**2*DSIN(w0*time(t))
  else
  F(t)=INTENS*DSIN(w0*time(t))
  end if
  !F(t)=INTENS*(Dsin(pi*time05/Ncyc/OPTCL))**2*DSIN(w0*time05)
  !F(t)=0.0d0 
  write(3,*)time(t)/OPTCL,F(t)
END DO  

DO t=1,T_STEP
  DO k=1,lmax
    DO j=1,N
      B(j,k,t)=exp(IM*F(t)*dt*r(j)*xLG(k))
    END DO
  END DO
END DO      

do j=1,N  
  if (r(j)<ro) then    
    absorb(j)=1.0d0
  else if(r(j)>=ro) then
    absorb(j)=(cos(pi*(r(j)-ro)/2.0d0/(Rmax-ro)))**0.25
  end if
  !if (r(j)<ro) then    
  !  absorb(j)=1.0d0
  !else if(r(j)>=ro) then
  !  absorb(j)=1.0d0/(1.0d0+exp(1.250d0*r(j)-ro))
  !end if
  write(4,*)r(j),real(absorb(j),kind=8)
end do

!   Set Hamiltonian
DO l=1,lmax
 DO i=1,N
   V(i,l)=-1.0d0/r(i)+real((l-1)*l,kind=8)/2.0d0/r(i)**2
 END DO   

 DO j=1,N
  DO i=1,N
   if(i==j) then
     D2(i,j)=real((N1+1)*(N1+2),kind=8)/(r_p(i)**2*6.0d0*(1.0d0-xLGL(i)**2))
   else 
     D2(i,j)=1.0d0/(r_p(i)*r_p(j)*(xLGL(i)-xLGL(j))**2)
   end if
  end do
 END DO

 DO j=1,N
   DO i=1,N   
       if (i==j) then
       H(i,j)=D2(i,j)+V(i,l) 
       else
       H(i,j)=D2(i,j)
       end if 
   END DO
 END DO
    


!     Solve eigenproblem.

  CALL DSYEV( 'Vectors', 'Upper', N, H ,N, W, WORK, LWORK, INFO )

!     Check for convergence.

  IF( INFO.GT.0 ) THEN
     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
     STOP  
  END IF

 
  DO i=1,n-l+1
    E(i+l-1,l)=W(i)
    
  END DO
  
  
  DO i=1,n-l+1  
    DO j=1,N
      Cnl(j,i+l-1,l)=H(j,i)
    end do
  end do
  
  !do i=1,N
  !  write(13,*)i,l,E(i,l)
  !end do
  
END DO

!  Get real eigenfunction

!call realwave(2,5)

!   Ground state
DO l=1,lmax
  DO j=1,N
    IF(l==1) THEN
      gl(j,l)=Cnl(j,1,1)
    ELSE
      gl(j,l)=0
    END IF
  END DO
END DO
!goto 1
!   Sij
DO l=1,lmax
  DO j=1,N
    DO i=1,j
      ZSUMM=0
      DO k=1,N
         IF(E(k,l).lt.30) &
         ZSUMM=ZSUMM+DCONJG(Cnl(i,k,l))*Cnl(j,k,l)*exp(-IM*E(k,l)*dt/2.0d0)
      END DO
      S(i,j,l) = ZSUMM
      S(j,i,l) = ZSUMM    !Symmetry of Sij!
    END DO
  END DO
END DO


!    Coordinate Transform matrix
DO l=1,lmax
  DO k=1,lmax
    T_l_th(k,l)=Legendre(l-1,xLG(k))*(-1.0d0)**(l-1)*DSQRT((2.0d0*l-1.0d0)/2.0d0)
    T_th_l(k,l)=Legendre(l-1,xLG(k))*(-1.0d0)**(l-1)*DSQRT((2.0d0*l-1.0d0)/2.0d0)*wLG(k)
  END DO
END DO



!   T-LOOP
DO t=1,T_STEP

  DO l=1,lmax
  CALL ZGEMV('N',N,N,one,S(1,1,l),N,gl(1,l),1,zero,gl2(1,l),1)
  END DO  

  CALL ZGEMM('N','T',N,lmax,lmax,one,gl2,N,T_l_th,lmax,zero,grth,N)

  do k=1,lmax
  do j=1,N 
  grth1(j,k)=grth(j,k)*B(j,k,t)
  end do
  end do

  CALL ZGEMM('N','N',N,lmax,lmax,one,grth1,N,T_th_l,lmax,zero,gl,N)

  DO l=1,lmax
  CALL ZGEMV('T',N,N,one,S(1,1,l),N,gl(1,l),1,zero,gl2(1,l),1)
  END DO
  
!  Absorb
  DO l=1,lmax
     DO j=1,n
        gl(j,l)=gl2(j,l)*absorb(j)
     END DO
  END DO

  call prob1
  write(8,*)time(t)/OPTCL,prob
!
  CALL ZGEMM('N','T',N,lmax,lmax,one,gl,N,T_l_th,lmax,zero,grth2,N)  


!   DIPOLE Length
  zsumm=0
  zsumm1=0
  do k=1,lmax    
  do j=1,n
  zsumm=zsumm+wLG(k)*Dconjg(grth2(j,k))*grth2(j,k)*r(j)*xLG(k)
  zsumm1=zsumm1+wLG(k)*Dconjg(grth2(j,k))*grth2(j,k)*(F(t)-xLG(k)/r(j)**2)
  end do
  end do
  dipl(t)=zsumm
  dipa(t)=zsumm1
  write(5,*)time(t)/optcl,real(dipl(t),kind=8),real(dipa(t),kind=8)
!   END OF TIME_LOOP
END DO

!   Comp. Spectrum
!call zfft1d(dipl,32768,0,wsave)
!call zfft1d(dipl,32768,-1,wsave)
!call ramp(25)
DO i=1,WP
  ZSUMM=0
  ZSUMM1=0
  DO t=1,T_STEP
    if(dipl(t)/=0.and.dipa(t)/=0) then
    !ZSUMM=ZSUMM+dipl(t)*cdexp(-IM*i*dwr*time(t))
    !ZSUMM1=ZSUMM1+dipa(t)*cdexp(-IM*i*dwr*time(t))/(i*dwr)**2
    ZSUMM=ZSUMM+dipl(t)*cdexp(-IM*i*dw*time(t))
    ZSUMM1=ZSUMM1+dipa(t)*cdexp(-IM*i*dw*time(t))/(i*dw)**2
    end if
  END DO
  !powr(i)=real(DCONJG(ZSUMM)*ZSUMM*dt*dt/tlengthr/tlengthr,kind=8)
  !powra(i)=real(DCONJG(ZSUMM1)*ZSUMM1*dt*dt/tlengthr/tlengthr,kind=8)
  powr(i)=real(DCONJG(ZSUMM)*ZSUMM*dt*dt/tlength/tlength,kind=8)
  powra(i)=real(DCONJG(ZSUMM1)*ZSUMM1*dt*dt/tlength/tlength,kind=8)
  !write(6,*)i*dwr/w0,Dlog10(powr(i)),Dlog10(powra(i))
  write(6,*)i*dw/w0,Dlog10(powr(i)),Dlog10(powra(i)) 
END DO  


call peak(3)
call peak(5)
call peak(7)
call peak(9)
call peak(11)
call peak(13)

!   Just to plot eigenfunction
  !do j=1,N
    !Rnl(j,i,l)=Cnl(j,i,l)*DSQRT((N1+1)*(N1+2)/2.0d0/r_p(j))*Legendre(N1+1,xLGL(j))  
  !  write(2,*)r(j),REAL(DCONJG(gl(j,1))*gl(j,1),kind=8)
  !end do
  !do l=1,lmax 
  ! call probl(l)
  ! write(7,*)l,prob
  !end do 
  


contains
subroutine prob2()
prob=0
do k=1,lmax
  do j=1,n
     prob=prob+Dconjg(grth2(j,k))*grth2(j,k)*wLG(k)
  end do
end do
write(*,*)'total probability of grtheta in the',T,'th step=',prob
end subroutine
subroutine prob1()
prob=0
do l=1,lmax
  do j=1,n
    prob=prob+Dconjg(gl(j,l))*gl(j,l)
  end do
end do
write(*,*)'total probability of gl in the',T,'th step=',prob
end subroutine
subroutine prob0()
do l=1,lmax
  do i=1,n
    prob=0
    do j=1,n
      prob=prob+DCONJG(cnl(j,i,l))*cnl(j,i,l)
    end do
    write(*,*) i,l,prob
  end do
end do
end subroutine
subroutine probl(m)
integer*8::m
prob=0
  do j=1,n
    prob=prob+Dconjg(gl(j,m))*gl(j,m)
  end do
write(*,*)'total probability of gl in l=',m,'is',prob
end subroutine


subroutine LGL_root()           
dx=(1.0d0/(N1+150))**2
!write(*,*)'dx=',dx
x1=-1
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
  r(i)=Lmap*(1.0d0+xLGL(i))/(1.0d0-xLGL(i)+alpha)
  r_p(i)=Lmap*(2.0d0+ALPHA)/((1.0d0-xLGL(i)+ALPHA)**2)
end do
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

subroutine realwave(i,l)
integer::i,l
do j=1,n
  Rnl(j,i,l)=Cnl(j,i,l)*DSQRT((N1+1)*(N1+2)/2.0d0/r_p(j))*Legendre(N1+1,xLGL(j)) 
  write(1,*)r(j),REAL(Cnl(j,i,l),KIND=8),REAL(Rnl(j,i,l),KIND=8)
end do 
end subroutine 

subroutine ramp(Ncyc1)
integer::Ncyc1
do t=1,T_STEP
  if(time(t)<Ncyc1*OPTCL) then
    dipl(t)=0
    dipa(t)=0
  end if 
  write(10,*)time(t)/OPTCL,real(dipl(t),kind=8),real(dipa(t),kind=8)
end do
i=0
do t=1,T_STEP
  if(time(t)>=Ncyc1*OPTCL) then
    i=i+1 
    time(t)=time(i)
  end if
  !write(11,*)t,time(t)
end do    
Tlengthr=(Ncyc-Ncyc1)*OPTCL
dwr=2*pi/Tlengthr
end subroutine   


subroutine peak(i)
integer::i,j,W1,W2
!W1=floor((i-1)*w0/dw)
!W2=floor((i+1)*w0/dw)
W1=floor((i-1)*w0/dwr)
W2=floor((i+1)*w0/dwr)
powr_peak=0
powra_peak=0
do j=W1,W2
powr_peak(j)=powr(j)
powra_peak(j)=powra(j)
!write(12,*)j,powr_peak(j),powra_peak(j)
end do
write(12,*)i,maxval(powr_peak),maxval(powra_peak)
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


END PROGRAM