implicit none
include 'mpif.h'
double precision cpu1,cpu2
!
!~~~ fundamental constants [all numbers are in SI units]~~~!
!
double precision, parameter :: pi=3.1415926535897932384626433832795D0,c=299792458.0D0
double precision, parameter :: mu0=4.0D-7*pi,eps0=1.0D0/(c*c*mu0)
double precision, parameter :: h=1.054571628D-34
double complex, parameter :: Im=(0.0,1.0)
double precision, parameter :: ev_to_radsec=2.0*pi*2.4180e14
!
!~~~ number of grid points & time steps ~~~!
!
integer, parameter :: Nt=200000,N_w=400
double precision, parameter :: omega_min=ev_to_radsec*1.5,omega_max=ev_to_radsec*4.0

integer, parameter :: Ny=1281,N_loc=40
double precision, parameter :: y0=-640.0D-9,yM=640.0D-9

integer, parameter :: Nx=3001
double precision, parameter :: x0=-1500.0e-9,xM=1500.0e-9

!
!~~~ CPML ~~~!
!
integer, parameter :: npml=19,m=3,ma=1 
double precision sigmaCPML,alphaCPML,kappaCPML
double precision psi_Hzy_1(Nx-1,npml-1),psi_Exy_1(Nx-1,npml)                              
double precision psi_Hzy_2(Nx-1,npml-1),psi_Exy_2(Nx-1,npml)
double precision be_y(npml),ce_y(npml),alphae_y(npml),sige_y(npml),kappae_y(npml)
double precision bh_y(npml-1),ch_y(npml-1),alphah_y(npml-1),sigh_y(npml-1),kappah_y(npml-1)
double precision den_ex(Nx),den_hx(Nx),den_ey(N_loc),den_hy(N_loc)

double precision psi_Hzx_1(npml-1,N_loc),psi_Eyx_1(npml,N_loc)
double precision psi_Hzx_2(npml-1,N_loc),psi_Eyx_2(npml,N_loc)
double precision be_x(npml),ce_x(npml),alphae_x(npml),sige_x(npml),kappae_x(npml)
double precision bh_x(npml-1),ch_x(npml-1),alphah_x(npml-1),sigh_x(npml-1),kappah_x(npml-1)

double precision psi_Hzy_1_inc(npml-1),psi_Exy_1_inc(npml)                              
double precision psi_Hzy_2_inc(npml-1),psi_Exy_2_inc(npml)

!
!~~~ scattered field zone ~~~!
!
integer, parameter :: i0 = 1 + (npml), i1 = Nx - (npml)
integer, parameter :: mj0=1,j0=11  !<--- - 590nm
integer, parameter :: mj1=28,j1=21 !<--- + 500nm

integer, parameter :: ms=29,js=21

double precision, parameter :: tau=0.36d-15,E0=1.0,omega=ev_to_radsec*3.0
double precision aBH(4)
double precision pulse(Nt)
double precision tmp1,tmp2,omega_P(N_w),SN(N_w,2)

!
!~~~ physical grid ~~~!
!
double precision, parameter :: dx=(xM-x0)/(Nx-1),dy=(yM-y0)/(Ny-1)
double precision x(Nx),xM2(Nx-1),y(N_loc),yM2(N_loc)
!
!~~~ dielectric constant of the host media (=1 for vacuum) ~~~!
!
double precision, parameter :: eps_delectric=1.0
!
!~~~ time step ~~~!
!
double precision, parameter :: dt=dy/(2.0*c)

!
!~~~ Drude model for Ag ~~~!
!
double precision, parameter :: eps_r=8.926,omegaD=ev_to_radsec*11.585,GammaD=ev_to_radsec*0.203
double precision, parameter :: A1=(2.0-GammaD*dt)/(2.0+GammaD*dt),A2=eps0*omegaD*omegaD*dt/(2.0+GammaD*dt)
double precision, parameter :: Ca=(eps_r*eps0/dt-0.5*A2)/(eps_r*eps0/dt+0.5*A2)
double precision, parameter :: Cb=1.0/(eps_r*eps0/dt+0.5*A2)
double precision, parameter :: Cc=0.5*(A1+1.0)/(eps_r*eps0/dt+0.5*A2)

double precision tmpE

double precision PDx(Nx-1,N_loc),PDy(Nx,N_loc)

logical FBx(Nx-1,N_loc),FBy(Nx,N_loc)

!
!~~~ Geometry ~~~!
!

double precision, parameter :: d_half =  100.0D-9
double precision, parameter :: Ag1 = -199.757575D-9, Ag2 = 200.252525D-9
double precision, parameter :: groove_depth = 100.0D-9, slit_depth = 400.0D-9
double precision, parameter :: groove_width = 100.0D-9, slit_width = 100.0D-9

double precision :: groove_d1, groove_d2, groove_w1, groove_w2
double precision :: slit_d1, slit_d2, slit_w1, slit_w2

!
!~~~ detection ~~~!
!
double precision tmp

integer, parameter :: iw1 = (Nx+1)/2+(d_half-100.0D-9)/dx
integer, parameter :: iw2 = iw1 + 200.0D-9/dx
double precision, parameter :: ywT = Ag2 - slit_depth - 400.0D-9 ! slit_d2 - 400nm
integer :: mwT, jwT

double precision Ex_temp(iw1:iw2,N_w,2),Hz_temp(iw1:iw2,N_w,2)
double precision Ex_temp_inc(iw1:iw2,N_w,2),Hz_temp_inc(iw1:iw2,N_w,2)
double precision Ex_w,Hz_w,av_x1,av_x2,av_y
double complex sum_int,cTMP1,cTMP2
double precision P_sct(N_w),P_inc(N_w)

!
!~~~ Grid Return ~~~!
!

integer :: Drude_Grid(Nx-1,Ny), FB_Grid(Nx-1,N_loc), a, mReturn

!
!~~~ EM field components ~~~!
!
double precision Ex(Nx-1,N_loc),Ey(Nx,N_loc),Hz(Nx-1,N_loc)
double precision Ex_inc(N_loc),Hz_inc(N_loc)
!
!~~~ other parameters and variables ~~~!
!
integer i,ii,j,jj,n,nn,k
double precision t
double precision, parameter :: dt_eps0=dt/eps0,dt_mu0=dt/mu0
!
!~~~ MPI part ~~~!
!
integer ierr,nprocs,myrank,j_glob,mfdtd,n1
integer :: istatus(MPI_STATUS_SIZE)
integer itag,ireq,itag1,itag2,itag3,itag4,itag5,itag6
double precision Ex_get(Nx-1),Ex_send(Nx-1)
double precision Hz_get(Nx-1),Hz_send(Nx-1)
double precision Ex_send_inc,Ex_get_inc
double precision Hz_send_inc,Hz_get_inc

!-------------------------!
 call MPI_INIT(ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

do nn=1,N_w
 omega_P(nn)=omega_min+(omega_max-omega_min)*(nn-1)/(N_w-1)
enddo

aBH(1)=0.353222222
aBH(2)=-0.488
aBH(3)=0.145
aBH(4)=-0.010222222

pulse=0.0
SN=0.0
do n=1,Nt
 t=dt*dble(n)
 if(t<=tau)then
   pulse(n)=E0*cos(omega*t)*( &
                 aBH(1)+ &
				 aBH(2)*cos(2.0*pi*t/tau)+ &
				 aBH(3)*cos(2.0*pi*2.0*t/tau)+ &
				 aBH(4)*cos(2.0*pi*3.0*t/tau))
  else
   pulse(n)=0.0
 endif

 do nn=1,N_w
  tmp1=sin(omega_P(nn)*t)
  tmp2=cos(omega_P(nn)*t)
  SN(nn,1)=SN(nn,1)+pulse(n)*tmp1
  SN(nn,2)=SN(nn,2)+pulse(n)*tmp2
 enddo
enddo

do nn=1,N_w
 tmp1=sqrt(SN(nn,1)**2+SN(nn,2)**2)
 SN(nn,2)=-atan2(SN(nn,1),SN(nn,2))
 SN(nn,1)=tmp1
enddo

!~~~ grid ~~~!
do i=1,Nx
 x(i)=x0+dx*(i-1)
enddo
do i=1,(Nx-1)
 xM2(i)=x0+dx*(i-1)+dx/2.0
enddo
do j=1,N_loc
 j_glob=myrank*N_loc+j
 y(j)=y0+dy*(j_glob-1)
enddo
do j=1,N_loc
 j_glob=myrank*N_loc+j
 yM2(j)=y0+dy*(j_glob-1)+dy/2.0
enddo

FBx=.false.
FBy=.false.

!~~~ structure ~~~!

groove_d1 = Ag2
groove_d2 = groove_d1 - groove_depth
groove_w1 = -1.0*d_half - groove_width/2.0
groove_w2 = -1.0*d_half + groove_width/2.0
slit_d1 = Ag2
slit_d2 = slit_d1 - slit_depth
slit_w1 = 1.0*d_half - slit_width/2.0
slit_w2 = 1.0*d_half + slit_width/2.0

do i=1,Nx-1
 do j=1,N_loc
 
  if( y(j) >= Ag1 .and. y(j) <= Ag2 )then ! True in the metal zone
   FBx(i,j) = .true.
  endif
  if( y(j) >= groove_d2 .and. xM2(i) >= groove_w1 .and. xM2(i) <= groove_w2)then ! False in the groove
   FBx(i,j) = .false.
  endif
  if( xM2(i) >= slit_w1 .and. xM2(i) <= slit_w2)then ! False in the slit
   FBx(i,j) = .false.
  endif 

 enddo
enddo

do i=1,Nx
 do j=1,N_loc

  if( yM2(j) >= Ag1 .and. yM2(j) <= Ag2 )then ! True in the metal zone
   FBy(i,j) = .true.
  endif
  if( yM2(j) >= groove_d2 .and. x(i) >= groove_w1 .and. x(i) <= groove_w2)then ! False in the groove
   FBy(i,j) = .false.
  endif
  if( x(i) >= slit_w1 .and. x(i) <= slit_w2)then ! False in the slit
   FBy(i,j) = .false.
  endif 

 enddo
enddo

do j = 1,N_loc
 do i = 1,Nx
  if(FBx(i,j))then
   FB_Grid(i,j) = 1
  else
   FB_Grid(i,j) = 0
  endif
 enddo
enddo

mwT = floor((ywT + (yM-y0)/2)/dx/N_loc) !floor(j_glob/N_loc)
jwT = 1 + (real((ywT + (yM-y0)/2)/dx/N_loc) - floor((ywT + (yM-y0)/2)/dx/N_loc))*N_loc !remainder(j_glob/N_loc)*N_loc

!-----------------------------------
!----------- Grid Return -----------
!-----------------------------------

!mReturn = nprocs/2
!
!itag = nprocs + 1
!do a = 0,nprocs-1
! itag = itag + 1
! if(myrank == a.and.a /= mReturn)then
!  call MPI_Send(FB_Grid(1,1), (Nx-1)*N_loc, MPI_Integer, mReturn, itag, MPI_COMM_WORLD,ierr)
! elseif(myrank == mReturn .and. a /= mReturn)then
!  call MPI_Recv(Drude_Grid(1,a*N_loc+1), (Nx-1)*N_loc, MPI_INTEGER, a, itag, MPI_COMM_WORLD,istatus,ierr)
! endif
!enddo
!
!if(myrank == mReturn)then
! do j = 1,N_loc
!  do i = 1,Nx
!   Drude_Grid(i,mReturn*N_loc+j) = FB_Grid(i,j)
!  enddo
! enddo
!endif
!
!if(myrank == mReturn)then
! open(file='Drude_Grid.dat',unit=42)
!  write(42,*) Drude_Grid
! close(unit=42)
!endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~ CPML vectors ~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sigmaCPML=0.8*(m+1)/(dx*(mu0/eps0*eps_delectric)**0.5)
alphaCPML=0.05
kappaCPML=5.0
!~~~ set CPML vectors ~~~!
do i=1,npml
 sige_x(i)=sigmaCPML*((npml-i)/(npml-1.0))**m
 alphae_x(i)=alphaCPML*((i-1.0)/(npml-1.0))**ma
 kappae_x(i)=1.0+(kappaCPML-1.0)*((npml-i)/(npml-1.0))**m
 be_x(i)=exp(-(sige_x(i)/kappae_x(i)+alphae_x(i))*dt/eps0)
 if( &
    (sige_x(i)==0.0).and. &
    (alphae_x(i)==0.0).and. & 
    (i==npml) &
   )then
   ce_x(i)=0.0
  else
   ce_x(i)=sige_x(i)*(be_x(i)-1.0)/(sige_x(i)+kappae_x(i)*alphae_x(i))/ kappae_x(i)
 endif
enddo

do i=1,npml-1
 sigh_x(i)=sigmaCPML*((npml-i-0.5)/(npml-1.0))**m
 alphah_x(i)=alphaCPML*((i-0.5)/(npml-1.0))**ma
 kappah_x(i)=1.0+(kappaCPML-1.0)*((npml-i-0.5)/(npml-1.0))**m
 bh_x(i)=exp(-(sigh_x(i)/kappah_x(i)+alphah_x(i))*dt/eps0)
 ch_x(i)=sigh_x(i)*(bh_x(i)-1.0)/(sigh_x(i)+kappah_x(i)*alphah_x(i))/kappah_x(i)
enddo

do j=1,npml
 sige_y(j)=sigmaCPML*((npml-j)/(npml-1.0))**m
 alphae_y(j)=alphaCPML*((j-1)/(npml-1.0))**ma
 kappae_y(j)=1.0+(kappaCPML-1.0)*((npml-j)/(npml-1.0))**m
 be_y(j)=exp(-(sige_y(j)/kappae_y(j)+alphae_y(j))*dt/eps0)
 if( &
    (sige_y(j)==0.0).and.&
    (alphae_y(j)==0.0).and. &
    (j==npml) &
   )then
   ce_y(j)=0.0
  else
   ce_y(j)=sige_y(j)*(be_y(j)-1.0)/(sige_y(j)+kappae_y(j)*alphae_y(j))/kappae_y(j)
 endif
enddo
   
do j=1,npml-1
 sigh_y(j)=sigmaCPML*((npml-j-0.5)/(npml-1.0))**m
 alphah_y(j)=alphaCPML*((j-0.5)/(npml-1.0))**ma
 kappah_y(j)=1.0+(kappaCPML-1.0)*((npml-j-0.5)/(npml-1.0))**m
 bh_y(j)=exp(-(sigh_y(j)/kappah_y(j)+alphah_y(j))*dt/eps0)
 ch_y(j)=sigh_y(j)*(bh_y(j)-1.0)/(sigh_y(j)+kappah_y(j)*alphah_y(j))/kappah_y(j)
enddo

den_hy=1.0/dy
if(myrank==0)then
  do j=1,N_loc
   if(j<=(npml-1))then
    den_hy(j)=1.0/(kappah_y(j)*dy)
   endif
  enddo
 elseif(myrank==(nprocs-1))then
  jj=npml-1
  do j=1,(N_loc-1)
   if(j>=((N_loc-1) - (npml-2)))then
     den_hy(j)=1.0/(kappah_y(jj)*dy)
     jj=jj-1
   endif
  enddo
endif

den_ey=1.0/dy
if(myrank==0)then
  do j=1,N_loc
   if(j<=npml)then
    den_ey(j)=1.0/(kappae_y(j)*dy)
   endif
  enddo
 elseif(myrank==(nprocs-1))then
  jj=npml
  do j=1,(N_loc-1)
   if(j>=((N_loc-1) - (npml-2)))then
     den_ey(j)=1.0/(kappae_y(jj)*dy)
     jj=jj-1
   endif
  enddo
endif

ii=npml-1
do i=1,Nx-1
 if(i<=(npml-1))then
   den_hx(i)=1.0/(kappah_x(i)*dx)
  elseif(i>=((Nx-1) - (npml-2)))then
   den_hx(i)=1.0/(kappah_x(ii)*dx)
   ii=ii-1
  else
   den_hx(i)=1.0/dx
 endif
enddo

ii=npml
do i=1,Nx-1
 if(i<=npml)then
   den_ex(i)=1.0/(kappae_x(i)*dx)
  elseif (i>=((Nx-1) - (npml-2)))then
   den_ex(i)=1.0/(kappae_x(ii)*dx)
   ii=ii-1
  else
   den_ex(i)=1.0/dx
 endif
enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~ end of CPML ~~~!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

Ex=0.0
Ey=0.0
Hz=0.0

PDx=0.0
PDy=0.0

psi_Hzy_1=0.0
psi_Exy_1=0.0
psi_Hzy_2=0.0
psi_Exy_2=0.0
psi_Hzx_1=0.0
psi_Eyx_1=0.0
psi_Hzx_2=0.0
psi_Eyx_2=0.0

Ex_get=0.0
Ex_send=0.0
Hz_get=0.0
Hz_send=0.0

Ex_inc=0.0
Hz_inc=0.0

Ex_get_inc=0.0
Ex_send_inc=0.0
Hz_get_inc=0.0
Hz_send_inc=0.0

psi_Hzy_1_inc=0.0
psi_Exy_1_inc=0.0
psi_Hzy_2_inc=0.0
psi_Exy_2_inc=0.0

Ex_temp=0.0
Hz_temp=0.0

Ex_temp_inc=0.0
Hz_temp_inc=0.0

if(myrank==mwT)then
 call cpu_time(cpu1)
 write(*,*) 'iw1 =', iw1
 write(*,*) 'iw2 =', iw2
 write(*,*) 'mwT =', mwT
 write(*,*) 'jwT =', jwT
endif

do n=1,Nt
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Hz ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
 do j=1,(nprocs-1) ! Send/Receive Ex
  itag=j
  itag1=j+1
  if(myrank==j)then
    do i=1,(Nx-1)
     Ex_send(i)=Ex(i,1)
    enddo    
    Ex_send_inc=Ex_inc(1)
    call MPI_SEND(Ex_send,(Nx-1),MPI_doUBLE_PRECISION,(j-1),itag,MPI_COMM_WORLD,ierr)
    call MPI_SEND(Ex_send_inc,1,MPI_doUBLE_PRECISION,(j-1),itag1,MPI_COMM_WORLD,ierr)
   elseif(myrank==(j-1))then
    call MPI_RECV(Ex_get,(Nx-1),MPI_doUBLE_PRECISION,j,itag,MPI_COMM_WORLD,istatus,ierr)
    call MPI_RECV(Ex_get_inc,1,MPI_doUBLE_PRECISION,j,itag1,MPI_COMM_WORLD,istatus,ierr)
  endif!send_recv
 enddo!nprocs

!-----------------------------
!------------- 0 -------------
!-----------------------------

if(myrank==0)then
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_hy(j))
  enddo
    
  j=N_loc
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex_get(i)-Ex(i,j))*den_hy(j))
 enddo

 do j=1,N_loc
! Left PML, Hz
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_x(i)*psi_Hzx_1(i,j)+ch_x(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
! Right PML, Hz
  ii=npml-1
  do i=(Nx-1) - (npml-2),Nx-1
   psi_Hzx_2(ii,j)=bh_x(ii)*psi_Hzx_2(ii,j)+ch_x(ii)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_2(ii,j)
   ii=ii-1
  enddo
 enddo

 do i=1,Nx-1
! Bottom PML, Hz
  do j=1,npml-1
   psi_Hzy_1(i,j)=bh_y(j)*psi_Hzy_1(i,j)+ch_y(j)*(Ex(i,j+1)-Ex(i,j))/dy
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzy_1(i,j)
  enddo
 enddo

!~~~~ incident ~~~~!
 do j=1,N_loc-1
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_inc(j+1)-Ex_inc(j))*den_hy(j)
 enddo
    
 j=N_loc
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_get_inc-Ex_inc(j))*den_hy(j)

! Bottom PML, Hz_inc
 do j=1,npml-1
  psi_Hzy_1_inc(j)=bh_y(j)*psi_Hzy_1_inc(j)+ch_y(j)*(Ex_inc(j+1)-Ex_inc(j))/dy
  Hz_inc(j)=Hz_inc(j)+dt_mu0*psi_Hzy_1_inc(j)
 enddo
endif

!----------------------------------
!------------ Interior ------------
!----------------------------------

if((myrank>0).and.(myrank<(nprocs-1)))then
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_hy(j))
  enddo
    
  j=N_loc
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex_get(i)-Ex(i,j))*den_hy(j))
 enddo

 do j=1,N_loc
! Left PML, Hz
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_x(i)*psi_Hzx_1(i,j)+ch_x(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
! Right PML, Hz
  ii=npml-1
  do i=(Nx-1) - (npml-2),Nx-1
   psi_Hzx_2(ii,j)=bh_x(ii)*psi_Hzx_2(ii,j)+ch_x(ii)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_2(ii,j)
   ii=ii-1
  enddo
 enddo

!~~~~ incident ~~~~!
 do j=1,N_loc-1
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_inc(j+1)-Ex_inc(j))*den_hy(j)
 enddo
    
 j=N_loc
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_get_inc-Ex_inc(j))*den_hy(j)

! scattered/total field updates
 if(myrank==mj0)then
  do i=i0,i1-1
   Hz(i,j0-1)=Hz(i,j0-1)-dt_mu0*Ex_inc(j0)/dy
  enddo
 endif
 
 if(myrank==mj1)then
  do i=i0,i1-1
   Hz(i,j1)=Hz(i,j1)+dt_mu0*Ex_inc(j1)/dy
  enddo
 endif
endif

!----------------------------------
!------------ nprocs-1 ------------
!----------------------------------

if(myrank==(nprocs-1))then !rank=(nprocs-1), here PML in y-direction is only applied to the top part
 do i=1,Nx-1
  do j=1,N_loc-1
   Hz(i,j)=Hz(i,j)+dt_mu0*((Ey(i,j)-Ey(i+1,j))*den_hx(i)+ &
			               (Ex(i,j+1)-Ex(i,j))*den_hy(j))
  enddo
 enddo
 
 do j=1,N_loc
! Left PML, Hz
  do i=1,npml-1
   psi_Hzx_1(i,j)=bh_x(i)*psi_Hzx_1(i,j)+ch_x(i)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_1(i,j)
  enddo
! Right PML, Hz
  ii=npml-1
  do i=(Nx-1) - (npml-2),Nx-1
   psi_Hzx_2(ii,j)=bh_x(ii)*psi_Hzx_2(ii,j)+ch_x(ii)*(Ey(i,j)-Ey(i+1,j))/dx
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzx_2(ii,j)
   ii=ii-1
  enddo
 enddo

! Top PML, Hz
 do i=1,Nx-1    
  jj=npml-1
  do j=(N_loc-1) - (npml-2),N_loc-1
   psi_Hzy_2(i,jj)=bh_y(jj)*psi_Hzy_2(i,jj)+ch_y(jj)*(Ex(i,j+1)-Ex(i,j))/dy
   Hz(i,j)=Hz(i,j)+dt_mu0*psi_Hzy_2(i,jj)
   jj=jj-1
  enddo
 enddo

!~~~~ incident ~~~~!
 do j=1,N_loc-1
  Hz_inc(j)=Hz_inc(j)+dt_mu0*(Ex_inc(j+1)-Ex_inc(j))*den_hy(j)
 enddo

! Top PML, Hz_inc
 jj=npml-1
 do j=(N_loc-1) - (npml-2),N_loc-1
  psi_Hzy_2_inc(jj)=bh_y(jj)*psi_Hzy_2_inc(jj)+ch_y(jj)*(Ex_inc(j+1)-Ex_inc(j))/dy
  Hz_inc(j)=Hz_inc(j)+dt_mu0*psi_Hzy_2_inc(jj)
  jj=jj-1
 enddo
endif

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
 do j=1,(nprocs-1) ! Send/Receive Hz
  itag=j
  itag1=j+1
  if(myrank==(j-1))then
    do i=1,Nx-1
     Hz_send(i)=Hz(i,N_loc)
    enddo
    Hz_send_inc=Hz_inc(N_loc)
    call MPI_SEND(Hz_send,(Nx-1),MPI_doUBLE_PRECISION,j,itag,MPI_COMM_WORLD,ierr)
    call MPI_SEND(Hz_send_inc,1,MPI_doUBLE_PRECISION,j,itag1,MPI_COMM_WORLD,ierr)
   elseif(myrank==j)then
    call MPI_RECV(Hz_get,(Nx-1),MPI_doUBLE_PRECISION,(j-1),itag,MPI_COMM_WORLD,istatus,ierr)
    call MPI_RECV(Hz_get_inc,1,MPI_doUBLE_PRECISION,(j-1),itag1,MPI_COMM_WORLD,istatus,ierr)
  endif!send_recv
 enddo!nprocs
 
!-----------------------------
!------------- 0 -------------
!-----------------------------

if(myrank==0)then
 do i=1,Nx-1
  do j=2,N_loc
   Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))*den_ey(j)
  enddo
 enddo
   
! Bottom PML, Ex
 do i=1,Nx-1
  do j=2,npml
   psi_Exy_1(i,j)=be_y(j)*psi_Exy_1(i,j)+ce_y(j)*(Hz(i,j)-Hz(i,j-1))/dy
   Ex(i,j)=Ex(i,j)+dt_eps0*psi_Exy_1(i,j)
  enddo
 enddo

!~~~ incident ~~~!
 do j=2,N_loc
  Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))*den_ey(j)
 enddo
   
! Bottom PML, Ex_inc
 do j=2,npml
  psi_Exy_1_inc(j)=be_y(j)*psi_Exy_1_inc(j)+ce_y(j)*(Hz_inc(j)-Hz_inc(j-1))/dy
  Ex_inc(j)=Ex_inc(j)+dt_eps0*psi_Exy_1_inc(j)
 enddo
endif

!----------------------------------
!------------ Interior ------------
!----------------------------------

if((myrank>0).and.(myrank<(nprocs-1)))then !no PML for y-direction here
 do i=1,Nx-1
  j=1
   if(FBx(i,j))then
     tmpE=Ca*Ex(i,j)+Cb*(Hz(i,j)-Hz_get(i))*den_ey(j)-Cc*PDx(i,j)
     PDx(i,j)=A1*PDx(i,j)+A2*(tmpE+Ex(i,j))
     Ex(i,j)=tmpE
    else
     Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz_get(i))*den_ey(j)
   endif
  
  do j=2,N_loc
   if(FBx(i,j))then
     tmpE=Ca*Ex(i,j)+Cb*(Hz(i,j)-Hz(i,j-1))*den_ey(j)-Cc*PDx(i,j)
     PDx(i,j)=A1*PDx(i,j)+A2*(tmpE+Ex(i,j))
     Ex(i,j)=tmpE
    else
     Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))*den_ey(j)
   endif
  enddo
 enddo

!~~~ incident ~~~!
 j=1
  Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_get_inc)*den_ey(j)
   
 do j=2,N_loc
  if((myrank==ms).and.(j==js))then !laser pulse
    Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))/dy+ &
 	                    pulse(n)
   else
	Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))*den_ey(j)
  endif
 enddo

! scattered/total field updates
 if(myrank==mj0)then
  do i=i0,i1-1
   Ex(i,j0)=Ex(i,j0)-dt_eps0*Hz_inc(j0-1)/dy
  enddo
 endif

 if(myrank==mj1)then
  do i=i0,i1-1
   Ex(i,j1)=Ex(i,j1)+dt_eps0*Hz_inc(j1)/dy
  enddo
 endif
endif

!----------------------------------
!------------ nprocs-1 ------------
!----------------------------------

if(myrank==(nprocs-1))then
 j=1
  do i=1,Nx-1
   Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz_get(i))*den_ey(j)
  enddo
   
 do i=1,Nx-1
  do j=2,(N_loc-1)
   Ex(i,j)=Ex(i,j)+dt_eps0*(Hz(i,j)-Hz(i,j-1))*den_ey(j)
  enddo
 enddo

! Top PML, Ex
 do i=1,Nx-1
  jj=npml
  do j=(N_loc-1) - (npml-2),N_loc-1
   psi_Exy_2(i,jj)=be_y(jj)*psi_Exy_2(i,jj)+ce_y(jj)*(Hz(i,j)-Hz(i,(j-1)))/dy
   Ex(i,j)=Ex(i,j)+dt_eps0*psi_Exy_2(i,jj)
   jj=jj-1
  enddo
 enddo

!~~~ incident ~~~!
 j=1
  Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_get_inc)*den_ey(j)
   
 do j=2,N_loc-1
  Ex_inc(j)=Ex_inc(j)+dt_eps0*(Hz_inc(j)-Hz_inc(j-1))*den_ey(j)
 enddo

! Top PML, Ex_inc
 jj=npml
 do j=(N_loc-1) - (npml-2),N_loc-1
  psi_Exy_2_inc(jj)=be_y(jj)*psi_Exy_2_inc(jj)+ce_y(jj)*(Hz_inc(j)-Hz_inc(j-1))/dy
  Ex_inc(j)=Ex_inc(j)+dt_eps0*psi_Exy_2_inc(jj)
  jj=jj-1
 enddo
endif

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Ey ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

!----------------------------------------
!---------- 0 through nprocs-2 ----------
!----------------------------------------

if((myrank>=0).and.(myrank<(nprocs-1)))then
 do i=2,Nx-1
  do j=1,N_loc
   if(FBy(i,j))then
     tmpE=Ca*Ey(i,j)+Cb*(Hz(i-1,j)-Hz(i,j))*den_ex(i)-Cc*PDy(i,j)
     PDy(i,j)=A1*PDy(i,j)+A2*(tmpE+Ey(i,j))
     Ey(i,j)=tmpE
	else
	 Ey(i,j)=Ey(i,j)+dt_eps0*(Hz(i-1,j)-Hz(i,j))*den_ex(i)
   endif
  enddo
 enddo

 do j=1,N_loc
! Left PML, Ey
  do i=2,npml
   if(FBy(i,j))then
    psi_Eyx_1(i,j)=be_x(i)*psi_Eyx_1(i,j)+ce_x(i)*(Hz(i-1,j)-Hz(i,j))/dx
    PDy(i,j) = PDy(i,j) + A2*(Cb*psi_Eyx_1(i,j))
    Ey(i,j) = Ey(i,j) + Cb*psi_Eyx_1(i,j)
   else
    psi_Eyx_1(i,j)=be_x(i)*psi_Eyx_1(i,j)+ce_x(i)*(Hz(i-1,j)-Hz(i,j))/dx
    Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_1(i,j)
   endif
  enddo
! Right PML, Ey
 ii=npml
 do i=(Nx-1) - (npml-2),Nx-1
  if(FBy(i,j))then
   psi_Eyx_2(ii,j)=be_x(ii)*psi_Eyx_2(ii,j)+ce_x(ii)*(Hz(i-1,j)-Hz(i,j))/dx
   PDy(i,j) = PDy(i,j) + A2*(Cb*psi_Eyx_2(ii,j))
   Ey(i,j) = Ey(i,j) + Cb*psi_Eyx_2(ii,j)
  else  
   psi_Eyx_2(ii,j)=be_x(ii)*psi_Eyx_2(ii,j)+ce_x(ii)*(Hz(i-1,j)-Hz(i,j))/dx
   Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_2(ii,j)
  endif
  
  ii=ii-1
  enddo
 enddo

! scattered/total field updates
 if(myrank==mj0)then
  do j=j0,N_loc
   Ey(i0,j)=Ey(i0,j)+dt_eps0*Hz_inc(j)/dy
   Ey(i1,j)=Ey(i1,j)-dt_eps0*Hz_inc(j)/dy
  enddo
 endif

 if((myrank>mj0).and.(myrank<mj1))then
  do j=1,N_loc
   Ey(i0,j)=Ey(i0,j)+dt_eps0*Hz_inc(j)/dy
   Ey(i1,j)=Ey(i1,j)-dt_eps0*Hz_inc(j)/dy
  enddo
 endif

 if(myrank==mj1)then
  do j=1,j1-1
   Ey(i0,j)=Ey(i0,j)+dt_eps0*Hz_inc(j)/dy
   Ey(i1,j)=Ey(i1,j)-dt_eps0*Hz_inc(j)/dy
  enddo
 endif
endif

!----------------------------------
!------------ nprocs-1 ------------
!----------------------------------

if(myrank==(nprocs-1))then
 do i=2,Nx-1
  do j=1,N_loc-1
   Ey(i,j)=Ey(i,j)+dt_eps0*((Hz(i-1,j)-Hz(i,j))*den_ex(i))
  enddo
 enddo

 do j=1,N_loc-1
! Left PML, Ey
  do i=2,npml
   psi_Eyx_1(i,j)=be_x(i)*psi_Eyx_1(i,j)+ce_x(i)*(Hz(i-1,j)-Hz(i,j))/dx
   Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_1(i,j)
  enddo
! Right PML, Ey
  ii=npml
  do i=(Nx-1) - (npml-2),Nx-1
   psi_Eyx_2(ii,j)=be_x(ii)*psi_Eyx_2(ii,j)+ce_x(ii)*(Hz(i-1,j)-Hz(i,j))/dx
   Ey(i,j)=Ey(i,j)+dt_eps0*psi_Eyx_2(ii,j)
   ii=ii-1
  enddo
 enddo
endif


!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~==========================~~~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~         detection          ~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!
!~~~~~~~~~~~~~~~~~~~~~~~~==========================~~~~~~~~~~~~~~~~~~~~~~~~!
!--------------------------------------------------------------------------!
!if(myrank==mwR)then
! j=jwR
! 
! do nn=1,N_w
!  tmp1=sin(omega_P(nn)*dt*dble(n))
!  tmp2=cos(omega_P(nn)*dt*dble(n))
!  do i=iw1,iw2
!   Ex_temp(i,nn,1)=Ex_temp(i,nn,1)+tmp1*(Ex(i-1,j)+Ex(i,j))/2.0
!   Ex_temp(i,nn,2)=Ex_temp(i,nn,2)+tmp2*(Ex(i-1,j)+Ex(i,j))/2.0
!  
!   Hz_temp(i,nn,1)=Hz_temp(i,nn,1)+tmp1*(Hz(i-1,j)+Hz(i,j)+Hz(i-1,j-1)+Hz(i,j-1))/4.0
!   Hz_temp(i,nn,2)=Hz_temp(i,nn,2)+tmp2*(Hz(i-1,j)+Hz(i,j)+Hz(i-1,j-1)+Hz(i,j-1))/4.0
!
!   Ex_temp_inc(i,nn,1)=Ex_temp_inc(i,nn,1)+tmp1*Ex_inc(j)
!   Ex_temp_inc(i,nn,2)=Ex_temp_inc(i,nn,2)+tmp2*Ex_inc(j)
!  
!   Hz_temp_inc(i,nn,1)=Hz_temp_inc(i,nn,1)+tmp1*(Hz_inc(j-1)+Hz_inc(j))/2.0
!   Hz_temp_inc(i,nn,2)=Hz_temp_inc(i,nn,2)+tmp2*(Hz_inc(j-1)+Hz_inc(j))/2.0
!  enddo
! enddo
!endif

if(myrank==mwT)then
 j=jwT

 if(jwT /= 1)then
  do nn=1,N_w
   tmp1=sin(omega_P(nn)*dt*dble(n))
   tmp2=cos(omega_P(nn)*dt*dble(n))
   do i=iw1,iw2
    Ex_temp(i,nn,1)=Ex_temp(i,nn,1)+tmp1*(Ex(i-1,j)+Ex(i,j))/2.0
    Ex_temp(i,nn,2)=Ex_temp(i,nn,2)+tmp2*(Ex(i-1,j)+Ex(i,j))/2.0
  
    Hz_temp(i,nn,1)=Hz_temp(i,nn,1)+tmp1*(Hz(i-1,j)+Hz(i,j)+Hz(i-1,j-1)+Hz(i,j-1))/4.0
    Hz_temp(i,nn,2)=Hz_temp(i,nn,2)+tmp2*(Hz(i-1,j)+Hz(i,j)+Hz(i-1,j-1)+Hz(i,j-1))/4.0

    Ex_temp_inc(i,nn,1)=Ex_temp_inc(i,nn,1)+tmp1*Ex_inc(j)
    Ex_temp_inc(i,nn,2)=Ex_temp_inc(i,nn,2)+tmp2*Ex_inc(j)
  
    Hz_temp_inc(i,nn,1)=Hz_temp_inc(i,nn,1)+tmp1*(Hz_inc(j-1)+Hz_inc(j))/2.0
    Hz_temp_inc(i,nn,2)=Hz_temp_inc(i,nn,2)+tmp2*(Hz_inc(j-1)+Hz_inc(j))/2.0
   enddo
  enddo
 else
  do nn=1,N_w
   tmp1=sin(omega_P(nn)*dt*dble(n))
   tmp2=cos(omega_P(nn)*dt*dble(n))
   do i=iw1,iw2
    Ex_temp(i,nn,1)=Ex_temp(i,nn,1)+tmp1*(Ex(i-1,j)+Ex(i,j))/2.0
    Ex_temp(i,nn,2)=Ex_temp(i,nn,2)+tmp2*(Ex(i-1,j)+Ex(i,j))/2.0
  
    Hz_temp(i,nn,1)=Hz_temp(i,nn,1)+tmp1*(Hz(i-1,j)+Hz(i,j)+Hz_get(i-1)+Hz_get(i))/4.0
    Hz_temp(i,nn,2)=Hz_temp(i,nn,2)+tmp2*(Hz(i-1,j)+Hz(i,j)+Hz_get(i-1)+Hz_get(i))/4.0

    Ex_temp_inc(i,nn,1)=Ex_temp_inc(i,nn,1)+tmp1*Ex_inc(j)
    Ex_temp_inc(i,nn,2)=Ex_temp_inc(i,nn,2)+tmp2*Ex_inc(j)
  
    Hz_temp_inc(i,nn,1)=Hz_temp_inc(i,nn,1)+tmp1*(Hz_get_inc+Hz_inc(j))/2.0
    Hz_temp_inc(i,nn,2)=Hz_temp_inc(i,nn,2)+tmp2*(Hz_get_inc+Hz_inc(j))/2.0
   enddo
  enddo
 endif
  
   if(mod(n,10000)==0)then
		call cpu_time(cpu2)
		write(*,*) 'cpu time [total hours]',(cpu2-cpu1)*(Nt/10000)/3600
		cpu1=cpu2
		write(*,*) n,Ex(Nx/2,N_loc/2),Ey(Nx/2,N_loc/2)
   endif
endif

enddo !Nt

!if(myrank==mwR)then
! do nn=1,N_w
!  do i=iw1,iw2
!   tmp2=sqrt(Ex_temp(i,nn,1)**2+Ex_temp(i,nn,2)**2)/SN(nn,1)
!   Ex_temp(i,nn,2)=-atan2(Ex_temp(i,nn,1),Ex_temp(i,nn,2))-SN(nn,2)
!   Ex_temp(i,nn,1)=tmp2
!    
!   tmp2=sqrt(Hz_temp(i,nn,1)**2+Hz_temp(i,nn,2)**2)/SN(nn,1)
!   Hz_temp(i,nn,2)=-atan2(Hz_temp(i,nn,1),Hz_temp(i,nn,2))-SN(nn,2)
!   Hz_temp(i,nn,1)=tmp2
!  enddo
!  
!  sum_int=(0.0D0,0.0D0)
!  do i=iw1,iw2
!   cTMP1=Ex_temp(i,nn,1)*exp(-Im*Ex_temp(i,nn,2))
!   cTMP2=Hz_temp(i,nn,1)*exp(-Im*Hz_temp(i,nn,2))
!   sum_int=sum_int+cTMP1*conjg(cTMP2)
!  enddo
!  P_sct(nn)=dreal(sum_int)
!
!
!  do i=iw1,iw2
!   tmp2=sqrt(Ex_temp_inc(i,nn,1)**2+Ex_temp_inc(i,nn,2)**2)/SN(nn,1)
!   Ex_temp_inc(i,nn,2)=-atan2(Ex_temp_inc(i,nn,1),Ex_temp_inc(i,nn,2))-SN(nn,2)
!   Ex_temp_inc(i,nn,1)=tmp2
!    
!   tmp2=sqrt(Hz_temp_inc(i,nn,1)**2+Hz_temp_inc(i,nn,2)**2)/SN(nn,1)
!   Hz_temp_inc(i,nn,2)=-atan2(Hz_temp_inc(i,nn,1),Hz_temp_inc(i,nn,2))-SN(nn,2)
!   Hz_temp_inc(i,nn,1)=tmp2
!  enddo
!  
!  sum_int=(0.0D0,0.0D0)
!  do i=iw1,iw2
!   cTMP1=Ex_temp_inc(i,nn,1)*exp(-Im*Ex_temp_inc(i,nn,2))
!   cTMP2=Hz_temp_inc(i,nn,1)*exp(-Im*Hz_temp_inc(i,nn,2))
!   sum_int=sum_int+cTMP1*conjg(cTMP2)
!  enddo
!  P_inc(nn)=dreal(sum_int)
! enddo
!
! open(file='R_2D_1slit-ADE_Drude-R_83nm-t_150nm-length_200nm.dat',unit=32)
! do nn=1,N_w
!  write(32,*) omega_P(nn)/ev_to_radsec,abs(P_sct(nn)/P_inc(nn))
! enddo
! close(unit=32)
!endif

if(myrank==mwT)then
 do nn=1,N_w
  do i=iw1,iw2
   tmp2=sqrt(Ex_temp(i,nn,1)**2+Ex_temp(i,nn,2)**2)/SN(nn,1)
   Ex_temp(i,nn,2)=-atan2(Ex_temp(i,nn,1),Ex_temp(i,nn,2))-SN(nn,2)
   Ex_temp(i,nn,1)=tmp2
    
   tmp2=sqrt(Hz_temp(i,nn,1)**2+Hz_temp(i,nn,2)**2)/SN(nn,1)
   Hz_temp(i,nn,2)=-atan2(Hz_temp(i,nn,1),Hz_temp(i,nn,2))-SN(nn,2)
   Hz_temp(i,nn,1)=tmp2
  enddo
  
  sum_int=(0.0D0,0.0D0)
  do i=iw1,iw2
   cTMP1=Ex_temp(i,nn,1)*exp(-Im*Ex_temp(i,nn,2))
   cTMP2=Hz_temp(i,nn,1)*exp(-Im*Hz_temp(i,nn,2))
   sum_int=sum_int+cTMP1*conjg(cTMP2)
  enddo
  P_sct(nn)=dreal(sum_int)


  do i=iw1,iw2
   tmp2=sqrt(Ex_temp_inc(i,nn,1)**2+Ex_temp_inc(i,nn,2)**2)/SN(nn,1)
   Ex_temp_inc(i,nn,2)=-atan2(Ex_temp_inc(i,nn,1),Ex_temp_inc(i,nn,2))-SN(nn,2)
   Ex_temp_inc(i,nn,1)=tmp2
    
   tmp2=sqrt(Hz_temp_inc(i,nn,1)**2+Hz_temp_inc(i,nn,2)**2)/SN(nn,1)
   Hz_temp_inc(i,nn,2)=-atan2(Hz_temp_inc(i,nn,1),Hz_temp_inc(i,nn,2))-SN(nn,2)
   Hz_temp_inc(i,nn,1)=tmp2
  enddo
  
  sum_int=(0.0D0,0.0D0)
  do i=iw1,iw2
   cTMP1=Ex_temp_inc(i,nn,1)*exp(-Im*Ex_temp_inc(i,nn,2))
   cTMP2=Hz_temp_inc(i,nn,1)*exp(-Im*Hz_temp_inc(i,nn,2))
   sum_int=sum_int+cTMP1*conjg(cTMP2)
  enddo
  P_inc(nn)=dreal(sum_int)
 enddo

 open(file='T_2D_1slit-ADE_Drude-R_83nm-t_150nm-length_200nm.dat',unit=32)
 do nn=1,N_w
  write(32,*) omega_P(nn)/ev_to_radsec,abs(P_sct(nn)/P_inc(nn))
 enddo
 close(unit=32)
endif

!---------------------------------------------------------------------!
 call MPI_FINALIZE(ierr)
!---------------------------------------------------------------------!
end
