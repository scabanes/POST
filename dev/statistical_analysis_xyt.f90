PROGRAM spectra_analysis


use netcdf

implicit none

integer :: iarg,narg,mt,istep,getfield,tstep,dimid_out(3),varid_out(4),latid,colatid,lonid,tid,idinfos,ninfos,idninfos,NetcdefStatData=1,hnid,fnid,fnrotid,hnrotid,fndiv1id,hndiv1id,fndiv2id,hndiv2id
character (len=100), dimension(:), allocatable :: arg
character (len=100) :: tmp_char
logical :: is_file,NotFirst_reading=.false.,First_reading=.true.
!----------------------------------------------------------------------
integer :: nlon=0,nlat=0,nalt,ntime,fieldsmap=0,R_sat = 58232000 ! R_sat est le rayon de saturne en m ! fieldsmap=1 j enregistre les cartes
real*8 :: Amp = 1.0,kf=62.,ptimestep=4756.50,radius = 58232000
real*8  :: factor, factor2, n=0
double precision, dimension (:), allocatable :: j_i_conversion,j_1d,i_1d,j_1d_tm,i_1d_tm,f_1d,h_1d,f_1d_tm,h_1d_tm,fn_rot,hn_rot,fn_div1,hn_div1,fn_div2,hn_div2
double precision, dimension(:,:), allocatable :: vort_r, vort_i, div_r, div_i,vortf_r, vortf_i, divf_r, divf_i, sf_r, sf_i, vp_r, vp_i
double precision, dimension(:,:), allocatable :: vort,term1_r,term1_i,term2_v,term2_w, term3,vort_uv_br,vort_uv_bi,vort_uv_cr,vort_uv_ci,vortv,vortu,uv_dot_grad_vort_r,uv_dot_grad_vort_i,rot_ugradu_i,rot_ugradu_r,div1_ugradu_i,div1_ugradu_r
double precision, dimension(:,:), allocatable :: j_2d,i_2d,j_2d_tm,i_2d_tm


integer :: mdab_v,mdab_s,mdab_smaller
!----------------------------------------------------------------------
double precision, dimension(:,:), allocatable ::  uoff,voff,ureg,vreg,uoff_f,voff_f,ureg_f,vreg_f,u_DYN,v_DYN,ug,vg,vm,wm,fvm,fwm,UvectVT_lon,UvectVT_lat,u_rot,v_rot,u_div,v_div,udotu,vgrad,wgrad,vgradN,wgradN
double precision, dimension(:,:), allocatable :: fu,fv,rotuwP
double precision, dimension (:), allocatable :: lat_DYN,lon_DYN,colat_DYN,lat_SP,lon_SP,colat_SP
double precision, dimension(:,:), allocatable ::  br,bi,cr,ci,fbr,fbi,fcr,fci,Fdivuwr,Fdivuwi,Frotuwr,Frotuwi,Jmn,uur,uui,bbr,bbi,ccr,cci,ar,ai,divgraduur,divgraduui,rotgraduur,rotgraduui,div2_ugradu_r,div2_ugradu_i
character (len=100) :: file_netcdf,file_spectra,nickname
integer :: idfile,ncid,ierror,idu,idv,idlat,idlon,idvw,jdvw,mdab,ndab,ivrt,jvrt,mdc,ndc,nt=1,ityp=0,isym=0,ioff,fgg,lon_dimid,lat_dimid,t_dimid,wmid,vmid,vortid, dimids(3),dimids_1d(2), idvti,idfile2,idtime,idnlon,idnlat,idSu,fuid,fvid,urotid,vrotid,udivid,vdivid
integer :: i,j,it,itt,in,jm,ig,mmax,l1,ll1,l2,vlwork,lsav,lGwork,lvhsec,lwwork,llwwork,sldwork
double precision, dimension (:), allocatable :: gmwork,work,ddwork,vwork,wsav,Gwork
double precision, dimension (:,:,:), allocatable :: vm_3D,wm_3D,vt_3D
real(16), parameter :: PI=4.D0*DATAN(1.D0)


!###################################################################################################################
!######################################################################################################################
!###################################################################################################################
!**********
!Initialisation
mt=1
istep=1
tstep=1
getfield=3
!**********
!Input reading
narg = command_argument_count()
allocate(arg(narg))
do iarg = 1, narg
  call get_command_argument(iarg,arg(iarg))
end do
iarg = 1
do while (iarg .le. narg)
  if ( arg(iarg) == '-mt' .or. arg(iarg) == '-istp' .or. arg(iarg) == '-tstp' .or. arg(iarg) == '-field') then
    select case(arg(iarg))
    case('-mt')
      read(arg(iarg+1),'(i10)' ) mt !extent of temporal length
    case('-istp')
      read(arg(iarg+1),'(i10)' ) istep !name of time axis
    case('-field')
      read(arg(iarg+1),'(i10)' ) getfield !name of time axis
    case('-tstp')
      read(arg(iarg+1),'(i10)' ) tstep !name of time axis
    end select
    iarg = iarg + 2
  elseif (arg(iarg) == '-h' .or. arg(iarg) == '--help') then
      print*,'Usage'
      print*,'spectra_analysis netcdfFile [option]'
      print*,'[-h or --help]	: brief help'
      print*,'[-mt int]	: final temporal file created (default: 0)'
      print*,'[-istp int]	: initial temporal step to create files (default: 0)'
      print*,'[-field int]	: select the field to decompose: 1.zonalV 2.meridV 3.totalV 4.stream 5.Vpoten (default: 3)'
      print*,'[-tstp int]	: temporal step to create files (default: 1)'
      print*,'[-mz int]	: extent of vertical mean (default: 0)'
      print*,''
      print*,'!!! Note that variables name are: for 1) latitude: lat, 2) longitude: lon, 3) zonal velocity: u, 4) meridional velocity: v'
      print*,'Compute the kinetic energy spectrum of a wind field on a longitude-latitude grid.'
      print*,'The grid must have a redundant point in longitude.'
      print*,'Ideally the analysis works better on a symetric grid (2N points in longitude and N points in latitude).'
      print*,'The output text file (called spectra by default) give the decomposition'
      print*,'of the velocity on the vectorial spherical harmonic basis.'
      stop 'End help'
  else
      file_netcdf = arg(iarg)
      iarg = iarg + 1
  end if
end do
print*,'create files using steps in time t =',istep,':',tstep,':',mt
!**********
!Check input/output files
inquire(file=trim(file_netcdf),exist=is_file)
if (.not. is_file) then
  print*,"no netcdf file: ",trim(file_netcdf)
  stop "Stopped"
end if 
!###################################################################################################################
!===================================================================================================================
!						 LOAD IN_PUT FILE 
!						    WIND FIELD
!===================================================================================================================
!###################################################################################################################
! ------ > u_DYN & v_DYN: matrice issues des champs vents DYNAMICO
! ici j obtiens deux matrices des composantes en longitude et latitude
! u_DYN et v_DYN definies sur lat=[-89.75,89.75] et lon=[0.25,359.75] dans DYN.
! avec nlon = 720 et nlat = 360. Les donnees sont donc organisees
! south to north suivant les indices croissant en latitude. 
! soit : 1) v_DYN[nlon,nlat] & v_DYN[nlon,nlat]
!	 2) lat = [-89.75,89.75] and lon = [0.25,359.75]
!        3) South --> North , increasing lat index
!	 4) Latitudinal and longitudinal component

! ------- > vm & wm: vm is a nlat by nlon array containing the
! the colatitudinal vector component.  wm is a nlat by nlon array
! containing the east longitudinal vector component.  Values in
!(vm,wm) are ordered from the northern to the southern hemisphere
! with increasing colatitude subscript.
! soit : 1) vm[nlat+1,nlon] & wm[nlat+1,nlon]
!	 2) colat = [0,180] and lon = [0,359.5]  --> poles included
!        3) North --> South , increasing colat index
!	 4) Colatitudinal and longitudinal component (i.e. vm = -v_DYN)
 call check( nf90_open(trim(file_netcdf),NF90_NOWRITE,idfile))

 call check( nf90_inq_dimid(idfile,trim('time_counter'),idtime) )
 call check( nf90_inq_dimid(idfile,trim('lon'),idnlon) )
 call check( nf90_inq_dimid(idfile,trim('lat'),idnlat) )

 call check( nf90_inquire_dimension(idfile,idtime,tmp_char,ntime))
 call check( nf90_inquire_dimension(idfile,idnlon,tmp_char,nlon))
 call check( nf90_inquire_dimension(idfile,idnlat,tmp_char,nlat))

 call check( nf90_inq_varid(idfile,trim('lat'),idlat))
 call check( nf90_inq_varid(idfile,trim('lon'),idlon))

 allocate(lat_DYN(nlat))
 allocate(lon_DYN(nlon))
 call check( nf90_get_var(idfile,idlat,lat_DYN,(/1/),(/nlat/)) ) !(lon,lat,iz,it)
 call check( nf90_get_var(idfile,idlon,lon_DYN,(/1/),(/nlon/)) ) !(lon,lat,iz,it)
!###################################################################################################################
!=================================================================================================================== 
!						     forcage
!===================================================================================================================
!##################################################################################################################
! Attention ici, l integrale se fait 
! Int_0^pi Int_0^2pi f(theta,phi) sin(theta) dtheta dphi avec 
! theta et phi la colatitude/latitude et la longitude. Dans physiq_
! mod.f90 la variable est la colatitude dans les calculs suivant ont 
! utilise egalement la colatitude (i.e. colat). Ceci est important
! puisque cos(lat) ~ sin(colat). 
! l expression du forcage de Taylor-green dans dynamico est:
! f = Amp * cos(lat) cos(kf*phi) sin(kf*colat)/ptimestep .etheta
!     Amp * cos(lat) cos(kf*phi) sin(kf*colat)/ptimestep .ephi
!     0 .ez
!			ou
! f = Amp * sin(colat) cos(kf*phi) sin(kf*colat)/ptimestep .etheta
!     Amp * sin(colat) cos(kf*phi) sin(kf*colat)/ptimestep .ephi
!     0 .ez
allocate(fu(nlon,nlat))
allocate(fv(nlon,nlat))
allocate(colat_DYN(nlat))
colat_DYN = (lat_DYN/(180./PI)) + (PI/2.) ! Comme dans la routine physiq_mod.f90 dans dynaico ou le forcage est impose.
do j=1,nlon
   do i=1,nlat
	fu(j,i) = Amp*sin(colat_DYN(i))*cos(kf*(lon_DYN(j)/(180./PI)))*sin(kf*colat_DYN(i))/ptimestep
	fv(j,i) = Amp*sin(colat_DYN(i))*sin(kf*(lon_DYN(j)/(180./PI)))*cos(kf*colat_DYN(i))/ptimestep
   end do
end do
!######################################################################################################################
!##################################### CREATE OUT_PUT FILE ############################################################
!######################################################################################################################
if(NetcdefStatData == 1) then
  print *, "Create file StatisticalData.nc! "
! ---------------------------DEFINITION OF DIMENSIONS
  call check( nf90_create("StatisticalData.nc", NF90_CLOBBER, ncid))
  call check( nf90_def_dim(ncid, "nlon", nlon, lon_dimid))
  call check( nf90_def_dim(ncid, "nlat", nlat+1, lat_dimid))  !!!!! to change
 ! call check( nf90_def_dim(ncid, "nlat", nlat, lat_dimid))
  call check( nf90_def_dim(ncid, "time", mt, t_dimid))
  dimids =  (/ lon_dimid, lat_dimid, t_dimid /)!
  dimids_1d =  (/ lat_dimid, t_dimid /)!
! ---------------------------DFINITION OF VARIBLES
  call check( nf90_def_var(ncid, "lon", NF90_DOUBLE, lon_dimid, lonid))
  call check( nf90_def_var(ncid, "lat", NF90_DOUBLE, lat_dimid, latid))
  call check( nf90_def_var(ncid, "time", NF90_REAL4, t_dimid, tid))
  call check( nf90_def_var(ncid, "colat", NF90_DOUBLE, lat_dimid, colatid))
! champ 1d
  call check( nf90_def_var(ncid, "fn_1d", NF90_DOUBLE, dimids_1d, fnid))
  call check( nf90_def_var(ncid, "hn_1d", NF90_DOUBLE, dimids_1d, hnid))
  !call check( nf90_def_var(ncid, "fn_rot", NF90_DOUBLE, dimids_1d, fnrotid))
  !call check( nf90_def_var(ncid, "hn_rot", NF90_DOUBLE, dimids_1d, hnrotid))
  !call check( nf90_def_var(ncid, "fn_div1", NF90_DOUBLE, dimids_1d, fndiv1id))
  !call check( nf90_def_var(ncid, "hn_div1", NF90_DOUBLE, dimids_1d, hndiv1id))
  !call check( nf90_def_var(ncid, "fn_div2", NF90_DOUBLE, dimids_1d, fndiv2id))
  !call check( nf90_def_var(ncid, "hn_div2", NF90_DOUBLE, dimids_1d, hndiv2id))
! champ 2d
  !call check( nf90_def_var(ncid, "wm", NF90_DOUBLE, dimids, wmid))
  !call check( nf90_def_var(ncid, "vm", NF90_DOUBLE, dimids, vmid))
  !call check( nf90_def_var(ncid, "vort", NF90_DOUBLE, dimids, vortid))
  !call check( nf90_def_var(ncid, "u_rot", NF90_DOUBLE, dimids, urotid))
  !call check( nf90_def_var(ncid, "v_rot", NF90_DOUBLE, dimids, vrotid))
  !call check( nf90_def_var(ncid, "u_div", NF90_DOUBLE, dimids, udivid))
  !call check( nf90_def_var(ncid, "v_div", NF90_DOUBLE, dimids, vdivid))
  
  !call check( nf90_def_var(ncid, "fu", NF90_DOUBLE, dimids, fuid))
  !call check( nf90_def_var(ncid, "fv", NF90_DOUBLE, dimids, fvid))
! ---------------------------VARIABLES FILL INSTRUCTIONS
! champ 1d
  call check( NF90_PUT_ATT  (ncid, fnid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncid, hnid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, fnrotid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, hnrotid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, hndiv1id, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, fndiv1id, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, hndiv2id, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, fndiv2id, "_FillValue", NF90_FILL_DOUBLE) )
! champ 2d
  call check( NF90_PUT_ATT  (ncid, tid, "_FillValue", NF90_FILL_REAL4) )
  !call check( NF90_PUT_ATT  (ncid, wmid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, vmid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, vortid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, urotid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, vrotid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, udivid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, vdivid, "_FillValue", NF90_FILL_DOUBLE) )

  !call check( NF90_PUT_ATT  (ncid, fuid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncid, fvid, "_FillValue", NF90_FILL_DOUBLE) )


  call check( nf90_enddef(ncid) )
end if
!######################################################################################################################
!######################################################################################################################
!###################################################################################################################
!=================================================================================================================== 
!						     START LOOP
!===================================================================================================================
!##################################################################################################################
!######################################################################################################################
!######################################################################################################################
do it = istep,mt,tstep !#----------------------- START TIME LOOP

if (NotFirst_reading) then
  nlat = nlat - 1
else
  print*,'First netcdf file reading...'
  NotFirst_reading = .true.
end if

allocate(u_DYN(nlon,nlat))
allocate(v_DYN(nlon,nlat))

print*, it
 call check( nf90_inq_varid(idfile,trim('u'),idu))
 call check( nf90_inq_varid(idfile,trim('v'),idv))

 call check( nf90_get_var(idfile,idu,u_DYN,(/1,1,1,it/),(/nlon,nlat,1,1/)) ) !(lon,lat,iz,it)
 call check( nf90_get_var(idfile,idv,v_DYN,(/1,1,1,it/),(/nlon,nlat,1,1/)) ) !(lon,lat,iz,it)

!print*, 'shape of u_DYN: ',shape(u_DYN)
!print*, 'shape of u_DYN: ',shape(v_DYN)
!print*, 'nlon = ',nlon
!print*, 'nlat = ',nlat
!print*, 'lat_DYN = ',lat_DYN
!print*, 'lon_DYN = ',lon_DYN
!#############################################################
!######################## BLOCK 0 ################################	REGULAR GRID
!##############################################################
! we initialy have a lat=[-89.75,89.75] et lon=[0.25,359.75]
! we shift to a grid lat=[-90,90] et lon=[0,359.5]
! It only a shift of 0.5* grid increment, so conserve South --> North
!*********************************************************vshifti function (initialisations for vshifte function)
lsav = 2*(2*nlat+nlon+16)
allocate(wsav(lsav)) ; wsav(:) = 0.0d0
ioff = 0
ierror=3
call vshifti(ioff,nlon,nlat,lsav,wsav,ierror)
select case (ierror)
  case(1) 
    print*,'ioff is not 0 or 1'
  case(2) 
    print*,'nlon < 4'
  case(3) 
    print*,'nlat < 3'
  case(4) 
    print*,'lsav < 2*(2*nlat+nlon+16)'
end select
!*********************************************************
! ------------WIND field
allocate(uoff(nlon,nlat)) ; uoff = u_DYN
allocate(voff(nlon,nlat)) ; voff = v_DYN
allocate(ureg(nlon,nlat+1)) ; ureg = 0.0d0
allocate(vreg(nlon,nlat+1)) ; vreg = 0.0d0
!-------------forcage
allocate(uoff_f(nlon,nlat)) ; uoff_f = fu
allocate(voff_f(nlon,nlat)) ; voff_f = fv
allocate(ureg_f(nlon,nlat+1)) ; ureg_f = 0.0d0
allocate(vreg_f(nlon,nlat+1)) ; vreg_f = 0.0d0
!-------------
ioff = 0
if (mod(nlon,2) == 0) then ! check if it is even
  lGwork = 2*nlon*(nlat+1)
else
  lGwork = nlon*(5*nlat+1)
end if
allocate(Gwork(lGwork)) ; Gwork(:)=0.0d0
ierror=3
call vshifte(ioff,nlon,nlat,uoff,voff,ureg,vreg,wsav,lsav,Gwork,lGwork,ierror)
call vshifte(ioff,nlon,nlat,uoff_f,voff_f,ureg_f,vreg_f,wsav,lsav,Gwork,lGwork,ierror)
select case (ierror)
  case(0) 
    print*,'no error'
  case(1) 
    print*,'ioff is not equal to 0 or 1'
  case(2) 
    print*,'nlon < 4'
  case(3) 
    print*,'nlat < 3'
  case(4) 
    print*,'lsave < 2*(nlon+2*nlat)+32'
  case(5) 
    print*,'lwork < 2*nlon*(nlat+1) for nlon even or lwork < nlon*(5*nlat+1) for nlon odd'
end select
nlat = nlat+1
print*,'nlat =',nlat-1,'---->',nlat
!#############################################################
!######################## BLOCK 0 ################################	GEO2MATH
!##############################################################
! 1) flip to North --> south (hemisphere ordered data)
! 2) from latitude to colatitued component, i.e. vg = -vm 
! 3) nlat become the first component and nlon the second.
! 4) The grid is now colat=[0,180] et lon=[0,359.5]
ig = 0 ! 0 if grid South->North & 1 if grid North->South : North (90) to south (-90)
!Southern->Northern hemisphere ordered data
allocate(ug(nlon,nlat)) ; ug = ureg
allocate(vg(nlon,nlat)) ; vg = vreg
!Northern->southern hemisphere ordered data
!------------WIND field
allocate(vm(nlat,nlon)) ; vm(:,:) = 0.0d0
allocate(wm(nlat,nlon)) ; wm(:,:) = 0.0d0
!------------forcage
allocate(fvm(nlat,nlon)) ; fvm(:,:) = 0.0d0
allocate(fwm(nlat,nlon)) ; fwm(:,:) = 0.0d0
!------------
allocate(gmwork(nlon*nlat)) ; gmwork(:) = 0.0d0
ierror=3
call geo2mathv(ig,nlon,nlat,ug,vg,vm,wm,gmwork)
call geo2mathv(ig,nlon,nlat,ureg_f,vreg_f,fvm,fwm,gmwork)
! #############################################     NOUVELLES COORDONNEES LAT-LON-COLAT
! vshifte et geo2mathv redefinissent la grille comme suit:
allocate(lon_SP(nlon))
allocate(lat_SP(nlat))
allocate(colat_SP(nlat))
do i=1,nlat
	lat_SP(i) = 90 - ((i-1)*PI/(nlat-1))*(180./PI)
	colat_SP(i) = ((i-1)*PI/(nlat-1))*(180./PI)
end do
do j=1,nlon
	lon_SP(j) = ((j-1)*2*PI/nlon)*(180./PI)
end do
!########################################################################################################################################
!========================================================================================================================================
!					   		  SPECTRAL ANALYSIS
!========================================================================================================================================
!########################################################################################################################################

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!										  Initialisation:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!**********Maximum value for m
if (mod(nlon,2) == 0) then
  mmax = min(nlat,nlon/2)
else
  mmax = min(nlat,(nlon+1)/2)
end if
!**********Maximum value for m in scalar (s) and vector (v) analysis
if (mod(nlon,2) == 0) then
  mdab_v = min(nlat,nlon/2)
else
  mdab_v = min(nlat,(nlon+1)/2)
end if
if (mod(nlon,2) == 0) then
  mdab_s = min(nlat,(nlon+2)/2)
else
  mdab_s = min(nlat,(nlon+1)/2)
end if
mdab_smaller = min(mdab_s,mdab_v)
!**********Maximum value for n
ndab=nlat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!										         Forcage:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allocate(fbr(mdab_v,ndab)) ; fbr(:,:) = 0.0d0
allocate(fbi(mdab_v,ndab)) ; fbi(:,:) = 0.0d0
allocate(fcr(mdab_v,ndab)) ; fcr(:,:) = 0.0d0
allocate(fci(mdab_v,ndab)) ; fci(:,:) = 0.0d0
 call VSA(nlon,nlat,nt,ityp,mdab_v,ndab,fvm,fwm,fbr,fbi,fcr,fci)
!
allocate(vortf_r(mdab_s,ndab)) ; vortf_r = 0.0d0
allocate(vortf_i(mdab_s,ndab)) ; vortf_i = 0.0d0
allocate(divf_r (mdab_s,ndab)) ; divf_r  = 0.0d0
allocate(divf_i (mdab_s,ndab)) ; divf_i  = 0.0d0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!										  Velocity field:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!******************************************************************
! Loop over total wavenumber
! (1) We omit n = 0 as factor = 0 for that wavenumber, and 1/factor
!     would break the streamfunction and velocity potential spectra,
!     which for n = 0 should be zero
! (2) We omit n = ndab-1 (i.e. the final total wavenumber filled is 
!     ndab-1) because this makes the result the same as if we ran 
!     vrtec on (uv_cr,uv_ci) then ran shaec to get the vorticity 
!     spectrum (equivalent for other quantities). This would be the 
!     proper (but longer) way to calculate these spectra. (NB this 
!     shortcut may not work if nlat no longer equals nlon/2+1)
!     If this is not done we must use vrtec rather than shsec in 
!     derived_fields.f90.
! (3) We omit n = ndab as ndab+1 is outside the array index range
!******************************************************************
! Spherical harmonic coeffs of vector fields
allocate(br(mdab_v,ndab)) ; br(:,:) = 0.0d0
allocate(bi(mdab_v,ndab)) ; bi(:,:) = 0.0d0
allocate(cr(mdab_v,ndab)) ; cr(:,:) = 0.0d0
allocate(ci(mdab_v,ndab)) ; ci(:,:) = 0.0d0
 call VSA(nlon,nlat,nt,ityp,mdab_v,ndab,vm,wm,br,bi,cr,ci)
! Spherical harmonic coeffs of scalar fields
allocate(vort_r(mdab_s,ndab)) ; vort_r = 0.0d0
allocate(vort_i(mdab_s,ndab)) ; vort_i = 0.0d0
allocate(div_r (mdab_s,ndab)) ; div_r  = 0.0d0
allocate(div_i (mdab_s,ndab)) ; div_i  = 0.0d0
allocate(sf_r  (mdab_s,ndab)) ; sf_r   = 0.0d0
allocate(sf_i  (mdab_s,ndab)) ; sf_i   = 0.0d0
allocate(vp_r  (mdab_s,ndab)) ; vp_r   = 0.0d0
allocate(vp_i  (mdab_s,ndab)) ; vp_i   = 0.0d0
!-----------------------------------------------------------------
do n = 1, ndab-2
  ! Scaling factor
  ! Need to divide through by radius as gradient operator is on unit sphere  
  factor = sqrt(1.0d0 * n * (n+1.0d0)) / radius
!*********************************************************************
!-----------------------  FORCAGE FIELD
  ! Vorticity spectrum: Spherepack 2.0 after Eq. 4.15
  vortf_r(1:mdab_smaller,n+1) =  fcr(1:mdab_smaller,n+1) * factor
  vortf_i(1:mdab_smaller,n+1) =  fci(1:mdab_smaller,n+1) * factor

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  divf_r (1:mdab_smaller,n+1) = -fbr(1:mdab_smaller,n+1) * factor
  divf_i (1:mdab_smaller,n+1) = -fbi(1:mdab_smaller,n+1) * factor
!*********************************************************************
!-----------------------  VELOCITY FIELD
  ! Vorticity spectrum: Spherepack 2.0 after Eq. 4.15
  vort_r(1:mdab_smaller,n+1) =  cr(1:mdab_smaller,n+1) * factor
  vort_i(1:mdab_smaller,n+1) =  ci(1:mdab_smaller,n+1) * factor

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  div_r (1:mdab_smaller,n+1) = -br(1:mdab_smaller,n+1) * factor
  div_i (1:mdab_smaller,n+1) = -bi(1:mdab_smaller,n+1) * factor

  ! Streamfunction spectrum: Spherepack 2.0 combination of Eqs. 4.1 
  ! for streamfunction, 4.6, and 4.8
  sf_r  (1:mdab_smaller,n+1) = -cr(1:mdab_smaller,n+1) / factor
  sf_i  (1:mdab_smaller,n+1) = -ci(1:mdab_smaller,n+1) / factor

  ! Velocity potential spectrum: Spherepack 2.0 combination of Eqs. 4.1 
  ! for velocity potential, 4.6, and 4.8
  vp_r  (1:mdab_smaller,n+1) =  br(1:mdab_smaller,n+1) / factor
  vp_i  (1:mdab_smaller,n+1) =  bi(1:mdab_smaller,n+1) / factor
!*********************************************************************
end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!									        Vorticity scalar:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!******************************************************************
! Synthesis of the vorticity coefficients to obtain the vorticity:
! vt(i,j) =  [-dv/dlambda + d(sint*w)/dtheta]/sint sur une sphere
! de rayon unite avec 1/raidus pour avoir la vorticite voir 
! coeffs vort_r et vort_i.
!******************************************************************
allocate(vort(nlat,nlon)) ; vort = 0.0d0
 call SHS(nlon,nlat,nt,isym,vort,mdab_s,ndab,vort_r,vort_i)
! rotational velocity synthesis
allocate(u_rot(nlat,nlon)) ; u_rot(:,:) = 0.0d0
allocate(v_rot(nlat,nlon)) ; v_rot(:,:) = 0.0d0
allocate(u_div(nlat,nlon)) ; u_div(:,:) = 0.0d0
allocate(v_div(nlat,nlon)) ; v_div(:,:) = 0.0d0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!							Rotational & divergent velocities vector:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !call ROTV(nlon,nlat,nt,isym,v_rot,u_rot,mdab_s,ndab,vort_r,vort_i)
 call DIVROTV(nlon,nlat,nt,isym,v_rot,u_rot,v_div,u_div,mdab_s,ndab,vort_r,vort_i,div_r,div_i)
v_rot = v_rot * radius
u_rot = u_rot * radius
v_div = v_div * radius
u_div = u_div * radius



!########################################################################################################################################
!========================================================================================================================================
!					     		SPECTRAL FLUXES --- div(u) = 0
!========================================================================================================================================
!########################################################################################################################################
! Initialization.
! Conversion factor from enstrophy to energy (BS83 Eq. 16)
allocate(j_i_conversion(ndab)) ; j_i_conversion = 0.0d0
do n = 1, ndab-1 
  j_i_conversion(n+1) = radius * radius / (n * (n + 1.0d0))
enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!							        Field: div . (vorticity * u_rot):
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==================================================================
! 1. Set up required data for this call
!==================================================================
allocate(term1_r(mdab_s,ndab))
allocate(term1_i(mdab_s,ndab))
allocate(term2_v(nlat,nlon))
allocate(term2_w(nlat,nlon))
allocate(term3  (nlat,nlon))
! Total transfers
term1_r = vort_r
term1_i = vort_i
term2_v = v_rot
term2_w = u_rot
term3   = vort
! Allocate and initialise the transfers and fluxes
allocate(j_2d   (mdab_s,ndab)) ; j_2d    = 0.0d0
allocate(i_2d   (mdab_s,ndab)) ; i_2d    = 0.0d0
allocate(j_2d_tm(mdab_s,ndab)) ; j_2d_tm = 0.0d0
allocate(i_2d_tm(mdab_s,ndab)) ; i_2d_tm = 0.0d0
allocate(j_1d   (ndab)) ; j_1d    = 0.0d0
allocate(i_1d   (ndab)) ; i_1d    = 0.0d0
allocate(j_1d_tm(ndab)) ; j_1d_tm = 0.0d0
allocate(i_1d_tm(ndab)) ; i_1d_tm = 0.0d0
allocate(f_1d   (ndab)) ; f_1d    = 0.0d0
allocate(h_1d   (ndab)) ; h_1d    = 0.0d0
allocate(f_1d_tm(ndab)) ; f_1d_tm = 0.0d0
allocate(h_1d_tm(ndab)) ; h_1d_tm = 0.0d0
!==================================================================
! 2. Calculate spherical harmonic decomposition of (v,w).grad vort 
!==================================================================
!******************************************************************
! From a standard vector calculus identity,
! (v,w).grad vort = grad.(vort(v,w)) - vort grad.(v,w) = grad.(vort(v,w))
! as grad.(v,w) = 0 for a divergence-free field
! Therefore we just need the spherical harmonic decomposition of 
! the divergence of vort(v,w) Allocate and initialise space.
!******************************************************************
allocate(vort_uv_br(mdab_v,ndab)) ; vort_uv_br = 0.0d0
allocate(vort_uv_bi(mdab_v,ndab)) ; vort_uv_bi = 0.0d0
allocate(vort_uv_cr(mdab_v,ndab)) ; vort_uv_cr = 0.0d0
allocate(vort_uv_ci(mdab_v,ndab)) ; vort_uv_ci = 0.0d0
! Field : vorticity * u_rot
allocate(vortv(nlat,nlon))
allocate(vortu(nlat,nlon))
vortv = term3*term2_v ! vort * v_rot
vortu = term3*term2_w ! vort * u_rot
! Vector Spectral Analysis
 call VSA(nlon,nlat,nt,ityp,mdab_v,ndab,vortv,vortu,vort_uv_br,vort_uv_bi,vort_uv_cr,vort_uv_ci)
!======================================================
! 3. Calculate divergence spectrum of vorticity * u_rot
!======================================================
! Allocate and initialise space
allocate(uv_dot_grad_vort_r(mdab_s,ndab)) ; uv_dot_grad_vort_r = 0.0d0
allocate(uv_dot_grad_vort_i(mdab_s,ndab)) ; uv_dot_grad_vort_i = 0.0d0
! factor to apply to divergent coefficients.
do n = 1, ndab-2

  ! Scaling factor
  ! Need to divide through by radius as gradient operator is on unit sphere  
  factor2 = sqrt(1.0d0 * n * (n+1.0d0)) / radius

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  uv_dot_grad_vort_r(1:mdab_smaller,n+1) = -vort_uv_br(1:mdab_smaller,n+1) * factor2
  uv_dot_grad_vort_i(1:mdab_smaller,n+1) = -vort_uv_bi(1:mdab_smaller,n+1) * factor2

enddo
!=====================================================
! 4. Calculate energy and enstrophy interaction terms
!=====================================================
! Allocate and initialise the transfers and fluxes
!
do n = 1, ndab-1
  j_2d(:,n+1) = -0.25d0 * 2.0d0 * (term1_r(:,n+1) * uv_dot_grad_vort_r(:,n+1) + term1_i(:,n+1) * uv_dot_grad_vort_i(:,n+1))
  i_2d(:,n+1) = j_2d(:,n+1) * j_i_conversion(n+1)
enddo
write(6,'("Day ",i3," net enstrophy transfer ",es20.13)'), 1, sum(j_2d(1,:)) + 2.0d0 * sum(j_2d(2:mdab_s,:))
write(6,'("Day ",i3," net energy    transfer ",es20.13)'), 1, sum(i_2d(1,:)) + 2.0d0 * sum(i_2d(2:mdab_s,:))

! n = 0 term
j_1d(1) = j_2d(1,1)
! n > 0 terms
do n = 1, ndab-1
  j_1d(n+1) = j_2d(1,n+1) + 2.0d0 * sum(j_2d(2:n+1,n+1))
enddo
! Convert to energy interaction term
i_1d(:) = j_1d(:) * j_i_conversion 

!====================================================================
! 5. Calculate energy and enstrophy spectral fluxes (BS13 Eqs 11-12)
!====================================================================
do n = 1, ndab-1 ! n > 0 terms (n = 0 term is zero)
  f_1d(n+1) = -sum(i_1d(1:n)) ! Energy flux
  h_1d(n+1) = -sum(j_1d(1:n)) ! Enstrophy flux
enddo
!f_1d    = 0.0d0
!h_1d    = 0.0d0
! call FLUXES(mdab_s,ndab,radius,vort_r,vort_i,uv_dot_grad_vort_r,uv_dot_grad_vort_i,f_1d,h_1d)


!########################################################################################################################################
!========================================================================================================================================
!					     		SPECTRAL FLUXES --- div(u) =/= 0
!========================================================================================================================================
!########################################################################################################################################



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!					             Field: - rot /\ (u /\ vorticity)  --> rotuw:
!   							    - div .  (u /\ vorticity)  --> divuw:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!******************************************************************
! Ici on calcul deux terme du produit scalaire (u,u.gradu) ou 
! u est la vitesse totale. d ailleurs je doit pouvoir suprimer
! le double appel a vhaec plus bas necessaire lorsque je m etais 
! d un cote la vitesse non div et de l autre la vitesse div. 
! c est la vitesse totale partout et donc vhaec pemet de trouver
! la divergence et le rotationel du terme u /\ vorticity.
!******************************************************************
!==================================================================
! 1. Set up required data for this call
!==================================================================
!--------------------------------------------->  - (u /\ vorticity)
allocate(UvectVT_lon(nlat,nlon)) ; UvectVT_lon(:,:) = 0.0d0
allocate(UvectVT_lat(nlat,nlon)) ; UvectVT_lat(:,:) = 0.0d0
do j=1,nlon
   do i=1,nlat
	UvectVT_lat(i,j) = -(wm(i,j)*vort(i,j))
	UvectVT_lon(i,j) = -(-vm(i,j)*vort(i,j))
   end do
end do
!do j=1,nlon
!   do i=1,nlat
!	UvectVT_lat(i,j) = -(u_div(i,j)*vort(i,j))
!	UvectVT_lon(i,j) = -(-v_div(i,j)*vort(i,j))
!   end do
!end do
!==================================================================
! 2. Calculate spherical harmonic decomposition of (v,w).grad vort 
!==================================================================
!----------->  - rot /\ (u /\ vorticity)  --> rotuw
!      	       - div .  (u /\ vorticity)  --> divuw
allocate(Fdivuwr(mdab_v,ndab)) ; Fdivuwr(:,:) = 0.0d0
allocate(Fdivuwi(mdab_v,ndab)) ; Fdivuwi(:,:) = 0.0d0
allocate(Frotuwr(mdab_v,ndab)) ; Frotuwr(:,:) = 0.0d0
allocate(Frotuwi(mdab_v,ndab)) ; Frotuwi(:,:) = 0.0d0
 call VSA(nlon,nlat,nt,ityp,mdab_v,ndab,UvectVT_lat,UvectVT_lon,Fdivuwr,Fdivuwi,Frotuwr,Frotuwi)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!					             		        Field: 1/2 Nabla ||u||^2:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!==================================================================
! 1. Champ scalaire u.u = ||u||^2
!==================================================================
allocate(udotu(nlat,nlon))
udotu(:,:) = 0.
do j=1,nlon
   do i=1,nlat
	udotu(i,j) = SQRT(wm(i,j)*wm(i,j) + vm(i,j)*vm(i,j))**2.
   end do
end do
!do j=1,nlon
!   do i=1,nlat
!	udotu(i,j) = SQRT(u_div(i,j)*u_div(i,j) + v_div(i,j)*v_div(i,j))**2.
!   end do
!end do
!==================================================================
! 2. Gradient of u.u
!==================================================================
!****************************************************************
! v is the colatitudinal and w is the east longitudinal component
! of the gradient: 
! v(i,j) = d(sf(i,j))/dtheta and 
! w(i,j) = 1/sint*d(sf(i,j))/dlambda
!****************************************************************
!psectral scalar analysis
allocate(uur(mdab_s,ndab)) ; uur(:,:) = 0.0d0
allocate(uui(mdab_s,ndab)) ; uui(:,:) = 0.0d0
 call SSA(nlon,nlat,nt,ityp,mdab_s,ndab,udotu,uur,uui)
! Grad
allocate(vgrad(nlat,nlon)) ; vgrad(:,:) = 0.0d0
allocate(wgrad(nlat,nlon)) ; wgrad(:,:) = 0.0d0
 call GRAD(nlon,nlat,nt,isym,vgrad,wgrad,uur,uui)
allocate(vgradN(nlat,nlon)) ; vgradN(:,:) = 0.0d0
allocate(wgradN(nlat,nlon)) ; wgradN(:,:) = 0.0d0
! Dimensionalisation
do j=1,nlon
   do i=1,nlat
	vgradN(i,j) = vgrad(i,j)/radius
	wgradN(i,j) = wgrad(i,j)/radius
   end do
end do
!==================================================================
! 2. div . grad ||u||^2
!==================================================================
allocate(divgraduur(mdab_v,ndab)) ; divgraduur(:,:) = 0.0d0 ! divergente
allocate(divgraduui(mdab_v,ndab)) ; divgraduui(:,:) = 0.0d0 ! divergente
allocate(rotgraduur(mdab_v,ndab)) ; rotgraduur(:,:) = 0.0d0 ! rotationnelle
allocate(rotgraduui(mdab_v,ndab)) ; rotgraduui(:,:) = 0.0d0 ! rotationnelle
 call VSA(nlon,nlat,nt,ityp,mdab_v,ndab,vgradN,wgradN,divgraduur,divgraduui,rotgraduur,rotgraduui)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!					             		              Spectra and Fluxes:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Allocate and initialise space
allocate(div1_ugradu_r(mdab_s,ndab)) ; div1_ugradu_r = 0.0d0
allocate(div1_ugradu_i(mdab_s,ndab)) ; div1_ugradu_i = 0.0d0
allocate(rot_ugradu_r(mdab_s,ndab)) ; rot_ugradu_r = 0.0d0
allocate(rot_ugradu_i(mdab_s,ndab)) ; rot_ugradu_i = 0.0d0
allocate(div2_ugradu_r(mdab_s,ndab)) ; div2_ugradu_r = 0.0d0
allocate(div2_ugradu_i(mdab_s,ndab)) ; div2_ugradu_i = 0.0d0
! factor to apply to divergent coefficients.
do n = 1, ndab-2

  ! Scaling factor
  ! Need to divide through by radius as gradient operator is on unit sphere  
  factor2 = sqrt(1.0d0 * n * (n+1.0d0)) / radius

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  div1_ugradu_r(1:mdab_smaller,n+1) = -Fdivuwr(1:mdab_smaller,n+1) * factor2
  div1_ugradu_i(1:mdab_smaller,n+1) = -Fdivuwi(1:mdab_smaller,n+1) * factor2

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  div2_ugradu_r(1:mdab_smaller,n+1) = -0.5*divgraduur(1:mdab_smaller,n+1) * factor2
  div2_ugradu_i(1:mdab_smaller,n+1) = -0.5*divgraduui(1:mdab_smaller,n+1) * factor2

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  rot_ugradu_r(1:mdab_smaller,n+1) = Frotuwr(1:mdab_smaller,n+1) * factor2
  rot_ugradu_i(1:mdab_smaller,n+1) = Frotuwi(1:mdab_smaller,n+1) * factor2

enddo

allocate(fn_rot(ndab)) ; fn_rot = 0.0d0
allocate(hn_rot(ndab)) ; hn_rot = 0.0d0
allocate(fn_div1(ndab)) ; fn_div1 = 0.0d0
allocate(hn_div1(ndab)) ; hn_div1 = 0.0d0
allocate(fn_div2(ndab)) ; fn_div2 = 0.0d0
allocate(hn_div2(ndab)) ; hn_div2 = 0.0d0

! rot(u) * - rot /\ (u /\ vorticity)
 call FLUXES(mdab_s,ndab,radius,vort_r,vort_i,rot_ugradu_r,rot_ugradu_i,fn_rot,hn_rot)
! div(u) * - div .  (u /\ vorticity)
 call FLUXES(mdab_s,ndab,radius,div_r,div_i,div1_ugradu_r,div1_ugradu_i,fn_div1,hn_div1)
! div(u) * 1/2 Laplacien ||u||^2
 call FLUXES(mdab_s,ndab,radius,div_r,div_i,div2_ugradu_r,div2_ugradu_i,fn_div2,hn_div2)



! -----------------------------------VORTICITY
!allocate(vt(nlat,nlon)) ; vt(:,:) = 0.0d0
! call curl(nlon,nlat,nt,isym,vt,cr,ci)
!######################################################################################################################################
!######################################################################################################################################
!######################################################################################################################################
!######################################################################################################################################
!######################################################################################################################################
!#####################################################################################################################################
!###################################################################################################################################
!##################################### SAVE OUT_PUT VARIABLES ####################################################################
!##############################################################################################################################
!##################################################################################################
! -- SAVE scalar & vector fields
if(NetcdefStatData == 1) then
  if (First_reading) then
    call check( NF90_PUT_VAR (ncid, lonid, lon_SP) )
    call check( NF90_PUT_VAR (ncid, latid, lat_SP) )
    call check( NF90_PUT_VAR (ncid, colatid, colat_SP) )
  else
    First_reading = .false.
  end if
!champ 2d
   call check( nf90_put_var( ncid, tid, it, start=(/it/) ) )
   !call check( nf90_put_var( ncid, wmid, transpose(wm),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncid, vmid, transpose(vm),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncid, vortid, transpose(vort),(/1,1,it/),(/nlon,nlat,1/) ) )

   !call check( nf90_put_var( ncid, urotid, transpose(u_rot),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncid, vrotid, transpose(v_rot),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncid, udivid, transpose(u_div),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncid, vdivid, transpose(v_div),(/1,1,it/),(/nlon,nlat,1/) ) )
!champ 1d
   call check( nf90_put_var( ncid, fnid, f_1d,(/1,it/),(/nlat,1/) ) )
   call check( nf90_put_var( ncid, hnid, h_1d,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncid, fnrotid, fn_rot,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncid, hnrotid, hn_rot,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncid, fndiv1id, fn_div1,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncid, hndiv1id, hn_div1,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncid, fndiv2id, fn_div2,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncid, hndiv2id, hn_div2,(/1,it/),(/nlat,1/) ) )

   !call check( nf90_put_var( ncid, fuid, transpose(fwm),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncid, fvid, transpose(fvm),(/1,1,it/),(/nlon,nlat,1/) ) )
end if
!##################################################################################
! -- SAVE coefficients
file_spectra = 'xspec_rotu_'!------------------------------------------- rot /\ U
nickname = 'vortic'
 call savecoeffs( it,file_netcdf,file_spectra,nickname,mdab_s,nlat,vort_r,vort_i )
!
file_spectra = 'xspec_divu_'!------------------------------------------- div . U
nickname = 'diverg'
 call savecoeffs( it,file_netcdf,file_spectra,nickname,mdab_s,nlat,div_r,div_i )
!
file_spectra = 'xspec_rotf_'!------------------------------------------- rot /\ forcage
nickname = 'rotatf'
 call savecoeffs( it,file_netcdf,file_spectra,nickname,mdab_s,nlat,vortf_r,vortf_i )
!
file_spectra = 'xspec_divf_'!------------------------------------------- div . forcage
nickname = 'diverf'
 call savecoeffs( it,file_netcdf,file_spectra,nickname,mdab_s,nlat,divf_r,divf_i )
!
!file_spectra = 'xspec_divuvorticity_'!---------------------------------- div . (w * u_rot)
!nickname = 'divuvt'
! call savecoeffs( it,file_netcdf,file_spectra,nickname,mmax,nlat,vort_uv_br,vort_uv_bi )
!
!file_spectra = 'xspec_rotuw_'!------------------------------------------- -rot /\ (U /\ w)
!nickname = 'rotuvt'
! call savecoeffs( it,file_netcdf,file_spectra,nickname,mmax,nlat,Frotuwr,Frotuwi )
!
!file_spectra = 'xspec_divuw_'!------------------------------------------- -div . (U /\ w)
!nickname = 'rotuvt'
! call savecoeffs( it,file_netcdf,file_spectra,nickname,mmax,nlat,Fdivuwr,Fdivuwi )
!
!file_spectra = 'xspec_udotu_'!------------------------------------------- || u . u ||
!nickname = 'udotu '
! call savecoeffs( it,file_netcdf,file_spectra,nickname,mmax,nlat,divgraduur,divgraduui )
!
!if (mod(nlon,2) == 0) then ! ---- Il se trouve que la decomposition champ scalaire avec shaec mdab/=mmax 
!  mdab=min(nlat,(nlon+2)/2)
!else
!  mdab=min(nlat,(nlon+1)/2)
!end if
!file_spectra = 'xspec_udotu_'!------------------------------------------- || u . u ||
!nickname = 'udotu '
! call savecoeffs( it,file_netcdf,file_spectra,nickname,mdab,nlat,uur,uui )


!file_spectra = 'xspec_W_'!------------------------------------------- || u . u ||
!nickname = 'udotu '
! call savecoeffs( it,file_netcdf,file_spectra,nickname,mdab_s,nlat,wr,wi )
!
!file_spectra = 'xspec_vort_'!------------------------------------------- 
!nickname = 'rotuvt'
! call savecoeffs( it,file_netcdf,file_spectra,nickname,mdab_s,nlat,vort_r,vort_i )

!##########################################################
!##################################### DEALLOCATE ############
!##########################################################
deallocate(u_DYN)
deallocate(v_DYN)
deallocate(wsav)
deallocate(uoff)
deallocate(voff)
deallocate(ureg)
deallocate(vreg)
deallocate(uoff_f)
deallocate(voff_f)
deallocate(ureg_f)
deallocate(vreg_f)
deallocate(Gwork)
deallocate(ug)
deallocate(vg)
deallocate(vm)
deallocate(wm)
deallocate(fvm)
deallocate(fwm)
deallocate(gmwork)
deallocate(br)
deallocate(bi)
deallocate(cr)
deallocate(ci)
deallocate(fbr)
deallocate(fbi)
deallocate(fcr)
deallocate(fci)
deallocate(uur)
deallocate(uui)
deallocate(v_rot)
deallocate(u_rot)
deallocate(lon_SP)
deallocate(lat_SP)
deallocate(colat_SP)
deallocate(UvectVT_lon)
deallocate(UvectVT_lat)
deallocate(Fdivuwr)
deallocate(Fdivuwi)
deallocate(Frotuwr)
deallocate(Frotuwi)
deallocate(udotu)
deallocate(vgrad)
deallocate(wgrad)
deallocate(vgradN)
deallocate(wgradN)
deallocate(divgraduur)
deallocate(divgraduui)
deallocate(rotgraduur) 
deallocate(rotgraduui)
deallocate(vort_r)
deallocate(vort_i)
deallocate(div_r)
deallocate(div_i) 
deallocate(vortf_r)
deallocate(vortf_i)
deallocate(divf_r)
deallocate(divf_i) 
deallocate(sf_r)  
deallocate(sf_i)  
deallocate(vp_r) 
deallocate(vp_i)
deallocate(vort)
deallocate(u_div)
deallocate(v_div)
deallocate(term1_r)
deallocate(term1_i)
deallocate(term2_v)
deallocate(term2_w)
deallocate(term3)
deallocate(j_1d)
deallocate(i_1d)
deallocate(j_1d_tm)
deallocate(i_1d_tm)
deallocate(f_1d)
deallocate(h_1d)
deallocate(f_1d_tm)
deallocate(h_1d_tm)
deallocate(j_2d)
deallocate(i_2d)
deallocate(j_2d_tm)
deallocate(i_2d_tm)
deallocate(vort_uv_br)
deallocate(vort_uv_bi)
deallocate(vort_uv_cr)
deallocate(vort_uv_ci)
deallocate(j_i_conversion)
deallocate(vortv)
deallocate(vortu)
deallocate(uv_dot_grad_vort_r)
deallocate(uv_dot_grad_vort_i)
deallocate(rot_ugradu_i)
deallocate(rot_ugradu_r)
deallocate(div1_ugradu_i)
deallocate(div1_ugradu_r)
deallocate(div2_ugradu_r)
deallocate(div2_ugradu_i)
deallocate(fn_rot)
deallocate(hn_rot)
deallocate(fn_div1)
deallocate(hn_div1)
deallocate(fn_div2)
deallocate(hn_div2)


end do
if(NetcdefStatData == 1) then
  call check( nf90_close(ncid) )
end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			  SUBROUTINES 		      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES I
!#							     check
!###################################################################################################################
!###################################################################################################################
contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check 
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES II 
!#							   savecoeffs
!###################################################################################################################
!###################################################################################################################
! ---------- SAVE DATA
subroutine savecoeffs(it,file_netcdf,file_spectra,nickname,mmax,nlat,coeffr,coeffi)
    implicit none
    integer :: it,pos1,pos2,pos3,mt1,mmax,nlat,i,j
    character (len=100) :: file_netcdf,suffix,file_name,file_spectra,nickname
    double precision, dimension(:,:) ::  coeffr,coeffi

	!###########################################################################################################
	!######################################  WRITE SPHERICAL HARMONIC COEFFICIENTS  #################################
	!****************************************************************************************************************
	!****************************************************************************************************************
	! 						Start Writing file
	pos1=scan(file_netcdf, "-", .FALSE.) ! ici on cherche la position des '-' qui encadre le premier chiffre=premiere iteration.
	pos2=scan(file_netcdf(pos1+1:), "-", .FALSE.)
	Read( file_netcdf(pos1+1:pos1+pos2-1), * )mt1 ! puis on extrait du nom du fichier en entree la valeur de la premiere iteration temp
	write(unit=suffix,fmt=*)mt1+(it-1)
	file_name = trim(file_spectra)//trim(adjustl(suffix))
	!**********Writing header of file_name
	open(unit=10,file=file_name,action="write",position="rewind")
	write(10,'(a10,a2)',advance='no') '#Spherical',' '
	write(10,'(a50)',advance='no') file_netcdf
	write(10,*)
	write(10,'(a11,a1)',advance='no') '#Coeffs-max',' '
	write(10,'(a4,i5,a5,i5,a31)',advance='no') 'nbM=',mmax,' nbN=',nlat,' '
	write(10,*)
	write(10,'(a1,a1)',advance='no') '#',' '
	write(10,'(a6,a54)',advance='no') nickname,' '
	!**********Writing coefficients
	do i=1,mmax
	write(10,*)
	write(10,'(a6,i3,a53)',advance='no') '# m = ',i-1,' '
	write(10,*)
	do j=1,nlat
	  write(10,'(i5,a6)',advance='no') j-1,' ' !cette ligne permet de mettre la colonne des indices n
	  write(10,'(e13.6E2,a10)',advance='no') coeffr(i,j),' '
	  write(10,'(e13.6E2,a10)',advance='no') coeffi(i,j),' '
	  if (j /= nlat) write(10,*)
	end do
	end do
	!**********END Writing
	close(10) 

	print *, "***SUCCESS writing ",trim(file_name)
	! 						End Writing file
	!****************************************************************************************************************
	!****************************************************************************************************************
end subroutine savecoeffs
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES III 
!#							       CURL
!###################################################################################################################
!###################################################################################################################
! ---------- curl
subroutine curl(nlon,nlat,nt,isym,rot,cr,ci)
    implicit none
    integer :: nlon,nlat,ll1,l1,l2,lshsec,lldwork,ierror,nt,isym,ivrt,jvrt,mdc,ndc,vlwork
    double precision, dimension (:), allocatable :: wshsec,ddwork,vwork
    double precision, dimension(:,:), allocatable :: rot,cr,ci
!************************************************************
! The vorticity field is obtained from cr and ci coefficients 
! of the rotational part of the wind field. It is important
! to note that the vorticity scalar field is computed as follow:
! vt(i,j) =  [-dv/dlambda + d(sint*w)/dtheta]/sint il manque donc
! 1/R pour avoir la vorticite de l ecoulement.
!*********************************************************shseci function (initialisations for Vhaec function)
if (mod(nlon,2) == 0) then ! check if it is even
  ll1 = min(nlat,(nlon+2)/2)
else
  ll1 = min(nlat,(nlon+1)/2)
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
lshsec=2*nlat*l2+3*((ll1-2)*(nlat+nlat-ll1-1))/2+nlon+15
allocate(wshsec(lshsec))
wshsec(:) = 0.
lldwork=nlat+1
allocate(ddwork(lldwork))
ddwork(:) = 0.
ierror=3
call shseci(nlat,nlon,wshsec,lshsec,ddwork,lldwork,ierror)
select case (ierror)
  case(0) 
    print*,'shseci: No ERROR in subroutine'
  case(1) 
    print*,'shseci: ERROR on nlat'
  case(2) 
    print*,'shseci: ERROR on nlong'
  case(3) 
    print*,'shseci: ERROR on lvhaec'
  case(4) 
    print*,'shseci: ERROR on ldwork'
end select
!**********Vrtec function (computes the vorticity)	
ivrt = nlat
jvrt = nlon
if (mod(nlon,2) == 0) then
  mdc=min(nlat,nlon/2)
else
  mdc=min(nlat,(nlon+1)/2)
end if
ndc = nlat
if (mod(nlon,2) == 0) then
  !l1 = min(nlat,nlon/2) 
  l1 = min(nlat,(nlon+2)/2) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! attention pas bon ici
else
  l1 = min(nlat,(nlon+1)/2) 
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
!
if (isym == 0) then
  vlwork = nlat*(nt*nlon+max(3*l2,nlon)+2*nt*l1+1)
else
  vlwork = l2*(nt*nlon+max(3*nlat,nlon)) + nlat*(2*nt*l1+1)
end if
allocate(vwork(vlwork))
vwork(:) = 0.
ierror = 3
call vrtec(nlat,nlon,isym,nt,rot,ivrt,jvrt,cr,ci,mdc,ndc,wshsec,lshsec,vwork,vlwork,ierror)
!**********End Vrtec
select case (ierror)
  case(0) 
    print*,'vrtec: No error'
  case(1) 
    print*,'vrtec: ERROR on nlat'
  case(2) 
    print*,'vrtec: ERROR on nlong'
  case(3) 
    print*,'vrtec: ERROR on ityp'
  case(4) 
    print*,'vrtec: ERROR on nt'
  case(5) 
    print*,'vrtec: ERROR on idvw'
  case(6) 
    print*,'vrtec: ERROR on jdvw'
  case(7) 
    print*,'vrtec: ERROR on mdab'
  case(8) 
    print*,'vrtec: ERROR on ndab'
  case(9) 
    print*,'vrtec: ERROR on lvhaec'
  case(10) 
    print*,'vrtec: ERROR on lwork'
end select
end subroutine curl
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES IV
!#							       VSA:
!#						    Vector Spectral Analysis
!###################################################################################################################
!###################################################################################################################
! ---------- VSA
subroutine VSA(nlon,nlat,nt,ityp,mdab,ndab,vm,wm,br,bi,cr,ci)
    implicit none
    integer :: nlon,nlat,l1,l2,lvhaec,ldwork,ierror,nt,ityp,idvw,jdvw,mdab,ndab,lwork
    double precision, dimension (:), allocatable :: wvhaec,dwork,work
    double precision, dimension(:,:), allocatable :: br,bi,cr,ci,vm,wm
!********************************************************************
! Here, we perform vector analysis of zonal_wind and merid_wind
! using vhaec and leading to gauss coefficients of spherical harmonic
! functions cr,ci,br,bi. Respectively of the rotational and divergent
!*********************************************************Vhaeci function (initialisations for Vhaec function)
if (mod(nlon,2) == 0) then ! check if it is even
  l1 = min(nlat,nlon/2)
else
  l1 = min(nlat,(nlon+1)/2) 
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
lvhaec=4*nlat*l2+3*max(l1-2,0)*(nlat+nlat-l1-1)+nlon+15
allocate(wvhaec(lvhaec))
wvhaec(:) = 0.
ldwork=  2*(nlat+2)
allocate(dwork(ldwork))
dwork(:) = 0.
ierror=3
call vhaeci(nlat,nlon,wvhaec,lvhaec,dwork,ldwork,ierror)
!********************************************************Vhaec function (calcul spectral coefficients cr, ci, br, bi. 
!**********Vhaec function (computes spectra coefficients)
idvw=nlat
jdvw=nlon
if (mod(nlon,2) == 0) then
  l1 = min(nlat,nlon/2)
else
  l1 = min(nlat,(nlon+1)/2) 
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
lvhaec=4*nlat*l2+3*max(l1-2,0)*(nlat+nlat-l1-1)+nlon+15
if (ityp .le. 2) then
  lwork=nlat*(2*nt*nlon+max(6*l2,nlon))
else
  lwork=l2*(2*nt*nlon+max(6*nlat,nlon))
end if
allocate(work(lwork)) ; work(:) = 0.0d0
ierror=3
call vhaec(nlat,nlon,ityp,nt,vm,wm,idvw,jdvw,br,bi,cr,ci,mdab,ndab,wvhaec,lvhaec,work,lwork,ierror)
select case (ierror)
  case(0) 
    print*,'no error in subroutine VSA'
end select
deallocate(work)
end subroutine VSA
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES V
!#							       SSA:
!#						    Scalar Spectral Analysis
!###################################################################################################################
!###################################################################################################################
! ---------- SSA
subroutine SSA(nlon,nlat,nt,ityp,mdab,ndab,sf,ar,ai)
    implicit none
    integer :: nlon,nlat,ll1,l2,nt,ierror,ityp,idvw,jdvw,mdab,ndab,lshaec,sldwork,sslwork
    double precision, dimension (:), allocatable :: wshaec,sdwork,sswork
    double precision, dimension(:,:), allocatable :: ar,ai,sf
!**********Shaeci function (initialisations for Shaec function)
if (mod(nlon,2) == 0) then
  ll1 = min(nlat,(nlon+2)/2)
else
  ll1 = min(nlat,(nlon+1)/2)
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
lshaec=2*nlat*l2+3*((ll1-2)*(nlat+nlat-ll1-1))/2+nlon+15
allocate(wshaec(lshaec))
wshaec(:) = 0.
sldwork=nlat+1
allocate(sdwork(sldwork))
sdwork(:) = 0.
ierror=3
call shaeci(nlat,nlon,wshaec,lshaec,sdwork,sldwork,ierror)
!**********Shaeci function result
select case (ierror)
 ! case(0) 
 !   print*,'Shaeci: no error in subroutine SSA'
  case(1) 
    print*,'Shaeci: ERROR on nlat'
  case(2) 
    print*,'Shaeci: ERROR on nlong'
  case(3) 
    print*,'Shaeci: ERROR on lshaec'
  case(4) 
    print*,'Shaeci: ERROR on sldwork'
end select
!********************************************************
!**********Shaec function (computes spectra coefficients)
if (ityp == 0) then
idvw=nlat
else
print*,'Check condition in doc'
end if
jdvw=nlon
if (ityp == 0) then
  sslwork=nlat*(nt*nlon+max(3*l2,nlon))
else
  sslwork=l2*(nt*nlon+max(3*nlat,nlon))
end if
allocate(sswork(sslwork))
sswork(:) = 0.
ierror=3
call shaec(nlat,nlon,ityp,nt,sf,idvw,jdvw,ar,ai,mdab,ndab,wshaec,lshaec,sswork,sslwork,ierror)
select case (ierror)
  case(0) 
    print*,'shaec: no error in subroutine SSA'
end select
end subroutine SSA
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES VI
!#							       GRAD:
!#						    	     GRADient
!###################################################################################################################
!###################################################################################################################
! ---------- SSA
subroutine GRAD(nlon,nlat,nt,isym,vgrad,wgrad,ar,ai)
    implicit none
    integer :: nlon,nlat,l1,l2,nt,ierror,isym,idvw,jdvw,mdab,ndab,glvhsec,glwwork,gllwwork
    double precision, dimension (:), allocatable :: gwvhsec,gwwork,gdwwork
    double precision, dimension(:,:), allocatable :: ar,ai,vgrad,wgrad

!********** function vhseci (initialisations for Shaec function)
if (mod(nlon,2) == 0) then
  l1 = min(nlat,nlon/2)!+1)
else
  l1 = min(nlat,(nlon+1)/2) 
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
glvhsec = 4*nlat*l2+3*max(l1-2,0)*(nlat+nlat-l1-1)+nlon+15
allocate(gwvhsec(glvhsec))
gwvhsec(:) = 0.
glwwork = 2*(nlat+2)!4*(nlat+1)
allocate(gwwork(glwwork))
gwwork(:) = 0.
ierror=3
call vhseci(nlat,nlon,gwvhsec,glvhsec,gwwork,glwwork,ierror)
select case (ierror)
  case(0) 
    print*,'vhseci : no error subroutine GRAD'
end select
!********************************************************
!********** Gradec function (computes spectra coefficients)
if (isym == 0) then
idvw=nlat
else
print*,'Check condition in doc'
end if
jdvw=nlon
if (mod(nlon,2) == 0) then
  mdab=min(nlat,(nlon+2)/2)
else
  mdab=min(nlat,(nlon+1)/2)
end if
ndab=nlat
if (isym == 0) then
  gllwwork = nlat*(2*nt*nlon+max(6*l2,nlon)) + nlat*(2*l1*nt+1) 
else
  gllwwork = l2*(2*nt*nlon+max0(6*nlat,nlon)) + nlat*(2*l1*nt+1)
end if
allocate(gdwwork(gllwwork))
gdwwork(:) = 0.
ierror=3
call gradec(nlat,nlon,isym,nt,vgrad,wgrad,idvw,jdvw,ar,ai,mdab,ndab,gwvhsec,glvhsec,gdwwork,gllwwork,ierror)
select case (ierror)
  case(0) 
    print*,'gradec: no error in subroutine GRAD'
end select
end subroutine GRAD
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES VII
!#							       SHS:
!#						   Scalar Harmonic Synthesis
!###################################################################################################################
!###################################################################################################################
! ---------- SHS
subroutine SHS(nlon,nlat,nt,isym,g,mdab_s,ndab,aar,aai)
    implicit none
    integer :: nlon,nlat,idg,jdg,ll1,l2,nt,isym,slshsec,slldwork,ierror,mdab_s,ndab,slwork	
    double precision, dimension (:), allocatable :: swshsec,sddwork,svwork
    double precision, dimension(:,:), allocatable :: aar,aai,g

!*********************************************************shseci function (initialisations for Vhaec function)
if (mod(nlon,2) == 0) then ! check if it is even
  ll1 = min(nlat,(nlon+2)/2)
else
  ll1 = min(nlat,(nlon+1)/2)
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
slshsec=2*nlat*l2+3*((ll1-2)*(nlat+nlat-ll1-1))/2+nlon+15
allocate(swshsec(slshsec)) ; swshsec(:) = 0.0d0
slldwork=nlat+1
allocate(sddwork(slldwork)) ; sddwork(:) = 0.0d0
ierror=3
call shseci(nlat,nlon,swshsec,slshsec,sddwork,slldwork,ierror)
deallocate(sddwork)
select case (ierror)
  case(0) 
    print*,'shseci: No ERROR in subroutine'
end select
!********************************************************
!********** Shsec function (computes spectra coefficients)
if (isym == 0) then
idg=nlat
else
print*,'Check condition in doc'
end if
jdg = nlon
if (isym == 0) then
  slwork = nlat*(nt*nlon+max(3*l2,nlon))
else
  slwork = l2*(nt*nlon+max(3*nlat,nlon))
end if
allocate(svwork(slwork)) ; svwork(:) = 0.0d0
ierror = 3
 call shsec(nlat,nlon,isym,nt,g,idg,jdg,aar,aai,mdab_s,ndab,swshsec,slshsec,svwork,slwork,ierror)
select case (ierror)
  case(0) 
    print*,'shsec: no error in subroutine SHS'
end select
deallocate(svwork)
end subroutine SHS
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES VIII
!#							       DIVROTV:
!#						   DIVergent ROTational Velocity
!###################################################################################################################
!###################################################################################################################
! ---------- DIVROTV
subroutine DIVROTV(nlon,nlat,nt,isym,v_rot,u_rot,v_div,u_div,mdab_s,ndab,wr,wi,dr,di)!(nlon,nlat,nt,isym,v_rot,u_rot,mdab_s,ndab,wr,wi)
    implicit none
    integer :: nlon,nlat,idvw,jdvw,l1,l2,nt,isym,slshsec,mdab_s,ndab,lvhsec,lwwork,ierror,llwwork
    double precision, dimension (:), allocatable :: wvhsec,wwork,dwwork
    double precision, dimension(:,:), allocatable :: wr,wi,dr,di,u_rot,v_rot,u_div,v_div,pertrb_rot,pertrb_div

!********** INITIALISATION
if (mod(nlon,2) == 0) then
  l1 = min(nlat,nlon/2)!+1)
else
  l1 = min(nlat,(nlon+1)/2) 
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
lvhsec = 4*nlat*l2+3*max(l1-2,0)*(nlat+nlat-l1-1)+nlon+15
allocate(wvhsec(lvhsec)) ; wvhsec(:) = 0.0d0
lwwork = 2*(nlat+2)
allocate(wwork(lwwork)) ; wwork(:) = 0.0d0
ierror=3
call vhseci(nlat,nlon,wvhsec,lvhsec,wwork,lwwork,ierror)
select case (ierror)
  case(0) 
    print*,'vhseci : no error subroutine'
end select
!********************************************************
!**********Ivrtec
idvw=nlat
jdvw=nlon
allocate(pertrb_rot(nlat,nlon)) ; pertrb_rot(:,:) = 0.0d0
if (mod(nlon,2) == 0) then
  l1 = min(nlat,nlon/2)
  !l1 = min(nlat,(nlon+1)/2)
else
  l1 = min(nlat,(nlon+1)/2) 
end if
if (mod(nlat,2) == 0) then
  l2 = nlat/2
else
  l2 = (nlat+1)/2
end if
if (isym == 0) then
  llwwork = nlat*(2*nt*nlon+max(6*l2,nlon) + 2*nt*l1 + 1)
else
  llwwork = l2*(2*nt*nlon+max(6*nlat,nlon)) + nlat*(2*nt*l1+1)
end if
allocate(dwwork(llwwork)) ; dwwork(:) = 0.0d0
ierror=3
call ivrtec(nlat,nlon,isym,nt,v_rot,u_rot,idvw,jdvw,wr,wi,mdab_s,ndab,wvhsec,lvhsec,dwwork,llwwork,pertrb_rot,ierror)
select case (ierror)
  case(0) 
    print*,'ivrtec : no error'
end select
!**************************************************************************************************************	
!**********idivec
allocate(pertrb_div(nlat,nlon)) ; pertrb_div(:,:) = 0.0d0
dwwork(:) = 0.0d0
ierror=3
call idivec(nlat,nlon,isym,nt,v_div,u_div,idvw,jdvw,dr,di,mdab_s,ndab,wvhsec,lvhsec,dwwork,llwwork,pertrb_div,ierror)
select case (ierror)
  case(0) 
    print*,'idivec : no error'
end select
end subroutine DIVROTV
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES IX
!#							     FLUXES:
!###################################################################################################################
!###################################################################################################################
! ---------- FLUXES
subroutine FLUXES(mdab_s,ndab,radius,Cmn1_r,Cmn1_i,Cmn2_r,Cmn2_i,fn,hn)
    implicit none
    integer :: mdab_s,ndab
    real*8 :: radius,n=0
    double precision, dimension (:), allocatable :: ji_conversion,jn,in,fn,hn
    double precision, dimension(:,:), allocatable :: jmn,imn,Cmn1_r,Cmn1_i,Cmn2_r,Cmn2_i

! Allocate and initialise the transfers and fluxes
allocate(jmn(mdab_s,ndab)) ; jmn    = 0.0d0
allocate(imn(mdab_s,ndab)) ; imn    = 0.0d0
allocate(jn(ndab)) ; jn    = 0.0d0
allocate(in(ndab)) ; in    = 0.0d0
!=====================================================
! Conversion factor
!=====================================================
allocate(ji_conversion(ndab)) ; ji_conversion = 0.0d0	
do n = 1, ndab-1 
  ji_conversion(n+1) = radius * radius / (n * (n + 1.0d0))
enddo
!=====================================================
! 4. Calculate energy and enstrophy interaction terms
!=====================================================
do n = 1, ndab-1
  jmn(:,n+1) = -0.25d0 * 2.0d0 * (Cmn1_r(:,n+1) * Cmn2_r(:,n+1) + Cmn1_i(:,n+1) * Cmn2_i(:,n+1))
  imn(:,n+1) = jmn(:,n+1) * ji_conversion(n+1)
enddo
write(6,'("Day ",i3," enstrophy transfer subb",es20.13)'), 1, sum(jmn(1,:)) + 2.0d0 * sum(jmn(2:mdab_s,:))
write(6,'("Day ",i3," net energy    transfer ",es20.13)'), 1, sum(imn(1,:)) + 2.0d0 * sum(imn(2:mdab_s,:))
! la somme etant de sum_(m = -n a n) on fait sum_(m=n)*2
! n = 0 term
jn(1) = jmn(1,1)
! n > 0 terms
do n = 1, ndab-1
  jn(n+1) = jmn(1,n+1) + 2.0d0 * sum(jmn(2:n+1,n+1))
enddo
! Convert to energy interaction term
in(:) = jn(:) * ji_conversion 
!====================================================================
! 5. Calculate energy and enstrophy spectral fluxes (BS13 Eqs 11-12)
!====================================================================
do n = 1, ndab-1 ! n > 0 terms (n = 0 term is zero)
  fn(n+1) = -sum(in(1:n)) ! Energy flux
  hn(n+1) = -sum(jn(1:n)) ! Enstrophy flux
enddo
end subroutine FLUXES
!###################################################################################################################
!######################################################################################################################
!###################################################################################################################
END PROGRAM spectra_analysis
