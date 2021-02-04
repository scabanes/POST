PROGRAM spectra_analysis


use netcdf

implicit none
! Ecrit par Simon Cabanes Septembre 2018...
!***************************************************************************************
!***************************************************************************************
!----------------------------------------------------------
!			PTS
!		Parameters to set:
!----------------------------------------------------------
! TGforcage =	1/ if a TG forage exist, 0/ if not
! Amp = 	est l'amplitude des perturbation dans le
! 		forcage TaylorGreen (TG)
! kf = 		est la frequence d injection des perturbat
!      		-ions du forcage TG.
! [omega]= s-1  Pour le calcul de la vorticite potentielle
! ptimestep = 	temps caracteristique d injection
!	      	des perturbation TG dans DYNAMICO
! NetcdefStatData = 1 to save netcedef data.
!----------------------------------------------------------
! Saving options:
! NetcdefStatData --> if 1, save data as a Netcdef file
! SavePhysicalFields --> if 1, save field maps
! SavePV --> if 1, save PV related data
! SaveSPA --> if 1, save Spectral related data.
!----------------------------------------------------------
! omega =  [1] -> 0.000165121 ;  [1/2] -> 0.00008256 ;  [1/4] -> 0.00004128
! ----------- Saturne:
!ptimestep=4756.50
real*8 :: Amp = 0.22,kf=56.,ptimestep= 3805.2,radius = 58232000, omega = 0.000165121, H = 375000!59000 ![radius] = m, depth: [H] = m, omega sat:0.000165121  ! [omega]= s-1 ! ptimestep = 1189.125 ; 3805.2 (itau_physics=10 ; 32) ;
! ----------- Jupiter:
!real*8 :: Amp = 0.22,kf=56.,ptimestep= 3805.2,radius = 69911000, omega = 0.0001750, H = 1000000 !omega Jupiter: 1.75 10-4 ! [omega]= s-1 !, depth: [H] = m from Nature Kaspi,  Rayon en metre Jupiter : radius = 71492000 (DYNAMICO) and 69911000 m pour les obs cassini de Jupiter.
!---------------------
integer :: TGforcage = 0,NetcdefStatData=1,SavePhysicalFields=1,SavePV=1,SaveSPA=1
! Latitude truncation for PV monotonization: here we supress polar values
!integer :: lat_PVmin = 2,lat_PVmax = 1 ! supress values at the pole only
integer :: lat_PVmin = 21,lat_PVmax = 20 ! supress values at the pole only
!integer :: lat_PVmin = 83,lat_PVmax = 82 ! supress values at 49 degres de latitude
!***************************************************************************************
!***************************************************************************************
integer :: iarg,narg,mt,istep,getfield,tstep,dimid_out(3),varid_out(4),latid,altid,colatid,lonid,tid,idinfos,ninfos,idninfos,Urmsid,hnid,PVzmid,PVsortedzmid,LMzmid,fnid,dstepsid,epsilonfid,ERid,ETid,EZid,hn_zmid,fn_zmid,Enid,En_divid,En_rotid,En_zmid,fnrotid,hnrotid,fndiv1id,hndiv1id,fndiv2id,hndiv2id
character (len=100), dimension(:), allocatable :: arg
character (len=100) :: tmp_char
logical :: is_file,NotFirst_reading=.false.,First_reading=.true.
!----------------------------------------------------------------------
integer :: nlon=0,nlat=0,nlat_mono=0,nalt,ntime,fieldsmap=0,ilat=0,ilon=0,ialt=0,typeS
real*8  :: factor, factor2, n=0,iz,ER=0.,EZ=0.,ET=0.,epsilonf_zm=0.,Idxtom=0.
double precision, dimension (:), allocatable :: fn_zm,hn_zm,fn_rot,hn_rot,fn_div1,hn_div1,fn_div2,hn_div2,Enmo_zm,Ene_zm,Urms,PVprof,ascend
double precision, dimension(:,:,:), allocatable :: vort_r, vort_i, div_r, div_i, sf_r, sf_i, vp_r, vp_i,PVsorted
double precision, dimension(:,:), allocatable :: vortf_r, vortf_i, divf_r, divf_i,Enmo,Ene,En_rot,En_div,PV_zm
double precision, dimension(:,:,:), allocatable :: vort,PV,term1_r,term1_i,term2_v,term2_w, term3,vort_uv_br,vort_uv_bi,vort_uv_cr,vort_uv_ci,vortv,vortu,uv_dot_grad_vort_r,uv_dot_grad_vort_i,rot_ugradu_i,rot_ugradu_r,div1_ugradu_i,div1_ugradu_r,fmn,hmn,Emn
double precision, dimension(:,:), allocatable :: fn,hn
integer :: mdab_v,mdab_s,mdab_smaller
!----------------------------------------------------------------------
double precision, dimension(:,:), allocatable :: uoff_f,voff_f,fvm,fwm,uoff,voff,vm,wm,fbr,fbi,fcr,fci,LM_zm,PVsorted_zm
double precision, dimension(:,:,:), allocatable :: vm3d,wm3d,u_DYN,v_DYN,UvectVT_lon,UvectVT_lat,u_rot,v_rot,u_div,v_div,udotu,vgrad,wgrad,vgradN,wgradN,L_M
double precision, dimension(:,:), allocatable :: fu,fv,rotuwP
double precision, dimension (:), allocatable :: time_counter_DYN,lat_DYN,lon_DYN,alt_DYN,colat_DYN,dsteps_DYN,lat_SP,lat_mono,lon_SP,colat_SP,epsilonf
double precision, dimension(:,:,:), allocatable ::  br,bi,cr,ci,Fdivuwr,Fdivuwi,Frotuwr,Frotuwi,Jmn,uur,uui,bbr,bbi,ccr,cci,ar,ai,divgraduur,divgraduui,rotgraduur,rotgraduui,div2_ugradu_r,div2_ugradu_i
character (len=100) :: file_netcdf,file_spectra,nickname
integer :: idfile,ncidF,ierror,idu,idv,idlat,idlon,idalt,idvw,jdvw,mdab,ndab,ivrt,jvrt,mdc,ndc,nt=1,ityp=0,isym=0,ioff,fgg,lon_dimid,lat_dimid,latmono_dimid,alt_dimid,t_dimid,wmid,vmid,vortid,Emnid,PVid,PVsortedid,deltaIdxdid, dimid3s(4),dimids(4),dimids_mono(4),dimids_1d(2),dimids_2d(3),dimidsmono_3d(3),dimids_3d(3), dimids_alttime(2), idvti,idfile2,idtime,idnlon,idnlat,idnalt,iddsteps,idSu,fuid,fvid,urotid,vrotid,udivid,vdivid
integer :: i,j,it,itt,in,jm,ig,l1,ll1,l2,vlwork,lvhsec,lwwork,llwwork,sldwork
double precision, dimension (:), allocatable :: work,ddwork,vwork
double precision, dimension (:,:,:), allocatable :: vm_3D,wm_3D,vt_3D
integer, dimension (:), allocatable :: L_Idx
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
! ------ > u_DYN & v_DYN: matrice issues des champs vents DYNAMICO ou OBSERVATIONS
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
 call check( nf90_inq_dimid(idfile,trim('altitude'),idnalt) )

 call check( nf90_inquire_dimension(idfile,idtime,tmp_char,ntime))
 call check( nf90_inquire_dimension(idfile,idnlon,tmp_char,nlon))
 call check( nf90_inquire_dimension(idfile,idnlat,tmp_char,nlat))
 call check( nf90_inquire_dimension(idfile,idnalt,tmp_char,nalt))

 call check( nf90_inq_varid(idfile,trim('time_counter'),idtime))
 call check( nf90_inq_varid(idfile,trim('lat'),idlat))
 call check( nf90_inq_varid(idfile,trim('lon'),idlon))
 call check( nf90_inq_varid(idfile,trim('altitude'),idalt))

 allocate(time_counter_DYN(ntime))
 allocate(lat_DYN(nlat))
 allocate(lon_DYN(nlon))
 allocate(alt_DYN(nalt))
 call check( nf90_get_var(idfile,idtime,time_counter_DYN,(/1/),(/ntime/)) ) !(lon,lat,iz,it)
 call check( nf90_get_var(idfile,idlat,lat_DYN,(/1/),(/nlat/)) ) !(lon,lat,iz,it)
 call check( nf90_get_var(idfile,idlon,lon_DYN,(/1/),(/nlon/)) ) !(lon,lat,iz,it)
 call check( nf90_get_var(idfile,idalt,alt_DYN,(/1/),(/nalt/)) ) !(lon,lat,iz,it)

 allocate(dsteps_DYN(ntime))
 call check( nf90_inq_varid(idfile,trim('dsteps'),iddsteps))
 call check( nf90_get_var(idfile,iddsteps,dsteps_DYN,(/1/),(/ntime/)) ) !(it)
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
  !call check( nf90_create("StatisticalData.nc", NF90_CLOBBER, ncid))
  call check( nf90_create(path="StatisticalData.nc", cmode=or(nf90_clobber,nf90_64bit_offset),ncid=ncidF))
  call check( nf90_def_dim(ncidF, "lon", nlon, lon_dimid))
  call check( nf90_def_dim(ncidF, "lat", nlat+1, lat_dimid))  !!!!! to change
  call check( nf90_def_dim(ncidF, "lat_mono", (nlat+1)-(lat_PVmin+lat_PVmax-1), latmono_dimid))  
  call check( nf90_def_dim(ncidF, "altitude", nalt, alt_dimid))
  call check( nf90_def_dim(ncidF, "time_counter", nf90_unlimited, t_dimid))
  dimids =  (/ lat_dimid, lon_dimid, alt_dimid, t_dimid /)!
  dimids_mono =  (/ latmono_dimid, lon_dimid, alt_dimid, t_dimid /)!
  dimid3s =  (/ lat_dimid, lat_dimid, alt_dimid, t_dimid /)!
  dimids_2d =  (/ lat_dimid, lon_dimid, t_dimid /)!
  dimids_1d =  (/ lat_dimid, t_dimid /)!
  dimids_3d =  (/ lat_dimid, alt_dimid, t_dimid /)!
  dimidsmono_3d =  (/ latmono_dimid, alt_dimid, t_dimid /)!
  dimids_alttime = (/ alt_dimid, t_dimid /)!
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! ---------------------------DFINITION OF VARIBLES
  call check( nf90_def_var(ncidF, "lon", NF90_DOUBLE, lon_dimid, lonid))
  call check( nf90_def_var(ncidF, "lat", NF90_DOUBLE, lat_dimid, latid))
  call check( nf90_def_var(ncidF, "altitude", NF90_DOUBLE, alt_dimid, altid))
  call check( nf90_def_var(ncidF, "time_counter", NF90_DOUBLE, t_dimid, tid))
  call check( nf90_def_var(ncidF, "colat", NF90_DOUBLE, lat_dimid, colatid))
!------------------------------------------------------------------------------------------
! from here time dependent varialbes
  call check( nf90_def_var(ncidF, "Urms", NF90_DOUBLE, dimids_alttime, Urmsid))

  call check( nf90_def_var(ncidF, "dsteps", NF90_DOUBLE, t_dimid, dstepsid))
  call check( nf90_def_var(ncidF, "ET", NF90_DOUBLE, t_dimid, ETid))
  call check( nf90_def_var(ncidF, "EZ", NF90_DOUBLE, t_dimid, EZid))
  call check( nf90_def_var(ncidF, "ER", NF90_DOUBLE, t_dimid, ERid))
  call check( nf90_def_var(ncidF, "epsilonf", NF90_DOUBLE, t_dimid, epsilonfid))
!------------------------------------------------------------------------------------------
if(SaveSPA == 1) then !--------------------------------------------------------------------- Spectral related
  call check( nf90_def_var(ncidF, "Enmo", NF90_DOUBLE, dimids_3d, En_zmid))
  call check( nf90_def_var(ncidF, "Ene", NF90_DOUBLE, dimids_3d, Enid))
  call check( nf90_def_var(ncidF, "En_rot", NF90_DOUBLE, dimids_3d, En_rotid))
  call check( nf90_def_var(ncidF, "En_div", NF90_DOUBLE, dimids_3d, En_divid))
  call check( nf90_def_var(ncidF, "fn_zm", NF90_DOUBLE, dimids_1d, fn_zmid))
  call check( nf90_def_var(ncidF, "hn_zm", NF90_DOUBLE, dimids_1d, hn_zmid))
  call check( nf90_def_var(ncidF, "fn", NF90_DOUBLE, dimids_3d, fnid))
  call check( nf90_def_var(ncidF, "hn", NF90_DOUBLE, dimids_3d, hnid))
  call check( nf90_def_var(ncidF, "Emn", NF90_DOUBLE, dimid3s, Emnid))

  !call check( nf90_def_var(ncidF, "fn_rot", NF90_DOUBLE, dimids_1d, fnrotid))
  !call check( nf90_def_var(ncidF, "hn_rot", NF90_DOUBLE, dimids_1d, hnrotid))
  !call check( nf90_def_var(ncidF, "fn_div1", NF90_DOUBLE, dimids_1d, fndiv1id))
  !call check( nf90_def_var(ncidF, "hn_div1", NF90_DOUBLE, dimids_1d, hndiv1id))
  !call check( nf90_def_var(ncidF, "fn_div2", NF90_DOUBLE, dimids_1d, fndiv2id))
  !call check( nf90_def_var(ncidF, "hn_div2", NF90_DOUBLE, dimids_1d, hndiv2id))
end if
if(SavePhysicalFields == 1) then !-------------------------------------------------------- 3D maps
    call check( nf90_def_var(ncidF, "wm", NF90_DOUBLE, dimids, wmid)) 
    call check( nf90_def_var(ncidF, "vm", NF90_DOUBLE, dimids, vmid)) 
    call check( nf90_def_var(ncidF, "vort", NF90_DOUBLE, dimids, vortid))
end if
if(SavePV == 1) then !-------------------------------------------------------------------- PV related
    ! 3D fields
    call check( nf90_def_var(ncidF, "PV", NF90_DOUBLE, dimids, PVid))
    call check( nf90_def_var(ncidF, "PVsorted", NF90_DOUBLE, dimids_mono, PVsortedid))
    call check( nf90_def_var(ncidF, "L_M", NF90_DOUBLE, dimids_mono, deltaIdxdid))
    ! 2D fields 
    call check( nf90_def_var(ncidF, "PV_zm", NF90_DOUBLE, dimids_3d, PVzmid))
    call check( nf90_def_var(ncidF, "PVsorted_zm", NF90_DOUBLE, dimidsmono_3d, PVsortedzmid))
    call check( nf90_def_var(ncidF, "LM_zm", NF90_DOUBLE, dimidsmono_3d, LMzmid))
end if

  !call check( nf90_def_var(ncidF, "u_rot", NF90_DOUBLE, dimids, urotid))
  !call check( nf90_def_var(ncidF, "v_rot", NF90_DOUBLE, dimids, vrotid))
  !call check( nf90_def_var(ncidF, "u_div", NF90_DOUBLE, dimids, udivid))
  !call check( nf90_def_var(ncidF, "v_div", NF90_DOUBLE, dimids, vdivid))
  
  !call check( nf90_def_var(ncidF, "fu", NF90_DOUBLE, dimids, fuid))
  !call check( nf90_def_var(ncidF, "fv", NF90_DOUBLE, dimids, fvid))


!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! ---------------------------VARIABLES FILL INSTRUCTIONS for time dependent variables

  call check( NF90_PUT_ATT  (ncidF, tid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, Urmsid, "_FillValue", NF90_FILL_DOUBLE) )

  call check( NF90_PUT_ATT  (ncidF, dstepsid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, ETid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, EZid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, ERid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, epsilonfid, "_FillValue", NF90_FILL_DOUBLE) )

!------------------------------------------------------------------------------------------
if(SaveSPA == 1) then !--------------------------------------------------------------------- Spectral related
  call check( NF90_PUT_ATT  (ncidF, En_zmid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, Enid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, En_rotid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, En_divid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, fn_zmid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, hn_zmid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, fnid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, hnid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, Emnid, "_FillValue", NF90_FILL_DOUBLE) )

  !call check( NF90_PUT_ATT  (ncidF, fnrotid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, hnrotid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, hndiv1id, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, fndiv1id, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, hndiv2id, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, fndiv2id, "_FillValue", NF90_FILL_DOUBLE) )
end if


if(SavePhysicalFields == 1) then !-------------------------------------------------------- 3D maps
  call check( NF90_PUT_ATT  (ncidF, wmid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, vmid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, vortid, "_FillValue", NF90_FILL_DOUBLE) ) 
end if
if(SavePV == 1) then !-------------------------------------------------------------------- PV related
    ! 3D fields 
  call check( NF90_PUT_ATT  (ncidF, PVid, "_FillValue", NF90_FILL_DOUBLE) ) 
  call check( NF90_PUT_ATT  (ncidF, PVsortedid, "_FillValue", NF90_FILL_DOUBLE) ) 
  call check( NF90_PUT_ATT  (ncidF, deltaIdxdid, "_FillValue", NF90_FILL_DOUBLE) ) 
    ! 2D fields 
  call check( NF90_PUT_ATT  (ncidF, PVzmid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, PVsortedzmid, "_FillValue", NF90_FILL_DOUBLE) )
  call check( NF90_PUT_ATT  (ncidF, LMzmid, "_FillValue", NF90_FILL_DOUBLE) )
end if


  !call check( NF90_PUT_ATT  (ncidF, urotid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, vrotid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, udivid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, vdivid, "_FillValue", NF90_FILL_DOUBLE) )

  !call check( NF90_PUT_ATT  (ncidF, fuid, "_FillValue", NF90_FILL_DOUBLE) )
  !call check( NF90_PUT_ATT  (ncidF, fvid, "_FillValue", NF90_FILL_DOUBLE) )


  call check( nf90_enddef(ncidF) )
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
!#--------------------------------------------------------------------------------------------------------------------------------------
do it = istep,mt,tstep !#----------------------------------------------------------------------------------------------- START TIME LOOP

if (NotFirst_reading) then
  nlat = nlat - 1
else
  print*,'First netcdf file reading...'
  NotFirst_reading = .true.
end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!										      LOAD FIELD:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allocate(u_DYN(nlon,nlat,nalt))
allocate(v_DYN(nlon,nlat,nalt))

print*, it
 call check( nf90_inq_varid(idfile,trim('u'),idu))
 call check( nf90_inq_varid(idfile,trim('v'),idv))

 call check( nf90_get_var(idfile,idu,u_DYN,(/1,1,1,it/),(/nlon,nlat,nalt,1/)) ) !(lon,lat,iz,it)
 call check( nf90_get_var(idfile,idv,v_DYN,(/1,1,1,it/),(/nlon,nlat,nalt,1/)) ) !(lon,lat,iz,it)

print*, 'shape of u_DYN: ',shape(u_DYN)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!										    Forcage (2D):
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!========================================================
! 1. Project forcage from DYNAMICO -- to --> Spherepack
!========================================================
! ------------WIND field grid South->North
allocate(uoff_f(nlon,nlat)) ; uoff_f = fu
allocate(voff_f(nlon,nlat)) ; voff_f = fv
!------------WIND Northern->southern hemisphere ordered data
allocate(fvm(nlat+1,nlon)) ; fvm(:,:) = 0.0d0
allocate(fwm(nlat+1,nlon)) ; fwm(:,:) = 0.0d0

 call GRIDMATHS(nlon,nlat,uoff_f,voff_f,fvm,fwm)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!								             Velocity field (3D):
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!========================================================
! 1. Project velocity from DYNAMICO -- to --> Spherepack
!========================================================
! WIND Northern->southern hemisphere ordered data
allocate(vm3d(nlat+1,nlon,nalt)) !; vm3d(:,:,:) = 0.0d0
allocate(wm3d(nlat+1,nlon,nalt)) !; wm3d(:,:,:) = 0.0d0
allocate(Urms(nalt))
!#------------------------------------------------------------------------------------------------------------------------------------
do iz = 1,nalt !#------------------------------------------------------------------------------------------------- START Altitude LOOP
! WIND field grid South->North
allocate(uoff(nlon,nlat)) ; uoff = u_DYN(:,:,iz)
allocate(voff(nlon,nlat)) ; voff = v_DYN(:,:,iz)
! WIND Northern->southern hemisphere ordered data
allocate(vm(nlat+1,nlon)) ; vm(:,:) = 0.0d0
allocate(wm(nlat+1,nlon)) ; wm(:,:) = 0.0d0

 call GRIDMATHS(nlon,nlat,uoff,voff,vm,wm)

vm3d(:,:,iz) = vm
wm3d(:,:,iz) = wm
!Urms(iz) = sqrt(SUM(vm*vm + wm*wm)/((nlat+1)*nlon))
Urms(iz) = sqrt(SUM(wm*wm)/((nlat+1)*nlon))
deallocate(uoff)
deallocate(voff)
deallocate(vm)
deallocate(wm)
end do!#------------------------------------------------------------------------------------------------------------- END Altitude LOOP
!#-------------------------------------------------------------------------------------------------------------------------------------
print*,'Urms = '
print*,Urms
!######################################################################################
! #############################################     NOUVELLES COORDONNEES LAT-LON-COLAT
! vshifte et geo2mathv redefinissent la grille comme suit:
!***********************************
nlat = nlat+1
print*,'nlat =',nlat-1,'---->',nlat
!***********************************
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
!***************************************************************
! Here we give Maximum value of m zonal spherical hamonic index.
! This is diffferent while analysing a vector or a scalar field.
! Then mdab_v and mdab_s are rescpetively maximum m for a vector
! and scalar field. Number of total index n remain the same i.e 
! ndab.
!***************************************************************
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
!										    Forcage (2D):
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!========================================================
! 2. Vector Spectral Analysis 2D
!========================================================
allocate(fbr(mdab_v,ndab)) ; fbr(:,:) = 0.0d0
allocate(fbi(mdab_v,ndab)) ; fbi(:,:) = 0.0d0
allocate(fcr(mdab_v,ndab)) ; fcr(:,:) = 0.0d0
allocate(fci(mdab_v,ndab)) ; fci(:,:) = 0.0d0
 call VSA2D(nlon,nlat,1,ityp,mdab_v,ndab,fvm,fwm,fbr,fbi,fcr,fci)

allocate(vortf_r(mdab_s,ndab)) ; vortf_r = 0.0d0
allocate(vortf_i(mdab_s,ndab)) ; vortf_i = 0.0d0
allocate(divf_r (mdab_s,ndab)) ; divf_r  = 0.0d0
allocate(divf_i (mdab_s,ndab)) ; divf_i  = 0.0d0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!								             Velocity field (3D):
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
allocate(br(mdab_v,ndab,nalt)) ; br(:,:,:) = 0.0d0
allocate(bi(mdab_v,ndab,nalt)) ; bi(:,:,:) = 0.0d0
allocate(cr(mdab_v,ndab,nalt)) ; cr(:,:,:) = 0.0d0
allocate(ci(mdab_v,ndab,nalt)) ; ci(:,:,:) = 0.0d0
 call VSA(nlon,nlat,nalt,ityp,mdab_v,ndab,vm3d,wm3d,br,bi,cr,ci)
! Spherical harmonic coeffs of scalar fields
allocate(vort_r(mdab_s,ndab,nalt)) ; vort_r = 0.0d0
allocate(vort_i(mdab_s,ndab,nalt)) ; vort_i = 0.0d0
allocate(div_r (mdab_s,ndab,nalt)) ; div_r  = 0.0d0
allocate(div_i (mdab_s,ndab,nalt)) ; div_i  = 0.0d0
allocate(sf_r  (mdab_s,ndab,nalt)) ; sf_r   = 0.0d0
allocate(sf_i  (mdab_s,ndab,nalt)) ; sf_i   = 0.0d0
allocate(vp_r  (mdab_s,ndab,nalt)) ; vp_r   = 0.0d0
allocate(vp_i  (mdab_s,ndab,nalt)) ; vp_i   = 0.0d0
!-----------------------------------------------------------------
do n = 1, ndab-2
!*********************************************************************
  ! Scaling factor
  ! Need to divide through by radius as gradient operator is on unit sphere  
  ! Here e.g. vort_r/i = n(n+1)^1/2 * cr/i /R and the streamfunction
  ! being vort = Laplacian sf, we have:
  !  n(n+1)^1/2 * cr/i * 1/R = -n(n+1)/R^2 sf_r/i and we can obtain the 
  ! streamfunction coefficients from sf_r/i = - R/n(n+1)^1/2. 
  factor = sqrt(1.0d0 * n * (n+1.0d0)) / radius
!*********************************************************************
!-----------------------  FORCAGE FIELD 2D
  ! Vorticity spectrum: Spherepack 2.0 after Eq. 4.15
  vortf_r(1:mdab_smaller,n+1) =  fcr(1:mdab_smaller,n+1) * factor
  vortf_i(1:mdab_smaller,n+1) =  fci(1:mdab_smaller,n+1) * factor

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  divf_r (1:mdab_smaller,n+1) = -fbr(1:mdab_smaller,n+1) * factor
  divf_i (1:mdab_smaller,n+1) = -fbi(1:mdab_smaller,n+1) * factor
!*********************************************************************
!-----------------------  VELOCITY FIELD 3D
  ! Vorticity spectrum: Spherepack 2.0 after Eq. 4.15
  vort_r(1:mdab_smaller,n+1,1:nalt) =  cr(1:mdab_smaller,n+1,1:nalt) * factor
  vort_i(1:mdab_smaller,n+1,1:nalt) =  ci(1:mdab_smaller,n+1,1:nalt) * factor

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  div_r (1:mdab_smaller,n+1,1:nalt) = -br(1:mdab_smaller,n+1,1:nalt) * factor
  div_i (1:mdab_smaller,n+1,1:nalt) = -bi(1:mdab_smaller,n+1,1:nalt) * factor

  ! Streamfunction spectrum: Spherepack 2.0 combination of Eqs. 4.1 
  ! for streamfunction, 4.6, and 4.8
  sf_r  (1:mdab_smaller,n+1,1:nalt) = -cr(1:mdab_smaller,n+1,1:nalt) / factor
  sf_i  (1:mdab_smaller,n+1,1:nalt) = -ci(1:mdab_smaller,n+1,1:nalt) / factor

  ! Velocity potential spectrum: Spherepack 2.0 combination of Eqs. 4.1 
  ! for velocity potential, 4.6, and 4.8
  vp_r  (1:mdab_smaller,n+1,1:nalt) =  br(1:mdab_smaller,n+1,1:nalt) / factor
  vp_i  (1:mdab_smaller,n+1,1:nalt) =  bi(1:mdab_smaller,n+1,1:nalt) / factor
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
!******** Vorticite relative **************************************
! truncation in latitude to cancel, at least, polar value. 
nlat_mono = nlat-(lat_PVmin+lat_PVmax-1)
allocate(lat_mono(nlat_mono)) ; lat_mono = 0.0d0
lat_mono = lat_SP(lat_PVmin:nlat-lat_PVmax)
!print*, lat_mono
!******************************************************************
!size nlat
allocate(vort(nlat,nlon,nalt)) ; vort = 0.0d0
allocate(PV(nlat,nlon,nalt)) ; PV = 0.0d0
allocate(PV_zm(nlat,nalt)) ; PV_zm = 0.0d0
! size nlat_mono
allocate(PVsorted(nlat_mono,nlon,nalt)) ; PVsorted = 0.0d0
allocate(PVsorted_zm(nlat_mono,nalt)) ; PVsorted_zm = 0.0d0
allocate(L_M(nlat_mono,nlon,nalt)) ; L_M = 0.0d0
allocate(LM_zm(nlat_mono,nalt)) ; LM_zm = 0.0d0
allocate(PVprof(nlat_mono)); PVprof = 0.0d0
allocate(ascend(nlat_mono)); ascend = 0.0d0 ! ascendent sorted data
allocate(L_Idx(nlat_mono)); L_Idx = 0 ! Index sorted 
 call SHS(nlon,nlat,nalt,isym,vort,mdab_s,ndab,vort_r,vort_i)
!******** Vorticite potentielle ***********************************
!                  &
!	     Monotonization
!******************************************************************
! Call a sorting function that: sorting(ndata,VtoSort,Vsorted,typeS,Idx)
! if typeS == 1 ascendent sorting & typeS == 2 descendent sorting
! The ascendent choice depend on the data order. The output are,
! ndata: is the length of VtoSort(ndata)
! Vtosort: is a single vector
! Vsorted: is a vector once Vtosort sorted it has the same length
! Idx: are index unsorted - index sorted --> typical sorted scale
!
! Tyical sorted scale is converted from index to meter via Idxtom
!******************************************************************
typeS=2
Idxtom = (PI*radius)/(nlat-1) ![(PI*radius)/(nlat-1)] = [Idxto m]= m et la conversion d indice en m en multipliant par le pas en m
!print*, Idxtom
do ialt = 1, nalt ! ----------- loop on the altitude
  do ilon = 1, nlon ! ----------- loop on the longitude
    PV(:,ilon,ialt) = ( vort(:,ilon,ialt) + 2.0d0 * omega * sin(lat_SP/(180./PI)) ) / H
    ! Sorting on instantaneous -- truncation in latitude --
    PVprof=PV(lat_PVmin:nlat-lat_PVmax,ilon,ialt)
    call sorting(nlat_mono,PVprof,ascend,typeS,L_Idx) 
    PVsorted(:,ilon,ialt) = ascend
    L_M(:,ilon,ialt) = abs(L_Idx)*Idxtom
  enddo ! ----------- end loop on the longitude
    ! Sorting on zonal mean PV profil
    PV_zm(:,ialt) = sum(PV(:,:,ialt),DIM=2)/nlon
    PVprof= PV_zm(lat_PVmin:nlat-lat_PVmax,ialt)
    call sorting(nlat_mono,PVprof,ascend,typeS,L_Idx) 
    PVsorted_zm(:,ialt) = ascend
    LM_zm(:,ialt) = abs(L_Idx)*Idxtom
enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!							Rotational & divergent velocities vector:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!***************************************************************
! Here we compute the rotational and divergent part of the
! velocity field, u_rot, v_rot and u_div, v_div.
!***************************************************************
! rotational velocity synthesis
allocate(u_rot(nlat,nlon,nalt)) ; u_rot(:,:,:) = 0.0d0
allocate(v_rot(nlat,nlon,nalt)) ; v_rot(:,:,:) = 0.0d0
allocate(u_div(nlat,nlon,nalt)) ; u_div(:,:,:) = 0.0d0
allocate(v_div(nlat,nlon,nalt)) ; v_div(:,:,:) = 0.0d0
 call DIVROTV(nlon,nlat,nalt,isym,v_rot,u_rot,v_div,u_div,mdab_s,ndab,vort_r,vort_i,div_r,div_i)
v_rot = v_rot * radius
u_rot = u_rot * radius
v_div = v_div * radius
u_div = u_div * radius

!########################################################################################################################################
!========================================================================================================================================
!					     		Kinetic Energy Spectra
!========================================================================================================================================
!########################################################################################################################################
allocate(Emn(mdab_s,ndab,nalt)) ; Emn    = 0.0d0
allocate(Enmo(ndab,nalt)) ; Enmo    = 0.0d0
allocate(Ene(ndab,nalt)) ; Ene    = 0.0d0
allocate(Enmo_zm(ndab)) ; Enmo_zm    = 0.0d0
allocate(Ene_zm(ndab)) ; Ene_zm    = 0.0d0
 call ES(mdab_s,ndab,nalt,radius,vort_r,vort_i,div_r,div_i,Emn,Enmo,Ene,Enmo_zm,Ene_zm,ER,EZ,ET,En_rot,En_div)
!########################################################################################################################################
!========================================================================================================================================
!					     		SPECTRAL FLUXES --- div(u) = 0
!========================================================================================================================================
!########################################################################################################################################


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!							        Field: div . (vorticity * u_rot):
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==================================================================
! 1. Set up required data for this call
!==================================================================
allocate(term1_r(mdab_s,ndab,nalt))
allocate(term1_i(mdab_s,ndab,nalt))
allocate(term2_v(nlat,nlon,nalt))
allocate(term2_w(nlat,nlon,nalt))
allocate(term3  (nlat,nlon,nalt))
! Total transfers
term1_r = vort_r
term1_i = vort_i
term2_v = v_rot
term2_w = u_rot
term3   = vort
! Allocate and initialise the transfers and fluxes

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
allocate(vort_uv_br(mdab_v,ndab,nalt)) ; vort_uv_br = 0.0d0
allocate(vort_uv_bi(mdab_v,ndab,nalt)) ; vort_uv_bi = 0.0d0
allocate(vort_uv_cr(mdab_v,ndab,nalt)) ; vort_uv_cr = 0.0d0
allocate(vort_uv_ci(mdab_v,ndab,nalt)) ; vort_uv_ci = 0.0d0
! Field : vorticity * u_rot
allocate(vortv(nlat,nlon,nalt))
allocate(vortu(nlat,nlon,nalt))
vortv = term3*term2_v ! vort * v_rot
vortu = term3*term2_w ! vort * u_rot
! Vector Spectral Analysis
 call VSA(nlon,nlat,nalt,ityp,mdab_v,ndab,vortv,vortu,vort_uv_br,vort_uv_bi,vort_uv_cr,vort_uv_ci)
!======================================================
! 3. Calculate divergence spectrum of vorticity * u_rot
!======================================================
! Allocate and initialise space
allocate(uv_dot_grad_vort_r(mdab_s,ndab,nalt)) ; uv_dot_grad_vort_r = 0.0d0
allocate(uv_dot_grad_vort_i(mdab_s,ndab,nalt)) ; uv_dot_grad_vort_i = 0.0d0
! factor to apply to divergent coefficients.
do n = 1, ndab-2

  ! Scaling factor
  ! Need to divide through by radius as gradient operator is on unit sphere  
  factor2 = sqrt(1.0d0 * n * (n+1.0d0)) / radius

  ! Divergence spectrum: Spherepack 2.0 after Eq. 4.16
  uv_dot_grad_vort_r(1:mdab_smaller,n+1,1:nalt) = -vort_uv_br(1:mdab_smaller,n+1,1:nalt) * factor2
  uv_dot_grad_vort_i(1:mdab_smaller,n+1,1:nalt) = -vort_uv_bi(1:mdab_smaller,n+1,1:nalt) * factor2

enddo
!=====================================================
! 4. Calculate energy and enstrophy interaction terms
!=====================================================
! Allocate and initialise the transfers and fluxes
allocate(fn(ndab,nalt)) ; fn = 0.0d0
allocate(hn(ndab,nalt)) ; hn = 0.0d0
allocate(fn_zm(ndab)) ; fn_zm = 0.0d0
allocate(hn_zm(ndab)) ; hn_zm = 0.0d0
!
allocate(fmn(mdab_s,ndab,nalt)) ; fmn = 0.0d0
allocate(hmn(mdab_s,ndab,nalt)) ; hmn = 0.0d0
 call FLUXES(mdab_s,ndab,nalt,radius,term1_r,term1_i,uv_dot_grad_vort_r,uv_dot_grad_vort_i,fn,hn,fn_zm,hn_zm,fmn,hmn)

!########################################################################################################################################
!========================================================================================================================================
!					     		Epsilon forcage = < u.f >
!========================================================================================================================================
!########################################################################################################################################
if (TGforcage==1) then
	allocate(epsilonf(nalt)) ; epsilonf		= 0.0d0
	call epsilonforcage(mdab_s,ndab,nalt,radius,vort_r,vort_i,div_r,div_i,vortf_r,vortf_i,divf_r,divf_i,epsilonf)
	epsilonf_zm = sum(epsilonf)/nalt
else
	allocate(epsilonf(nalt)) ; epsilonf		= 0.0d0
end if






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
    call check( NF90_PUT_VAR (ncidF, lonid, lon_SP) )
    call check( NF90_PUT_VAR (ncidF, latid, lat_SP) )
    call check( NF90_PUT_VAR (ncidF, altid, alt_DYN) )
    call check( NF90_PUT_VAR (ncidF, colatid, colat_SP) )
  else
    First_reading = .false.
  end if
!------------------------------------------------------------------------------------------
   call check( nf90_put_var( ncidF, tid, time_counter_DYN(it), start=(/it/) ) )
   call check( nf90_put_var( ncidF, Urmsid, Urms,(/1,it/),(/nalt,1/) ) )

   call check( nf90_put_var( ncidF, dstepsid, dsteps_DYN(it), start=(/it/) ) )
   call check( nf90_put_var( ncidF, ETid, ET, start=(/it/) ) )
   call check( nf90_put_var( ncidF, EZid, EZ, start=(/it/) ) )
   call check( nf90_put_var( ncidF, ERid, ER, start=(/it/) ) )
   call check( nf90_put_var( ncidF, epsilonfid, epsilonf_zm, start=(/it/) ) )

!------------------------------------------------------------------------------------------
if(SaveSPA == 1) then !--------------------------------------------------------------------- Spectral related
   call check( nf90_put_var( ncidF, En_zmid, Enmo,(/1,1,it/),(/nlat,nalt,1/) ) )
   call check( nf90_put_var( ncidF, Enid, Ene,(/1,1,it/),(/nlat,nalt,1/) ) )
   call check( nf90_put_var( ncidF, En_rotid, En_rot,(/1,1,it/),(/nlat,nalt,1/) ) )
   call check( nf90_put_var( ncidF, En_divid, En_div,(/1,1,it/),(/nlat,nalt,1/) ) )
   call check( nf90_put_var( ncidF, fn_zmid, fn_zm,(/1,it/),(/nlat,1/) ) )
   call check( nf90_put_var( ncidF, hn_zmid, hn_zm,(/1,it/),(/nlat,1/) ) )
   call check( nf90_put_var( ncidF, fnid, fn,(/1,1,it/),(/nlat,nalt,1/) ) )
   call check( nf90_put_var( ncidF, hnid, hn,(/1,1,it/),(/nlat,nalt,1/) ) )
   call check( nf90_put_var( ncidF, Emnid, Emn,(/1,1,1,it/),(/nlat,nlat,nalt,1/) ) )

 !call check( nf90_put_var( ncidF, fnrotid, fn_rot,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncidF, hnrotid, hn_rot,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncidF, fndiv1id, fn_div1,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncidF, hndiv1id, hn_div1,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncidF, fndiv2id, fn_div2,(/1,it/),(/nlat,1/) ) )
   !call check( nf90_put_var( ncidF, hndiv2id, hn_div2,(/1,it/),(/nlat,1/) ) )

end if
if(SavePhysicalFields == 1) then !----------------------------------------------------------- 3D maps
     call check( nf90_put_var( ncidF, wmid, wm3d,(/1,1,1,it/),(/nlat,nlon,nalt,1/) ) ) 
     call check( nf90_put_var( ncidF, vmid, vm3d,(/1,1,1,it/),(/nlat,nlon,nalt,1/) ) ) 
     call check( nf90_put_var( ncidF, vortid, vort,(/1,1,1,it/),(/nlat,nlon,nalt,1/) ) )
end if
if(SavePV == 1) then !---------------------------------------------------------------------- PV related
    ! 3D fields 
   call check( nf90_put_var( ncidF, PVid, PV,(/1,1,1,it/),(/nlat,nlon,nalt,1/) ) )
   call check( nf90_put_var( ncidF, PVsortedid, PVsorted,(/1,1,1,it/),(/nlat_mono,nlon,nalt,1/) ) ) 
   call check( nf90_put_var( ncidF, deltaIdxdid, L_M,(/1,1,1,it/),(/nlat_mono,nlon,nalt,1/) ) ) 
    ! 2D fields 
   call check( nf90_put_var( ncidF, PVzmid, PV_zm,(/1,1,it/),(/nlat,nalt,1/) ) )
   call check( nf90_put_var( ncidF, PVsortedzmid, PVsorted_zm,(/1,1,it/),(/nlat_mono,nalt,1/) ) )
   call check( nf90_put_var( ncidF, LMzmid, LM_zm,(/1,1,it/),(/nlat_mono,nalt,1/) ) )
end if

   !call check( nf90_put_var( ncidF, urotid, transpose(u_rot),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncidF, vrotid, transpose(v_rot),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncidF, udivid, transpose(u_div),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncidF, vdivid, transpose(v_div),(/1,1,it/),(/nlon,nlat,1/) ) )

   !call check( nf90_put_var( ncidF, ERid, ER,(/it/),(/1/) ) )
   !call check( nf90_put_var( ncidF, EZid, EZ,(/it/),(/1/) ) )

   !call check( nf90_put_var( ncidF, fuid, transpose(fwm),(/1,1,it/),(/nlon,nlat,1/) ) )
   !call check( nf90_put_var( ncidF, fvid, transpose(fvm),(/1,1,it/),(/nlon,nlat,1/) ) )
end if

!##########################################################
!##################################### DEALLOCATE ############
!##########################################################
deallocate(u_DYN)
deallocate(v_DYN)
deallocate(uoff_f)
deallocate(voff_f)
deallocate(vm3d)
deallocate(wm3d)
deallocate(Urms)
deallocate(fvm)
deallocate(fwm)
deallocate(br)
deallocate(bi)
deallocate(cr)
deallocate(ci)
deallocate(fbr)
deallocate(fbi)
deallocate(fcr)
deallocate(fci)
deallocate(vortf_r)
deallocate(vortf_i)
deallocate(divf_r)
deallocate(divf_i)
deallocate(v_rot)
deallocate(u_rot)
deallocate(lon_SP)
deallocate(lat_SP)
deallocate(colat_SP)
deallocate(vort_r)
deallocate(vort_i)
deallocate(div_r)
deallocate(div_i) 
deallocate(sf_r)  
deallocate(sf_i)  
deallocate(vp_r) 
deallocate(vp_i)
deallocate(vortv)
deallocate(vortu)
deallocate(vort)
deallocate(lat_mono)
deallocate(PV)
deallocate(PV_zm)
deallocate(PVprof)
deallocate(PVsorted)
deallocate(PVsorted_zm)
deallocate(ascend)
deallocate(L_M)
deallocate(LM_zm)
deallocate(L_Idx)
deallocate(u_div)
deallocate(v_div)
deallocate(term1_r)
deallocate(term1_i)
deallocate(term2_v)
deallocate(term2_w)
deallocate(term3)
deallocate(hn)
deallocate(fn)
deallocate(hmn)
deallocate(fmn)
deallocate(hn_zm)
deallocate(fn_zm)
deallocate(vort_uv_br)
deallocate(vort_uv_bi)
deallocate(vort_uv_cr)
deallocate(vort_uv_ci)
deallocate(uv_dot_grad_vort_r)
deallocate(uv_dot_grad_vort_i)
deallocate(Emn)
deallocate(Enmo)
deallocate(Ene)
deallocate(Enmo_zm)
deallocate(Ene_zm)
deallocate(epsilonf)
deallocate(En_rot)
deallocate(En_div)

end do!#---------------------------------------------------------------------------------------------------------------- ENS TIME LOOP
!#------------------------------------------------------------------------------------------------------------------------------------
if(NetcdefStatData == 1) then
  call check( nf90_close(ncidF) )
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
!#							 SUBROUTINES III 
!#							    GRIDMATHS
!###################################################################################################################
!###################################################################################################################
! ---------- GRIDMATHS
subroutine GRIDMATHS(nlon,nlat,Guoff,Gvoff,Gvm,Gwm)
    implicit none
    integer :: nlon,nlat,nlatp1,lsav,ioff,ierror,lGwork,ig
    double precision, dimension (:), allocatable :: wsav,Gwork,gmwork
    double precision, dimension(:,:), allocatable :: Guoff,Gvoff,ureg,vreg,ug,vg,Gvm,Gwm

!#############################################################
!######################## BLOCK 0 ################################	REGULAR GRID
!##############################################################
!*************************************************************
! we initialy have a lat=[-89.75,89.75] et lon=[0.25,359.75]
! we shift to a grid lat=[-90,90] et lon=[0,359.5]
! It only a shift of 0.5* grid increment, so conserve South --> North
!*************************************************************
! --------------------vshifti function (initialisations for vshifte function)
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
! ------------WIND field
allocate(ureg(nlon,nlat+1)) ; ureg = 0.0d0
allocate(vreg(nlon,nlat+1)) ; vreg = 0.0d0
ioff = 0
if (mod(nlon,2) == 0) then ! check if it is even
  lGwork = 2*nlon*(nlat+1)
else
  lGwork = nlon*(5*nlat+1)
end if
allocate(Gwork(lGwork)) ; Gwork(:)=0.0d0
ierror=3
call vshifte(ioff,nlon,nlat,Guoff,Gvoff,ureg,vreg,wsav,lsav,Gwork,lGwork,ierror)
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
!#############################################################
!######################## BLOCK 0 ################################	GEO2MATH
!##############################################################
!*************************************************************
! 1) flip to North --> south (hemisphere ordered data)
! 2) from latitude to colatitued component, i.e. vg = -vm 
! 3) nlat become the first component and nlon the second.
! 4) The grid is now colat=[0,180] et lon=[0,359.5]
ig = 0 ! 0 if grid South->North & 1 if grid North->South : North (90) to south (-90)
!*************************************************************
!Southern->Northern hemisphere ordered data
allocate(ug(nlon,nlat+1)) ; ug = ureg
allocate(vg(nlon,nlat+1)) ; vg = vreg
!------------
allocate(gmwork(nlon*(nlat+1))) ; gmwork(:) = 0.0d0
ierror=3
nlatp1 = nlat+1
call geo2mathv(ig,nlon,nlatp1,ug,vg,Gvm,Gwm,gmwork)
end subroutine GRIDMATHS
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
subroutine VSA(nlon,nlat,nt,ityp,mdab,ndab,vm,wm,bbr,bbi,ccr,cci)
    implicit none
    integer :: nlon,nlat,l1,l2,lvhaec,ldwork,ierror,nt,ityp,idvw,jdvw,mdab,ndab,lwork
    double precision, dimension (:), allocatable :: wvhaec,dwork,work
    double precision, dimension(:,:,:), allocatable :: bbr,bbi,ccr,cci,vm,wm

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
call vhaec(nlat,nlon,ityp,nt,vm,wm,idvw,jdvw,bbr,bbi,ccr,cci,mdab,ndab,wvhaec,lvhaec,work,lwork,ierror)
select case (ierror)
  case(0) 
    print*,'no error in subroutine VSA'
end select
deallocate(work)
end subroutine VSA
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES IV
!#							       VSA2D:
!#						    Vector Spectral Analysis 2D
!###################################################################################################################
!###################################################################################################################
! ---------- VSA2D
subroutine VSA2D(nlon,nlat,nt,ityp,mdab,ndab,vm,wm,bbr,bbi,ccr,cci)
    implicit none
    integer :: nlon,nlat,l1,l2,lvhaec,ldwork,ierror,nt,ityp,idvw,jdvw,mdab,ndab,lwork
    double precision, dimension (:), allocatable :: wvhaec,dwork,work
    double precision, dimension(:,:), allocatable :: bbr,bbi,ccr,cci,vm,wm ! only this line change compare to subroutine VSA

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
call vhaec(nlat,nlon,ityp,nt,vm,wm,idvw,jdvw,bbr,bbi,ccr,cci,mdab,ndab,wvhaec,lvhaec,work,lwork,ierror)
select case (ierror)
  case(0) 
    print*,'no error in subroutine VSA2D'
end select
deallocate(work)
end subroutine VSA2D
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
    double precision, dimension(:,:,:), allocatable :: aar,aai,g

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
    double precision, dimension (:), allocatable :: wvhsec,wwork,dwwork,pertrb_rot,pertrb_div
    double precision, dimension(:,:,:), allocatable :: wr,wi,dr,di,u_rot,v_rot,u_div,v_div

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
!allocate(pertrb_rot(nlat,nlon)) ; pertrb_rot(:,:) = 0.0d0
allocate(pertrb_rot(nt)) ; pertrb_rot(:) = 0.0d0
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
!allocate(pertrb_div(nlat,nlon)) ; pertrb_div(:,:) = 0.0d0
allocate(pertrb_div(nt)) ; pertrb_div(:) = 0.0d0
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
subroutine FLUXES(mdab_s,ndab,nalt,radius,Cmn1_r,Cmn1_i,Cmn2_r,Cmn2_i,fn,hn,fn_zm,hn_zm,fmn,hmn)
    implicit none
    integer :: mdab_s,ndab,nalt,id,im
    real*8 :: radius,n=0,m=0
    double precision, dimension (:), allocatable :: ji_conversion,fn_zm,hn_zm,jn_zm,in_zm
    double precision, dimension(:,:), allocatable :: jmn_zm,imn_zm,jn,in,fn,hn
    double precision, dimension(:,:,:), allocatable :: jmn,imn,Cmn1_r,Cmn1_i,Cmn2_r,Cmn2_i,fmn,hmn

! Allocate and initialise the transfers and fluxes
allocate(jmn(mdab_s,ndab,nalt)) ; jmn    = 0.0d0
allocate(imn(mdab_s,ndab,nalt)) ; imn    = 0.0d0
allocate(jmn_zm(mdab_s,ndab)) ; jmn_zm = 0.0d0
allocate(imn_zm(mdab_s,ndab)) ; imn_zm = 0.0d0
allocate(jn(ndab,nalt)) ; jn    = 0.0d0
allocate(in(ndab,nalt)) ; in    = 0.0d0
allocate(jn_zm(ndab)) ; jn_zm = 0.0d0
allocate(in_zm(ndab)) ; in_zm = 0.0d0
!
!allocate(jmn(mdab_s,ndab,nalt)) ; fmn    = 0.0d0			----------ici commente je sais pas pourquoi
!allocate(imn(mdab_s,ndab,nalt)) ; hmn    = 0.0d0			----------ici commente je sais pas pourquoi

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
!*************************************************************************
! 2D Enstrophy and energy interaction terms (BS83 Eqs 15-16 / BS13 Eqs 8-9)
! Factor of two is required as there are two terms in the square brackets
!*************************************************************************
! Jmn(m,n,altitude)
write(6,'(a)') "2D enstrophy and energy interaction terms for each day ..."
do id = 1, nalt
  do n = 1, ndab-1
    !jmn(:,n+1,id) = -0.25d0 * 2.0d0 * (Cmn1_r(:,n+1,id) * Cmn2_r(:,n+1,id) + Cmn1_i(:,n+1,id) * Cmn2_i(:,n+1,id))
    jmn(:,n+1,id) = -0.25d0 * (Cmn1_r(:,n+1,id) * Cmn2_r(:,n+1,id) + Cmn1_i(:,n+1,id) * Cmn2_i(:,n+1,id))
    imn(:,n+1,id) = jmn(:,n+1,id) * ji_conversion(n+1)
  enddo
  write(6,'("Lev ",i3," net enstrophy transfer ",es20.13)'), id, sum(jmn(1,:,id)) + 2.0d0 * sum(jmn(2:mdab_s,:,id))
  write(6,'("Lev ",i3," net energy    transfer ",es20.13)'), id, sum(imn(1,:,id)) + 2.0d0 * sum(imn(2:mdab_s,:,id))
enddo
! Jmn(m,n) -> Mean in altitude
write(6,'(a)') "2D enstrophy and energy interaction terms, time average ..."
do n = 1, ndab-1
  do m = 0, mdab_s-1
    !jmn_zm(m+1,n+1) = -0.25d0 * 2.0d0 * sum(Cmn1_r(m+1,n+1,:) * Cmn2_r(m+1,n+1,:) + Cmn1_i(m+1,n+1,:) * Cmn2_i(m+1,n+1,:)) / nalt
    jmn_zm(m+1,n+1) = -0.25d0 * sum(Cmn1_r(m+1,n+1,:) * Cmn2_r(m+1,n+1,:) + Cmn1_i(m+1,n+1,:) * Cmn2_i(m+1,n+1,:)) / nalt
  enddo
  imn_zm(:,n+1) = jmn_zm(:,n+1) * ji_conversion(n+1)
enddo
!*************************************************************************
! 1D enstrophy and energy interaction terms (BS83 Eqs 15-16 / BS13 Eqs 8-9)
! Factor of two is required for m>0 to account for +/- m values
!*************************************************************************
write(6,'(a)') "1D enstrophy and energy interaction terms for each day ..."
do id = 1, nalt
  ! n = 0 term
  jn(1,id) = jmn(1,1,id)
  ! n > 0 terms
  do n = 1, ndab-1
    jn(n+1,id) = jmn(1,n+1,id) + 2.0d0 * sum(jmn(2:n+1,n+1,id))
  enddo
  ! Convert to energy interaction term
  in(:,id) = jn(:,id) * ji_conversion 
enddo

write(6,'(a)') "1D enstrophy and energy interaction terms, time average ..."
! n = 0 term
jn_zm(1) = jmn_zm(1,1)
! n > 0 terms
do n = 1, ndab-1
  jn_zm(n+1) = jmn_zm(1,n+1) + 2.0d0 * sum(jmn_zm(2:n+1,n+1))
enddo
! Convert to energy interaction term
in_zm(:) = jn_zm(:) * ji_conversion

!====================================================================
! 5. Calculate energy and enstrophy spectral fluxes (BS13 Eqs 11-12)
!====================================================================
write(6,'(a)') "1D enstrophy and energy spectral fluxes for each day ..."
do id = 1, nalt
do im = 1, mdab_s
  do n = 1, ndab-1 ! n > 0 terms (n = 0 term is zero)
    fmn(im,n+1,id) = -sum(imn(im,1:n,id)) ! Energy flux
    hmn(im,n+1,id) = -sum(jmn(im,1:n,id)) ! Enstrophy flux
  enddo
enddo
enddo

write(6,'(a)') "1D enstrophy and energy spectral fluxes for each day ..."
do id = 1, nalt
  do n = 1, ndab-1 ! n > 0 terms (n = 0 term is zero)
    fn(n+1,id) = -sum(in(1:n,id)) ! Energy flux
    hn(n+1,id) = -sum(jn(1:n,id)) ! Enstrophy flux
  enddo
enddo

write(6,'(a)') "1D enstrophy and energy spectral fluxes, time average ..."
 do n = 1, ndab-1 ! n > 0 terms (n = 0 term is zero)
   fn_zm(n+1) = -sum(in_zm(1:n)) ! Energy flux
   hn_zm(n+1) = -sum(jn_zm(1:n)) ! Enstrophy flux
enddo

end subroutine FLUXES
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES X
!#							       ES:
!#							 Energy Spectra
!###################################################################################################################
!###################################################################################################################
! ---------- ES
subroutine ES(mdab_s,ndab,nalt,radius,cmn_r,cmn_i,bmn_r,bmn_i,Emn,Enmo,Ene,Enmo_zm,Ene_zm,ER,EZ,ET,En_rot,En_div)
    implicit none
    integer :: mdab_s,ndab,nalt,id
    real*8 :: radius,n=0,m=0,ER,EZ,ET
    double precision, dimension (:), allocatable :: ji_conversion,En_zm,Enmo_zm,Ene_zm
    double precision, dimension(:,:), allocatable :: Emn_zm,En,Enmo,Ene,En_rot,En_div
    double precision, dimension(:,:,:), allocatable :: cmn_r,cmn_i,bmn_r,bmn_i,Emn,Emn_rot,Emn_div

!allocate(Emn(mdab_s,ndab,nalt)) ; Emn    = 0.0d0
allocate(Emn_zm(mdab_s,ndab)) ; Emn_zm    = 0.0d0
allocate(En(ndab,nalt)) ; En    = 0.0d0
allocate(En_zm(ndab)) ; En_zm    = 0.0d0
!div-rot
allocate(Emn_rot(mdab_s,ndab,nalt)) ; Emn_rot    = 0.0d0
allocate(Emn_div(mdab_s,ndab,nalt)) ; Emn_div    = 0.0d0
allocate(En_rot(ndab,nalt)) ; En_rot    = 0.0d0
allocate(En_div(ndab,nalt)) ; En_div    = 0.0d0
!allocate(Enmo(ndab,nalt)) ; Enmo    = 0.0d0
!allocate(Ene(ndab-1,nalt)) ; Ene    = 0.0d0
!allocate(Enmo_zm(ndab)) ; Enmo_zm    = 0.0d0
!allocate(Ene_zm(ndab-1)) ; Ene_zm    = 0.0d0
!=====================================================
! 1. Conversion factor
!=====================================================
allocate(ji_conversion(ndab)) ; ji_conversion = 0.0d0	
do n = 1, ndab-1 
  ji_conversion(n+1) = radius * radius / (n * (n + 1.0d0))
enddo
!=====================================================
! 2. Calculate energy spectra
!=====================================================
! (m,n,altitude)
write(6,'(a)') "2D enstrophy and energy interaction terms for each day ..."
do id = 1, nalt
  do n = 1, ndab-1
    Emn(:,n+1,id) = 0.25d0*ji_conversion(n+1)*(cmn_r(:,n+1,id)*cmn_r(:,n+1,id) + cmn_i(:,n+1,id)*cmn_i(:,n+1,id) + bmn_r(:,n+1,id)*bmn_r(:,n+1,id) + bmn_i(:,n+1,id)*bmn_i(:,n+1,id))
    Emn_rot(:,n+1,id) = 0.25d0*ji_conversion(n+1)*(cmn_r(:,n+1,id)*cmn_r(:,n+1,id) + cmn_i(:,n+1,id)*cmn_i(:,n+1,id))
    Emn_div(:,n+1,id) = 0.25d0*ji_conversion(n+1)*(bmn_r(:,n+1,id)*bmn_r(:,n+1,id) + bmn_i(:,n+1,id)*bmn_i(:,n+1,id))
  enddo
enddo
! Emn(m,n) -> Mean in altitude
write(6,'(a)') "2D enstrophy and energy interaction terms, time average ..."
do n = 1, ndab-1
  do m = 0, mdab_s-1
    Emn_zm(m+1,n+1) = 0.25d0*ji_conversion(n+1)*sum(cmn_r(m+1,n+1,:)*cmn_r(m+1,n+1,:) + cmn_i(m+1,n+1,:)*cmn_i(m+1,n+1,:) + bmn_r(m+1,n+1,:)*bmn_r(m+1,n+1,:) + bmn_i(m+1,n+1,:)*bmn_i(m+1,n+1,:)) / nalt
  enddo
enddo 
!*************************************************************************
! Factor of two is required for m>0 to account for +/- m values
!*************************************************************************
! (m,n,altitude)
write(6,'(a)') "1D enstrophy and energy interaction terms for each day ..."
do id = 1, nalt
  ! n = 0 term
  En(1,id) = Emn(1,1,id)
  Enmo(1,id) = Emn(1,1,id)
  En_rot(1,id) = Emn_rot(1,1,id)
  En_div(1,id) = Emn_div(1,1,id)
  ! n > 0 terms
  do n = 1, ndab-1
    En(n+1,id) = Emn(1,n+1,id) + 2.0d0 * sum(Emn(2:n+1,n+1,id))
    En_rot(n+1,id) = Emn_rot(1,n+1,id) + 2.0d0 * sum(Emn_rot(2:n+1,n+1,id))
    En_div(n+1,id) = Emn_div(1,n+1,id) + 2.0d0 * sum(Emn_div(2:n+1,n+1,id))
    !En_rot(n+1,id) = 2.0d0 * sum(Emn_rot(2:n+1,n+1,id)) ! here I look at all modes except m=0 for the rot component
    !En_div(n+1,id) = 2.0d0 * sum(Emn_div(2:n+1,n+1,id)) ! here I look at all modes except m=0 for the div component
    Enmo(n+1,id) = Emn(1,n+1,id)
    Ene(n+1,id) = 2.0d0 * sum(Emn(2:n+1,n+1,id))
  enddo
enddo
! Emn(n) -> Mean in altitude
write(6,'(a)') "1D enstrophy and energy interaction terms, time average ..."
! n = 0 term
En_zm(1) = Emn_zm(1,1)
Enmo_zm(1) = Emn_zm(1,1)
! n > 0 terms
do n = 1, ndab-1
  En_zm(n+1) = Emn_zm(1,n+1) + 2.0d0 * sum(Emn_zm(2:n+1,n+1))
  Enmo_zm(n+1) = Emn_zm(1,n+1)
  Ene_zm(n+1) = 2.0d0 * sum(Emn_zm(2:n+1,n+1))
enddo
!=====================================================
! 3. Integral de l energie
!=====================================================
ET = sum(Enmo_zm + Ene_zm)
EZ = sum(Enmo_zm)
ER = sum(Ene_zm)
end subroutine ES
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES XI
!#							 epsilonforcage:
!#							 epsilon forcage
!###################################################################################################################
!###################################################################################################################
! ---------- epsilonf
subroutine epsilonforcage(mdab_s,ndab,nalt,radius,vort_r,vort_i,div_r,div_i,vortf_r,vortf_i,divf_r,divf_i,epsilonf)
    implicit none
    integer :: mdab_s,ndab,nalt,id
    real*8 :: radius,n=0,m=0
    double precision, dimension (:), allocatable :: ji_conversion,epsilonf
    double precision, dimension(:,:,:), allocatable :: epsilmn,vort_r,vort_i,div_r,div_i
    double precision, dimension(:,:), allocatable :: epsiln,epsilne,epsilnmo,vortf_r,vortf_i,divf_r,divf_i

allocate(epsilmn(mdab_s,ndab,nalt)) ; epsilmn	= 0.0d0
allocate(epsilnmo(ndab,nalt)) ; epsilnmo	= 0.0d0
allocate(epsilne(ndab-1,nalt)) ; epsilne	= 0.0d0
allocate(epsiln(ndab,nalt)) ; epsiln		= 0.0d0
!allocate(epsilonf(nalt)) ; epsilonf		= 0.0d0
!=====================================================
! 1. Conversion factor
!=====================================================
allocate(ji_conversion(ndab)) ; ji_conversion = 0.0d0	
do n = 1, ndab-1 
  ji_conversion(n+1) = radius * radius / (n * (n + 1.0d0))
enddo
!=====================================================
! 2. Calculate u.f
!=====================================================
!***********************************************************
! u.f = a*a/(n(n+1)) * R{conj[rot(u)]*rot(f) + conj[div(u)]*div(f)}
! or rot(u) = vort_r + i*vort_i et rot(f) = vortf_r + i*vortf_i,
! on a conj[rot(u)]*rot(f) = 0.5*[vort_r*vortf_r + vort_i*vortf_i]
! de meme pour conj[div(u)]*div(f)
! A factor *2 account on m/=0 when m=0 coeff is 1/2
! (m,n,altitude)
!***********************************************************
write(6,'(a)') "2D < u.f > for each day ..."
do id = 1, nalt
  do n = 1, ndab-1
    epsilmn(:,n+1,id) = 0.25d0*ji_conversion(n+1)*( vort_r(:,n+1,id)*vortf_r(:,n+1) + vort_i(:,n+1,id)*vortf_i(:,n+1) + div_r(:,n+1,id)*divf_r(:,n+1) + div_i(:,n+1,id)*divf_i(:,n+1) )
  enddo
enddo
!*************************************************************************
! Factor of two is required for m>0 to account for +/- m values
!*************************************************************************
! (m,n,altitude)
write(6,'(a)') "1D enstrophy and energy interaction terms for each day ..."
do id = 1, nalt
  ! n = 0 term
  epsilne(1,id) = epsilmn(1,1,id)
  epsilnmo(1,id) = epsilmn(1,1,id)
  ! n > 0 terms
  do n = 1, ndab-1
    epsiln(n+1,id) = epsilmn(1,n+1,id) + 2.0d0 * sum(epsilmn(2:n+1,n+1,id))
    epsilnmo(n+1,id) = epsilmn(1,n+1,id)
    epsilne(n+1,id) = 2.0d0 * sum(epsilmn(2:n+1,n+1,id))
  enddo
   epsilonf(id) = sum(epsiln(:,id))!*1./2.
enddo
print*,'epsilonf(iz) = ', epsilonf

end subroutine epsilonforcage
!###################################################################################################################
!###################################################################################################################
!#							 SUBROUTINES XII
!#							     sorting:
!#					               profil monotonization
!###################################################################################################################
!###################################################################################################################
!******** Vorticite potentielle ***********************************
!                  &
!	     Monotonization
!******************************************************************
! Call a sorting function that: sorting(ndata,VtoSort,Vsorted,typeS,Idx)
! if typeS == 1 ascendent sorting & typeS == 2 descendent sorting
! The ascendent choice depend on the data order. The output are,
! ndata: is the length of VtoSort(ndata)
! Vtosort: is a single vector
! Vsorted: is a vector once Vtosort sorted it has the same length
! Idx: are index unsorted - index sorted --> typical sorted scale
!******************************************************************
subroutine sorting(nlat,VtoSort,Vsorted,typeS,L_Idx)
    double precision, dimension (:), allocatable :: VtoSort,Vsorted
    double precision :: t,tt
    integer, dimension (nlat) :: Idx,IdxSorted,L_Idx
    integer :: nlat,i,j,typeS

! if typeS == 1 ascendent sorting & typeS == 2 descendent sorting
!allocate(aa(nlat),dd(nlat)) ; 
Vsorted = VtoSort 
! generate indices
do i  =1,nlat
Idx(i)=i
IdxSorted(i)=i
enddo
! arranging in increasing order
if(typeS == 1)then
do i=1,nlat
   do j=i,nlat
      if(Vsorted(i)>Vsorted(j))then
	! sorted vector
	t=Vsorted(j)
	Vsorted(j)=Vsorted(i)
	Vsorted(i)=t
	! sorted indices
	tt=IdxSorted(j)
	IdxSorted(j)=IdxSorted(i)
	IdxSorted(i)=tt
      endif
   enddo
enddo
! arranging in decreasing order
elseif(typeS == 2)then
do i=1,nlat
   do j=i,nlat
      if(Vsorted(i)<Vsorted(j))then
	! sorted vector
	t=Vsorted(j)
	Vsorted(j)=Vsorted(i)
	Vsorted(i)=t
	! sorted indices
	tt=IdxSorted(j)
	IdxSorted(j)=IdxSorted(i)
	IdxSorted(i)=tt
      endif
   enddo
enddo
endif
L_Idx = Idx-IdxSorted
end subroutine sorting
!###################################################################################################################
!######################################################################################################################
!###################################################################################################################
END PROGRAM spectra_analysis
