PROGRAM spectra_analysis


use netcdf

implicit none
!**********************************************************
integer :: iarg,narg,mt,istep,getfield,estep,dimid_out(3),varid_out(4),latid,altid,colatid,lonid,tid,idinfos,ninfos,idninfos,Urmsid,hnid,fnid,dstepsid,epsilonfid,ERid,ETid,EZid,hn_zmid,fn_zmid,Enid,En_divid,En_rotid,En_zmid,fnrotid,hnrotid,fndiv1id,hndiv1id,fndiv2id,hndiv2id
character (len=100), dimension(:), allocatable :: arg
character (len=100) :: tmp_char
logical :: is_file,NotFirst_reading=.false.,First_reading=.true.
!----------------------------------------------------------------------
character (len=100) :: file_netcdf
integer :: idfile,idtime,idnlon,idnlat,idnalt,idnu,idnv
integer :: nlon=0,nlat=0,nalt,ntime,iddsteps,nv
integer :: idlat,idlon,idalt,idu,idv,latlen,ierror
integer :: ndims, nvars, natts, unlimdimid, nformat
double precision, dimension (:), allocatable :: lat_DYN,lon_DYN,alt_DYN,dsteps_DYN
double precision, dimension(:,:,:), allocatable :: u_DYN
character(len=*), parameter :: latname = 'lat'
!###################################################################################################################
!######################################################################################################################
!###################################################################################################################
!**********
!Initialisation
mt=1
istep=1
estep=1
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
  if ( arg(iarg) == '-mt' .or. arg(iarg) == '-istp' .or. arg(iarg) == '-estp' .or. arg(iarg) == '-field') then
    select case(arg(iarg))
    case('-mt')
      read(arg(iarg+1),'(i10)' ) mt !extent of temporal length
    case('-istp')
      read(arg(iarg+1),'(i10)' ) istep !name of time axis
    case('-field')
      read(arg(iarg+1),'(i10)' ) getfield !name of time axis
    case('-estp')
      read(arg(iarg+1),'(i10)' ) estep !name of time axis
    end select
    iarg = iarg + 2
  elseif (arg(iarg) == '-h' .or. arg(iarg) == '--help') then
      print*,'Usage'
      print*,'spectra_analysis netcdfFile [option]'
      print*,'[-h or --help]	: brief help'
      print*,'[-mt int]	: final temporal file created (default: 0)'
      print*,'[-istp int]	: initial temporal step to create files (default: 0)'
      print*,'[-field int]	: select the field to decompose: 1.zonalV 2.meridV 3.totalV 4.stream 5.Vpoten (default: 3)'
      print*,'[-estp int]	: temporal step to create files (default: 1)'
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
print*,'create files using steps in time t =',istep,':',estep,':',mt
print*,file_netcdf
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
 call check(  nf90_open(trim(file_netcdf),NF90_NOWRITE,idfile))
!ierror=-99999
!ierror=nf90_inquire(idfile, ndims, nvars, natts, unlimdimid, nformat)

!allocate(u_DYN(720,360,1))
!ierror = nf90_inq_varid(idfile,trim('u'),idu)
!ierror = nf90_get_var(idfile,idu,u_DYN,(/1,1,1,1/),(/720,360,1,1/))

 !call check( nf90_open(trim(file_netcdf),NF90_NOWRITE,idfile))

 !call check( nf90_inquire(idfile, ndims, nvars, natts, unlimdimid, nformat))

 !call check( nf90_inq_dimid(idfile,trim('u'),idtime) )
 !call check( nf90_inquire_dimension(idfile,idtime,tmp_char,ntime))

 !call check( nf90_inq_dimid(idfile,latname,idnlat) )
 !call check( nf90_inq_dimid(idfile,latname,idnlon) )
 !call check( nf90_inquire_dimension(idfile, idnlat, len=latlen) )

 !call check( nf90_inquire_dimension(idfile,idnu,tmp_char,nu))

 !call check( nf90_inq_varid(idfile,trim('u'),idu))


print*,'rien'
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
END PROGRAM spectra_analysis
