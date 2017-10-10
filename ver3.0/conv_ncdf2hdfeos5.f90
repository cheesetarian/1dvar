      module conv_ncdf2hdfeos5

!------ contains routine to read in netCDF output from 1DVAR retrieval,
!------  and output to a HDF-EOS5 format file
!------  D Duncan, CSU, 7/26/16, last updated 9/28/16

      use netcdf
      use define_csu1dvar

      implicit none
      include 'hdfeos5.inc'

      contains

      subroutine conv(infile,outfile)

      logical :: fexists
      character(len=256) :: infile,outfile

      integer  he5_swopen,he5_swcreate,he5_swdefdim,he5_swdefgfld,he5_swdefdfld,he5_swwrfld,he5_swdetach,he5_swclose,he5_swwrdmeta
      external he5_swopen,he5_swcreate,he5_swdefdim,he5_swdefgfld,he5_swdefdfld,he5_swwrfld,he5_swdetach,he5_swclose,he5_swwrdmeta

      !integer, parameter :: CMPRS = 2 !level of compression
      integer :: nscans,npix,nch,nvar,nz=16
      integer :: time_varid!, sclat_varid, sclon_varid, scalt_varid
      integer :: lat_varid, lon_varid!, press_varid, tbsim_varid, tbobs_varid
      !integer :: tbdif_varid, sal_varid, slp_varid, 
      integer :: tpw_varid, wsp_varid, lwp_varid, wdir_varid, sst_varid, reysst_varid
      integer :: chi_varid, sun_varid, land_varid, qual_varid
      integer :: ncid, post_varid, tai93_varid
      integer :: pix_dimid, scan_dimid, chan_dimid
      integer :: time_dim, scan_dim, chan_dim
      integer :: post_dim,var_dim
      !character(len=40)  :: sdate, sprior
      !character(len=40)  :: sfreq
      character(len=10)  :: sgran,siter,sthresh,chi1,chith,version,smiss
      character(len=10)  :: satellite,sensor,rtmodel,emodel,convolve
      character(len=20)  :: ppversion
      character(len=100) :: cfile, syfile

      integer*2,   allocatable :: scantime(:,:)
      real*8,    allocatable :: tai93time(:)
      integer*1, allocatable :: qflag(:,:)
      integer*1, allocatable :: landout(:,:)
      integer*2, allocatable :: sunout(:,:)

      real, allocatable :: latout(:,:)
      real, allocatable :: lonout(:,:)
      real, allocatable :: tpwout(:,:)
      real, allocatable :: wspout(:,:)
      real, allocatable :: lwpout(:,:)
      !real, allocatable :: sstout(:,:)
      real, allocatable :: chiout(:,:)

      real, allocatable :: reyout(:,:)
      !real, allocatable :: slpout(:,:)
      real, allocatable :: post(:,:,:)
      !real, allocatable :: tbsim(:,:,:)
      !real, allocatable :: tbobs(:,:,:)
      !real, allocatable :: tbdif(:,:,:)

      INTEGER(4), PARAMETER   :: FLOAT64 = HE5T_NATIVE_DOUBLE!real8/double
      INTEGER(4), PARAMETER   :: FLOAT32 = HE5T_NATIVE_FLOAT !real
      INTEGER(4), PARAMETER   :: INT32   = HE5T_NATIVE_INT   !int4/long
      INTEGER(4), PARAMETER   :: INT16   = HE5T_NATIVE_INT16 !int2
      INTEGER(4), PARAMETER   :: INT8    = HE5T_NATIVE_INT8  !byte
      character(1000) :: fieldlist
      integer,parameter :: numparams = 15 ! # output parameters
      !integer(8) :: count
      integer(8) :: start1, stride1, edge1
      integer(8) :: start(2), stride(2), edge(2)
      !integer(8) :: start3(3), stride3(3), edge3(3)
      integer :: fieldtype(numparams)!,fielddim(numparams),fieldarray(numparams)
      integer :: swid, swfid, status

! read in input/output filenames, then get npix/nscans from input file
        !print*,trim(infile)
        !print*,trim(outfile)
      inquire(file=infile,exist = fexists)
      if (fexists .eq. .false.) then
        write(6,*) 'Input file ',trim(infile),' not found!'
      endif

!--- start opening and reading netcdf file
      call check(nf90_open(trim(infile), nf90_NoWrite, ncid))

!--- Read nscan dimension value
      call check(nf90_inq_dimid(ncid, "nscan",  scan_dimid))
      call check(nf90_inquire_dimension(ncid, scan_dimid, len = nscans))
      !write(*,*) 'nscans = ',nscans
      call check(nf90_inq_dimid(ncid, "npix",  pix_dimid))
      call check(nf90_inquire_dimension(ncid, pix_dimid, len = npix))
      !write(*,*) 'npix = ',npix
      call check(nf90_inq_dimid(ncid, "nchan",  chan_dimid))
      call check(nf90_inquire_dimension(ncid, chan_dimid, len = nch))
      !write(*,*) 'nchan = ',nch

      chan_dim  = nch
      time_dim  = 6
      var_dim   = nvar
      post_dim  = 4 !nvar+1

      !---Allocate arrays---

      allocate (scantime(nscans,6))
      allocate (tai93time(nscans))
      allocate (qflag(nscans,npix))
      allocate (landout(nscans,npix))
      allocate (sunout(nscans,npix))

      allocate (latout(nscans,npix))
      allocate (lonout(nscans,npix))
      allocate (tpwout(nscans,npix))
      allocate (wspout(nscans,npix))
      allocate (lwpout(nscans,npix))
      allocate (reyout(nscans,npix))
      !allocate (sstout(nscans,npix))
      allocate (chiout(nscans,npix))

      !allocate (slpout(nscans,npix))
      allocate (post(nscans,npix,post_dim))
      !allocate (tbsim(nscans,npix,nch))
      !allocate (tbobs(nscans,npix,nch))
      !allocate (tbdif(nscans,npix,nch))
      !allocate (wvp_prof(nscans,npix,nz))
      !allocate (tmp_prof(nscans,npix,nz))
      !allocate (lwp_prof(nscans,npix,nz))

!---  with vars allocated, read in global attributes
       call check(nf90_get_att(ncid,nf90_global,"Satellite",satellite))
       call check(nf90_get_att(ncid,nf90_global,"Sensor",sensor))
       call check(nf90_get_att(ncid,nf90_global,"Version",version))
       call check(nf90_get_att(ncid,nf90_global,"RT_Model",rtmodel))
       call check(nf90_get_att(ncid,nf90_global,"Emis_Model",emodel))
       call check(nf90_get_att(ncid,nf90_global,"Preprocessor_Version",ppversion))
       call check(nf90_get_att(ncid,nf90_global,"Convolution_Resolution",convolve))
       call check(nf90_get_att(ncid,nf90_global,"Calibration_File",cfile))
       call check(nf90_get_att(ncid,nf90_global,"SY_File",syfile))
       !call check(nf90_get_att(ncid,nf90_global,"Missing_Value",smiss))
       call check(nf90_get_att(ncid,nf90_global, &
         "ChiSquare_Threshold",chith))
       call check(nf90_get_att(ncid,nf90_global, &
         "ChiSquare_Threshold_High_Quality",chi1))
        !print*,satellite,sensor,version,rtmodel,emodel
        !print*,cfile,syfile,smiss,chith,chi1
        
!---  get variable IDs from input file
      call check(nf90_inq_varid(ncid,"latitude",lat_varid))
      call check(nf90_inq_varid(ncid,"longitude",lon_varid))
      call check(nf90_inq_varid(ncid,"sc_time",time_varid))
      call check(nf90_inq_varid(ncid,"time",tai93_varid))
      call check(nf90_inq_varid(ncid,"atmosphere_mass_content_of_water_vapor",tpw_varid))
      call check(nf90_inq_varid(ncid,"wind_speed",wsp_varid))
      call check(nf90_inq_varid(ncid,"atmosphere_mass_content_of_cloud_liquid_water",lwp_varid))
      call check(nf90_inq_varid(ncid,"chi_squared",chi_varid))
      !call check(nf90_inq_varid(ncid,"sea_surface_subskin_temperature",sst_varid))
      call check(nf90_inq_varid(ncid,"sunglint_angle",sun_varid))
      call check(nf90_inq_varid(ncid,"land_area_fraction",land_varid))
      call check(nf90_inq_varid(ncid,"quality_flag",qual_varid))
      call check(nf90_inq_varid(ncid,"posterior_error",post_varid))
      !call check(nf90_inq_varid(ncid,"tb_difference",tbdif_varid))
      call check(nf90_inq_varid(ncid,"reynolds_sst",reysst_varid))
!---  read variables
      call check(nf90_get_var(ncid,lat_varid,latout))
      call check(nf90_get_var(ncid,lon_varid,lonout))
      call check(nf90_get_var(ncid,time_varid,scantime))
      call check(nf90_get_var(ncid,tai93_varid,tai93time))
      call check(nf90_get_var(ncid,tpw_varid,tpwout))
      call check(nf90_get_var(ncid,wsp_varid,wspout))
      call check(nf90_get_var(ncid,lwp_varid,lwpout))
      call check(nf90_get_var(ncid,chi_varid,chiout))
      !call check(nf90_get_var(ncid,sst_varid,sstout))
      call check(nf90_get_var(ncid,sun_varid,sunout))
      call check(nf90_get_var(ncid,land_varid,landout))
      call check(nf90_get_var(ncid,qual_varid,qflag))
      call check(nf90_get_var(ncid,post_varid,post))
      !call check(nf90_get_var(ncid,tbdif_varid,tbdif))
      call check(nf90_get_var(ncid,reysst_varid,reyout))

!---  close input file
      call check(nf90_close(ncid))

!---- prepare fields needed for HDF-EOS5 output file
      fieldlist='Latitude,Longitude,Time,TimeHR,'//&
        'TotalPrecipitableWater,WindSpeed,LiquidWaterPath,'//&
        'ReynoldsSST,ChiSquared,QualityFlag,'//&
        'ErrorTPW,ErrorWind,ErrorLWP,'//& !ErrorSST,'//&
        'SunGlintAngle,LandPercentage' !,TbDifference'

      fieldtype(1)  = FLOAT32    ! Latitude  (real)
      fieldtype(2)  = FLOAT32    ! Longitude (real)
      fieldtype(3)  = FLOAT64    ! Time, TAI93 (double)
      fieldtype(4)  = INT16      ! Time, ymdhms (int2)
      fieldtype(5)  = FLOAT32    ! TPW (real)
      fieldtype(6)  = FLOAT32    ! Wind (real)
      fieldtype(7)  = FLOAT32    ! LWP (real)
      fieldtype(8)  = FLOAT32    ! ReynoldsSST (real)
      fieldtype(9)  = FLOAT32    ! ChiSq (real)
      fieldtype(10) = INT8       ! Quality Flag (byte)
      fieldtype(11) = FLOAT32    ! Posterior error, TPW (real)
      fieldtype(12) = FLOAT32    ! Posterior error, Wind (real)
      fieldtype(13) = FLOAT32    ! Posterior error, LWP (real)
      fieldtype(14) = INT16      ! Sunglint Angle (int2)
      fieldtype(15) = INT8       ! Land Percent in FOV (byte)
      !fieldtype(14) = FLOAT32    ! Posterior error, SST (real)
      !fieldtype(17) = FLOAT32    ! Tb Residual (real)

!----- open HDF-EOS5 file
      swfid = he5_swopen(trim(outfile),he5f_acc_trunc)
!----- create swath (as opposed to point or level)
      swid  = he5_swcreate(swfid,'AMSR2_Level2_Ocean_Suite')

!----- define dimensions of output arrays
      status = he5_swdefdim(swid,'nscans',nscans)  
      status = he5_swdefdim(swid,'npix',npix)  
      !status = he5_swdefdim(swid,'nchannels',nch)  
      status = he5_swdefdim(swid,'ntime',time_dim)  
!----- create geolocation fields  
      status = he5_swdefgfld(swid,'Latitude','nscans,npix','',&
        fieldtype(1),he5_hdfe_nomerge)
      status = he5_swdefgfld(swid,'Longitude','nscans,npix','',&
        fieldtype(2),he5_hdfe_nomerge)
      status = he5_swdefgfld(swid,'Time','nscans','',&
        fieldtype(3),he5_hdfe_nomerge)
!----- define data fields
      status = he5_swdefdfld(swid,'TimeHR','nscans','',&
        fieldtype(4),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'TotalPrecipitableWater','nscans,npix','',&
        fieldtype(5),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'WindSpeed','nscans,npix','',&
        fieldtype(6),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'LiquidWaterPath','nscans,npix','',&
        fieldtype(7),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'ReynoldsSST','nscans,npix','',&
        fieldtype(8),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'ChiSquared','nscans,npix','',&
        fieldtype(9),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'QualityFlag','nscans,npix','',&
        fieldtype(10),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'ErrorTPW','nscans,npix','',&
        fieldtype(11),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'ErrorWind','nscans,npix','',&
        fieldtype(12),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'ErrorLWP','nscans,npix','',&
        fieldtype(13),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'SunGlintAngle','nscans,npix','',&
        fieldtype(14),he5_hdfe_nomerge)
      status = he5_swdefdfld(swid,'LandPercentage','nscans,npix','',&
        fieldtype(15),he5_hdfe_nomerge)
      !status = he5_swdefdfld(swid,'ErrorSST','nscans,npix','',&
      !  fieldtype(14),he5_hdfe_nomerge)
      !status = he5_swdefdfld(swid,'TbDifference','nscans,npix,nchannels','',&
      !  fieldtype(17),he5_hdfe_nomerge)
       
!----- write data fields,
      !  all output arrays of size [nscans,npix]
      start(1) = 0
      start(2) = 0
      stride(1) = 1
      stride(2) = 1
      edge(1) = nscans
      edge(2) = npix 

      status = he5_swwrfld(swid,'Latitude',start,stride,edge,latout)
      call chkstatus(status,'swwrfld lat')
      status = he5_swwrfld(swid,'Longitude',start,stride,edge,lonout)
      call chkstatus(status,'swwrfld lon')
      status = he5_swwrfld(swid,'TotalPrecipitableWater',start,stride,edge,tpwout)
      call chkstatus(status,'swwrfld tpw')
      status = he5_swwrfld(swid,'WindSpeed',start,stride,edge,wspout)
      call chkstatus(status,'swwrfld wind')
      status = he5_swwrfld(swid,'LiquidWaterPath',start,stride,edge,lwpout)
      call chkstatus(status,'swwrfld lwp')
      status = he5_swwrfld(swid,'ReynoldsSST',start,stride,edge,reyout)
      call chkstatus(status,'swwrfld reysst')
      !status = he5_swwrfld(swid,'SeaSurfaceTemperature',start,stride,edge,sstout)
      !call chkstatus(status,'swwrfld sst')
      status = he5_swwrfld(swid,'ChiSquared',start,stride,edge,chiout)
      call chkstatus(status,'swwrfld chisq')
      status = he5_swwrfld(swid,'QualityFlag',start,stride,edge,qflag)
      call chkstatus(status,'swwrfld qual')
      status = he5_swwrfld(swid,'ErrorTPW',start,stride,edge,post(:,:,1))
      call chkstatus(status,'swwrfld errortpw')
      status = he5_swwrfld(swid,'ErrorWind',start,stride,edge,post(:,:,2))
      call chkstatus(status,'swwrfld errorwind')
      status = he5_swwrfld(swid,'ErrorLWP',start,stride,edge,post(:,:,3))
      call chkstatus(status,'swwrfld errorlwp')
      !status = he5_swwrfld(swid,'ErrorSST',start,stride,edge,post(:,:,4))
      !call chkstatus(status,'swwrfld errorsst')
      status = he5_swwrfld(swid,'SunGlintAngle',start,stride,edge,sunout)
      call chkstatus(status,'swwrfld sunglint')
      status = he5_swwrfld(swid,'LandPercentage',start,stride,edge,landout)
      call chkstatus(status,'swwrfld landpct')

      edge(2) = 6 ! year/mon/day/hour/min/sec
      status = he5_swwrfld(swid,'TimeHR',start,stride,edge,scantime)
      call chkstatus(status,'swwrfld timehr')

      start1 = 0
      stride1= 1
      edge1 = nscans
      status = he5_swwrfld(swid,'Time',start1,stride1,edge1,tai93time)
      call chkstatus(status,'swwrfld time')

      ! doesnt actually work with 3 dimensions? if wanting to output Tb
      !  residuals, split up into 2d arrays, one for each channel...
      !start3(:) = 0
      !stride3(:)= 1
      !edge3(1) = nscans
      !edge3(2) = npix
      !edge3(3) = nch
      !status = he5_swwrfld(swid,'TbDifference',start3,stride3,edge3,tbdif)
      !call chkstatus(status,'swwrfld tbdif')
        
      ! detach from swath interface
      status = he5_swdetach(swid)
      call chkstatus(status,'swdetach')

      ! close hdfeos5 file
      status = he5_swclose(swfid)
      call chkstatus(status,'swclose')


      end subroutine conv

!______________________________________________________________________

      subroutine check(status)
      integer, intent ( in) :: status
      if (status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if
      end subroutine check  

      subroutine chkStatus(status, prtMsg)
!---  this routine simply checks the status code and writes an Error
!---  message if status < 0,  These are fatal errors so process stops.      
      integer           :: status
      character(len=*)  :: prtMsg
!      write(6,*) status, prtMsg
      if(status .lt. 0) then
           write(*,*) 'FATAL error from: ', trim(prtMsg)
           stop
      endif
      return
      end subroutine chkStatus

      end module conv_ncdf2hdfeos5
