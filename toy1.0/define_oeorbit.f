      MODULE define_oeorbit

       USE CRTM_MODULE
       implicit none
       save

!--- these are shared variables, used by more than one 
!---  subroutine and therefore commonly defined here

!......................................................................
! DEFINE EVERYTHING FOR THE TOY 1DVAR HERE:
       character(10) :: sensor    = 'GMI'    ! GMI, AMSR2 supported
       integer,parameter :: nch   = 13       ! #channels (certain combos only!)
       real          :: latitude  =  -55.13  ! deg N
       real          :: longitude = -165.70  ! deg E

       integer :: year = 2014
       integer :: month= 4
       integer :: day  = 1
       real    :: hour = 10
       real    :: minu = 10
       real    :: seco = 26
       character(len=8) :: alldate = '20140401' !make these match!
       
       logical :: useanalysis     = .true. ! if true, reads in analysis data
       logical :: merra           = .false. ! merra first guess (not EC)
        integer,parameter :: maxchans = 15
       real :: Tb_in(maxchans) = 
     >          (/161.46,  89.52, 183.48, 116.67, 204.44,
     >            -999.9, 211.24, 155.98, 252.22, 224.02,
     >            266.23, 264.97, -999.9, 250.44, 260.90/)
       real :: EIAs(2) = (/53.0,49.31/)
       real :: azz     = 117.2 ! satellite azimuth angle
       ! add more here later? overrides for SST, cloud height, wind
       ! direction, first guesses, which vars are retrieved, etc.
       real :: wdir = 123.2 ! wind direction (deg clockwise from N=0)
       ! emissivity model choice:
       logical :: run_rss = .false. ! if false, FASTEM6 runs
       real :: psst

!......................................................................


       character(10) :: code_version = 'toy1.0' ! version of code

!______________________________________________________________________
       logical :: dontoutput     = .false. ! if true, code stops after OE step
       logical :: runseaice      = .false. ! if true, force SST=271 for pix w/ reynolds seaice
       logical :: includeicecloud= .false. ! should match IWP being retrieved or not!
       !logical :: outrt,run_rss,fastem4 ! set by script calling code (see oeorbit.f)
       logical :: fastem4

!       VARS TO BE MODIFIED FOR TESTING/RUNNING 
       !integer  :: scstart=660, pxstart=1, scend=660, pxend=221
       !integer  :: scstart=1, pxstart=1, scend=3000, pxend=221
       integer  :: scstart=1, pxstart=1, scend=4000, pxend=243
       ! if end values are too high, nscans/npix will be end values
       
       integer,parameter    :: nz = 16  !# layers in atmosphere -- not negotiable now
       integer,parameter    :: nbins = 33 , npc = 3 
         ! SST bins, # PCs in EOF LUT files, SST bins for Sy matrices/offsets
       integer,parameter    :: nn = 9 ! max # of iterations to find convergence
       integer,parameter    :: nvar= 7 !# of parameters with potential to be retrieved
       integer,parameter    :: nretvar = 2+npc ! # parameters being retrieved
       real                 :: conv_factor = 50.0 !value to divide nchan
                               !by to determine convergence (Rodgers Eq5.33)
       real                 :: varmult = 1.8 ! a priori variance multiplier (wsp only)
       real                 :: cld_effrad = 12.2 ! cld water droplet eff radius [um]
       real                 :: ice_effrad = 60.0 ! cld ice eff radius [um]
       real                 :: chis1 = 1.5 ! chisq threshold for 'High Quality'
! reminder of variable order: magn EOFs 1-3, WIND, log(LWP), IWP, SST.
!     Bounds on retrieved elements -- all that are capable of being retrieved!
      real :: x_min(nvar) = (/-4.0,-4.0,-4.0,
     >                         0.5,-5.0, 0.0, 270.0/)
      real :: x_max(nvar) = (/ 4.0, 4.0, 4.0,
     >                        35.0,0.12, 0.3, 307.0/) ! increased to 307K

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real :: final_tpw ! pass tpw between subroutines with this?
        integer,parameter :: knd1 = 1
        integer,parameter :: knd2 = 2
        integer,parameter :: knd8 = 8
        real              :: c_sigma(maxchans)    !channel errors
! --------
       real,allocatable   :: oe_tbs(:), test_out(:) ! just for testing sim vs obs tbs
       character(len=256) :: orig_file
       character(len=100) :: input_file, output_file, spec_file 
       character(len=100) :: cal_file, lmskfile, sy_file
       character(len=10)  :: satellite, nchannels, nfreqs
       !character(len=10)  :: satellite, sensor, nchannels, nfreqs
       character(len=20)  :: ppversion
       character(len=3)   :: satcode
       character(len=1)   :: ascdesc ! A or D, (N if neither) from filename
       integer(2)         :: npix, nscans, oepix, oelin
       integer            :: first_good_date  !long integer string of YYYYMMDD'      
       integer(2)         :: avail(maxchans) ! whether ch exists or not
       character(6)       :: freqs(maxchans) ! whether ch exists or not

!--- allocatable arrays used by multiple subroutines
       real, allocatable    :: Tbb(:,:,:), Tb_diff(:,:,:), Tb_sim(:,:,:)! (pix,scan,nchan)
       real, allocatable    :: poster(:,:,:)  ! posteriori error, stddev
       real, allocatable    :: sat_eia(:,:,:)  ! earth incidence angle (pix,scan,chan)
       real, allocatable    :: eia_out(:,:,:)  ! earth incidence angle (pix,scan,2)
       real, allocatable    :: save_tprof(:,:,:),save_mrprof(:,:,:) ! temperature/MR profile output
       real, allocatable    :: save_slp(:,:),save_wdir(:,:) ! SLP, Wind dir
       real, allocatable    :: save_clwc(:,:,:),save_ciwc(:,:,:) ! CLW/CIW
       real, allocatable    :: lat(:,:), lon(:,:)!, sfc_type(:,:) ! (pix,scan)
       real, allocatable    :: sclat(:), sclon(:), scalt(:) !spacecraft lat/lon/altitude
       real, allocatable    :: scorient(:) ! spacecraft forward/backward
       integer, allocatable :: scad(:) !spacecraft orbital direction, asc/descending
       integer, allocatable :: stdtime(:,:) ! time
       integer(1),allocatable ::lsmask(:,:),sfc_type(:,:),lo_flag(:,:) !landmask array
       !integer(2), allocatable :: savessdex(:,:)
       integer(1), allocatable :: save_iter(:,:)
       integer(1), allocatable :: sglint(:,:), qualflag(:,:) ! new for BCMX
       real,allocatable     :: sst(:,:),sstgrid(:,:) !ice and sst grids
       logical, allocatable :: icegrid(:,:)
       real,allocatable     :: freeqs(:) !frequencies (only needed if using RSS model)
       integer    :: nfreeqs

       real                 :: miss_flt = -9999.9
       integer, parameter   :: miss_int = -9999    !missing integer value
       integer, parameter   :: miss_byt = -99      !missing byte value

!--- OE variables passed back to retrieval
       real                 :: xr(nvar)  ! retrieval variables in OE
       real                 :: chisq  ! Chi Squared
       real, allocatable    :: oe_output(:,:,:),screen(:,:)
       real                 :: last_oeout(nvar) ! used for a closer first guess 
       integer              :: loeop=1,loeos=1 ! save pix/scan for fg
       integer              :: run_count=0
       real                 :: avg_iter ! averge # iterations required for retrieval

!--- values from spec table passed back to retrieval
      real :: chisq_out_thresh ! threshold above which oe_output set to missing

!--- variables for ERA-derived LUT table of means/sigmas
      character(len=90) erafile
      integer*2,allocatable :: era_wm(:,:),era_ws(:,:)
      integer*2 :: ssdex
      real, allocatable :: mprof(:,:), eofs(:,:,:), peofs(:,:)
      real :: mrmp(nz)
      integer,parameter :: gm=4 !LUT grid multiplier -- set in LUT creation!
      integer,parameter :: nlon=512*gm,nlat=256*gm!,nmo=12
      real :: losize=.703125, lasize=.701760 ! need to match LUT creation!
      ! gsize is exact for lons, lats start at 89.463 and go by ~.701760
      real,allocatable :: sy(:,:), sy_i(:,:)
      real,allocatable :: toffsets(:),s_toffsets(:,:),s_sy(:,:,:)

!--- CRTM-specific shared variables
        INTEGER :: n_channels  ! dimension 'L'
        INTEGER :: n_profiles  ! dimension 'M'
        INTEGER :: n_sensors   ! dimension 'N'
        integer :: n_layers
        TYPE(CRTM_ChannelInfo_type), ALLOCATABLE :: chinfo(:) !N
        TYPE(CRTM_Geometry_type),ALLOCATABLE :: geo(:) !M
        TYPE(CRTM_Options_type), ALLOCATABLE :: opt(:) !M
        ! Forward declarations
        TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm(:) !M
        TYPE(CRTM_Surface_type), ALLOCATABLE :: sfc(:) !M
        TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:) !LxM      
       
! -- P-level vars, in meters: deltaz between layers, layer t/b height
        real :: lpress(nz+1)= (/100,  200, 300,  400,  500,  550,  600,
     >      650,  700,  750,  800,  850, 900, 925, 950, 975, 1000/)
        real :: dp(nz) =     (/ 100,  100, 100,  100,   50,   50,   50,
     >       50,   50,   50,   50,   50,  25,  25,  25,  25/)
        real :: dz(nz) =     (/2965, 2669, 1997, 1629,  717,  665, 621,
     >      582,  549,  518,  491,  466, 225, 219, 215, 210/)
        real :: lz(nz+1) =   (/14822,11858,9189, 7192, 5563, 4846,4181,
     >     3560, 2978, 2429, 1911, 1420, 954, 729, 509, 295, 85/)
        ! altitudes in [m], derived from mean geopotential at P levels
        real :: pressave(nz)

        ! for new (in BCM3) bounds on EOF coefficients
        real,allocatable :: mr_sigmas(:,:),mr_mins(:,:),mr_maxes(:,:)
        real :: eof_lw_cov(3,nbins)

        ! calculating satellite azimuthal angle
        real,allocatable :: sataz(:,:)
        ! variables read in from binary EC file passed from IDL routine
        real :: winddir(512,256), pwinddir
        real :: windmag(512,256), pwindmag
        real :: eclwp(512,256), peclwp
        real :: eciwp(512,256), peciwp
        real :: ecslp(512,256), pecslp
        real :: ect(512,256,nz+1), pect(nz+1), apect(nz)
        real :: ecmr(512,256,nz+1), pecmr(nz+1)
        real :: pplev(nz+1)
        
        character(6) :: granny ! granule #
        integer(4) :: granule
        integer :: nEIA

      END MODULE define_oeorbit
