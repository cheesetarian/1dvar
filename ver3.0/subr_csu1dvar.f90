      module subr_csu1dvar

!------ contains subroutines for reading in pre-processor files,
!------  reading and processing ancillary data
! added geos5_pp read routine and updated to f90 syntax 12/14/16

       use define_csu1dvar

       implicit none

       contains

      subroutine read_l1c_pp(input_file) ! for TMI L1C!!!
       implicit none
     
!---  routine to read the GPM pre-processor file
       character(len=100) :: input_file
       !real    :: chan_freq(15)

       integer :: iscan,ipix,reccnt=0,i,j,ios,icnt(maxchans)=0,rlun
       integer :: ic, idate
       character(len=128) :: blank = ' '      
       real,parameter  :: Tbmin=50.0, Tbmax=325.0
       logical         :: igood
       
!---  input orbit header structure
       
       type :: OrbitHdr              
           character(len=12) :: satellite   !orbit header 440bytes
           character(len=12) :: sensor
           character(len=12) :: pp_version  
           character(len=128):: radfile
           character(len=128):: pdbfile
           character(len=128):: calfile
           integer           :: granule
           integer           :: nscans
           integer           :: npixels
           integer           :: nchans
           real              :: chan_freqs(maxchans)
           character(len=40) :: comment
       end type OrbitHdr
          
!---  input time structure

       type :: Date         
           integer(kind=knd2):: year
           integer(kind=knd2):: month
           integer(kind=knd2):: day
           integer(kind=knd2):: hour
           integer(kind=knd2):: minute
           integer(kind=knd2):: second
       end type Date    

!---  input scan header structure
       
       type :: ScanHdr
           type(Date)         :: ScanDate
           real               :: Sclat
           real               :: Sclon
           real               :: Scalt
           real               :: Scorient
       end type ScanHdr

!---  input pixel structure

       type :: DataRec     
           real     :: latitude
           real     :: longitude      
           real     :: Tbs(maxchans)       !Tb channels       
           real     :: eia(maxchans)       !earth incident angle
           real     :: Twb                 !Wet Bulb Temperature
           real     :: LR        !Lapse rate in lowest 500m             
           real     :: tcwv 
           real     :: skint               !skin temperature
           real     :: T2m       !2 meter temperature index
           integer  :: qualflag  ! L1C quality flag (different from L1CR read!)
           integer(kind=knd1) :: sunglint_angle
           integer(kind=knd1) :: surface_type_index
           integer(kind=knd1) :: snow_cover_index
           integer(kind=knd1) :: orolift_index  
       end type DataRec
                
!---  assign structure names
       
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: pix

!---  open input preprocessor file

       call gprof_lun(rlun)
       open(unit=rlun,file=trim(input_file),access='stream',status='old',iostat = ios)
        if(ios .ne. 0) then
          write(*,*),' Error opening pp file'
          stop
        endif

!---  read orbit header

       read(rlun, iostat=ios)  orbhdr       !read out orbit header
       !if(ios .ne. 0) call GPM_reprt_err(11,rlun,blank)
!reading preproc orbhdr

!---  write orbit header information

       !write(log_lun,'(a,a)')' ORBITAL HEADER INFO'
       !write(*,*)' Radiometer file= ',trim(orbhdr%radfile)
       !write(log_lun,'(a,a)')' profile DB file= ',trim(orbhdr%pdbfile)
       !write(*,*)' Calibrtion file= ',trim(orbhdr%calfile)
       !write(log_lun,'(a,a)')' comment       = ',trim(orbhdr%comment)       
       !write(*,*)' PP version    = ',trim(orbhdr%pp_version)
       !write(*,*)' satellite     = ',trim(orbhdr%satellite)
       !write(*,*)' sensor        = ',trim(orbhdr%sensor)
       !write(log_lun,'(a,i7)')' granule num  = ',orbhdr%granule
       !write(*,*)' number scans = ',orbhdr%nscans
       !write(*,*)' pixels/scan  = ',orbhdr%npixels       
       !write(*,*)' number chans = ',orbhdr%nchans
       !write(*,*) ' chan freqs = ',orbhdr%chan_freqs
!--- assign these for global variables

       satellite   = orbhdr%satellite
       sensor      = orbhdr%sensor
       ppversion   = orbhdr%pp_version
       convolve    = trim("None")
       orig_file   = orbhdr%radfile
!       pdbfile     = orbhdr%pdbfile
       cal_file    = orbhdr%calfile
       npix        = orbhdr%npixels          
       nscans      = orbhdr%nscans
       !write(granny,'(i06)') (orbhdr.granule+200000)
       granule     = orbhdr%granule
!       nchannels   = orbhdr%nchans
!       chan_freq   = orbhdr%chan_freqs
              
!--- define channel availability array
  
       chan_avail = 0
       do i = 1,maxchans
         if(orbhdr%chan_freqs(i) .gt. 0) chan_avail(i) = 1
       enddo 
       if(nch .ne. sum(chan_avail)) then
         if(nch.eq.7) chan_avail=(/0,0,1,1,1,0,1,1,1,1,0,0,0,0,0/)
         if(nch.eq.5) chan_avail=(/0,0,1,1,1,0,1,1,0,0,0,0,0,0,0/)
       endif
       !write(*,*)'  chan_avail= ',chan_avail(:)

!--- allocate memory for pre-processor scan and pixel data
       allocate (sclat(nscans),sclon(nscans),scalt(nscans))
       allocate (tai93time(nscans))
        tai93time(:) = miss_flt ! not yet calculated for TMI
       allocate (scorient(nscans),scad(nscans))
       allocate (stdtime(nscans,6))

       allocate (lat(npix,nscans))
       allocate (lon(npix,nscans))
       allocate (Tbb(npix,nscans,maxchans))
       allocate (sfc_type(npix,nscans),lo_flag(npix,nscans))
       allocate (sst(npix,nscans))
       allocate (sat_eia(npix,nscans,maxchans))
       allocate (save_tprof(npix,nscans,nz),save_mrprof(npix,nscans,nz))
       allocate (save_clwc(npix,nscans,nz),save_ciwc(npix,nscans,nz))
       allocate (save_slp(npix,nscans),save_wdir(npix,nscans))
       allocate (Tb_diff(npix,nscans,nch),Tb_sim(npix,nscans,nch)) ! all
       allocate (poster(npix,nscans,nvar+1)) ! all variables (nvar)
       allocate (eia_out(npix,nscans,nEIA)) ! 9 for TMI
       allocate (save_iter(npix,nscans)) ! 
       allocate (sataz(npix,nscans)) ! azimuthal angle, rel to N
       allocate (sglint(npix,nscans))!,qualflag(npix,nscans))
       allocate (save_sal(npix,nscans))

!--- loop over and read all scans

       igood = .false.
       do iscan = 1, nscans
       
         read(rlun, iostat=ios)  scnhdr !read in scan header
         if(ios .ne. 0) print*,'scnhdr read promblem'
        
         stdtime(iscan,1) = scnhdr%scandate%year
         stdtime(iscan,2) = scnhdr%scandate%month 
         stdtime(iscan,3) = scnhdr%scandate%day
         stdtime(iscan,4) = scnhdr%scandate%hour
         stdtime(iscan,5) = scnhdr%scandate%minute
         stdtime(iscan,6) = scnhdr%scandate%second

         first_good_date = 0
         if(first_good_date .eq. 0) then
            idate = stdtime(iscan,1)*10000 + stdtime(iscan,2)*100 + stdtime(iscan,3)
            if(idate .gt. 19870101 .and. idate.lt.20500101) then
               first_good_date = idate
            endif
         endif

         sclat(iscan) = scnhdr%sclat
         sclon(iscan) = scnhdr%sclon
         scalt(iscan) = scnhdr%scalt    
         scorient(iscan) = scnhdr%scorient

!--     read all pixel data in each scan

         do ipix = 1, npix
         
           read(rlun, iostat=ios)  pix                         !read in
!pixel header
           !if(ios .ne. 0) call GPM_reprt_err(14,ipix,blank)
          
!           if(pix%latitude  .le. latmin .or.     !lat/lon subset chk 
!     >        pix%latitude  .ge. latmax .or.       
!     >        pix%longitude .le. lonmin .or. 
!     >        pix%longitude .ge. lonmax) then
!                lat(ipix,iscan)       = pix%latitude
!                lon(ipix,iscan)       = pix%longitude 
!                tbs(ipix,iscan,:)     = miss_flt
!                eia(ipix,iscan,:)     = miss_flt
!                Twb(ipix,iscan)       = miss_flt
!                LR(ipix,iscan)        = miss_flt
!                ppskint(ipix,iscan)   = miss_flt
!                ppT2m(ipix,iscan)     = miss_flt
!                pptcwv(ipix,iscan)    = miss_flt
!                sglint(ipix,iscan)   = miss_byt
!                sfccode(ipix,iscan)   = miss_byt
!                orolifti(ipix,iscan)  = miss_byt
!                snowci(ipix,iscan)    = miss_byt
!                pixel_status(ipix,iscan) = 1                 !pixel
!outside of lat/long bounds
           !else            
                reccnt = reccnt + 1                     
                lat(ipix,iscan)       = pix%latitude
                lon(ipix,iscan)       = pix%longitude
                Tbb(ipix,iscan,:)     = pix%tbs
                sat_eia(ipix,iscan,:) = pix%eia
!                Twb(ipix,iscan)       = pix%Twb
!                LR(ipix,iscan)        = pix%LR
!                ppskint(ipix,iscan)   = pix%skint
!                ppT2m(ipix,iscan)     = pix%T2m
!                pptcwv(ipix,iscan)    = pix%tcwv                
                sglint(ipix,iscan)   = pix%sunglint_angle
!                sfccode(ipix,iscan)   = pix%surface_type_index
!                orolifti(ipix,iscan)  = pix%orolift_index
!                snowci(ipix,iscan)    = pix%snow_cover_index

!---           quality control the Tbs - set pixel status if bad
                
!           endif                   
         enddo   !ipix
       enddo  !iscan
       
       !write(*,*)' number of sat pixels read  = ',reccnt
       !write(log_lun,*)' number of sat pixels read  = ',reccnt
       !if(reccnt .eq. 0) call GPM_reprt_err(15,reccnt,blank)
       if(reccnt .eq. 0) print*,'Oh, no!'
        
       
       return      
      end subroutine read_l1c_pp

!-----------------------------------------------------------------
      subroutine read_l1cr_pp(input_file)
       implicit none
     
!---  routine to read the (modified) GPM pre-processor file
       character(len=100) :: input_file
       !real    :: chan_freq(15)

       integer :: iscan,ipix,reccnt=0,i,j,ios,icnt(maxchans)=0,rlun
       integer :: ic, idate
       character(len=128) :: blank = ' '      
       real,parameter  :: Tbmin=50.0, Tbmax=325.0
       logical         :: igood
       
!---  input orbit header structure
       
       type :: OrbitHdr              
           character(len=12) :: satellite   !orbit header 440bytes
           character(len=12) :: sensor
           character(len=12) :: pp_version  
           character(len=128):: radfile
           character(len=128):: pdbfile
           character(len=128):: calfile
           integer           :: granule
           integer           :: nscans
           integer           :: npixels
           integer           :: nchans
           real              :: chan_freqs(maxchans)
           character(len=40) :: comment
       end type OrbitHdr
          
!---  input time structure

       type :: Date         
           integer(kind=knd2):: year
           integer(kind=knd2):: month
           integer(kind=knd2):: day
           integer(kind=knd2):: hour
           integer(kind=knd2):: minute
           integer(kind=knd2):: second
       end type Date    

!---  input scan header structure
       
       type :: ScanHdr
           type(Date)         :: ScanDate
           real               :: Sclat
           real               :: Sclon
           real               :: Scalt
           real               :: Scorient
       end type ScanHdr

!---  input pixel structure

       type :: DataRec     
           real     :: latitude
           real     :: longitude      
           real     :: Tbs(maxchans)       !Tb channels       
           real     :: eia(maxchans)       !earth incident angle
           real     :: Twb                 !Wet Bulb Temperature
           real     :: LR        !Lapse rate in lowest 500m             
           real     :: tcwv 
           real     :: skint               !skin temperature
           real     :: T2m       !2 meter temperature index
           integer(kind=knd1) :: sunglint_angle
           integer(kind=knd1) :: surface_type_index
           integer(kind=knd1) :: snow_cover_index
           integer(kind=knd1) :: orolift_index  
       end type DataRec
                
!---  assign structure names
       
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: pix

!---  open input preprocessor file

       call gprof_lun(rlun)
       open(unit=rlun,file=trim(input_file),access='stream',status='old',iostat = ios)
        if(ios .ne. 0) then
          write(*,*),' Error opening pp file'
          stop
        endif

!---  read orbit header

       read(rlun, iostat=ios)  orbhdr       !read out orbit header
       !if(ios .ne. 0) call GPM_reprt_err(11,rlun,blank)
!reading preproc orbhdr

!---  write orbit header information

       !write(log_lun,'(a,a)')' ORBITAL HEADER INFO'
       !write(*,*)' Radiometer file= ',trim(orbhdr%radfile)
       !write(log_lun,'(a,a)')' profile DB file= ',trim(orbhdr%pdbfile)
       !write(*,*)' Calibrtion file= ',trim(orbhdr%calfile)
       !write(log_lun,'(a,a)')' comment       = ',trim(orbhdr%comment)       
       !write(*,*)' PP version    = ',trim(orbhdr%pp_version)
       !write(*,*)' satellite     = ',trim(orbhdr%satellite)
       !write(*,*)' sensor        = ',trim(orbhdr%sensor)
       !write(log_lun,'(a,i7)')' granule num  = ',orbhdr%granule
       !write(*,*)' number scans = ',orbhdr%nscans
       !write(*,*)' pixels/scan  = ',orbhdr%npixels       
       !write(*,*)' number chans = ',orbhdr%nchans
       !write(*,*) ' chan freqs = ',orbhdr%chan_freqs
!--- assign these for global variables

       satellite   = orbhdr%satellite
       sensor      = orbhdr%sensor
       ppversion   = orbhdr%pp_version
       convolve    = trim("None")
       orig_file   = orbhdr%radfile
!       pdbfile     = orbhdr%pdbfile
       cal_file    = orbhdr%calfile
       npix        = orbhdr%npixels          
       nscans      = orbhdr%nscans
       !write(granny,'(i06)') (orbhdr.granule+100000)
       granule     = orbhdr%granule
!       nchannels   = orbhdr%nchans
!       chan_freq   = orbhdr%chan_freqs
              
!--- define channel availability array
  
       chan_avail = 0
       do i = 1,maxchans
         if(orbhdr%chan_freqs(i) .gt. 0) chan_avail(i) = 1
       enddo 
       if(nch .ne. sum(chan_avail)) then
         if(nch.eq.9) chan_avail=(/1,1,1,1,1,0,1,1,1,1,0,0,0,0,0/)
         if(nch.eq.7) chan_avail=(/0,0,1,1,1,0,1,1,1,1,0,0,0,0,0/)
         if(nch.eq.5) chan_avail=(/0,0,1,1,1,0,1,1,0,0,0,0,0,0,0/)
       endif
       !write(*,*)'  chan_avail= ',chan_avail(:)

!--- allocate memory for pre-processor scan and pixel data

       allocate (sclat(nscans),sclon(nscans),scalt(nscans))
       allocate (tai93time(nscans))
        tai93time(:) = miss_flt
       allocate (scorient(nscans),scad(nscans))
       allocate (stdtime(nscans,6))

       allocate (lat(npix,nscans))
       allocate (lon(npix,nscans))       
       allocate (Tbb(npix,nscans,maxchans))    
       allocate (sfc_type(npix,nscans),lo_flag(npix,nscans))
       allocate (sst(npix,nscans))
       allocate (sat_eia(npix,nscans,maxchans))
       allocate (save_tprof(npix,nscans,nz),save_mrprof(npix,nscans,nz))
       allocate (save_clwc(npix,nscans,nz),save_ciwc(npix,nscans,nz))
       allocate (save_slp(npix,nscans),save_wdir(npix,nscans))
       allocate (Tb_diff(npix,nscans,nch),Tb_sim(npix,nscans,nch)) ! all GMI chans used
       allocate (poster(npix,nscans,nvar+1)) ! all variables (nvar)
       allocate (eia_out(npix,nscans,nEIA)) ! 2 for GMI
       allocate (save_iter(npix,nscans)) ! 
       allocate (sataz(npix,nscans)) ! azimuthal angle, rel to N
       allocate (sglint(npix,nscans))!,qualflag(npix,nscans)) 
       allocate (save_sal(npix,nscans))

!--- loop over and read all scans

       igood = .false.
       do iscan = 1, nscans
       
         read(rlun, iostat=ios)  scnhdr !read in scan header
         if(ios .ne. 0) print*,'scnhdr read promblem'
        
         stdtime(iscan,1) = scnhdr%scandate%year
         stdtime(iscan,2) = scnhdr%scandate%month 
         stdtime(iscan,3) = scnhdr%scandate%day
         stdtime(iscan,4) = scnhdr%scandate%hour
         stdtime(iscan,5) = scnhdr%scandate%minute
         stdtime(iscan,6) = scnhdr%scandate%second

         first_good_date = 0
         if(first_good_date .eq. 0) then
            idate = stdtime(iscan,1)*10000 + stdtime(iscan,2)*100 + stdtime(iscan,3)
            if(idate .gt. 19870101 .and. idate.lt.20500101) then
               first_good_date = idate
            endif
         endif

         sclat(iscan) = scnhdr%sclat
         sclon(iscan) = scnhdr%sclon
         scalt(iscan) = scnhdr%scalt    
         scorient(iscan) = scnhdr%scorient

!--     read all pixel data in each scan

         do ipix = 1, npix
         
           read(rlun, iostat=ios)  pix                         !read in
!pixel header
           !if(ios .ne. 0) call GPM_reprt_err(14,ipix,blank)
          
!           if(pix%latitude  .le. latmin .or.     !lat/lon subset chk 
!     >        pix%latitude  .ge. latmax .or.       
!     >        pix%longitude .le. lonmin .or. 
!     >        pix%longitude .ge. lonmax) then
!                lat(ipix,iscan)       = pix%latitude
!                lon(ipix,iscan)       = pix%longitude 
!                tbs(ipix,iscan,:)     = miss_flt
!                eia(ipix,iscan,:)     = miss_flt
!                Twb(ipix,iscan)       = miss_flt
!                LR(ipix,iscan)        = miss_flt
!                ppskint(ipix,iscan)   = miss_flt
!                ppT2m(ipix,iscan)     = miss_flt
!                pptcwv(ipix,iscan)    = miss_flt
!                sglinta(ipix,iscan)   = miss_byt
!                sfccode(ipix,iscan)   = miss_byt
!                orolifti(ipix,iscan)  = miss_byt
!                snowci(ipix,iscan)    = miss_byt
!                pixel_status(ipix,iscan) = 1                 !pixel
!outside of lat/long bounds
           !else            
                reccnt = reccnt + 1                     
                lat(ipix,iscan)       = pix%latitude
                lon(ipix,iscan)       = pix%longitude
                Tbb(ipix,iscan,:)     = pix%tbs
                sat_eia(ipix,iscan,:) = pix%eia
!                Twb(ipix,iscan)       = pix%Twb
!                LR(ipix,iscan)        = pix%LR
!                ppskint(ipix,iscan)   = pix%skint
!                ppT2m(ipix,iscan)     = pix%T2m
!                pptcwv(ipix,iscan)    = pix%tcwv                
                sglint(ipix,iscan)    = pix%sunglint_angle
!                sfccode(ipix,iscan)   = pix%surface_type_index
!                orolifti(ipix,iscan)  = pix%orolift_index
!                snowci(ipix,iscan)    = pix%snow_cover_index

!---           quality control the Tbs - set pixel status if bad
                
!           endif                   
         enddo   !ipix
       enddo  !iscan
       
       !write(*,*)' number of sat pixels read  = ',reccnt
       !write(log_lun,*)' number of sat pixels read  = ',reccnt
       !if(reccnt .eq. 0) call GPM_reprt_err(15,reccnt,blank)
       if(reccnt .eq. 0) print*,'Oh, no!'
        
       
       return      
      end subroutine read_l1cr_pp

!-------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_geos5_pp(input_file) 

       ! new geos5 preprocessor input (12/14/16)

	! DD 12/14 -- started work here but plenty more to do

!---  routine to read the (modified) GPM pre-processor file
       character(len=100) :: input_file

       integer :: iscan,ipix,reccnt=0,i,j,ios,icnt(maxchans)=0,rlun
       integer :: ic, idate
       character(len=128) :: blank = ' '      
       real,parameter  :: Tbmin=50.0, Tbmax=325.0
       logical         :: igood
       
!---  input orbit header structure
       
       type :: OrbitHdr              
           character(len=12) :: satellite   !orbit header 440bytes
           character(len=12) :: sensor
           character(len=12) :: pp_version  
           character(len=128):: radfile
           character(len=128):: pdbfile
           character(len=128):: calfile
           integer           :: granule
           integer           :: nscans
           integer           :: npixels
           integer           :: nchans
           real              :: chan_freqs(maxchans)
           character(len=40) :: comment
       end type OrbitHdr
          
!---  input time structure

       type :: Date         
           integer(kind=knd2):: year
           integer(kind=knd2):: month
           integer(kind=knd2):: day
           integer(kind=knd2):: hour
           integer(kind=knd2):: minute
           integer(kind=knd2):: second
       end type Date    

!---  input scan header structure
       
       type :: ScanHdr
           type(Date)         :: ScanDate
           real               :: Sclat
           real               :: Sclon
           real               :: Scalt
           !real               :: Scorient
           real*8             :: TAI93 !time in s since 1993
       end type ScanHdr

!---  input pixel structure

       type :: DataRec     
           real     :: latitude
           real     :: longitude      
           real     :: Tbs(maxchans)       !Tb channels       
           real     :: eia(maxchans)       !earth incident angle
           real     :: Twb                 !Wet Bulb Temperature
           real     :: tcwv 
           real     :: skint               !skin temperature
           real     :: T2m                 !2 meter temperature index
           real     :: azimuth             !azimuth angle
           real     :: windmag             !magnitude of wind speed at 10m
           real     :: slp                 !sea level pressure
	   integer*2          :: tprof(nz)        !t profile (K*10)
           integer*2          :: wvprof(nz)       !wv mixratio (g/kg*100)
           integer*2          :: qualflag         !quality flag from L1R
           integer*2          :: winddir          !direction, deg from N
           integer*1          :: lo_pct(4) ! land pct at each FOV size
           integer(kind=knd1) :: sunglint_angle
           integer(kind=knd1) :: surface_type_index
           integer(kind=knd1) :: snow_cover_index
           integer(kind=knd1) :: orolift_index  
       end type DataRec
                
!---  assign structure names
       
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: pix

!---  open input preprocessor file

       call gprof_lun(rlun)
       open(unit=rlun,file=trim(input_file),access='stream',status='old',iostat = ios)
        if(ios .ne. 0) then
          write(*,*),' Error opening pp file'
          stop
        endif

	print*,'czech'

!---  read orbit header

       read(rlun, iostat=ios)  orbhdr       !read out orbit header
       !if(ios .ne. 0) call GPM_reprt_err(11,rlun,blank)
!reading preproc orbhdr

!---  write orbit header information

       !write(log_lun,'(a,a)')' ORBITAL HEADER INFO'
       !write(*,*)' Radiometer file= ',trim(orbhdr%radfile)
       !write(log_lun,'(a,a)')' profile DB file= ',trim(orbhdr%pdbfile)
       !write(*,*)' Calibrtion file= ',trim(orbhdr%calfile)
       !write(log_lun,'(a,a)')' comment       = ',trim(orbhdr%comment)       
       !write(*,*)' PP version    = ',trim(orbhdr%pp_version)
       write(*,*)' satellite     = ',trim(orbhdr%satellite)
       !write(*,*)' sensor        = ',trim(orbhdr%sensor)
       !write(log_lun,'(a,i7)')' granule num  = ',orbhdr%granule
       write(*,*)' number scans = ',orbhdr%nscans
       write(*,*)' pixels/scan  = ',orbhdr%npixels       
       write(*,*)' number chans = ',orbhdr%nchans
       write(*,*) ' chan freqs = ',orbhdr%chan_freqs
	print*,'czech'
!--- assign these for global variables

       satellite   = orbhdr%satellite
       sensor      = orbhdr%sensor
       ppversion   = orbhdr%pp_version
       convolve    = trim("None")
       orig_file   = orbhdr%radfile
!       pdbfile     = orbhdr%pdbfile
       cal_file    = orbhdr%calfile
       npix        = orbhdr%npixels          
       nscans      = orbhdr%nscans
       granule     = orbhdr%granule
!       nchannels   = orbhdr%nchans
              
!--- define channel availability array
  
       chan_avail = 0
       do i = 1,maxchans
         if(orbhdr%chan_freqs(i) .gt. 0) chan_avail(i) = 1
       enddo 
       if(nch .ne. sum(chan_avail)) then
         if(nch.eq.9) chan_avail=(/1,1,1,1,1,0,1,1,1,1,0,0,0,0,0/)
         if(nch.eq.7) chan_avail=(/0,0,1,1,1,0,1,1,1,1,0,0,0,0,0/)
         if(nch.eq.5) chan_avail=(/0,0,1,1,1,0,1,1,0,0,0,0,0,0,0/)
       endif
       write(*,*)'  chan_avail= ',chan_avail(:)

!--- allocate memory for pre-processor scan and pixel data

       allocate (sclat(nscans),sclon(nscans),scalt(nscans))
       allocate (tai93time(nscans))
       ! tai93time(:) = miss_flt
       allocate (scorient(nscans),scad(nscans))
       allocate (stdtime(nscans,6))

       allocate (lat(npix,nscans))
       allocate (lon(npix,nscans))       
       allocate (Tbb(npix,nscans,maxchans))    
       allocate (sfc_type(npix,nscans),lo_flag(npix,nscans))
       allocate (sst(npix,nscans))
       allocate (sat_eia(npix,nscans,maxchans))
       allocate (save_tprof(npix,nscans,nz),save_mrprof(npix,nscans,nz))
       allocate (save_clwc(npix,nscans,nz),save_ciwc(npix,nscans,nz))
       allocate (save_slp(npix,nscans),save_wdir(npix,nscans))
       allocate (Tb_diff(npix,nscans,nch),Tb_sim(npix,nscans,nch)) ! all GMI chans used
       allocate (poster(npix,nscans,nvar+1)) ! all variables (nvar)
       allocate (eia_out(npix,nscans,nEIA)) ! 2 for GMI
       allocate (save_iter(npix,nscans)) ! 
       allocate (sataz(npix,nscans)) ! azimuthal angle, rel to N
       allocate (sglint(npix,nscans))!,qualflag(npix,nscans)) 
       allocate (save_sal(npix,nscans))

!--- loop over and read all scans

       igood = .false.
       do iscan = 1, nscans
       
         read(rlun, iostat=ios)  scnhdr !read in scan header
         if(ios .ne. 0) print*,'scnhdr read promblem'
        
         stdtime(iscan,1) = scnhdr%scandate%year
         stdtime(iscan,2) = scnhdr%scandate%month 
         stdtime(iscan,3) = scnhdr%scandate%day
         stdtime(iscan,4) = scnhdr%scandate%hour
         stdtime(iscan,5) = scnhdr%scandate%minute
         stdtime(iscan,6) = scnhdr%scandate%second

         first_good_date = 0
         if(first_good_date .eq. 0) then
            idate = stdtime(iscan,1)*10000 + stdtime(iscan,2)*100 + stdtime(iscan,3)
            if(idate .gt. 19870101 .and. idate.lt.20500101) then
               first_good_date = idate
            endif
         endif

         sclat(iscan) = scnhdr%sclat
         sclon(iscan) = scnhdr%sclon
         scalt(iscan) = scnhdr%scalt    
         scorient(iscan) = scnhdr%scorient

!--     read all pixel data in each scan

         do ipix = 1, npix
         
           read(rlun, iostat=ios)  pix                         !read in
!pixel header
           !if(ios .ne. 0) call GPM_reprt_err(14,ipix,blank)
          
!           if(pix%latitude  .le. latmin .or.     !lat/lon subset chk 
!     >        pix%latitude  .ge. latmax .or.       
!     >        pix%longitude .le. lonmin .or. 
!     >        pix%longitude .ge. lonmax) then
!                lat(ipix,iscan)       = pix%latitude
!                lon(ipix,iscan)       = pix%longitude 
!                tbs(ipix,iscan,:)     = miss_flt
!                eia(ipix,iscan,:)     = miss_flt
!                Twb(ipix,iscan)       = miss_flt
!                LR(ipix,iscan)        = miss_flt
!                ppskint(ipix,iscan)   = miss_flt
!                ppT2m(ipix,iscan)     = miss_flt
!                pptcwv(ipix,iscan)    = miss_flt
!                sglinta(ipix,iscan)   = miss_byt
!                sfccode(ipix,iscan)   = miss_byt
!                orolifti(ipix,iscan)  = miss_byt
!                snowci(ipix,iscan)    = miss_byt
!                pixel_status(ipix,iscan) = 1                 !pixel
!outside of lat/long bounds
           !else            
                reccnt = reccnt + 1                     
                lat(ipix,iscan)       = pix%latitude
                lon(ipix,iscan)       = pix%longitude
                Tbb(ipix,iscan,:)     = pix%tbs
                sat_eia(ipix,iscan,:) = pix%eia
                Twb(ipix,iscan)       = pix%Twb
                ppskint(ipix,iscan)   = pix%skint
                ppT2m(ipix,iscan)     = pix%T2m
                pptcwv(ipix,iscan)    = pix%tcwv                
                sglint(ipix,iscan)    = pix%sunglint_angle
                sfccode(ipix,iscan)   = pix%surface_type_index
                orolifti(ipix,iscan)  = pix%orolift_index
                snowci(ipix,iscan)    = pix%snow_cover_index

!---           quality control the Tbs - set pixel status if bad
                
!           endif                   
         enddo   !ipix
       enddo  !iscan
       
       !write(*,*)' number of sat pixels read  = ',reccnt
       !write(log_lun,*)' number of sat pixels read  = ',reccnt
       !if(reccnt .eq. 0) call GPM_reprt_err(15,reccnt,blank)
       if(reccnt .eq. 0) print*,'Oh, no!'
        
      return
      end subroutine read_geos5_pp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_l1r_pp(input_file) 

       ! updated to work with 1605_OE preprocessor (5/2/16)

       character(len=100) :: input_file
       integer            :: iscan, ipix, err, reccnt=0, icnt(6)=0
       integer            :: lun,ios, idate, i
       integer            :: n_arguments

!--- input orbit header structure
       
       type :: OrbitHdr
           character(len=256):: origfile
           character(len=90):: comment 
           character(len=20) :: pp_version        
           character(len=10) :: convolved
           character(len=10) :: satellite
           character(len=10) :: sensor
           integer           :: granule
           integer(2)        :: nscans
           integer(2)        :: npixels
           integer(2)        :: nchans
           real              :: chan_freqs(12)
           real              :: missing
           character(len=46):: spares  !!orbit hdr total =500 bytes
       end type OrbitHdr
          
!--- input time structure

       type :: Date         
           integer(2)        :: year
           integer(2)        :: month
           integer(2)        :: day
           integer(2)        :: hour
           integer(2)        :: minute
           integer(2)        :: second
       end type Date    

!--- input scan header structure
       
       type :: ScanHdr
           type(Date) :: ScanDate
           real*8     :: TAI93
           real       :: Sclat
           real       :: Sclon
           real       :: Scalt
       end type ScanHdr

!--- input pixel structure

       type :: DataRec
           real       :: latitude
           real       :: longitude
           real       :: Tb(12) !tchan=12 set in pp
           real       :: EIA(12)
           real       :: azimuth
           real       :: LOFlag 
           real       :: SunglintAM2 ! new!
       end type DataRec

!--- assign structure names
       
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: pix

!--- open input file

       call gprof_lun(lun)
       open(unit=lun,file=trim(input_file),access='stream',status='old',iostat = ios)
        if(ios .ne. 0) then
          write(*,*),' Error opening pp file'
          stop
        endif

!--- read orbit header

       read(lun, iostat=ios)  orbhdr               !read out orbit header

!--- move orbit header information to program variables

       orig_file   = orbhdr.origfile
       ppversion   = orbhdr.pp_version
       convolve    = trim(orbhdr.convolved)
       satellite   = orbhdr.satellite ! get satellite, sensor, pp_ver
       sensor      = orbhdr.sensor    !   from pp file
!       nchans      = orbhdr.nchans
       granule     = orbhdr.granule
       nscans      = orbhdr.nscans
       npix        = orbhdr.npixels
       !chan_freq   = orbhdr.chan_freqs
       !write(granny,'(i06)') (orbhdr.granule+100000)
       !miss_flt    = orbhdr.missing

!--- set channel availability flag
!       do i = 1, 12 ! -- pp now uses 12 channels
!         if(orbhdr.chan_freqs(i) .eq. miss_flt) then
!              chan_avail(i) = 0 !.false.
!         else
!              chan_avail(i) = 1 !.true.
!         endif
!       enddo
!       if(nch .ne. sum(chan_avail)) then
         if(nch.eq.12)chan_avail=(/1,1,1,1,1,1,1,1,1,1,1,1,0,0,0/)
         if(nch.eq.11)chan_avail=(/1,0,1,1,1,1,1,1,1,1,1,1,0,0,0/)
         if(nch.eq.10)chan_avail=(/0,0,1,1,1,1,1,1,1,1,1,1,0,0,0/)
         if(nch.eq.9) chan_avail=(/0,0,1,1,1,1,1,0,1,1,1,1,0,0,0/)
         if(nch.eq.7) chan_avail=(/0,0,0,0,1,1,1,0,1,1,1,1,0,0,0/)
         if(nch.eq.5) chan_avail=(/0,0,0,0,1,1,1,0,1,1,0,0,0,0,0/)
!       endif


!--- check if there are any scans to process - if not then exit

       if(nscans .eq. 0) write(*,*)'nscans is 0!'

!--- allocate memory for pixel level data

       allocate (stdtime(nscans,6))
       allocate (sclat(nscans), sclon(nscans), scalt(nscans))       
       allocate (tai93time(nscans))
       allocate (scad(nscans), scorient(nscans))
       allocate (lat(npix,nscans))
       allocate (lon(npix,nscans))
       allocate (sfc_type(npix,nscans),lo_flag(npix,nscans))
       allocate (sst(npix,nscans))
       allocate (Tbb(npix,nscans,maxchans))    
       allocate (sat_eia(npix,nscans,maxchans)) 
       allocate (save_tprof(npix,nscans,nz),save_mrprof(npix,nscans,nz))
       allocate (save_clwc(npix,nscans,nz),save_ciwc(npix,nscans,nz))
       allocate (save_slp(npix,nscans),save_wdir(npix,nscans))
       allocate (Tb_diff(npix,nscans,nch),Tb_sim(npix,nscans,nch)) ! all chans used
       allocate (poster(npix,nscans,nvar+1)) ! all variables (nvar)
       allocate (eia_out(npix,nscans,nEIA)) ! 2 for AMSR2 (in case using b-scan later)
       allocate (save_iter(npix,nscans)) 
       allocate (sataz(npix,nscans)) ! azimuthal angle, rel to N
       allocate (sglint(npix,nscans))!,qualflag(npix,nscans)) 
       allocate (save_sal(npix,nscans))
              
!--- loop over and read all scans
       first_good_date = 0
       do iscan = 1, orbhdr.nscans
       
         read(lun, iostat=ios)  scnhdr          !write out scan header
               
         stdtime(iscan,1) = scnhdr.scandate.year
         stdtime(iscan,2) = scnhdr.scandate.month 
         stdtime(iscan,3) = scnhdr.scandate.day
         stdtime(iscan,4) = scnhdr.scandate.hour
         stdtime(iscan,5) = scnhdr.scandate.minute
         stdtime(iscan,6) = scnhdr.scandate.second
          
         if(first_good_date .eq. 0) then       
            idate = stdtime(iscan,1)*10000 + stdtime(iscan,2)*100 + stdtime(iscan,3)     
            if(idate .gt. 19870101 .and. idate.lt.20500101) then
               first_good_date = idate
            endif       
         endif

         sclat(iscan)     = scnhdr.sclat ! currently all miss_flt from pp!
         sclon(iscan)     = scnhdr.sclon
         scalt(iscan)     = scnhdr.scalt
         tai93time(iscan) = scnhdr.TAI93
         scorient(iscan)  = 0 ! or 180 for backwards, but AM2 always fwd!
         
!--      read all pixel data
         do ipix = 1, npix
         
           read(lun, iostat=ios)  pix                !read pixel header
           reccnt = reccnt + 1
                         
           lat(ipix,iscan)       = pix.latitude 
           lon(ipix,iscan)       = pix.longitude
                           
           Tbb(ipix,iscan,1:12)     = pix.Tb(:)         !12channels
           sat_eia(ipix,iscan,1:12) = pix.EIA(:)         !12channels
           sataz(ipix,iscan)        = pix.azimuth  !
           lo_flag(ipix,iscan)      = int(pix.LOFlag)   ! NEW!
           sglint(ipix,iscan)       = int(pix.sunglintAM2)
           !^ calculated from sun azimuth, sun elevation, 
           ! and EIA in L1R pp
         enddo   !ipix
       enddo  !iscan
       !print*,'scnhdr date:',scnhdr.scandate.year,scnhdr.scandate.month
       
!--- close input file   
       close(lun)
       !write(*,*)'   num pixel read     = ',reccnt
       !print*,' random Tbs from pixel: ',Tbb(179,800,:)

      return
      end subroutine read_l1r_pp

      subroutine read_landmask
       implicit none         

       integer(2) :: ixsize, iysize, irecl, lmaskres 
       character(len=2)   :: cres, sres
       integer            :: db_lun, ios, i,j
 
!--- get Mask grid resolution from filename
       read(lmskfile(15:16),'(i2)',iostat=ios) lmaskres
!      write(*,*)'   landmask resolution          : ',lmaskres
       
!--- set size of landmask file       
       ixsize = lmaskres * 360 
       iysize = lmaskres * 180
       irecl = ixsize / 4
       
!--- allocate landmask variable
       allocate (lsmask(ixsize,iysize))!gets deallocated later

!--- open landmask file
!      write(*,*)'   reading landmask             : ',
!    >                   trim(lmskfile)       
       call gprof_lun(db_lun)
       !open(unit=db_lun,file=trim(lmskfile),
       open(unit=db_lun,file='binary/landmask56-32_16.bin',&
	form='unformatted', status='old',recl = irecl,access='direct',iostat=ios)
       if(ios.ne.0) write(*,*)'open error on landmask'

!--- read in land-sea mask file
       do j = 1,iysize
         read(db_lun,rec=j,iostat=ios)(lsmask(i,j),i=1,ixsize)
         if(ios.ne.0) write(*,*)'read error on landmask'
       enddo
       close(db_lun)            
       return

      end subroutine read_landmask


      subroutine read_reynolds_sst
       implicit none

!       stealing the gprof sst reader!
!-----------------------------------------------------------------
! This routine will read in a Reynolds 0.25X0.25 daily SST and Sea 
! ice file selected from the time as specified in stdtime.  If the 
! individual SST file is not available, the most recent file is used.
!-----------------------------------------------------------------
    
       integer(2),parameter :: xsize=1440, ysize=720
       logical :: fexists
       character(len=256):: tfile
       character(len=100) :: dir_sst 
       character(len=8)  :: cdate
       character(len=4)  :: cyear
       integer :: i, j, k, m, read_date
       integer :: iyear,imon,iday, iyr,imo,ida       
       
       integer :: rlun, ios
       real    :: sbuff(1440)
       integer(2) :: ibuff(1440), err, anom            
       integer :: numday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

       integer(2),allocatable :: fice(:,:)
       integer(2),allocatable :: fsst(:,:)

       dir_sst = '/xdata/drandel/gpm/ppingest/sstdata_0.25/' !hard-coded!
       allocate (icegrid(xsize,ysize) )   
       allocate (sstgrid(xsize,ysize) )  
       allocate (fsst(xsize,ysize) )         !temp sst from file
       allocate (fice(xsize,ysize) )         !temp sea-ice from file
                
!--- Find the Reynolds SST filename for this date 
       read_date = first_good_date       
       !write(*,*) '   read_date = ', read_date
        
       do i = 1,30                !search back 30 days for file
       
         write(cdate,'(i8)') read_date
         read(cdate(1:4),'(i4)') iyear
         read(cdate(5:6),'(i2)') imon
         read(cdate(7:8),'(i2)') iday
             !print*,cdate,iyear,imon,iday 
         if(read_date.lt.20020601 .or. read_date.gt.20111004) then
!AVHRR+AMSR starts this date               
               tfile = trim(dir_sst) // cdate(1:4) // '/avhrr-only-v2.' // cdate
               inquire(file = tfile,exist = fexists)     
         else               
               tfile = trim(dir_sst) // cdate(1:4) // '/amsr-avhrr-v2.' // cdate
               inquire(file = tfile,exist = fexists)
         endif 
         if(i .eq. 5) then
            write(*,*)' WARNING: sst more than 5 days old'
         endif  
         if(fexists) exit              !file for requested date found

         if(iday-1 .gt. 0) then        !if no file, look back in time
            read_date = read_date - 1  !still within month
         else            
            imon = imon - 1               !look in previous month
            if(imon .lt. 1) then          !look in previous year
                imon = 12
                iyear = iyear - 1
            endif
            if(mod(iyear,4) .eq. 0) numday(2) = 29
            iday = numday(imon)
            read_date = iyear *10000 + imon*100 + iday           
         endif   
       enddo
       if(.not. fexists) then       
            write(*,*)'no SST file : ',trim(tfile)
       endif 

       call gprof_lun(rlun)
!      write(*,*)'   reading sst/ice, filename    : ',
!    >                  trim(tfile)
       open(rlun,file=trim(tfile), form='unformatted',status='old',&
	 access='sequential',convert='big_endian',iostat=ios)
       if(ios.ne.0) write(*,*)'error opening SST file'

!--- read in SST and ICE grids
       
       read(rlun,iostat=ios)iyr,imo,ida,((fsst(i,j),i=1,xsize),j=ysize,1,-1) !flip SST for NP up
       !write(*,*)'writing out sst grid!: ',fsst
       read(rlun,iostat=ios)iyr,imo,ida,((anom,i=1,xsize),j=1,ysize)  !don't use anom field
       read(rlun,iostat=ios)iyr,imo,ida,((err,i=1,xsize),j=1,ysize)  !don't use err field
       if(ios.ne.0) write(*,*)'issue 3' 
       read(rlun,iostat=ios)iyr,imo,ida,((fice(i,j),i=1,xsize),j=ysize,1,-1) !flip ICE for NP up
       if(ios.ne.0) write(*,*)'issue 4'
       
!--- close file
       close(rlun)
      
!--- shift and scale sst field to center on 0 degrees

       do j = 1,ysize
         do i = 1,720
           if(fsst(i+720,j) .ne. -999) then
               sstgrid(i,j) = fsst(i+720,j) / 100.
           else
               sstgrid(i,j) = miss_flt
           endif
           
           if(fsst(i,j) .ne. -999) then    
               sstgrid(i+720,j)= fsst(i,j) / 100.
           else
               sstgrid(i+720,j) = miss_flt
           endif
         enddo
       enddo
       
!--- shift ice field to center on 0 degrees
     
       do j = 1,ysize
         do i = 1,720
           ibuff(i)     = fice(i,j)
           fice(i,j)    = fice(i+720,j)
           fice(i+720,j)= ibuff(i)
         enddo
       enddo
                         
!--- assign ice to icegrid logical grid       
     
       do i = 1,xsize
         do j = 1, ysize
            if(fice(i,j) .gt. 0 .and. fice(i,j).le. 100) then
                icegrid(i,j) = .true.
            else
                icegrid(i,j) = .false.
            endif
         enddo
       enddo

c--- call routine to fill in AMSRE derived Caspian Sea Ice

       !call read_seaice_amsre_csea(imon,iday)       

       deallocate (fsst, fice)

      return
      end subroutine read_reynolds_sst

      subroutine assign_sfc_type

!------------------------------------------------------
!          routine to assign the sensor specific 
!          resolution landmask and sea-ice climatology 
!          to the lat-lon scan points
!
!        outputs in sfc_type are: 1 = 100% water (0-2% land)
!                                 0 = EVERYTHING ELSE
!                                 2 = sea ice - (1-100% ice)
!-------------------------------------------------------
                                        
      implicit none

      real               :: alat, alon
      real               :: rinc, rinci
      integer            :: lin, pix, isea, ir, iland, isst!, icnt(9)=0
      integer(2) :: ilat, ilon, islat,islon, imaskres
      integer(2) :: iclat, iclon, lmaskres
      integer(2) :: i,j,ixsize,iysize

!--- decode the landmask resolution from the landmask filename

      read(lmskfile(15:16),'(i2)') lmaskres

        !write(*,*)' landmask file: ',trim(lmskfile)
        !write(*,*)' landmask res: ',lmaskres

!--- set input landmask grid size

      ixsize = lmaskres * 360                   !landmask size
      iysize = lmaskres * 180 
      rinc     = (1. / float(lmaskres))  * 0.5  !1/2 of landmask grid increment

      imaskres = 4                             !ice - 0.25 x 0.25 degree 
      rinci    = (1. / float(imaskres)) * 0.5  !1/2 of Sea-ice grid increment

      isea = 0
      iland = 0
      do lin = 1, nscans                   !nscans = number of scans
        do pix = 1, npix           !NPIX = number of pixels per scan
          
          alat = lat(pix,lin)
          alon = lon(pix,lin)            !input centered on 0 deg
          ilat = nint(lmaskres *(90.+rinc-alat))    !assign land array element
          if(ilat .gt. iysize) ilat = iysize        !chk for northpole
          ilon = nint(lmaskres * (180+rinc+alon))   !for landmask
          !write(*,*)' ilat/ilon: ',ilat,ilon
          if(ilon .gt. ixsize) ilon = 1
          if ((ilat .lt. 1) .or. (ilon .lt. 1) )  then ! boundary check
            sfc_type(pix,lin) = 0
            if(satcode.eq.'GMI'.or.satcode.eq.'TMI') then
              lo_flag(pix,lin) = -1 ! NEW
            endif
          else   
 
            if(satcode.eq.'GMI'.or.satcode.eq.'TMI') then 
              lo_flag(pix,lin)=lsmask(ilon,ilat) ! for non-AMSR

              if(lsmask(ilon,ilat) .ge. 0 .and. & !0-2% land = ocean -- CHANGE!
                lsmask(ilon,ilat) .le. landlimit) then   ! changed from 2% cutoff!!
                 sfc_type(pix,lin) = 1 
              elseif(lsmask(ilon,ilat) .gt. landlimit) then
                 sfc_type(pix,lin) = 0 
              endif
            endif                                
            if(satcode.eq.'AM2') then ! 0-2% land = ocean for AM2!
              if(lo_flag(pix,lin).ge.0 .and.lo_flag(pix,lin).le.landlimit) then
                sfc_type(pix,lin)=1

        ! adding in new expansion of screening, d duncan, 9/28/16
               if((lin<=(nscans-3) .and.(lo_flag(pix,lin+1)>landlimit &
               .or.lo_flag(pix,lin+2)>landlimit & 
               .or.lo_flag(pix,lin+3)>landlimit))   .or. &
               (lin<=(nscans-2) .and.(lo_flag(pix,lin+1)>landlimit &
               .or.lo_flag(pix,lin+2)>landlimit))   .or. &
               (lin<=(nscans-1).and.(lo_flag(pix,lin+1)>landlimit)).or. &
               (lin>=4 .and.(lo_flag(pix,lin-1)>landlimit &
               .or.lo_flag(pix,lin-2)>landlimit &
               .or.lo_flag(pix,lin-3)>landlimit))   .or. &
               (lin>=3 .and.(lo_flag(pix,lin-1)>landlimit &
               .or.lo_flag(pix,lin-2)>landlimit))   .or. &
               (lin>=2 .and.(lo_flag(pix,lin-1)>landlimit))       .or. &
               (pix<=(npix-3) .and.(lo_flag(pix+1,lin)>landlimit &
               .or.lo_flag(pix+2,lin)>landlimit &
               .or.lo_flag(pix+3,lin)>landlimit))   .or. &
               (pix<=(npix-2) .and.(lo_flag(pix+1,lin)>landlimit &
               .or.lo_flag(pix+2,lin)>landlimit))   .or. &
               (pix<=(npix-1) .and. lo_flag(pix+1,lin)>landlimit) .or. &
               (pix>=4 .and.(lo_flag(pix-1,lin)>landlimit &
               .or.lo_flag(pix-2,lin)>landlimit &
               .or.lo_flag(pix-3,lin)>landlimit))   .or. &
               (pix>=3 .and.(lo_flag(pix-1,lin)>landlimit &
               .or.lo_flag(pix-2,lin)>landlimit))   .or. &
               (pix>=2 .and. lo_flag(pix-1,lin)>landlimit) .or. &
               (pix>=2.and.lin>=2 .and. & ! +/- 1 diagonal 
                lo_flag(pix-1,lin-1)>landlimit)     .or. &
               (pix<=(npix-1).and.lin>=2 &
               .and.lo_flag(pix+1,lin-1)>landlimit) .or. &
               (pix>=2.and.lin<=(nscans-1) &
               .and.lo_flag(pix-1,lin+1)>landlimit) .or. &
               (pix<=(npix-1).and.lin<=(nscans-1) &
               .and.lo_flag(pix+1,lin+1)>landlimit) .or. &
               (pix>=3.and.lin<=(nscans-2) & ! +/- 2 diagonal 
               .and.lo_flag(pix-2,lin+2)>landlimit) .or. &
               (pix>=3.and.lin>=3 &
               .and.lo_flag(pix-2,lin-2)>landlimit) .or. &
               (pix<=(npix-2).and.lin<=(nscans-2) &
               .and.lo_flag(pix+2,lin+2)>landlimit) .or. &
               (pix<=(npix-2).and.lin>=3 &
               .and.lo_flag(pix+2,lin-2)>landlimit)) &
              sfc_type(pix,lin)=0

              elseif(lo_flag(pix,lin) > landlimit) then 
                sfc_type(pix,lin)=0
              endif
            endif
          endif
           
          if(sfc_type(pix,lin) .eq. 0) cycle !dont chk seaice over
                                             !non-ocean
                  
          iclat = nint(float(imaskres)*(90.+rinci-alat)) !assign ice array elements
          if(iclat .gt. 180*imaskres) iclat = 180*imaskres  !chk for southpole
          iclon = nint(imaskres * (180+rinci+alon))         
          if(iclon .gt. (imaskres*360)) iclon = 1
          if ((iclat .lt. 1) .or. (iclat .gt. imaskres*180) .or. & !boundary check 
             (iclon .lt. 1) .or. (iclon .gt. imaskres*360))  then
              isea = isea + 1
              sfc_type(pix,lin) = 0 ! just set to no retrieval
              cycle
          endif

          if(icegrid(iclon,iclat)) then
            !if(runseaice) then
              sfc_type(pix,lin) = 2 ! doesn't matter if ice over ocean or coast... 
            !else
            !  sfc_type(pix,lin) = 0
            !endif
          endif
          
        enddo  !pix
      enddo   !lin
            
      !write(*,*)' land mask values(one lon,all lats): ',lsmask(1,:)  
      !write(*,*)' sfc_type values(one pixel,all scans): ',sfc_type(1,:)  

      deallocate (lsmask)       !free up the landmask memory
      deallocate (icegrid)      !free up the icegrid

      return
      end subroutine assign_sfc_type

!------------------------------------------------------------------------------

      subroutine assign_sst

      implicit none

      real               :: alat, alon, rincs
      real               :: scl, sst1, sst2
      real               :: slat, slat1, slat2
      real               :: slon, slon1, slon2
      integer            :: lin, pix, isst1=0, isst2=0
      integer(2)         :: ilat, ilon, islat, islon, smaskres
      integer(2)         :: islat1, islat2, islon1, islon2
      integer(2)         :: i,j,ixsize,iysize
      

      smaskres = 4           ! ppd  (4 = 0.25 X 0.25) 
      rincs    = (1. / float(smaskres))  * 0.5 !1/2 of grid increment

      do lin = 1, nscans
        do pix = 1, npix
                !write(*,*)'sfc_type(pix,lin)=',sfc_type(pix,lin)
          if(sfc_type(pix,lin) .ne. 1) then
                 sst(pix,lin) = miss_flt   ! 20=land, 30=coast, ice handled below
                 !write(*,*)'sst(pix,lin)=',sst(pix,lin)
                 cycle      
          endif
          alat = lat(pix,lin)            !input lat
          alon = lon(pix,lin)            !input lon centered on 0 deg
          alon=alon + 180.0              !shift for grid counter 

          islat = nint(smaskres*(90.+rincs-alat)) !assign sst array elements
          slat = smaskres*(90.+rincs-alat)
          !if(islat .eq. 181*smaskres) islat = 180*smaskres !chk forsouthpole
          islon = nint(smaskres*(rincs+alon))
          slon = smaskres*(rincs+alon)

          if(islon .eq. (smaskres*360+1)) islon = 1
          if ((islat .lt. 1) .or. (islat .gt. smaskres*180) .or. &!boundary check 
             (islon .lt. 1) .or. (islon .gt. smaskres*360))  then
              sst(pix,lin) = miss_flt
              cycle
          else
              if(sstgrid(islon,islat) .lt. -1.8) then ! why this cutoff?
                  sst(pix,lin) = miss_flt
                  sfc_type(pix,lin) = 0 !31  !assign
                  cycle
              else
                  slat1= slat - islat ! lat index fraction -.5 to .5
                  slon1= slon - islon ! lon index fraction -.5 to .5
                  if(slat1.eq.0.0.and.slon1.eq.0.0) then
                    sst(pix,lin) = sstgrid(islon,islat) + 273.15 
                    cycle
                  endif
                  if(islat.eq.1.or.islat.eq.720.or.islon.eq.1.or.islon.eq.1440) then 
                    sst(pix,lin) = sstgrid(islon,islat) + 273.15 
                    cycle
                  endif
                  if(minval(sstgrid(islon-1:islon+1,islat-1:islat+1)).lt.-1.8) then
                    sst(pix,lin) = sstgrid(islon,islat) + 273.15 
                    cycle
                  endif

                  if(slat1.gt.0) slat2=sstgrid(islon,islat+1)-sstgrid(islon,islat)
                  if(slat1.lt.0) slat2=sstgrid(islon,islat-1)-sstgrid(islon,islat)
                  if(slon1.gt.0) slon2=sstgrid(islon+1,islat)-sstgrid(islon,islat)
                  if(slon1.lt.0) slon2=sstgrid(islon-1,islat)-sstgrid(islon,islat)

                  sst(pix,lin) = 273.15 + sstgrid(islon,islat) + slat2*abs(slat1) + slon2*abs(slon1)
              endif
          endif
              
        enddo   !pix
      enddo     !lin

      deallocate (sstgrid)       !free up the sstgrid

      return
      end subroutine assign_sst


      subroutine prep_oe
        implicit none

        real :: freq,temp,pres,rho
        integer :: bin,era_lun1,era_lun2,era_lun3,rsslun,ios,syl,i,j,k
        character(len=4) :: yr
        character(len=2) :: mos(12) = (/'1','2','3','4','5','6','7','8','9','10','11','12'/)
        character(len=2) :: das(31) = (/'1','2','3','4','5','6',&
         '7','8','9','10','11','12','13','14','15','16','17',&
        '18','19','20','21','22','23','24','25','26','27','28','29','30','31'/)
        character(len=2) :: mo, da
        character(len=3) :: emis
        character(len=1) :: spc
        character(len=100):: sydir
        sydir = '/tdata1/dduncan/oe/sy/'

      allocate(oe_output(npix,nscans,9),screen(npix,nscans)) !last dimension hard-coded!
      allocate(toffsets(nch),s_toffsets(nbins,nch),s_sy(nch,nch,nbins))
      allocate(sy(nch,nch),sy_i(nch,nch)) 

      write(spc,'(I1)') npc
        !spc = '3' ! for testing
      call gprof_lun(syl)

      emis = 'F' ! default
      if(run_rss) emis = 'R'
      if(useanalysis) emis = trim(emis) // 'a'
      if(useanalysis.ne..true.) emis = trim(emis) // 'c' !climatology
      if(nretvar.gt.5) emis=trim(emis)//'s'
        !emis='Fc' !OVERRIDE FOR TESTING!
        
      sy_file = trim(sydir)//satcode// '.SY.'//trim(emis)// '.'//trim(spc)//'.'//trim(n_ch)//'.v2.1.bin'
      open(unit=syl,file=trim(sy_file),access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open sy_file',sy_file
        read(syl) sy
        read(syl) toffsets
        read(syl) s_sy
        read(syl) s_toffsets
      close(syl)

      mr_sigmas(:,1) = (/1.33,1.27,1.37,1.45,1.47,1.44,1.40,1.34,1.30,&
       1.25,1.20,1.15,1.15,1.14,1.11,1.11,1.12,1.12,1.10,1.07,1.08,&
       1.05,1.02,1.00,0.98,1.03,1.05,1.09,1.20,1.20,1.20,1.20,1.20/)
      mr_sigmas(:,2) = (/0.76,0.87,0.95,0.94,0.93,0.87,0.88,0.91,0.91,&
       0.89,0.86,0.85,0.88,0.92,0.93,0.95,0.95,0.96,0.94,0.91,0.86,&
       0.79,0.74,0.68,0.63,0.60,0.60,0.66,0.66,0.66,0.66,0.66,0.66/)
      mr_sigmas(:,3) = 0.40

      call gprof_lun(syl)
      open(unit=syl,file= 'binary/SA-EOF_LW.3.v2.0c.bin',access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open eof-lw cov file'
        read(syl) eof_lw_cov ! 3 x nbins
      close(syl)
     
!--- allocate and read in EC LUT data
      allocate(era_wm(nlon,nlat), era_ws(nlon,nlat))
      allocate(mprof(nbins,nz),eofs(nbins,nz,6))!npc))
      allocate(peofs(nz,6))

      write(yr,'(I4)') stdtime(5,1)
      mo = mos(stdtime(5,2)) !use as index
      da = das(stdtime(5,3)) !use as index

      erafile = 'binary/eof_mr.03.v3.bin' ! v3 has consistency across bins
      call gprof_lun(era_lun1)
      open(unit=era_lun1,file=erafile,access='stream',status='old',form='unformatted',iostat = ios)
       if(ios .ne. 0) write(*,*)' cant open EOF LUT file'
        read(era_lun1) mprof
        read(era_lun1) eofs
      close(era_lun1)

      erafile ='/tdata1/dduncan/oe/LUT/ws_n128x4_vC.'//trim(mo)//'.bin'
       call gprof_lun(era_lun1)
      open(unit=era_lun1,file=erafile,access='stream',status='old',form='unformatted',iostat = ios)
       if(ios .ne. 0) write(*,*)' cant open EC winds file'
        read(era_lun1) era_wm
        read(era_lun1) era_ws
      close(era_lun1)

      end subroutine prep_oe

      subroutine qual_czech
       implicit none
!----- new 11/17/15, D Duncan
!      this will screen out anomalously high TPW
!       pixels in which it is likely raining lightly or LWP
!       is being traded off for more WV.

       integer :: lin,pix,ct,s1,p1
       integer,parameter :: sd=3,pd=5
       real    :: sigmatp,meantp
       real    :: meth=5.0,sith=6.5,tottp,eachtp((2*pd+1)*(2*sd+1))
       real    :: temptp(2*pd+1,2*sd+1), tempch(2*pd+1,2*sd+1)
       
       screen(:,:) = oe_output(:,:,1)
       do lin = 1+sd,nscans-sd
         do pix = 1+pd,npix-pd
           if(oe_output(pix,lin,1).lt.0) cycle
           temptp(:,:) = oe_output((pix-pd):(pix+pd),(lin-sd):(lin+sd),1)
           tempch(:,:) = oe_output((pix-pd):(pix+pd),(lin-sd):(lin+sd),5)
           ct=0
           tottp=0.0
           do s1=1,2*sd+1
             do p1=1,2*pd+1
               if(temptp(p1,s1).gt.0.and.tempch(p1,s1).lt.chis1.and.temptp(p1,s1).ne.oe_output(pix,lin,1)) then
                 ct=ct+1
                 tottp=tottp+temptp(p1,s1)
                 eachtp(ct) = temptp(p1,s1)
               endif
             enddo
           enddo
           if(ct.lt.15) cycle ! need enough to get real background sense
           meantp = tottp/real(ct)
           sigmatp = SQRT(SUM((eachtp(1:ct)-meantp)**2)/real(ct-1))
           if(oe_output(pix,lin,1).gt.(meantp+meth) .or. sigmatp.gt.sith) then
             screen(pix,lin) = -996
           endif 
         enddo
       enddo
       oe_output(:,:,1) = screen(:,:)

      end subroutine qual_czech

      subroutine dealloc
! --- Deallocate variables used in output procedure, now ncdf
       deallocate(sst,lat,lon)
       deallocate(oe_output,screen)
       deallocate(stdtime,tai93time)
       deallocate(sclat,sclon,scalt,scorient,scad)
       deallocate(Tbb,sat_eia,Tb_diff,Tb_sim,eia_out,poster)
       deallocate(save_tprof,save_mrprof,save_slp,save_wdir)
       deallocate(save_clwc,save_ciwc,save_sal)
       deallocate(sfc_type,lo_flag)
       deallocate(save_iter,sataz,save_inc)
       deallocate(sglint)
       deallocate(era_wm,era_ws)
       deallocate(eofs,mprof,peofs)
       deallocate(oe_tbs,test_out)
       deallocate(sy,sy_i,s_sy,toffsets,s_toffsets)
       end subroutine dealloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine gprof_lun(ilun)
!    
!** This routine gets an open logical unit number starting with 100
!
       integer  :: ilun
       logical  :: llun
!       
       do ilun = 100,201
         inquire(unit=ilun, opened=llun)
         if(.not. llun) exit
       enddo
       return
      end subroutine gprof_lun

!--------------------------------------------------------------------
      subroutine calc_az(slat,ipi,ascdes,sc_orient,rot_ang)

!       calculating azimuthal angle of pixel, 
!        relative to due North, for GMI/TMI's footprint

!      slat = spacecraft lat
!      ipi  = position in scan (pixel number, 1-221)   
!      ascdes= ascending/descending pass (1=asc, -1=desc)    
!      sc_orient = spacecraft orientation (180=backward, 0=forward)

      real :: slat,sc_orient
      integer :: nfovhigh, iatt, ascdes, ipi!, badorient=0
      real*8 :: satincl != 65.00000000 ! satellite inclination, relative to E
      real*8 :: scanangle != 140.00
      real*8 :: dpi, d2r, theta, gtrack
      real   :: rot_ang

      if(satcode.eq.'GMI') then
        nfovhigh  = 221 ! # FOVs in high freq chans, thus # pix in L1C/L1CR
        satincl   = 65.000000
        scanangle = 140.00
      endif
      if(satcode.eq.'TMI') then
        nfovhigh  = 104 ! # FOVs in high freq chans, thus # pix in L1C/L1CR
        satincl   = 35.000000
        scanangle = 130.00
      endif

      iatt = 0
      dpi = 4.0d+0*atan(1.0) ! define Pi
      d2r = dpi / 180.0d+0   ! degrees to radians

        ! take into account direction of spacecraft -- rear/forward
      if(sc_orient.eq.180.0) iatt=-1
      if(sc_orient.eq.0.0)   iatt=1
      !if(iatt.eq.0) then
      !  badorient = badorient+1 ! counter for pixels with orientation
        ! not forward or backward, if wanting to output later to file?
        ! algorithm will skip this
        ! pixel due to missing Tbs anyway, so let it go...
      !endif

        ! make sure sat_lat works in below code
      if(slat.gt.satincl) slat = satincl
      if(slat.lt.(-1.0*satincl)) slat = -1.0*satincl

        ! contribution to azimuthal angle by pixel position
      theta = ((-0.5d+0*SCANANGLE) + (SCANANGLE*(ipi-1))/(NFOVHIGH-1.0d+0))
        ! which is different if spacecraft moving backwards!
      if(iatt .eq. -1) then
        theta=((-0.5d+0*scanangle -180.0) +(scanangle*(ipi-1))/(nfovhigh-1.0d+0))
      endif

        ! contribution of lat on azimuthal angle
      gtrack = ascdes * (90.0d+0 - asin(cos(satincl*d2r)/cos(abs(slat)*d2r))/d2r)

        ! now find angle between footprint long axis (from sat dir)
        !  and North, clockwise. gotta rotate 90 degrees, too.
      rot_ang = (90.0d+0 - gtrack - theta)
      if(rot_ang .ge. 360.0) rot_ang = rot_ang-360.0

      end subroutine calc_az

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!-----------------------------------------------------------------------
      subroutine import_sal
      ! opens/reads binary file of Aquarius monthly SSS data from 2014, 
      !  (smoothed and formatted offline) to 720x360 grid, starting at
      !  0E,90N with mis_val=-999.9, type real
      ! D Duncan, CSU, 2/18/16

        integer :: lun,ios
        character(len=80) sfile
        character(len=180) dir
        character(len=2) :: money
        character(len=2) :: mons(12) = (/'01','02','03','04','05','06','07','08','09','10','11','12'/)
        
        money = mons(stdtime(500,2)) !use as index
        !dir = 'binary/' ! kept local, for now
        dir = '/tdata1/dduncan/Aquarius_SSS/'
        sfile = trim(dir) // 'sss.'//money//'.v1.bin'

        call gprof_lun(lun)
        open(unit=lun,file=trim(sfile),access='stream',status='old',form='unformatted',iostat = ios)
         if(ios .ne. 0) write(*,*)' cant open sss file'
          read(lun) sss
        close(lun)

      end subroutine import_sal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module subr_csu1dvar
