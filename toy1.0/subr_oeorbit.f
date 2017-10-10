      module subr_oeorbit

       use define_oeorbit

       implicit none

       contains


!-----------------------------------------------------------------
      subroutine read_l1cr_pp(input_file)
       implicit none
     
!---  routine to read the (modified) GPM pre-processor file
       character(len=100) :: input_file
       real    :: chan_freq(15)

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
       open(unit=rlun,file=trim(input_file),access='stream',
     >      status='old',iostat = ios)
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
       orig_file   = orbhdr%radfile
!       pdbfile     = orbhdr%pdbfile
       !cal_file    = orbhdr%calfile
       npix        = orbhdr%npixels          
       nscans      = orbhdr%nscans
       write(granny,'(i06)') (orbhdr.granule+100000)
       granule     = orbhdr%granule
!       nchannels   = orbhdr%nchans
!       chan_freq   = orbhdr%chan_freqs
              
!--- define channel availability array
  
       !chan_avail = 0
       !do i = 1,maxchans
       !  if(orbhdr%chan_freqs(i) .gt. 0) chan_avail(i) = 1
       !enddo 
       !write(*,*)'  chan_avail= ',chan_avail(:)

       !avail = chan_avail

!--- allocate memory for pre-processor scan and pixel data

       allocate (sclat(nscans),sclon(nscans),scalt(nscans))
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
       allocate (eia_out(npix,nscans,2)) ! just low/high freqs
       allocate (save_iter(npix,nscans)) ! 
       allocate (sataz(npix,nscans)) ! azimuthal angle, rel to N
       allocate (sglint(npix,nscans),qualflag(npix,nscans)) 

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
            idate = stdtime(iscan,1)*10000 + stdtime(iscan,2)*100 +
     >              stdtime(iscan,3)
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
      subroutine read_l1r_pp(input_file) 

       ! updated to work with 1511_OE preprocessor (11/9/15)

       character(len=3)   :: satcode 
       character(len=100) :: input_file
       integer            :: iscan, ipix, err, reccnt=0, icnt(6)=0
       integer            :: lun,ios, idate, i
       integer            :: n_arguments

       real    :: chan_freq(12)
       logical :: chan_avail(12)

!--- input file header structure

       type :: fileHdr       
           character(len=80) :: sstdir
           character(len=80) :: ancdir       
           character(len=80) :: lmskfile
           character(len=80) :: kfile       
           character(len=80) :: pdbfile
           character(len=80) :: lpdbfile
           character(len=80) :: cpdbfile
           character(len=5)  :: landretsat !file hdr total =565 bytes
       end type fileHdr              

!--- input orbit header structure
       
       type :: OrbitHdr
           character(len=256):: origfile
           character(len=100):: comment 
           character(len=20) :: pp_version        
           character(len=10) :: satellite
           character(len=10) :: sensor
           integer           :: granule
           integer(2)        :: nscans
           integer(2)        :: npixels
           integer(2)        :: nchans
           real              :: chan_freqs(12)
           real              :: chan_err(12)
           real              :: lwp_cutoff
           real              :: missing
           character(len=90):: spares  !CHANGE!!   !orbit hdr total =600 bytes
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
           real       :: Sclat
           real       :: Sclon
           real       :: Scalt
           ! need to add scorient field to use AMSR2 as well!!!!
       end type ScanHdr

!--- input pixel structure

       type :: DataRec
           real       :: latitude
           real       :: longitude
           real       :: Tb(12)
           real       :: EIA(12)
           real       :: azimuth
           real       :: LOFlag ! new!
       end type DataRec

!--- assign structure names
       
       type(fileHdr)  :: fhdr
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: pix

!--- open input file

       call gprof_lun(lun)
       open(unit=lun,file=trim(input_file),access='stream',
     >      status='old',iostat = ios)
        if(ios .ne. 0) then
          write(*,*),' Error opening pp file'
          stop
        endif
!--- read file header
       read(lun, iostat=ios)  fhdr           !read out file header
        if(ios.ne.0) write(*,*),'error reading file header'


!--- read orbit header

       read(lun, iostat=ios)  orbhdr               !read out orbit header

!--- move orbit header information to program variables

       orig_file    = orbhdr.origfile
       ppversion  = orbhdr.pp_version
       satellite   = orbhdr.satellite ! get satellite, sensor, pp_ver
       sensor      = orbhdr.sensor    !   from pp file
!       nchans      = orbhdr.nchans
       granule     = orbhdr.granule
       nscans      = orbhdr.nscans
       npix        = orbhdr.npixels
       chan_freq   = orbhdr.chan_freqs
       write(granny,'(i06)') (orbhdr.granule+100000)
!       chan_err    = orbhdr.chan_err
!       lwp_cutoff  = orbhdr.lwp_cutoff
       !miss_flt    = orbhdr.missing

!--- set channel availability flag

       do i = 1, 12 ! -- pp now uses 12 channels
         if(chan_freq(i) .eq. miss_flt) then
              chan_avail(i) = .false.
         else
              chan_avail(i) = .true.
         endif
       enddo


!--- check if there are any scans to process - if not then exit

       if(nscans .eq. 0) write(*,*)'nscans is 0!'

!--- allocate memory for pixel level data

       allocate (stdtime(nscans,6))
       allocate (sclat(nscans), sclon(nscans), scalt(nscans))       
        ! need to get scorient, scad somehow for AMSR2!!!!
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
       allocate (eia_out(npix,nscans,2)) ! just one for all AM2 freqs
       allocate (save_iter(npix,nscans)) 
       allocate (sataz(npix,nscans)) ! azimuthal angle, rel to N
       allocate (sglint(npix,nscans),qualflag(npix,nscans)) 
              
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
            idate = stdtime(iscan,1)*10000 + stdtime(iscan,2)*100 + 
     >              stdtime(iscan,3)     
            if(idate .gt. 19870101 .and. idate.lt.20500101) then
               first_good_date = idate
            endif       
         endif

         sclat(iscan) = scnhdr.sclat ! currently all miss_flt from pp!
         sclon(iscan) = scnhdr.sclon
         scalt(iscan) = scnhdr.scalt
         scorient(iscan) = 0 ! or 180 for backwards, but AM2 always fwd!
         
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
           sglint(ipix,iscan)       = miss_byt ! no sunglint info in L1R!
           ! sun azimuth & elevation are stored, but no explicit
           ! sunglint information. so, set to missing here.
                                 
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
       open(unit=db_lun,file='binary/landmask56-32_16.bin',
     >      form='unformatted', status='old',recl = irecl, 
     >      access='direct',iostat=ios)
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
!       read_date = first_good_date       
       !write(*,*) '   read_date = ', read_date
        
       do i = 1,30                !search back 30 days for file
       
         !write(cdate,'(i8)') read_date
         !read(cdate(1:4),'(i4)') iyear
         !read(cdate(5:6),'(i2)') imon
         !read(cdate(7:8),'(i2)') iday
             !print*,cdate,iyear,imon,iday 
!         if(read_date.lt.20020601 .or. read_date.gt.20111004) then
!AVHRR+AMSR starts this date               
!               tfile = trim(dir_sst) // cdate(1:4) // 
!     >                 '/avhrr-only-v2.' // cdate
               tfile = trim(dir_sst) // alldate(1:4) // 
     >                 '/avhrr-only-v2.' // alldate
               inquire(file = tfile,exist = fexists)     
!         else               
!               tfile = trim(dir_sst) // cdate(1:4) // 
!     >                 '/amsr-avhrr-v2.' // cdate
!               inquire(file = tfile,exist = fexists)
!         endif 
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
       open(rlun,file=trim(tfile), form='unformatted',
     >      status='old', access='sequential',convert='big_endian',
     >      iostat=ios)
       if(ios.ne.0) write(*,*)'error opening SST file'

!--- read in SST and ICE grids
       
       read(rlun,iostat=ios)iyr,imo,ida,
     >           ((fsst(i,j),i=1,xsize),j=ysize,1,-1) !flip SST for NP up
       !write(*,*)'writing out sst grid!: ',fsst
       read(rlun,iostat=ios)iyr,imo,ida,
     >           ((anom,i=1,xsize),j=1,ysize)  !don't use anom field
       read(rlun,iostat=ios)iyr,imo,ida,
     >           ((err,i=1,xsize),j=1,ysize)  !don't use err field
       if(ios.ne.0) write(*,*)'issue 3' 
       read(rlun,iostat=ios)iyr,imo,ida,
     >           ((fice(i,j),i=1,xsize),j=ysize,1,-1) !flip ICE for NP up
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
            if(trim(sensor).eq.'GMI') lo_flag(pix,lin) = -1 ! NEW
          else   
 
            if(trim(sensor).eq.'GMI') then 
              lo_flag(pix,lin)=lsmask(ilon,ilat) ! for non-AMSR

              if(lsmask(ilon,ilat) .ge. 0 .and.   !0-5% land = ocean -- CHANGE!
     >           lsmask(ilon,ilat) .le. 5) then   ! changed from 2% cutoff!!
                 sfc_type(pix,lin) = 1 
              elseif(lsmask(ilon,ilat) .gt. 5) then
                 sfc_type(pix,lin) = 0 
              endif
            endif                                
            if(trim(sensor).eq.'AMSR2') then ! 0-2% land = ocean for AM2!
              if(lo_flag(pix,lin).ge.0 .and. lo_flag(pix,lin).le.3) then
                sfc_type(pix,lin)=1
              elseif(lo_flag(pix,lin).gt.3) then 
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
          if ((iclat .lt. 1) .or. (iclat .gt. imaskres*180) .or. !boundary check 
     >        (iclon .lt. 1) .or. (iclon .gt. imaskres*360))  then
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

!      do lin = 1, nscans
!        do pix = 1, npix
                !write(*,*)'sfc_type(pix,lin)=',sfc_type(pix,lin)
!          if(sfc_type(pix,lin) .ne. 1) then
!                 sst(pix,lin) = miss_flt   ! 20=land, 30=coast, ice handled below
                 !write(*,*)'sst(pix,lin)=',sst(pix,lin)
!                 cycle      
!          endif
!          alat = lat(pix,lin)            !input lat
!          alon = lon(pix,lin)            !input lon centered on 0 deg
          alat = latitude
          alon = longitude
          alon=alon + 180.0              !shift for grid counter 

          islat = nint(smaskres*(90.+rincs-alat)) !assign sst array elements
          if(islat .eq. 181*smaskres) islat = 180*smaskres !chk forsouthpole
          islon = nint(smaskres*(rincs+alon))

!          if(islon .eq. (smaskres*360+1)) islon = 1
!          if ((islat .lt. 1) .or. (islat .gt. smaskres*180) .or. !boundary check 
!     >        (islon .lt. 1) .or. (islon .gt. smaskres*360))  then
!              sst(pix,lin) = miss_flt
!              cycle
!          else
!              if(sstgrid(islon,islat) .lt. -1.8) then ! why this cutoff?
!                  sst(pix,lin) = miss_flt
!                  sfc_type(pix,lin) = 0 !31  !assign
!                  cycle
!              else
!                  sst(pix,lin) = sstgrid(islon,islat) + 273.15 
                  psst = sstgrid(islon,islat) + 273.15 
!              endif
!          endif
              
!        enddo   !pix
!      enddo     !lin

      deallocate (sstgrid)       !free up the sstgrid

      return
      end subroutine assign_sst


      subroutine prep_oe
        implicit none

        real :: freq,temp,pres,rho
        integer :: bin,era_lun1,era_lun2,era_lun3,rsslun,ios,syl,i,j,k
        character(len=2) :: charnch 
        character(len=3) :: ver='vW' ! LUT version
        character(len=4) :: yr
        character(len=2) :: mos(12) = (/'1','2','3','4','5','6',
     >                                '7','8','9','10','11','12'/)
        character(len=2) :: das(31) = (/'1','2','3','4','5','6',
     >     '7','8','9','10','11','12','13','14','15','16','17',
     >   '18','19','20','21','22','23','24','25','26','27','28',
     >   '29','30','31'/)
        character(len=2) :: mo, da
        character(len=3) :: emis
        character(len=1) :: spc
        character(len=100):: sydir
        logical :: fexists
        sydir = '/tdata1/dduncan/oe/sy/'

      !allocate(oe_output(npix,nscans,9),screen(npix,nscans)) !last dimension hard-coded!
      allocate(toffsets(nch),s_toffsets(nbins,nch),s_sy(nch,nch,nbins))
      allocate(sy(nch,nch),sy_i(nch,nch)) ! local currently, with nch in name.

      if(nch.ge.10) write(charnch,'(I2)'), nch !
      if(nch.le.9)  write(charnch,'(I1)'), nch !
      write(spc,'(I1)') npc
        !spc = '3' ! for testing
      call gprof_lun(syl)
!      changed Sy reading 7/21/15 DD, again 8/27, again 10/04 (just
!      directory call instead of within current dir)

      emis = 'F' ! default
      if(run_rss) emis = 'R'
      if(fastem4) emis = 'F4'
      if(useanalysis) emis = trim(emis) // 'a'
      if(useanalysis.ne..true.) emis = trim(emis) // 'c' !climatology
        !emis='Fc' !OVERRIDE FOR TESTING!
        
      !sy_file = trim(sydir)//trim(sensor)// '.SY.v2d.'
      sy_file = trim(sydir)//trim(sensor)// '.SY.'//trim(emis)// 
     >  '.'//trim(spc)//'.'//trim(charnch)//'.v5.bin'
      if(trim(sensor).eq.'AMSR2') sy_file=trim(sydir)//trim(sensor)//
     > '.SY.'//trim(emis)//'.'//trim(spc)//'.'//trim(charnch)//'.v4.bin'
        inquire(file = sy_file,exist = fexists)
       if(.not.fexists) then
         print*,'SY file doesnt exist, sorry!'
         stop
       endif 
        !print*,sy_file
      open(unit=syl,file=trim(sy_file),
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open sy_file',sy_file
        read(syl) sy
        read(syl) toffsets
        read(syl) s_sy
        read(syl) s_toffsets
      close(syl)

      allocate(mr_sigmas(nbins,3))!npc)) ! simpler!
      allocate(mr_mins(nbins,3),mr_maxes(nbins,3))

      call gprof_lun(syl)
      open(unit=syl,file='binary/mreof_limits.' //
     >  trim(spc) // '.bin',
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open MREOF analysis file'
        read(syl) mr_sigmas
        read(syl) mr_mins
        read(syl) mr_maxes
      close(syl)

      call gprof_lun(syl)
      open(unit=syl,file= 'binary/SA-EOF_LW.' // trim(spc) // '.bin',
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open eof-lw cov file'
        read(syl) eof_lw_cov ! 3 x nbins
      close(syl)
     
!--- allocate and read in EC LUT data
      allocate(era_wm(nlon,nlat), era_ws(nlon,nlat))
      allocate(mprof(nbins,nz),eofs(nbins,nz,6))!npc))
      allocate(peofs(nz,6))

      !write(yr,'(I4)') stdtime(5,1)
      !mo = mos(stdtime(5,2)) !use as index
      !da = das(stdtime(5,3)) !use as index
      write(yr,'(I4)') year
      mo = mos(month) !use as index
      da = das(day) !use as index

      erafile = 'binary/eof_mr.03.v3.bin' ! v3 has consistency across bins
      call gprof_lun(era_lun1)
      open(unit=era_lun1,file=erafile,access='stream',
     >     status='old',form='unformatted',iostat = ios)
       if(ios .ne. 0) write(*,*)' cant open EOF LUT file'
        read(era_lun1) mprof
        read(era_lun1) eofs
      close(era_lun1)

      erafile ='/tdata1/dduncan/oe/LUT/ws_n128x4_' // trim(ver) //
     >    '_' // trim(mo) // '.' // trim(yr) // '.bin'
       call gprof_lun(era_lun1)
      open(unit=era_lun1,file=erafile,access='stream',
     >     status='old',form='unformatted',iostat = ios)
       if(ios .ne. 0) write(*,*)' cant open EC winds LUT file'
        read(era_lun1) era_wm
        read(era_lun1) era_ws
      close(era_lun1)

      end subroutine prep_oe

      subroutine qual_czech
       use define_oeorbit
       implicit none
!----- new 11/17/15, D Duncan
!      this will screen out anomalously high TPW
!       pixels in which it is likely raining lightly or LWP
!       is being traded off for more WV.

       integer :: lin,pix,ct,s1,p1
       integer,parameter :: sd=3,pd=5
       real    :: sigmatp,meantp
       real    :: meth=4.0,sith=5.5,tottp,eachtp((2*pd+1)*(2*sd+1))
       real    :: temptp(2*pd+1,2*sd+1), tempch(2*pd+1,2*sd+1)
       
       screen(:,:) = oe_output(:,:,1)
       do lin = 1+sd,nscans-sd
         do pix = 1+pd,npix-pd
           if(oe_output(pix,lin,1).lt.0) cycle
           temptp(:,:) = oe_output((pix-pd):(pix+pd),
     >       (lin-sd):(lin+sd),1)
           tempch(:,:) = oe_output((pix-pd):(pix+pd),
     >       (lin-sd):(lin+sd),5)
           ct=0
           tottp=0.0
           do s1=1,2*sd+1
             do p1=1,2*pd+1
               if(temptp(p1,s1).gt.0.and.tempch(p1,s1).lt.chis1
     >            .and.temptp(p1,s1).ne.oe_output(pix,lin,1)) then
                 ct=ct+1
                 tottp=tottp+temptp(p1,s1)
                 eachtp(ct) = temptp(p1,s1)
               endif
             enddo
           enddo
           if(ct.lt.15) cycle ! need enough to get real background sense
           meantp = tottp/real(ct)
           sigmatp = SQRT(SUM((eachtp(1:ct)-meantp)**2)/real(ct-1))
           if(oe_output(pix,lin,1).gt.(meantp+meth) .or.
     >        sigmatp.gt.sith) then
             screen(pix,lin) = -996
           endif 
         enddo
       enddo
       oe_output(:,:,1) = screen(:,:)

      end subroutine qual_czech

      subroutine output_oe
       use define_oeorbit
       implicit none

       integer   :: ios, ou_lun, lin, pix, a=1, b=1,c, cha, cc
       integer   :: temptime(8)
!--- create ouput structures
       type :: Date         
          integer(2):: year
          integer(2):: month
          integer(2):: day
          integer(2):: hour
          integer(2):: minute
          integer(2):: second
       end type Date    ! 12 bytes
       
       type :: OrbitHdr         ! same header for this AND RT version
          character(len=10) :: Satellite
          character(len=10) :: Sensor
          character(len=10) :: Algver
          character(len=20) :: Ppver
          character(len=100):: OrigFile    !150
          character(len=100):: CalFile     !250
          character(len=100):: SpecFile     !350
          character(len=100):: SyFile     !450
          character(len=6)  :: Freqs(maxchans)   !+90 540
          integer(2)        :: Avail(maxchans)   !+30 570
          integer(1)        :: Nchannels         !+1  571
          integer(1)        :: NPCs              !+1  572
          integer(4)        :: Granule           !+4  576
          integer(2)        :: NEIA              !+2  578
          integer(2)        :: NLayers           !+2  580
          real              :: PLevels(nz+1)     !+4*17 648

          integer(2)        :: NScans           !+2 650
          integer(2)        :: NPixels          !+2 652

          real              :: Avgiter          !+4 656
          real              :: MissingData      !+4 660
          real              :: ChiSqThresh      !+4 664
          type(Date)        :: CreationDate     !+12 676
       end type OrbitHdr
       
       type :: ScanHdr
          type(Date) :: ScanDate
          real       :: Sclat 
          real       :: Sclon
          real       :: Scalt
       end type ScanHdr

       ! reminder of variable order: EOF coeff 1-3, WIND, LWP, IWP, SST.
       !  output order of vars: TPW, LWP, WIND, IWP, ChiSq, SST, lapse
       !  EOF coeff 1-3
       type :: DataRec
          real       :: latitude
          real       :: longitude
          
          real       :: TPW
          real       :: LWP
          real       :: WIND
          real       :: IWP
          real       :: CHISQ
          real       :: SST
          real       :: ReySST

          integer(2) :: EOFC(3) ! EOF coefficients -- 100*real
          integer(2) :: WindDir ! Deg CW from N (0-360)

          real       :: Tbdiff(maxchans) ! extras are left missing/blank
          real       :: posteriori(nvar+1) ! now 8!

          integer(1) :: NIterations ! new for BCMX
          integer(1) :: SunGlintAngle
          integer(1) :: LandOceanPct
          integer(1) :: QualityFlag
       end type DataRec
       
!--- assign structure names
       
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: rec
      
!---open output file

       call gprof_lun(ou_lun)
       open(unit=ou_lun,file=trim(output_file),access='stream',
     >      status='unknown',iostat = ios)
       if(ios .ne. 0) write(*,*)' problem opening output file'

       orbhdr%satellite           = trim(satellite)
       orbhdr%sensor              = trim(sensor)
       orbhdr%algver              = trim(code_version)
       orbhdr%ppver               = trim(ppversion)
       orbhdr%origfile            = trim(orig_file)
       orbhdr%calfile             = trim(cal_file)
       orbhdr%specfile            = trim(spec_file)
       orbhdr%syfile              = trim(sy_file)
       orbhdr%freqs               = freqs(:) 
       orbhdr%avail               = avail(:)
       orbhdr%nchannels           = nch
       orbhdr%npcs                = npc
       orbhdr%granule             = granule
       orbhdr%neia                = nEIA
       orbhdr%nlayers             = nz
       orbhdr%plevels             = lpress(:)
       orbhdr%nscans              = nscans
       orbhdr%npixels             = npix
       orbhdr%avgiter             = avg_iter
       orbhdr%missingdata         = miss_flt !-9999.9
       orbhdr%chisqthresh         = chisq_out_thresh
       call date_and_time(values=temptime)
       orbhdr%CreationDate%year   = temptime(1)
       orbhdr%CreationDate%month  = temptime(2)
       orbhdr%CreationDate%day    = temptime(3)
       orbhdr%CreationDate%hour   = temptime(5)
       orbhdr%CreationDate%minute = temptime(6)
       orbhdr%CreationDate%second = temptime(7)

!--- write out orbit header info and cluster info

       write(ou_lun, iostat=ios)  orbhdr
       if(ios .ne. 0) write(*,*)' error writing orb hdr'

!--- loop over all scans

       do lin = 1, nscans         !outputs all scans in file

         scnhdr%scandate%year   = stdtime(lin,1)
         scnhdr%scandate%month  = stdtime(lin,2) 
         scnhdr%scandate%day    = stdtime(lin,3)
         scnhdr%scandate%hour   = stdtime(lin,4)
         scnhdr%scandate%minute = stdtime(lin,5)
         scnhdr%scandate%second = stdtime(lin,6)
        
         scnhdr%sclat = sclat(lin) ! these are all miss_flt now for AMSR!
         scnhdr%sclon = sclon(lin) ! -- issue with the preprocessor.
         scnhdr%scalt = scalt(lin)

         write(ou_lun, iostat=ios)  scnhdr    !write out scan header
         if(ios .ne. 0) write(*,*)' error writing scan header'
       
       ! reminder of variable order: EOF coeff 1-3, WIND, LWP, IWP, SST.
       !  output order of vars: TPW, LWP, WIND, IWP, ChiSq, SST, 
       !  EOF coeff 1-3

         do pix = 1, npix            !loop over all pixels
                
           rec%latitude           = lat(pix,lin)
           rec%longitude          = lon(pix,lin)
           
           rec%TPW                = oe_output(pix,lin,1) 
           rec%LWP                = oe_output(pix,lin,2)
           ! if OE's LWP is damn near zero, call it 0!
           if(oe_output(pix,lin,2) .gt. 0 .and. 
     >        oe_output(pix,lin,2) .lt. 0.001) rec%LWP = 0.0
           rec%WIND               = oe_output(pix,lin,3)
           rec%IWP                = oe_output(pix,lin,4)
           rec%CHISQ              = oe_output(pix,lin,5)
           rec%SST                = oe_output(pix,lin,6)
           rec%ReySST             = sst(pix,lin)
           if(rec%TPW.gt.0) then
             rec%EOFC(1:npc)      = 
     >          nint(100*oe_output(pix,lin,7:(6+npc)))
             if(npc.lt.3) rec%EOFC((npc+1):3) = miss_int
           else
             rec%EOFC(:) = miss_int
           endif
           rec%WindDir       = nint(save_wdir(pix,lin)) ! Deg CW from N (0-360)

           cc = 1
           do cha = 1, maxchans !nch
             if(avail(cha).eq.1) then
               rec%Tbdiff(cha) = Tb_diff(pix,lin,cc)
               cc = cc+1
             else
               rec%Tbdiff(cha) = miss_flt
             endif
           enddo

           rec%posteriori(:)     = poster(pix,lin,:)

           rec%NIterations       = save_iter(pix,lin)
           rec%SunGlintAngle     = sglint(pix,lin)
           rec%LandOceanPct      = lo_flag(pix,lin)

           rec%QualityFlag  = 3 ! default -- no convergence,other
           ! Not run or no convergence = 3
           if(rec%CHISQ.le.chis1 .and. rec%TPW.gt.0) rec%QualityFlag=0
           ! High Quality = 0
           if(rec%CHISQ.gt.chis1 .and. rec%TPW.gt.0) rec%QualityFlag=1
           ! Lower Quality = 1
           if(rec%TPW.eq.-996) then
           ! Screened out for TPW irregularity (LW/WV issue)
             rec%QualityFlag = 4
             rec%TPW         = miss_flt
           endif
           if(rec%TPW.eq.-998) rec%QualityFlag=5 ! Bad pixel, Land/Ice, not run
           if(sglint(pix,lin).lt.20 .and. sglint(pix,lin).ge.0 !-99=no sun in L1C
     >        .and. trim(sensor).eq.'GMI') then
           ! Sun Glint Issue = 2
             rec%QualityFlag   = 2
             rec%TPW           = miss_flt        
             rec%LWP           = miss_flt        
             rec%WIND          = miss_flt        
             rec%IWP           = miss_flt        
             rec%SST           = miss_flt        
             rec%CHISQ         = miss_flt        
           endif

           write(ou_lun, iostat=ios)  rec    !write output data structure
           if(ios .ne. 0) write(*,*) ' error writing data out'
         enddo  !pix
       enddo   !lin
                 
       close(ou_lun) 
       
! --- Deallocate variables used in output procedure
       deallocate(sst,lat,lon)
       deallocate(oe_output,screen)
       deallocate(stdtime)
       deallocate(sclat,sclon,scalt,scorient,scad)
       deallocate(Tbb,sat_eia,Tb_diff,Tb_sim,eia_out,poster)
       deallocate(save_tprof,save_mrprof,save_slp,save_wdir)
       deallocate(save_clwc,save_ciwc)
       deallocate(sfc_type,lo_flag)
       deallocate(save_iter,sataz)
       deallocate(sglint)
       deallocate(era_wm,era_ws)
       deallocate(eofs,mprof,peofs)
       deallocate(oe_tbs,test_out)
       deallocate(sy,sy_i,s_sy,toffsets,s_toffsets)

       return
      end subroutine output_oe


      subroutine output_oe_rt
       use define_oeorbit
       implicit none

       integer   :: ios, ou_lun, lin, pix, a=1, b=1,c, cha, cc
       integer   :: beg, beg2, beg3
       integer   :: temptime(8)
       character(120) :: rt_outfile,inbetween1,inbetween2
       logical           :: backsurch=.true.

!--- create ouput structures
       type :: Date         
          integer(2):: year
          integer(2):: month
          integer(2):: day
          integer(2):: hour
          integer(2):: minute
          integer(2):: second
       end type Date    ! 12 bytes
       
       type :: OrbitHdr         
          character(len=10) :: Satellite
          character(len=10) :: Sensor
          character(len=10) :: Algver
          character(len=20) :: Ppver
          character(len=100):: OrigFile    !150
          character(len=100):: CalFile     !250
          character(len=100):: SpecFile     !350
          character(len=100):: SyFile     !450
          character(len=6)  :: Freqs(maxchans)   !+90 540
          integer(2)        :: Avail(maxchans)   !+30 570
          integer(1)        :: Nchannels         !+1  571
          integer(1)        :: NPCs              !+1  572
          integer(4)        :: Granule           !+4  576
          integer(2)        :: NEIA              !+2  578
          integer(2)        :: NLayers           !+2  580
          real              :: PLevels(nz+1)     !+4*17 648

          integer(2)        :: NScans           !+2 650
          integer(2)        :: NPixels          !+2 652

          real              :: Avgiter          !+4 656
          real              :: MissingData      !+4 660
          real              :: ChiSqThresh      !+4 664
          type(Date)        :: CreationDate     !+12 676
       end type OrbitHdr
       
       type :: ScanHdr
          type(Date) :: ScanDate ! given as orbit start time for AM2 (pp issue)
          real       :: Sclat ! missing for AM2
          real       :: Sclon ! missing for AM2
          real       :: Scalt
          real       :: Scorient ! meaningless for AMSR2, but eh
       end type ScanHdr

       ! reminder of variable order: EOF coeff 1-3, WIND, LWP, IWP, SST.
       !  output order of vars: TPW, LWP, WIND, IWP, ChiSq, SST, 
       !  EOF coeff 1-3
       type :: DataRec
          real       :: latitude
          real       :: longitude

          real       :: SLP ! in mbar
          real       :: WindDir ! 0-360 
          
          real       :: TPW
          real       :: LWP
          real       :: WIND  
          real       :: IWP
          real       :: CHISQ
          real       :: SST
          real       :: REYSST

          ! TbDiff, L1CR, TbSim
          real       :: TbObs(maxchans) ! extras are left missing/blank
          real       :: TbSim(maxchans) ! extras are left missing/blank
          real       :: TbDiff(maxchans) ! extras are left missing/blank
          real       :: eiaout(2)
          real       :: MRProf(nz)
          real       :: TProf(nz)
          real       :: CLWC(nz)
          real       :: CIWC(nz)

          integer(1) :: NIterations ! new for BCMX
          integer(1) :: SunGlintAngle
          integer(1) :: LandOceanPct
          integer(1) :: QualityFlag

       end type DataRec
       
!--- assign structure names
       
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: rec
      

       beg = scan(output_file,'/',backsurch) !+ 1
       inbetween1 = trim(output_file(1:beg))//'RT/'
       beg2 = scan(output_file,'BIN',backsurch) !+ 1
       inbetween2 = trim(output_file(beg+1:beg2+3))
       rt_outfile = trim(inbetween1)//trim(inbetween2)
       beg3 = scan(rt_outfile,'.BIN',backsurch) -3 !+ 1
       rt_outfile = trim(rt_outfile(1:beg3)) // trim('RT.BIN')

!---open output file
       call gprof_lun(ou_lun)
       open(unit=ou_lun,file=trim(rt_outfile),
     >      access='stream',status='unknown',iostat = ios)
       if(ios .ne. 0) write(*,*)' problem opening output file'

       orbhdr%satellite           = trim(satellite)
       orbhdr%sensor              = trim(sensor)
       orbhdr%algver              = trim(code_version)
       orbhdr%ppver               = trim(ppversion)
       orbhdr%origfile            = trim(orig_file)
       orbhdr%calfile             = trim(cal_file)
       orbhdr%specfile            = trim(spec_file)
       orbhdr%syfile              = trim(sy_file)
       orbhdr%freqs               = freqs(:) ! trim()?
       orbhdr%avail               = avail(:)
       orbhdr%nchannels           = nch ! byte
       orbhdr%npcs                = npc ! byte
       orbhdr%granule             = granule
       orbhdr%neia                = nEIA
       orbhdr%nlayers             = nz
       orbhdr%plevels             = lpress(:)
       orbhdr%nscans              = nscans
       orbhdr%npixels             = npix
       orbhdr%avgiter             = avg_iter
       orbhdr%missingdata         = miss_flt !-9999.9
       orbhdr%chisqthresh         = chisq_out_thresh
       call date_and_time(values=temptime) ! others are UTC hours, ms
       orbhdr%CreationDate%year   = temptime(1)
       orbhdr%CreationDate%month  = temptime(2)
       orbhdr%CreationDate%day    = temptime(3)
       orbhdr%CreationDate%hour   = temptime(5)
       orbhdr%CreationDate%minute = temptime(6)
       orbhdr%CreationDate%second = temptime(7)


!--- write out orbit header info and cluster info

       write(ou_lun, iostat=ios)  orbhdr
       if(ios .ne. 0) write(*,*)' error writing orb hdr'

!--- loop over all scans

       do lin = 1, nscans         !outputs all scans in file

         scnhdr%scandate%year   = stdtime(lin,1)
         scnhdr%scandate%month  = stdtime(lin,2) 
         scnhdr%scandate%day    = stdtime(lin,3)
         scnhdr%scandate%hour   = stdtime(lin,4)
         scnhdr%scandate%minute = stdtime(lin,5)
         scnhdr%scandate%second = stdtime(lin,6)
        
         scnhdr%sclat = sclat(lin) 
         scnhdr%sclon = sclon(lin)
         scnhdr%scalt = scalt(lin)
         scnhdr%scorient = scorient(lin)

         write(ou_lun, iostat=ios)  scnhdr    !write out scan header
         if(ios .ne. 0) write(*,*)' error writing scan header'
       
       ! reminder of variable order: EOF coeff 1-3, WIND, LWP, IWP, SST.
       !  output order of vars: TPW, LWP, WIND, IWP, ChiSq, SST,
       !  EOF coeff 1-3

         do pix = 1, npix            !loop over all pixels
                
           rec%latitude           = lat(pix,lin)
           rec%longitude          = lon(pix,lin)
           
           rec%SLP                = save_slp(pix,lin)
           rec%WindDir            = save_wdir(pix,lin) ! Deg CW from N (0-360)

           rec%TPW                = oe_output(pix,lin,1) 
           rec%LWP                = oe_output(pix,lin,2)
           ! if OE's LWP is damn near zero, call it 0!
           if(oe_output(pix,lin,2) .gt. 0 .and. 
     >        oe_output(pix,lin,2) .lt. 0.001) rec%LWP = 0.0
           rec%WIND               = oe_output(pix,lin,3)
           rec%IWP                = oe_output(pix,lin,4)
           rec%CHISQ              = oe_output(pix,lin,5)
           rec%SST                = oe_output(pix,lin,6)
           rec%REYSST             = sst(pix,lin)

           cc = 1
           do cha = 1, maxchans !nch
             if(avail(cha).eq.1) then
               rec%TbObs(cha)  = Tbb(pix,lin,cha) ! level 1Tbs! (no offsets applied)
               rec%TbSim(cha)  = Tb_sim(pix,lin,cc) !
               rec%TbDiff(cha) = Tb_diff(pix,lin,cc) ! new for BCMV version
               cc = cc+1
             else
               rec%TbObs(cha)  = miss_flt
               rec%TbSim(cha)  = miss_flt
               rec%TbDiff(cha) = miss_flt
             endif
           enddo

           rec%eiaout(1:neia)    = eia_out(pix,lin,1:neia)
           !rec%posteriori(:)     = poster(pix,lin,:)

           rec%MRProf(:)    = save_mrprof(pix,lin,:)
           rec%TProf(:)     = save_tprof(pix,lin,:)
           rec%CLWC(:)      = save_clwc(pix,lin,:)
           rec%CIWC(:)      = save_ciwc(pix,lin,:)

           rec%NIterations  = save_iter(pix,lin)
           rec%SunGlintAngle= sglint(pix,lin)
           rec%LandOceanPct = lo_flag(pix,lin)

           rec%QualityFlag  = 3 ! default
           if(rec%CHISQ.le.chis1 .and. rec%TPW.gt.0) rec%QualityFlag= 0
           if(rec%CHISQ.gt.chis1 .and. rec%TPW.gt.0) rec%QualityFlag= 1
           if(rec%TPW.eq.-996) then
             rec%QualityFlag = 4
             rec%TPW         = miss_flt
           endif
           ! Screened out for TPW irregularity (LW/WV issue)
           if(sglint(pix,lin).le.20 .and. sglint(pix,lin).ge.0 !-99 = no sun in L1C
     >        .and. trim(sensor).eq.'GMI') then
             rec%QualityFlag  = 2
             rec%TPW   = miss_flt        
             rec%LWP   = miss_flt        
             rec%WIND  = miss_flt        
             rec%IWP   = miss_flt        
             rec%SST   = miss_flt        
             rec%CHISQ = miss_flt        
           endif
  
           write(ou_lun, iostat=ios)  rec    !write output data structure
           if(ios .ne. 0) write(*,*) ' error writing data out'
         enddo  !pix
       enddo   !lin
                 
       close(ou_lun) 
        !print*,'finished outrt'
       
      return

      end subroutine output_oe_rt

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

!----------------------------------------------------------------------
      subroutine read_specs !(spec_file)
      
       implicit none
       
       integer nf, np, ipos
       !real,allocatable    :: fov_xt(:),fov_dt(:)
       !character(len=100):: spec_file
             
       !open(unit=19, file=spec_file, form='formatted', status='old')

       !read(19,*) nch 
       !read(19,'(a3)') satcode
       !read(19,*)  freqs
       !read(19,*)  avail ! channels available 
       !read(19,*)  c_sigma
       !read(19,*)  chisq_out_thresh
       if(trim(sensor).eq.'GMI') then
        freqs(:) = (/'10.65V','10.65H',
     >  '18.70V','18.70H','23.80V','-999.9',
     >  '36.64V','36.64H','89.00V','89.00H','166.0V','166.0H', 
     >  '-999.9','183V+3','183V+7'/)
        if(nch.eq.13) avail(:)=(/1,1,1,1,1,0,1,1,1,1,1,1,0,1,1/)
        if(nch.eq.9)  avail(:)=(/1,1,1,1,1,0,1,1,1,1,0,0,0,0,0/) !no 166/183
        if(nch.eq.7)  avail(:)=(/0,0,1,1,1,0,1,1,1,1,0,0,0,0,0/) !no 10
        if(nch.eq.5)  avail(:)=(/0,0,1,1,1,0,1,1,0,0,0,0,0,0,0/) !no 89
       endif
       if(trim(sensor).eq.'AMSR2') then
        freqs(:) = (/'6.925V','6.925H','10.65V','10.65H',
     >  '18.70V','18.70H','23.80V','23.80H',
     >  '36.50V','36.50H','89.00V','89.00H',
     >  '-999.9','-999.9','-999.9'/) 
        if(nch.eq.12) avail(:)=(/1,1,1,1,1,1,1,1,1,1,1,1,0,0,0/)
        if(nch.eq.10) avail(:)=(/0,0,1,1,1,1,1,1,1,1,1,1,0,0,0/) !no 6
        if(nch.eq.9)  avail(:)=(/0,0,1,1,1,1,1,0,1,1,1,1,0,0,0/) !no 23H
        if(nch.eq.7)  avail(:)=(/0,0,0,0,1,1,1,0,1,1,1,1,0,0,0/) !no 10
        if(nch.eq.5)  avail(:)=(/0,0,0,0,1,1,1,0,1,1,0,0,0,0,0/) !no 10
       endif
       chisq_out_thresh = 4.0
      end subroutine read_specs
!--------------------------------------------------------------------
      subroutine calc_az(slat,ipi,ascdes,sc_orient,rot_ang)

!       calculating azimuthal angle of pixel, 
!        relative to due North, for GMI's footprint

!      slat = spacecraft lat
!      ipi  = position in scan (pixel number, 1-221)   
!      ascdes= ascending/descending pass (1=asc, -1=desc)    
!      sc_orient = spacecraft orientation (180=backward, 0=forward)

      real :: slat,sc_orient
      integer :: nfovhigh, iatt, ascdes, ipi!, badorient=0
      real*8 :: satincl = 65.00000000 ! satellite inclination, relative to E
      real*8 :: scanangle = 140.00
      real*8 :: dpi, d2r, theta, gtrack
      real   :: rot_ang

      nfovhigh= 221 ! # FOVs in high freq chans, thus # pix in L1C/L1CR
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
      theta = ((-0.5d+0*SCANANGLE) +
     >  (SCANANGLE*(ipi-1))/(NFOVHIGH-1.0d+0))
        ! which is different if spacecraft moving backwards!
      if(iatt .eq. -1) then
        theta=((-0.5d+0*scanangle -180.0)
     >   +(scanangle*(ipi-1))/(nfovhigh-1.0d+0))
      endif

        ! contribution of lat on azimuthal angle
      gtrack = ascdes * (90.0d+0 -
     >  asin(cos(satincl*d2r)/cos(abs(slat)*d2r))/d2r)

        ! now find angle between footprint long axis (from sat dir)
        !  and North, clockwise. gotta rotate 90 degrees, too.
      rot_ang = (90.0d+0 - gtrack - theta)
      if(rot_ang .ge. 360.0) rot_ang = rot_ang-360.0

      end subroutine calc_az

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine import_anal
      ! opens correct binary file of EC/MERRA 6-hrly analysis data 
      !  (computed offline) for orbit time

        character(len=40) wdfile, wdfile2
        character(len=180) first, second, third, ecdir
        integer :: bin,era_lun,ios
        integer :: ec_lun
        integer :: hodex, dadex, midscan
        character(len=4) :: yr
        character(len=2) :: mons(12) = (/'01','02','03','04','05','06',
     >                                '07','08','09','10','11','12'/)
        character(len=2) :: mos(12) = (/'1','2','3','4','5','6',
     >                                '7','8','9','10','11','12'/)
        character(len=2) :: das(31) = (/'01','02','03','04','05','06',
     >     '07','08','09','10','11','12','13','14','15','16','17',
     >   '18','19','20','21','22','23','24','25','26','27','28',
     >   '29','30','31'/)
        character(len=2) :: hos(4) = (/'00','06','12','18'/)
        character(len=2) :: money, mo, da, ho
        
        ecdir = '/tdata1/dduncan/oe/ec_bin/'
        if(merra) ecdir='/tdata1/dduncan/oe/merra_bin/'

        midscan = floor(float(nscans)/2.0)+1
      !write(yr,'(i4)') stdtime(midscan,1)
      !mo = mos(stdtime(midscan,2)) !use as index
      !money = mons(stdtime(midscan,2)) !use as index
      !dadex = stdtime(midscan,3)
      !hodex = nint((stdtime(midscan,4)+3.0+stdtime(midscan,5)/60.0)/6.0)
      write(yr,'(i4)') year
      mo = mos(month) !use as index
      money = mons(month) !use as index
      dadex = day
      hodex = nint((hour+3.0+minu/60.0)/6.0)
        if(dadex.le.0.or.hodex.le.0) then
          print*,'EC index value problem!'
          stop
        endif
      if(hodex.eq.5) then
        hodex=1
        dadex = dadex+1
      endif
      if(hodex.eq.0) then
        hodex=4
        dadex = dadex-1
      endif
      ho = hos(hodex)
      da = das(dadex) !use as index
      if(dadex.eq.0 .or. dadex.gt.31) then
        print*,'bad day index val!'
        stop
      endif

      call gprof_lun(ec_lun)
      open(unit=ec_lun,file=trim(ecdir)//trim(yr)//trim(money)//
     >     trim(da)//trim(ho)//'.bin',access='stream',
     >     status='old',form='unformatted',iostat = ios)
       if(ios .ne. 0) write(*,*)' cant open temp_ec file'
        read(ec_lun) winddir
        read(ec_lun) windmag
        read(ec_lun) ecslp
        if(merra.eq..false.)read(ec_lun) eclwp
        if(merra.eq..false.)read(ec_lun) eciwp
        read(ec_lun) ecmr
        read(ec_lun) ect
      close(ec_lun)

      end subroutine import_anal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module subr_oeorbit
