      program csu1dvar
!
!-------------------------------------------------------------------
!       program to run an orbit of Tbs through the OE
!
!   D Duncan, CSU  July 2014
!       last updated March 2016     
!-------------------------------------------------------------------
     
      USE define_csu1dvar
      USE subr_csu1dvar
      USE optest
      USE output
      USE conv_ncdf2hdfeos5
 
      implicit none

      integer           :: n_arguments, iargc, ios 
      integer           :: nfr 
      real		:: cputime, taz
      integer 		:: cpu1,cpu2!,begi
      integer           :: s,p
      character(len=6)  :: extra, extra2
      logical           :: backsearch=.true.
      character(len=256):: outhe5

      lmskfile = 'landmask56-32_16.bin' ! in binary dir

!--- Get command line inputs - <infile> <outfile> <caltable> 
!                              <emiss> <runrt>
      n_arguments = iargc() 
      if(n_arguments .eq. 7) then
        call getarg(1, satcode)     !3char sensor name
        call getarg(2, n_ch)         !# channels
        call getarg(3, input_file)  !pp binary output file
        call getarg(4, output_file) !binary output file
        call getarg(5, cal_file)    ! calibration file used in pp
        call getarg(6, extra)       ! run rss emissivity model?
        call getarg(7, extra2)      ! extra output for RT model?
      else    ! not correct number of arguments
         write(*,*)' error in command line arguments '
         write(*,*)' <satcode> <nchannels> 
     >               <input_file> <output_file> <cal_file> 
     >               <emis_model> <runrt>'
         stop
      endif

      read(n_ch,'(I2)') nch
!        print*,n_ch,nch

!--- echo startup information
!      write(*,*)' satellite/sensor    : ', satcode     
!      write(*,*)' number channels     : ', nch     
!      write(*,*)' input  file name    : ', trim(input_file)
!      write(*,*)' output file name    : ', trim(output_file)     
      
        ! FASTEM6 is default, but F4 or RSS/MW can be selected via call
        if(trim(extra) .eq. 'R') then 
          run_rss=.true.
!          write(*,*)'using RSS emissivities'
        else
          run_rss=.false.
!          write(*,*)'using FASTEM6 emissivities'
        endif
        if(trim(extra2) .eq. 'RT') then
          outrt=.true.
          !write(*,*)'outputting full RT suite! (too)'
        else
          outrt=.false.
        endif

      nEIA = 2 ! default
      if(satcode.eq.'AM2'.or.satcode.eq.'AME') then
        nEIA = 2
        freqs(:) = (/'6.925V','6.925H','10.65V','10.65H','18.70V',
     >               '18.70H','23.80V','23.80H','36.50V','36.50H',
     >               '89.00V','89.00H','-999.9','-999.9','-999.9'/)
      endif
      if(satcode.eq.'TMI') then
        nEIA = 9
        freqs(:) = (/'10.65V','10.65H','19.35V','19.35H','21.30V',
     >               '-999.9','37.00V','37.00H','85.50V','85.50H',
     >               '-999.9','-999.9','-999.9','-999.9','-999.9'/)
      endif
      if(satcode.eq.'GMI') then
        nEIA = 2
        freqs(:) = (/'10.65V','10.65H','18.70V','18.70H','23.80V',
     >               '-999.9','36.64V','36.64H','89.00V','89.00H',
     >               '166.0V','166.0H','-999.9','183+3V','183+7V'/)
      endif
!--- read pp input file -
!      write(*,*),' starting read pp file'
!        print*,'satcode: ',satcode
      if(satcode.eq.'GMI') call read_l1cr_pp(input_file) 
      if(satcode.eq.'TMI') call read_l1c_pp(input_file) 
      !if(satcode.eq.'AM2') call read_l1r_pp(input_file) 
      if(satcode.eq.'AM2') call read_geos5_pp(input_file) 
!      write(*,*)' done reading pp file'


      if(satcode.eq.'GMI' .or. satcode.eq.'TMI') then ! AM2 has azimuth in L1 files!!
        ! new section for winds read-in...
        scad(1) = 1 ! starts ascending, say
        scad(nscans) = -1 ! ends descending, say
        do s = 2, nscans-1
          if((sclat(s+1)-sclat(s-1)).gt.0) scad(s)=1
          if((sclat(s+1)-sclat(s-1)).lt.0) scad(s)=-1
          if((sclat(s+1)-sclat(s-1)).eq.0) scad(s)=1 ! ?
        enddo
        do s = 1, nscans
          do p = 1, npix
            call calc_az(sclat(s),p,scad(s),scorient(s),taz)
            sataz(p,s) = taz
          enddo
        enddo      
      endif

      !if(useanalysis) then
!       import analysis data for wind direction, etc.
!        write(*,*)' importing analysis data'
        call import_anal
!        write(*,*)' analysis import finished'
      !endif        

!      write(*,*)' reading in aquarius salinity climatology grid'
      call import_sal !
!      write(*,*)' salinity read in'

!      write(*,*)' reading in reynolds sst'
      ! should just need to define data directory, give a date
      call read_reynolds_sst ! gets date from pp read-in
!      write(*,*)' sst read in'

!      write(*,*)' reading in landmask'
      call read_landmask  
!      write(*,*)' landmask read in' 

!      write(*,*)' assigning sfc type'
      call assign_sfc_type 
!      write(*,*)' sfc types assigned' ! 1=100% ocean, 0=others

!      write(*,*)' assigning sst'
      call assign_sst ! sst dir defined in subroutine!
!      write(*,*)' sst assigned'
        
      !write(*,*)'nch,nfreq: ', nch, nfreq

      !write(*,*)' initializing CRTM'
      call crtm_innit
      !write(*,*)' CRTM initialized'
 
      !write(*,*)' calling OE'
      call prep_oe
      !write(*,*)' done prep_oe, starting retrieval'
      call retrieval_nr
      !write(*,*)' OE run'
      call qual_czech

        call cpu_time(cputime)
        cputime = cputime/60.
        cpu1 = int(cputime)
        cpu2 = nint((cputime-cpu1)*60)
        !write(*,*)'CPU time (m:s) = ',cpu1,':',cpu2

          if(dontoutput) then 
            print*,'not outputting! (set in definition file)'
            stop
          endif

      !write(*,*)' writing out netCDF output file'
        outdebug = .false.
      call output_nc(outdebug,outrt)
      !write(*,*)' data have been output'


      outhe5 = trim(output_file)//'.he5'
      if(am2oper) then
       !write(*,*)' converting netCDF to hdfeos5 file'
       call conv(output_file,outhe5)
       !write(*,*)' done conversion'
      endif

      call dealloc
!--- all done!
      !write(*,*) 'successful completion'

      end
