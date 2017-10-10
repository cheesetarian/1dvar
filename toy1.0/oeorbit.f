      program oeorbit
!
!-------------------------------------------------------------------
!       program to run an orbit of Tbs through the OE
!
!   D Duncan, CSU  July 2014
!       last updated January 2016     
!-------------------------------------------------------------------
     
      USE define_oeorbit
      USE subr_oeorbit
      USE nr_subrs 
 
      implicit none

      integer           :: n_arguments, iargc, ios 
      integer           :: nfr 
      real		:: cputime, taz
      integer 		:: cpu1,cpu2, begi
      integer           :: s,p
      character(len=6)  :: extra, extra2
      logical           :: backsearch=.true.

      !lmskfile = 'landmask56-32_16.bin' ! in binary dir

!--- Get command line inputs - <specfile> <infile> <outfile> <caltable> 
!                              <emiss> <runrt>
      n_arguments = iargc() 
      if(n_arguments .gt. 0) then
      !if(n_arguments .eq. 6) then
        print*,'too many arguments!'
        stop
        !call getarg(1, spec_file)   !specification file
        !call getarg(2, input_file)  !pp binary output file
        !call getarg(3, output_file) !binary output file
        !call getarg(4, cal_file)    ! calibration file used in pp
        !call getarg(5, extra)       ! run rss emissivity model?
        !call getarg(6, extra2)      ! extra output for RT model?
      !else    ! not correct number of arguments
!         write(*,*)' error in command line arguments '
!         write(*,*)' <spec_file> <input_file>
!     >               <output_file> <cal_file> <emis_model> <runrt>'
!         stop
      endif

!--- read in the spec file
      !write(*,*)' reading in spec file ',trim(spec_file)
      call read_specs !(spec_file)

!--- echo startup information
!      write(*,*)' satellite/sensor    : ', satcode     
!      write(*,*)' number channels     : ', trim(nchannels)     
!      write(*,*)' input  file name    : ', trim(input_file)
!      write(*,*)' output file name    : ', trim(output_file)     
      
        ! FASTEM6 is default, but F4 or RSS/MW can be selected via call
!        if(trim(extra) .eq. 'R') then 
!          run_rss=.true.
!!          write(*,*)'using RSS emissivities'
!        else
!          run_rss=.false.
!!          write(*,*)'using FASTEM6 emissivities'
!        endif
!        if(trim(extra) .eq. '4') then 
!          fastem4=.true.
!!          write(*,*)'using FASTEM4 emissivities'
!        else
!          fastem4=.false.
!!          write(*,*)'using FASTEM6 emissivities'
!        endif
!        if(trim(extra2) .eq. 'RT') then
!          outrt=.true.
!          !write(*,*)'outputting full RT suite! (too)'
!        else
!          outrt=.false.
!        endif


!--- read pp input file -
      !write(*,*),' starting read pp file'
!!        print*,'satcode: ',satcode
!      if(satcode.eq.'GMI') call read_l1cr_pp(input_file) 
!      if(satcode.eq.'AM2') call read_l1r_pp(input_file) 
!      !write(*,*)' done reading pp file'
!
!      nEIA = 2
!      ascdesc = 'N' ! neither, initialized
!      if(satcode.eq.'AM2'.or.satcode.eq.'AME') then
!        nEIA = 2
!      ! can calculate asc/desc, but dont need to get 'scad' for AMSR
!      ! satellite azimuth angle done in read_l1r_pp
!        begi = scan(output_file,'.BIN',backsearch) -4 !+ 1
!        ascdesc = trim(output_file(begi:begi))
!        !print*,'A/D (should be): ',ascdesc
!      endif

!      if(satcode.eq.'GMI') then ! AM2 has azimuth in L1 files!!
!        ! new section for winds read-in...
!        scad(1) = 1 ! starts ascending, say
!        scad(nscans) = -1 ! ends descending, say
!        do s = 2, nscans-1
!          if((sclat(s+1)-sclat(s-1)).gt.0) scad(s)=1
!          if((sclat(s+1)-sclat(s-1)).lt.0) scad(s)=-1
!          if((sclat(s+1)-sclat(s-1)).eq.0) scad(s)=1 ! ?
!        enddo
!        do s = 1, nscans
!          do p = 1, npix
!            call calc_az(sclat(s),p,scad(s),scorient(s),taz)
!            sataz(p,s) = taz
!          enddo
!        enddo      
!      endif

      !if(useanalysis) then
!       import analysis data for wind direction, etc.
!        write(*,*)' importing analysis data'
        call import_anal
!        write(*,*)' analysis import finished'
      !endif        

!      write(*,*)' reading in reynolds sst'
      ! should just need to define data directory, give a date
      call read_reynolds_sst ! gets date from pp read-in
!      write(*,*)' sst read in'

!      write(*,*)' reading in landmask'
!      call read_landmask  
!      write(*,*)' landmask read in' 

!      write(*,*)' assigning sfc type'
!      call assign_sfc_type 
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
!      call qual_czech

!      if(outrt) then
!        !print*,'writing out full RT file'
!        call output_oe_rt
!      endif 
!      !write(*,*)' outputting data to file'
!      call output_oe
!      !write(*,*)' data have been output'

!        call cpu_time(cputime)
!        cputime = cputime/60.
!        cpu1 = int(cputime)
!        cpu2 = nint((cputime-cpu1)*60)
        !write(*,*)'CPU time (m:s) = ',cpu1,':',cpu2

!--- all done!
      !write(*,*) 'successful completion'

      end
