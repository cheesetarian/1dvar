      module pp_commoncode
      use pp_definition

      contains

!--------------------------------------------------------------------------

      subroutine GPM_calibration
        implicit none

!---  this routine reads in the TB calibration coefficents from the calibration
!---  table and calculates and applies the calibration offset in Kelvin to the TBs.
!--   The calibration table is a text file

        character(len=400)  :: line, cal_file_spec
        integer             :: ios, j, iscan,ipix,ichan,clun
        integer             :: m,m1,npts,idum,ndum
        real                :: tmp1,tmp2,tmp
        real                :: cal1,cal2,delta_cal

        type :: struct_cal
           integer          :: npts(10)      ! Number of Calibration tie points 
           real             :: temp(10,50)   ! Temperature values for tie points
           real             :: tcal(10,50)   ! Calibration values for tie points
        end type struct_cal

!--- assign structure names

       type(struct_cal)  :: cal

!--- get LUN for cal file

       call gprof_lun(clun)
       
!--- open inter-calibration table, a text file 

       cal_file_spec = trim(cal_file)
       !cal_file_spec = trim(dir_ancil) // trim(cal_file)
       !cal_file_spec = trim(dir_ancil) // '/' // trim(cal_file)
       
       write(*,*)'  opening cal table : ',trim(cal_file_spec) 
       open(unit=clun,file=trim(cal_file_spec),action="read",
     >        status='old',readonly,iostat = ios)
       if(ios .ne. 0) call reprt_error(141)

!-- Read in TB Calibration coefficents

!-- read all comment lines and test

       ReadComments: do
         read (clun, '(A)',iostat=ios) line
	 if(ios .ne. 0) call reprt_error(142)
         if (line (1:1) /= "#") exit ReadComments
       enddo ReadComments

!-- Need to backup a line since no more comment lines        
       backspace (clun)

!-- Read in the channel number and number of tie points for that channel
       ReadData: do
         read (clun, *,iostat=ios) ichan, npts
        print*,ichan,npts
	 if(ios .ne. 0) call reprt_error(142)
         cal.npts(ichan)=npts
         if (cal.npts(ichan).gt.0) then
             backspace (clun)
             read (clun, *,iostat=ios) idum,ndum,
     >                     (cal.temp(ichan, j), j=1,cal.npts(ichan))
             if(ios .ne. 0) call reprt_error(142)
!             write(*,*) "temp values at ichan : ",ichan,
!     >                   (cal.temp(ichan,j), j=1,cal.npts(ichan))
             read (clun, *,iostat=ios) idum,ndum,
     >                     (cal.tcal(ichan, j), j=1,cal.npts(ichan))
             if(ios .ne. 0) call reprt_error(142)
!             write(*,*) "tcal values at ichan : ",ichan,
!     >                   (cal.tcal(ichan,j), j=1,cal.npts(ichan))
         else
             read (clun, *,iostat=ios) ichan, npts
             if(ios .ne. 0) call reprt_error(142)
         endif
         if (ichan .eq. maxchans) exit ReadData
       enddo ReadData
c       write(*,*) " exit ReadData "
       close (clun)
       
!-- Adjust environmental scene sensor pixels    */

       do iscan = 1, nscans
         do ipix = 1, npix
           do ichan = 1, maxchans
             if(Tb(ipix,iscan,ichan) .ne. miss_flt)  then
               if (cal.npts(ichan) .gt. 0) then

                 tmp  = Tb(ipix,iscan,ichan)
                 tmp1 = cal.temp(ichan,1)               ! setting up first temp
                 tmp2 = cal.temp(ichan,cal.npts(ichan)) ! setting up last  temp
                 cal1 = cal.tcal(ichan,1)               ! setting up first cal
                 cal2 = cal.tcal(ichan,cal.npts(ichan)) ! setting up last  cal

                 if (tmp .le. tmp1) then
                   delta_cal = cal1
                 elseif (tmp .ge. tmp2) then
                   delta_cal = cal2
                 else
                   do m=1,cal.npts(ichan)
                     if (tmp .gt. cal.temp(ichan,m)) then
                       tmp1 = cal.temp(ichan,m)
                       cal1 = cal.tcal(ichan,m)
                       m1 = m
                     endif
                   enddo
                   if ((m1+1) .le. cal.npts(ichan)) then
                     tmp2 = cal.temp(ichan,m1+1)
                     cal2 = cal.tcal(ichan,m1+1)
                   else
                     tmp2 = cal.temp(ichan,m1)
                     cal2 = cal.tcal(ichan,m1)
                   endif
                   delta_cal=cal1+(cal2-cal1)*(tmp-tmp1)/(tmp2-tmp1)
                 endif
                 Tb(ipix,iscan,ichan)=Tb(ipix,iscan,ichan)+delta_cal
               endif
             endif
           enddo
         enddo
       enddo

      return
      end subroutine GPM_calibration

!------------------------------------------------------------------------------------

      subroutine reprt_error(err)
       implicit none
      
      
!-- Subroutine to report errors and terminate execution if necessary.  Upon error,
!-- the output_pp_error routine is called.
                    
      integer  :: err
      character(len=30) :: fe=' FATAL ERR Termination code = '
      character(len=40) :: errmess
      
 8    format(a40,i3,a)         !output format for termination error lines     

!--- select error message depending on err code

      select case(err)

      case(4)
        errmess = ' Unexpected missing channel Tbs'
        write(*,8)fe,err,errmess
      case(5)
        errmess = ' Unable to assign LUN'
        write(*,8)fe,err,errmess	
      case(6)
        errmess = ' Input file does not exist'
        write(*,8)fe,err,errmess	
      case(7)
        errmess = ' Number of scans in granule = 0'
        write(*,8)fe,err,errmess  
      case(8)
        errmess = ' Error in TKopen '      
        write(*,8)fe,err,errmess
      case(9)
        errmess = ' Error in Tkgetmetaint NumberScansGranule'
        write(*,8)fe,err,errmess
      case(10)
        errmess = ' Error in Tkgetmetaint NumberPixels'
        write(*,8)fe,err,errmess
      case(11)
        errmess = ' Error in Tkgetmetaint GranuleNumber'
        write(*,8)fe,err,errmess
      case(12)
        errmess = ' Error in Tkgetmetastring sat info'
        write(*,8)fe,err,errmess
      case(13)
        errmess = ' Error in Tkreadscan'
      case(14)
        errmess = ' Error in Tkclose'
        write(*,8)fe,err,errmess 

      case(141)
        errmess = ' Error opening calibration table'
        write(*,8)fe,err,errmess 
      case(142)
        errmess = ' Error reading cal table header'
        write(*,8)fe,err,errmess 
      case(143)
        errmess = ' Error reading cal table data'
        write(*,8)fe,err,errmess 
             
      case(15)
        errmess = ' Error opening model prep file'
        write(*,8)fe,err,errmess  
      case(16)
        errmess = ' Error retrieving model parameters'
        write(*,8)fe,err,errmess  
      case(17)
        errmess = ' Error retrieving correct model date'
        write(*,8)fe,err,errmess 
      case(171)
        errmess = ' Error setting correct model hour'
        write(*,8)fe,err,errmess 
      case(172)
        errmess = ' Error beyond 27 model hours'
        write(*,8)fe,err,errmess 
      
      case(18)
        errmess = ' Error L1C file has no good lat/lons'
        write(*,8)fe,err,errmess
      case(19)
        errmess = ' Error Model Data is missing'
        write(*,8)fe,err,errmess 
      
      case(21)
        errmess = ' Error reading LandMask in make_surfgrid'
        write(*,8)fe,err,errmess  
      case(22)
        errmess = ' Error reading emissgrid in make_surfgrid'
        write(*,8)fe,err,errmess  
      case(23)
        errmess = ' Error reading autosnow in make_surfgrid'
        write(*,8)fe,err,errmess  
      case(24)
        errmess = ' Error making snow global grid'
        write(*,8)fe,err,errmess
      case(241)
        errmess = ' Error finding NH or SH snow file'
        write(*,8)fe,err,errmess   
      
      case(25)
        errmess = ' Error Lat/Lon OB for sfccode assignment'
        write(*,8)fe,err,errmess
      case(26)
        errmess = ' Error Lat/Lon OB for snowci assignment'
        write(*,8)fe,err,errmess 

      case(28)
        errmess = ' Error calculating orographic lifting'
        write(*,8)fe,err,errmess 

      case(30)
        errmess = ' Error opening output file'
        write(*,8)fe,err,errmess  
      case(31)
        errmess = ' Error writing orbit header'
        write(*,8)fe,err,errmess  
      case(32)
        errmess = ' Error writing scan header'
        write(*,8)fe,err,errmess  
      case(33)
        errmess = ' Error writing pixel data'
        write(*,8)fe,err,errmess  
      
      case(40)
        errmess = ' Error opening error output file'
        write(*,8)fe,err,errmess  
      case(41)
        errmess = ' Error writing error orbit header'
        write(*,8)fe,err,errmess 

      end select
      
      call output_pp_error(errmess)
      
      end subroutine reprt_error

!------------------------------------------------------------------------------
      
      subroutine output_pp_error(errmess)
       implicit none

!--- this output routine gets called from the reprt_error routine if a fatal 
!--- error is identified.  The file header and orbit header (with nscans = 0) 
!--- are written into the output file, and pp program terminates
       
        character(len=40)    :: errmess
        integer              :: ios, begi, ic, olun
        logical              :: backsearch=.true.
    
        type :: OrbitHdr                          !Orbit header = 500 bytes          
           character(len=12) :: satellite
           character(len=12) :: sensor
           character(len=12) :: pp_version
           character(len=128):: radfile
           character(len=128):: prfdbfile
	   character(len=128):: calfile
           integer           :: granule
           integer           :: nscans
           integer           :: npixels              !chans must be in following order:
           integer           :: nchans               !10v,10h,19v,19h,22v,22h,37v,37h,
           real              :: chan_freq(maxchans)  !91v,91h,165v,165h,183_7h,183_3h,183_1h  
           character(len=40) :: comment
        end type OrbitHdr              
 
!--- assign structure names
  
       type(OrbitHdr) :: orbhdr

       write(*,*)' starting output_pp_error'
         
!--- open output file

       call gprof_lun(olun)
       open(unit=olun,file=trim(output_file),access='stream',
     >      status='unknown',iostat = ios)
       if(ios .ne. 0)  then
          write(6,*)' error opening output file in error write'
	  stop 1
       endif

!--- fill up values and write orbit header

       orbhdr%satellite   = sat_name
       orbhdr%sensor      = sensor_name
       orbhdr%pp_version  = pp_version
       begi = scan(radfile,'/',backsearch) + 1
       orbhdr%radfile     = trim(radfile(begi:128))
       orbhdr%prfdbfile   = prfdbfile
       orbhdr%calfile     = cal_file
       orbhdr%granule     = granule
       orbhdr%nscans      = 0
       orbhdr%npixels     = 0
       orbhdr%nchans      = nchans
       orbhdr%chan_freq   = chan_freq
       orbhdr%comment     = trim(errmess)

       write(olun, iostat=ios)  orbhdr                            !write out orbit header
       if(ios .ne. 0)  then
          write(6,*)' error writing orbhdr in error write'
	  stop 1
       endif

!--- close output file   
    
       close(olun)

       write(6,*)'PreProcessor unsuccessful'
       stop 1
       
      end subroutine output_pp_error

!-----------------------------------------------------------------------
      
      subroutine output_pp
       implicit none

!---  this routine will create the 'pre-processor' output file IF there was no
!---  error reading the input file.  If there was, the reprt_error routine
!---  is called and a different output routine:   output_pp_error  is called to
!---  write out the headers only.
       
        integer              :: ios, iscan, ichan,ipix, reccnt=0, begi
        integer              :: ic, olun
        logical              :: backsearch=.true.
       
        type :: OrbitHdr                          !Orbit header = 400 bytes          
           character(len=12) :: satellite
           character(len=12) :: sensor
           character(len=12) :: pp_version
           character(len=128):: radfile
           character(len=128):: prfdbfile
           character(len=128):: calfile
           integer           :: granule
           integer           :: nscans               
           integer           :: npixels              
           integer           :: nchans                   !actual number of channels with data
           real              :: chan_freq(maxchans)   
           character(len=40) :: comment
        end type OrbitHdr              
          
        type :: Date         
           integer(kind=knd2):: year
           integer(kind=knd2):: month
           integer(kind=knd2):: day
           integer(kind=knd2):: hour
           integer(kind=knd2):: minute
           integer(kind=knd2):: second
        end type Date    
       
        type :: ScanHdr
           type(Date) :: ScanDate
           real       :: Sclat
           real       :: Sclon
           real       :: Scalt
           real*8     :: TAI93 ! time (s) since 1/1/1993
        end type ScanHdr

        type :: DataRec
           real               :: latitude
           real               :: longitude
           real               :: Tb(maxchans)     !Tb (K)
           real               :: eia(maxchans)    !earth incidence angle	   
           real               :: Twb              !wet bulb temp (K)
           real               :: tcwv             !model input TPW (mm)
           real               :: skin_temp        !T of skin (K)
           real               :: T2m              !2 meter temperature (K)
           real               :: azimuth          !satellite azimuth angle
           real               :: windmag          !wind magnitude (m/s)
           real               :: slp              !SLP in hPa
           integer*2          :: tprof(nz)        !t profile (K*10)
           integer*2          :: wvprof(nz)       !wv mixratio (g/kg*100)
           integer*2          :: qualflag         !quality flag from L1R	   	   	   
           integer*2          :: winddir          !direction, deg from N  	   	   
           integer*1          :: lo_pct(4)        !L1R land/ocean pct
           integer*1          :: sunglint_angle
           integer*1          :: surface_type_index
           integer*1          :: snow_cover_index
           integer*1          :: orolift_index
        end type DataRec
        
!--- assign structure names
  
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: pix
         
!--- open output file

       call gprof_lun(olun)
       open(unit=olun,file=trim(output_file),access='stream',
     >      status='unknown',iostat = ios)
       if(ios .ne. 0) call reprt_error(30)                             
       
!--- fill up values and write orbit header----------------------------12 

       orbhdr%satellite   = sat_name
       orbhdr%sensor      = sensor_name
       orbhdr%pp_version  = pp_version
       begi = scan(radfile,'/',backsearch) + 1
       orbhdr%radfile     = trim(radfile(begi:128))
       orbhdr%prfdbfile   = prfdbfile
       orbhdr%calfile     = cal_file
       orbhdr%granule     = granule
       orbhdr%nscans      = nscans
       orbhdr%npixels     = npix
       orbhdr%nchans      = nchans
       orbhdr%chan_freq   = chan_freq
       orbhdr%comment = 'preprocessor ok'
 
       write(olun, iostat=ios)  orbhdr                            !write out orbit header
       if(ios .ne. 0) call reprt_error(31)

!--- loop over and write out scans

       do iscan = 1, nscans
       
!--      fill up all scan specific values and write scan header
       
         scnhdr%scandate%year   = stdtime(iscan,1)
         scnhdr%scandate%month  = stdtime(iscan,2) 
         scnhdr%scandate%day    = stdtime(iscan,3)
         scnhdr%scandate%hour   = stdtime(iscan,4)
         scnhdr%scandate%minute = stdtime(iscan,5)
         scnhdr%scandate%second = stdtime(iscan,6)
 
         !scnhdr%TAI93 = tai93time(iscan) ! new 6/9/16
         scnhdr%TAI93 = tai93time(iscan+oscans) ! new 6/9/16
         scnhdr%sclat = sclat(iscan)
         scnhdr%sclon = sclon(iscan)
         scnhdr%scalt = scalt(iscan)
 
         write(olun, iostat=ios)  scnhdr                        !write out scan header
         if(ios .ne. 0) call reprt_error(32)

!--      fill the pixel structure and write out
  
         do ipix = 1, npix

           pix%latitude  = lat(ipix,iscan)          !fill lat/lon
           pix%longitude = lon(ipix,iscan)
           pix%tb        = Tb(ipix,iscan,:)           !Tbs (K)
           pix%eia       = eia(ipix,iscan,:)          !earth incidence angle (deg)
           pix%Twb       = Twb(ipix,iscan)            !Wet Bulb Temp (K)
           pix%tcwv      = tcwv(ipix,iscan)           !TPW (mm)
           pix%skin_temp = skint(ipix,iscan)          !skin Temp (K)
           pix%T2m       = T2m(ipix,iscan)            !2m Temp (K)
           pix%azimuth   = satAZ(ipix,iscan)          !azimuth angle (deg from N)
           pix%windmag   = wind(ipix,iscan)           !wind magnitude (m/s) 
           pix%slp       = slp(ipix,iscan)            !SLP (hPa)
           pix%tprof(:)  = nint(100.0*tprof(ipix,iscan,:))     !T profile (K *100) 
           pix%wvprof(:) = nint(1000.0*wvmr(ipix,iscan,:))   !WV profile (g/kg *1000) 
           pix%qualflag  = 0 !qualflag(ipix,iscan)      !L1C/L1R quality flag
           pix%winddir   = winddir(ipix,iscan)        !wind direction (deg from N) 
           pix%lo_pct(1:4)=int(lo_flag(ipix,iscan,1:4)) ! NEW for L1R version
           pix%sunglint_angle     = sglinta(ipix,iscan) ! now calculated!
           pix%surface_type_index = sfccode(ipix,iscan)
           pix%snow_cover_index   = snowci(ipix,iscan)  
           pix%orolift_index      = 0 !orolifti(ipix,iscan)

           reccnt = reccnt + 1
           write(olun, iostat=ios) pix                   !write out pixel data
           if(ios .ne. 0) call reprt_error(33)

         enddo   !ipix
       enddo  !iscan
       
!--- close output file   
    
       close(olun)
       
D      write(*,*)'   number of records written    : ',reccnt

      return
      end subroutine output_pp
    
!-----------------------------------------------------------------------
	
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
       if(ilun .eq. 201) call reprt_error(5)
       return
      end subroutine gprof_lun
	      
!-----------------------------------------------------------------------

      subroutine julian(iyr,imo,ida,ijul)
       implicit none
!*******************************************************************

! this routine inputs year month and day and calculates the julian
! day of the year. inputs are in integer (ie. 2003, 12, 04)

!  written by dave randel cira/csu
!          in august  1987

*******************************************************************

       integer :: iyr,imo,ida,ijul, mon,irem,irem2,irem3
       integer :: numday(12) =(/31,28,31,30,31,30,31,31,30,31,30,31/)

!**check if leap year

       numday(2) = 28
       irem = mod(iyr,4)
       irem2 = mod(iyr,100)
       irem3 = mod(iyr,400)
       if((irem.eq.0 .and. irem2.ne.0) .or. irem3.eq.0) then
          numday(2) = 29
       endif
       
!**calculate julian day from month and day

       ijul = 0
       if(imo .ne. 1) then
          do mon = 1,imo-1
             ijul = ijul + numday(mon)
          enddo
       endif
       ijul = ijul + ida

      return
      end subroutine julian

      subroutine calc_sunglint

        real*8 :: deg2rad =  0.017453292520
        real*8 :: rad2deg = 57.295779513082

        integer :: i,j
        real*8 :: sza,saz,ei,x1,x2,y1,y2,z1,z2
        real*8 :: sga

        do i = 1, nscans
          do j = 1, npix!l
           ! check for flag values in any of the inputs:
           if(EIA(j,i,1)<-180 .or. sunAZ(j,i)<-180 
     >        .or.sunEL(j,i)<-180) then
            sga = -99.0
           else
            ei = EIA(j,i,1)*deg2rad
            saz = sunAZ(j,i)*deg2rad
            sza = sunEL(j,i) + EIA(j,i,1)
            ! check the solar zenith angle (sza) for the condition of 
            !the apparent visible solar disk fully below the horizon  
            !0.833 is added to 90.0 degrees to allow for the effect
            !of atmospheric refraction and solar disk size.  Reference
            !http://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html  
            if(sza .gt. 90.833) then
              sga = -88.0 ! below horizon -- aka nighttime / okay!
            else
              sza = sza*deg2rad
              x1 = sin(ei)
              y1 = 0.0
              z1 = cos(ei)
              x2 = sin(sza) * cos(saz)
              y2 = sin(sza) * sin(saz)
              z2 = cos(sza)

              ! sun glint angle is defined as the angular separation 
              !  between the reflected view vector and sun direction
              sga = rad2deg*acos(x1*x2 + y1*y2 + z1*z2)
            endif
           endif
           sunglint(j,i) = sga
           sglinta(j,i) = int(sunglint(j,i))
          enddo
        enddo

       end subroutine calc_sunglint

!----------------------------------------------------------------------------------

       end module pp_commoncode
       
