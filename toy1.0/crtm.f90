
        MODULE CRTM

          USE CRTM_MODULE 
          USE define_oeorbit
          USE RSS_RTM

          implicit none

          CONTAINS

        SUBROUTINE CRTM_INNIT

        ! Initialize arrays/elements needed to run CRTM. No need to 
        !  do this more than once! (common vars stored in define subr)

        CHARACTER(20), ALLOCATABLE :: sensor_id(:) !N

        INTEGER :: n_absorbers, n_clouds, n_aerosols
        INTEGER :: alloc_stat, err_stat, n 

        n_sensors = 1
        n_profiles = 1 ! code is set for this to be fixed at 1

        allocate( sensor_id(n_sensors), chinfo(n_sensors), &
                  STAT = alloc_stat )
        IF ( alloc_stat /= 0 ) print*,'problem1'

        ! add lines here for supporting other sensors! (will need to
        !  modify some other areas too though...)
        sensor_id(1) = 'EMPTY' !sensor codes are in appendix B of doc
        if(trim(sensor).eq.'AMSR2') sensor_id = (/'amsr2_gcom-w1'/) 
        if(trim(sensor).eq.'AMSR-E') sensor_id = (/'amsre_aqua'/) 
        if(trim(sensor).eq.'GMI') sensor_id = (/'gmi_gpm'/) 
        if(trim(sensor_id(1)) .eq.'EMPTY') then
          print*,'sensor select error in crtm_innit!'
        endif
          !print*,' sensor_id check!: ',sensor_id

        ! next, call CRTM initialisation function
        if(fastem4) then 
          err_stat = CRTM_Init(sensor_id,chinfo, &
                      File_Path = 'binary/ODAS/', & ! set ODAS/ODPS here!
                      Quiet = .TRUE., & ! suppress output from screen
                      Load_CloudCoeff   = .TRUE., & ! only for scattering!
                      Load_AerosolCoeff = .FALSE., &
                    MWwaterCoeff_File = 'FASTEM4.MWwater.EmisCoeff.bin')
        else ! FASTEM6 is default
          err_stat = CRTM_Init(sensor_id,chinfo, &
                      File_Path = 'binary/ODAS/', & !set ODAS/ODPS here!
                      Quiet = .TRUE., & ! suppress output from screen
                      Load_CloudCoeff   = .TRUE., & ! only for scattering!
                      Load_AerosolCoeff = .FALSE.,& !)
                    MWwaterCoeff_File = 'FASTEM6.MWwater.EmisCoeff.bin')
        endif
        if(err_stat /= SUCCESS) print*,'PROBLEM2!'
       
        ! subsetting channels of sensor.
        if(trim(sensor).eq.'AMSR2' .or. trim(sensor).eq.'AMSR-E') then
         ! dont forget-- AM2 has 6, 7, 10, 18, 23, 37, 89 chans. 14 indices!
         !  AME has 12 indices (no 7.3GHz chan)
         if(nch.eq.12) then 
           nfreeqs = 6
           allocate(freeqs(nfreeqs))
           freeqs=(/6.925,10.65,18.7,23.8,36.5,89.0/)
           if(trim(sensor).eq.'AMSR2') then
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                     Channel_Subset = (/1,2,5,6,7,8,9,10,11,12,13,14/) )
           else 
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1),  &
                     Channel_Subset = (/1,2,3,4,5,6,7,8,9,10,11,12/) )
           endif
         endif
         if(nch.eq.10) then 
           nfreeqs = 5
           allocate(freeqs(nfreeqs))
           freeqs=(/10.65,18.7,23.8,36.5,89.0/)
           if(trim(sensor).eq.'AMSR2') then
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                     Channel_Subset = (/5,6,7,8,9,10,11,12,13,14/) )
           else 
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1),  &
                     Channel_Subset = (/3,4,5,6,7,8,9,10,11,12/) )
           endif
         endif
         if(nch.eq.9) then
           nfreeqs = 5
           allocate(freeqs(nfreeqs))
           freeqs=(/10.65,18.7,23.8,36.5,89.0/)
           if(trim(sensor).eq.'AMSR2') then
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                     Channel_Subset = (/5,6,7,8,9,11,12,13,14/) )
           else 
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1),  &
                     Channel_Subset = (/3,4,5,6,7,9,10,11,12/) )
           endif
         endif
         if(nch.eq.7) then
           nfreeqs = 4
           allocate(freeqs(nfreeqs))
           freeqs=(/10.65,18.7,23.8,36.5/)
           if(trim(sensor).eq.'AMSR2') then
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                   Channel_Subset = (/5,6,7,8,9,11,12/) ) ! no 89s
           else 
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1),  &
                   Channel_Subset = (/3,4,5,6,7, 9,10/) )
           endif
         endif
         if(nch.eq.5) then
           nfreeqs = 3
           allocate(freeqs(nfreeqs))
           freeqs=(/18.7,23.8,36.5/)
           if(trim(sensor).eq.'AMSR2') then
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1),  &
               Channel_Subset = (/7,8,9,11,12/) )
           else 
             err_stat = CRTM_ChannelInfo_Subset( chinfo(1),  &
               Channel_Subset = (/5,6,7, 9,10/) )
           endif
         endif
        endif

        if(trim(sensor).eq.'GMI') then 
         ! 10V/H, 19V/H, 23V, 37V/H, 89V/H, 166V/H, 183V+/-3, 183V+/-7
         if(nch.ge.9) then
                nfreeqs = 5
                allocate(freeqs(nfreeqs))
                freeqs=(/10.65,18.7,23.8,36.64,89.0/)
                err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                         Channel_Subset = (/1,2,3,4,5,6,7,8,9/))
         endif
         if(nch.eq.7) then
                nfreeqs = 4
                allocate(freeqs(nfreeqs))
                freeqs = (/10.65,18.7,23.8,36.64/)
                err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                  Channel_Subset = (/1,2,3,4,5,6,7/) ) ! 10-37!!! no 89!
                  !Channel_Subset = (/3,4,5,6,7,8,9/) )
         endif
         if(nch.eq.5) then 
                nfreeqs = 3
                allocate(freeqs(nfreeqs))
                freeqs = (/18.7,23.8,36.64/)
                err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                  Channel_Subset = (/3,4,5,6,7/) )
         endif
        endif
        if(err_stat /= SUCCESS) print*,'Ch subset fuckup!'

        allocate( geo(n_profiles), opt(n_profiles), &
                  atm(n_profiles), sfc(n_profiles), &
                  STAT = alloc_stat )
        if(alloc_stat /= SUCCESS) print*,'alloc error 1'
           

          ! get the number of channels to process for current sensor
          n_channels = CRTM_ChannelInfo_n_Channels( chinfo(1) )
          allocate(rts(n_channels,n_profiles),stat=alloc_stat)
          if(alloc_stat /= 0) print*,'alloc error 2'

          ! Create atmosphere:
          n_layers = nz ! set in definition module
          n_absorbers = 2 ! ozone is default number 2-- 2 is minimum!
          n_clouds = 1 ! default
          if(includeicecloud) n_clouds=2
          n_aerosols = 0
          Call CRTM_Atmosphere_Create( atm,  &
                                     n_layers,n_absorbers, &
                                     n_clouds,n_aerosols )
          if( any(.not. CRTM_Atmosphere_Associated( atm )) ) then
            print*,'Problem with atm initialization'
          endif
          Call CRTM_RTSolution_Create( rts, n_layers )
          if( any(.not. CRTM_RTSolution_Associated( rts )) ) then
            print*,'Problem with RTS init'
          endif
          ! Allocate the options structures
          Call CRTM_Options_Create( opt,    &
                                    n_channels )
          if(any(.not.CRTM_Options_Associated(opt)))print*,'Bad opt!'

          RETURN
        END SUBROUTINE CRTM_INNIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE RUN_CRTM(tave,wvmratio,lwc,iwc,esst,wsp,tbs_out)

          ! since structures were initialized prior (see above subr),
          !  here we simply read in necessary physical attributes and
          !  call the rt model

          INTEGER :: err_stat,bindex,tindex,i, bicedex,ticedex,fc,cha
          REAL :: tave(nz)
          REAL :: wvmratio(nz), lwp_part, iwp_part, lwc(nz), iwc(nz)
          REAL :: tbs_out(nch)
          REAL :: esst, wsp
          REAL(4) :: etot(2),ezero(2),ewind(2),edirstokes(4),eharm(2,4), rangle ! for RSS model

          ! start by defining elements of profile:
          ! read these in from profile subr (mostly)
          atm(1)%Climatology        = US_STANDARD_ATMOSPHERE 
          ! table 4.5 in documentation--
          ! 'tropical','midlatitude_summer', etc.
          atm(1)%Absorber_Id(1:2)   = (/ H2O_ID , O3_ID /)
          atm(1)%Absorber_Units(1:2)= (/ MASS_MIXING_RATIO_UNITS, &
                                         VOLUME_MIXING_RATIO_UNITS /)
          ! level_pressure gives P at discrete points, i.e. K+1, whereas
          ! normal pressure gives P [hPa] at midpoints, K points
          !atm(1)%Temperature = (/ toa,...,btm_layer /) ! etc.
          atm(1)%Level_Pressure = pplev(:) !ec_plev(:) ! level pressure array
          atm(1)%Pressure = pressave(:) !pave(:) ! average pressure in layer
          atm(1)%Temperature = tave(:)
          ! WV given in mass mixing ratio (g/kg)
          atm(1)%Absorber(:,1) = wvmratio(:)  !wv is absorber 1
          atm(1)%Absorber(:,2) = 0.0 !Ozone is absorber 2-- completely
                                     !unimportant at microwave frequencies

          ! only doing one cloud layer, hence 1 subscript...
          atm(1)%Cloud(1)%Type = Water_Cloud !possibilities in table 4.8
          atm(1)%Cloud(1)%Effective_Radius(:) = cld_effrad 
          ! units microns -- value comes from Han et al (1994)
          atm(1)%Cloud(1)%Water_Content(:) = lwc(:) ! [kg/m^2]
          !atm(1)%Cloud(1)%Water_Content(tindex:bindex) = lwp_part ! [kg/m^2]
        !print*,'clwp: ',atm(1)%Cloud(1)%Water_Content(:)
          if(includeicecloud) then
            atm(1)%Cloud(2)%Type = Ice_Cloud
            atm(1)%Cloud(2)%Water_Content(:) = iwc(:)
            atm(1)%Cloud(2)%Effective_Radius(:)=ice_effrad
          endif 

          ! assign sfc properties -- fractional coverage
          !  these have default values... just define ones that matter!
          !sfc(1)%Land_Coverage = 0.0
          sfc(1)%Water_Coverage = 1.0_fp
          !sfc(1)%Snow_Coverage = 0.0
          !sfc(1)%Ice_Coverage = 0.0
          sfc(1)%Water_Type = 1 ! doesnt matter for MW sensors!
          sfc(1)%Water_Temperature = esst ! [K]
          sfc(1)%Wind_Speed = wsp ! [m/s]
          !if(useanalysis) sfc(1)%Wind_Direction = pwinddir 
          sfc(1)%Wind_Direction = pwinddir 
          sfc(1)%Salinity = 35.0 ! [%o]

          ! assign 'geo' structure properties -- (satellite angles, etc)
          !  read directly from pp file. not assuming eia is static!
          geo(1)%Sensor_Zenith_Angle = EIAs(1) !10V eia
          !geo(1)%Sensor_Zenith_Angle = sat_eia(oepix,oelin,1) !10V eia
          geo(1)%Sensor_Azimuth_Angle = azz
          !geo(1)%Sensor_Azimuth_Angle = sataz(oepix,oelin)

          ! See Section 4.6.4 for help with these optional arguments
          !opt(1)%Use_n_Streams = .TRUE.
          !opt(1)%n_Streams = 4 ! valid values are 2,4,6,8,16 if ^ True
          if(includeicecloud) then
            opt(1)%Include_Scattering = .TRUE.
          else
            opt(1)%Include_Scattering = .FALSE.
          endif

        ! now include ability to run RSS emissivity model !
          if(run_rss) then
            rangle = pwinddir - azz ! NEW-- including dir
            !rangle = pwinddir - sataz(oepix,oelin) ! NEW-- including dir
            if(rangle.lt.0) rangle=rangle+360.0    ! for RSS model
            fc = 1
            opt(1)%Use_Emissivity = .TRUE. ! use user-defined emissivity
            do cha = 1, nfreeqs
              call find_surface_tb( freq=freeqs(cha), &
                        tht=EIAs(1), phir=rangle, &
                        surtep=esst, sal=35.0, ssws=wsp, & ! NEW!
                        e0=ezero, ewind=ewind, edirstokes=edirstokes)
                        !tht=sat_eia(oepix,oelin,1), phir=rangle, &
              etot(1)=ezero(1)+ewind(1)+edirstokes(1)
              etot(2)=ezero(2)+ewind(2)+edirstokes(2)
                if(maxval(etot(:)).gt.1.0.or.minval(etot(:)).lt.0.0) then
                  print*,'BAD EMISSIVITIES!'
                  print*,'emis V/H: ',etot(1:2)
                  print*,'ezeros: ',ezero(1:2)
                  print*,'ewind: ',ewind(1:2)
                  print*,'ewdir: ',edirstokes(1:2)
                  print*,'inputs: ',freeqs(cha),EIAs(1)
                  !print*,'inputs: ',freeqs(cha),sat_eia(oepix,oelin,1)
                  print*,'-- ',rangle,esst,wsp
                  stop
                endif
              if(trim(sensor).eq.'GMI' .and. freeqs(cha).eq.23.8) then
                opt(1)%Emissivity(fc) = etot(1) ! V only! ! wont work for 9ch AMSR
                fc = fc+1
              else if(trim(sensor).eq.'AMSR2'.and.freeqs(cha).eq.23.8 .and.nch.eq.9) then
                opt(1)%Emissivity(fc) = etot(1) ! V only! ! for 9ch AMSR
                fc = fc+1
              else
                opt(1)%Emissivity(fc:fc+1) = etot(1:2)
                fc = fc+2
              endif
            enddo
            !print*,'RSS emissivities: ',opt(1)%Emissivity(:)
          endif

          ! Finally... Call the Fwd model (rts holds all output info)
          !err_stat = CRTM_Forward( atm, sfc, geo, chinfo, rts, &
          !            options=opt)

          if(nch.le.9.or.trim(sensor).eq.'AMSR2'.or.trim(sensor).eq.'AMSRE') then
            err_stat = CRTM_Forward( atm, sfc, geo, chinfo, rts, options=opt)
            if(err_stat /= SUCCESS) print*,'Fwd model error!'
            do i = 1, nch
              tbs_out(i) = real(rts(i,1)%Brightness_Temperature)
            enddo
          endif

          if(nch.gt.9 .and. trim(sensor).eq.'GMI') then
            err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                  Channel_Subset =(/1,2,3,4,5,6,7,8,9/))
            geo(1)%Sensor_Zenith_Angle = EIAs(1)   !10V eia
            !geo(1)%Sensor_Zenith_Angle = sat_eia(oepix,oelin,1)!10V eia
            n_channels = CRTM_ChannelInfo_n_Channels( chinfo(1) )
            deallocate(rts) !new n_ch, so de then re-allocate
            allocate(rts(n_channels,n_profiles),stat=err_stat)
            if(err_stat /= 0) print*,'alloc error, low freq'
            Call CRTM_RTSolution_Create( rts, n_layers )
            if( any(.not. CRTM_RTSolution_Associated( rts )) ) then
              print*,'Problem with RTS init, high freq!'
            endif
            if(run_rss.eq..false.) then 
              err_stat= CRTM_Forward(atm,sfc,geo,chinfo,rts,options=opt)
            endif
            if(run_rss .and. includeicecloud) then 
              !opt(1)%Use_Emissivity = .FALSE. !RSS doesn't work at hi freqs!
              err_stat= CRTM_Forward(atm,sfc,geo,chinfo,rts,options=opt)
            endif
            if(run_rss .and. includeicecloud.eq..false.) then 
              err_stat= CRTM_Forward(atm,sfc,geo,chinfo,rts)
            endif
            if(err_stat /= SUCCESS) print*,'Fwd model error!'

             do i = 1, 9 !set so that if nch>9, all lower freqs are used!
              tbs_out(i) = real(rts(i,1)%Brightness_Temperature)
             enddo
             !need to re-initialize some things !
             if(nch.eq.13) then
              err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                  Channel_Subset =(/10,11,12,13/)) !all
             endif
             if(nch.eq.12) then
              err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                  Channel_Subset =(/10,11,12/)) ! no 183+/-7!
             endif
             if(nch.eq.11) then
              err_stat = CRTM_ChannelInfo_Subset( chinfo(1), &
                  Channel_Subset =(/11,13/)) !166H,183+/-7 !!
             endif
             geo(1)%Sensor_Zenith_Angle = EIAs(2) !hi freq eia
             !geo(1)%Sensor_Zenith_Angle = sat_eia(oepix,oelin,14)!183 eia
             n_channels = CRTM_ChannelInfo_n_Channels( chinfo(1) )
             deallocate(rts) !new n_ch, so de- then re-allocate
             allocate(rts(n_channels,n_profiles),stat=err_stat)
             if(err_stat /= 0) print*,'alloc error, high freq'
             Call CRTM_RTSolution_Create( rts, n_layers )
             if( any(.not. CRTM_RTSolution_Associated( rts )) ) then
               print*,'Problem with RTS init, high freq!'
             endif

             err_stat = CRTM_Forward( atm, sfc, geo, chinfo, rts,options=opt)
             if(err_stat /= SUCCESS) print*,'Fwd model error!'
             do i = 1,nch-9
               tbs_out(i+9) = real(rts(i,1)%Brightness_Temperature)
             enddo

          endif ! gt 9ch,GMI if 
         test_out = tbs_out 
          !print*,'obs tbs: ',oe_tbs(:)
          !print*,'tbs_out: ',tbs_out(:)
          !print*,'emiss: ',rts(:,1)%Surface_Emissivity

        RETURN
        END SUBROUTINE RUN_CRTM

        SUBROUTINE CLEANUP_CRTM

          ! need to destroy and deallocate 
          integer :: err_stat 

          err_stat = CRTM_Destroy(chinfo)
          CALL CRTM_Atmosphere_Destroy(atm)
          deallocate(rts,geo,opt,chinfo,sfc)
          deallocate(freeqs)

        END SUBROUTINE CLEANUP_CRTM

        END module CRTM
