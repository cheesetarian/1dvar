      module optest

      USE define_csu1dvar
      USE CRTM 

      implicit  none

!       non-raining 1DVAR retrieval algorithm for 
!        microwave imagers
!      
!       Last updated 
!       D Duncan, CSU, July 2016
!       
!----------------------------------------------------------------------

      contains

      subroutine retrieval_nr
       
      real         ::  oe_lat, oe_lon, oe_sst
      real         ::  worst
      integer*2    :: lodex,ladex, eclo,ecla,pma,pmo !index for lon,lat
      real         :: eclod,eclad ! diff from grid midpoint and lat/lon geo
      integer      :: lev,z,flag

!     Arrays for Apriori Parameters and associated errors
      real  :: xa(nvar), sa(nvar,nvar), noretrieval 
      real  :: xmin(nretvar), xmax(nretvar)
      
!     Retrieval/Error diagnostics
      real         Amatrix(nvar), sigma(nvar)

!     Counters, etc.
      integer      niter, a, b, c, d, e, mc, i
      integer :: tot_iter=0,noret_count=0,notrun_count=0 ! counters for oe
      real :: lapse
      real :: dmr(nz), tpwerr_lyr(nz,3), tpw_posterr
      real :: extra_ssd, extra_mr(nz), extra_eo(nz,6)
      real :: apecmr(nz)
      real :: offs(maxchans),fssd

      allocate(oe_tbs(nch),test_out(nch))
      allocate(save_inc(npix,nscans,3))
      save_inc(:,:,:) = -997

      ! initialize Sa matrix - set all to zero.
      sa(:,:) = 0.0
      
      !invert Sy matrix here (if static, no need to iterate over it!)
      call inverse(sy,sy_i,nch)

            mc = maxchans
            if(satcode.eq.'AM2') mc = 12 ! change to 14 if using 7GHz too
        
        !initialize oe_output fields, tb diff to -997 (not run!)
            oe_output(:,:,:) = -997
            Tb_diff(:,:,:) = -997
            Tb_sim(:,:,:) = -997
            save_slp(:,:) = -997
            save_wdir(:,:) = -997
            save_tprof(:,:,:) = -997
            save_mrprof(:,:,:) = -997
            save_clwc(:,:,:) = -997
            save_ciwc(:,:,:) = -997
            save_iter(:,:) = 0 
            eia_out(:,:,:) = -997
            save_sal(:,:) = -997

            ! for outputting purposes, save low/hi freq EIA
            if(satcode.eq.'GMI') then
              eia_out(:,:,1) = sat_eia(:,:,1) !low freq
              eia_out(:,:,2) = sat_eia(:,:,14)!hi freq
            else if(satcode.eq.'TMI') then
              do c = 1, 9 
                d = c
                if(c.ge.6) d = c+1 !since 23H skipped
                eia_out(:,:,c) = sat_eia(:,:,d) ! TMI frequencies have different EIAs!
              enddo
            else
              eia_out(:,:,1) = sat_eia(:,:,1) ! same for all AM2 chans?
              eia_out(:,:,2) = sat_eia(:,:,1) ! if using b-scan, change this
            endif

        if(scend.gt.nscans) scend=nscans
        if(pxend.gt.npix)   pxend=npix


!!!!!!!! ---- BEGIN LOOP OVER PIX/SCANS ---- !!!!!!!

            do oelin = scstart, scend 
              !write(*,*),'--   scan = ',oelin
              do oepix = pxstart, pxend
               !print*,'---   pix = ',oepix

               ! define lat,sst,tbs for this pixel before calling OE
               oe_lat = lat(oepix,oelin) 
               oe_lon = lon(oepix,oelin)
               if(oe_lon.lt.lomi.or.oe_lon.gt.loma.or. ! run selected lat/lon box only
     >            oe_lat.lt.lami.or.oe_lat.gt.lama) then ! new in ver1.2
                 !oe_output(oepix,oelin,:) = -998
                 notrun_count= notrun_count+1
                 cycle
               endif
               if(oe_lon.lt.0) oe_lon=lon(oepix,oelin)+360.0
               oe_sst = sst(oepix,oelin)
               ! fix SST, force retrieval for sea ice pixels
               if(sfc_type(oepix,oelin).eq.2) then
                 if(runseaice) then
                   oe_sst=272.5
                 else
                   oe_output(oepix,oelin,:) = -997 ! sea ice marker!
                   notrun_count= notrun_count+1
                   cycle 
                 endif
               endif ! assumes that sfc_type=1 is only other type!

               ! check for 100% ocean, sun glint, SST
               if(sfc_type(oepix,oelin) .eq. 0 .or.
     >            sat_eia(oepix,oelin,1) .lt. 0 .or.
     >            oe_sst .lt. 271 .or.
     >            (sglint(oepix,oelin).lt.20 !don't run sunglint pixels
     >            .and.sglint(oepix,oelin).ge.0) ) then 
                  oe_output(oepix,oelin,:) = -998 ! land/bad pixel marker!
                  notrun_count= notrun_count+1
                  cycle ! high quality oceanic Tbs only!
               endif
               
               ssdex = floor((oe_sst-271.0)/1.0)+1 ! rounds down. 271.9 yields 1, etc.
               if(ssdex.lt.1) ssdex = 1         ! first SST bin if <271K...
               if(ssdex.gt.nbins) ssdex = nbins ! max SST bin

               c = 1
               fssd = ((oe_sst-271.5)+1.0)-real(ssdex)
               if(sum(chan_avail).ne.nch) then
                 print*,'nch NE channels available!'
                 stop
               endif
               do d = 1, mc !max # channels
                 if(chan_avail(d).eq.1) then
                   if((fssd.lt.0.and.ssdex.eq.1) .or. fssd.eq.0 .or.
     >                (fssd.gt.0.and.ssdex.eq.nbins)) then
                     offs(c)=s_toffsets(ssdex,c)
                   endif
                   if(fssd.lt.0) offs(c)=s_toffsets(ssdex,c) +
     >             (s_toffsets(ssdex-1,c)-s_toffsets(ssdex,c))*abs(fssd)
                   if(fssd.gt.0) offs(c)=s_toffsets(ssdex,c) +
     >             (s_toffsets(ssdex+1,c)-s_toffsets(ssdex,c))*abs(fssd)

                  ! all sensors have fwd-model tb offsets applied
                   oe_tbs(c) = Tbb(oepix,oelin,d) + offs(c) 

                   c = c+1
                 endif
               enddo 
               if(c.ne.nch+1) print*,'Tb assignment problem!'
               worst = minval(oe_tbs)
               if(worst .eq. miss_flt) then
                 oe_output(oepix,oelin,:) = -998
                 cycle
               endif
        
        ! Grid starts now at 0E, 90N !!
        lodex = nint(oe_lon/(losize/gm))+1 ! should be exact!!
        ladex = nint(abs(89.4628-oe_lat)/(lasize/gm))+1 ! more complicated...
        if(ladex.le.0 .or. ladex.ge.nlat+1) ladex=1 !just in case...
        !will work up to NPol due to abs, won't quite work at SPol but no matter
          !print*,'ssdex, sst: ',ssdex,oe_sst
          !print*,'lo/la/mo: ',oe_lon,oe_lat,stdtime(oelin,2)    
          !print*,'dexes: ',lodex,ladex    

      ! define pixel wind direction from input EC grid
      eclo = floor(oe_lon/losize)+1  !EC grid index, lon
      ecla = floor(abs(89.4628-oe_lat)/lasize)+1 !EC grid index, lat
      eclod = ( oe_lon - (eclo*losize-losize/2.0) )/losize
      eclad = ( oe_lat - (89.4628-ecla*lasize+lasize/2.0) )/lasize
      pmo = nint(abs(eclod)/eclod)
      pma = -1*nint(abs(eclad)/eclad)  ! just sign
      if((eclo+pmo) .le.0 .or. (eclo+pmo).gt.nlon) pmo=0
      if((ecla+pma) .le.0 .or. (ecla+pma).gt.nlat) pma=0
      pwindmag=windmag(eclo,ecla) +
     >  abs(eclod)*(windmag(eclo+pmo,ecla)-windmag(eclo,ecla)) +
     >  abs(eclad)*(windmag(eclo,ecla+pma)-windmag(eclo,ecla))
      pwinddir=winddir(eclo,ecla) ! no interpolation, can't average degrees
      pecslp= 0.01 * (ecslp(eclo,ecla) +
     >  abs(eclod)*(ecslp(eclo+pmo,ecla)-ecslp(eclo,ecla)) +
     >  abs(eclad)*(ecslp(eclo,ecla+pma)-ecslp(eclo,ecla)))
      pplev(1:nz) = lpress(1:nz)
      pplev(nz+1) = pecslp
      if(pecslp.le.lpress(nz)) pplev(nz+1)=lpress(nz)+1.0
      dp(nz) = pplev(nz+1)-pplev(nz)
      pect(:) = ect(eclo,ecla,:) +
     >  abs(eclod)*(ect(eclo+pmo,ecla,:)-ect(eclo,ecla,:)) +
     >  abs(eclad)*(ect(eclo,ecla+pma,:)-ect(eclo,ecla,:))
      pecmr(:) = ecmr(eclo,ecla,:) +
     >  abs(eclod)*(ecmr(eclo+pmo,ecla,:)-ecmr(eclo,ecla,:)) +
     >  abs(eclad)*(ecmr(eclo,ecla+pma,:)-ecmr(eclo,ecla,:))
      do i = 1,nz ! get layer averages of T, MR, P
        apect(i) = (pect(i+1)+pect(i))/2.0
        apecmr(i)= (pecmr(i+1)+pecmr(i))/2.0
        pressave(i) = (pplev(i+1)+pplev(i))/2.0
        if(apecmr(i).lt.0.0) apecmr(i)=0.0
      enddo

      if(sss(floor(oe_lon*2.0)+1,floor((90-oe_lat)*2.0)+1).lt.20) then
        if(psss.lt.20) psss= 35.0 ! default if clim grid point missing
        ! if psss already defined, just use that again (mis_val=-999.9)
      else
        psss = sss(floor(oe_lon*2.0)+1,floor((90-oe_lat)*2.0)+1) ! new, ver1.2
      endif

      save_slp(oepix,oelin)  = pplev(nz+1)
      save_wdir(oepix,oelin) = pwinddir
      save_sal(oepix,oelin)  = psss

      save_inc(oepix,oelin,1)= sum(1.0/(1000.0*9.81) *
     >             apecmr(:)*(dp(:))*100.0)
      save_inc(oepix,oelin,2)= pwindmag
      save_inc(oepix,oelin,3)= eclwp(eclo,ecla)+
     >  abs(eclod)*(eclwp(eclo+pmo,ecla)-eclwp(eclo,ecla)) +
     >  abs(eclad)*(eclwp(eclo,ecla+pma)-eclwp(eclo,ecla))
    
      ! Apriori error variances 
      ! Values of 0.00001 mean the parameters are fixed in the retrieval, ie not retrieved.
      ! Must not set any values to 0.0 since internal matrix 
      ! inversion will fail in the retrieval algorithm.

      ! reminder of variable order: EOF1-3 coeffs, WIND, LWP, IWP, SST.
      noretrieval = 0.00001
      ! only works if lwp is 5th state variable!
      do e = 1, npc ! take stddev values from EC analysis
        sa(e,e) = mr_sigmas(ssdex,e)
        sa(e,5) = eof_lw_cov(e,ssdex)
        sa(5,e) = eof_lw_cov(e,ssdex) ! OFFDIAGONAL LWP/EOF COVARIANCES
      enddo

      sa(4,4) = (1.8*real(era_ws(lodex,ladex))*.01)**2 !10m WIND SIGMA SQUARED (m^2/s^2)
      sa(5,5) = 2.0**2 ! LOG10(LWP) SIGMA SQUARED
      sa(6,6) = noretrieval
      if(includeicecloud) sa(6,6) = 0.1**2 ! who knows!?
      sa(7,7) = 0.62**2 ! -- SST SIGMA SQUARED (K^2)
      !sa(7,7) = noretrieval ! -- SST SIGMA SQUARED (K^2)
      if(npc.lt.3) then
        do e = npc+1, 3 
          sa(e,e) = noretrieval
        enddo
      endif

      d = 0
      do e = 1, nvar 
        if(sa(e,e).ne.noretrieval) then
          d=d+1
          xmin(d) = x_min(e)
          xmax(d) = x_max(e)
        endif
        if(e.eq.7.and.sa(e,e).ne.noretrieval) then
          xmin(d) = oe_sst - 1.5 ! set special bounds on SST retrieval
          xmax(d) = oe_sst + 1.5
        endif
      enddo
      if(d.ne.nretvar) then 
        print*,'nretvar doesnt match given variances!'
        stop 
      endif
        

!  Apriori state vector
! NOTE! -- Prior values for EOF coefficients can matter a lot!
!       -- Have to be non-zero, and it's better if 1 or 2 are
!       -- negative. Big values will bias the answer, and 
!       -- values too near zero cause more avg iterations needed!
       ! reminder of variable order: EOF1-3, WIND, LWP, IWP, SST.
      xa(1) =  0.02  ! EOF1 - something nominal, small, non-zero
      xa(2) = -0.01  ! EOF2 - something nominal, small, non-zero
      xa(3) =  0.01  ! EOF3 - something nominal, small, non-zero
      xa(4) = real(era_wm(lodex,ladex))*.01 !10 METER SFC WIND IN M/S
      xa(5) = -3.0 ! (0.001mm)         ![LOG10(LWP)] 
      xa(6) = 0.0001 ! whether retrieving or not, set here
      xa(7) = oe_sst        !SST IN K 
      mrmp(:) = mprof(ssdex,:) ! 
      peofs(:,:) = eofs(ssdex,:,:)
      lapse = 6.26 ! [K/km] from Wilheit 2013

          peofs(:,:) = eofs(ssdex,:,:) ! initialize
          ! smoothing of mrmp/eofs to lessen grid artifacts...
          extra_ssd = (oe_sst-271.5+1) - ssdex
          if(oe_sst .gt. 271.5 .and. oe_sst .lt. 303.5) then
            do i = 1, nz 
              if(extra_ssd.gt.0) then 
                extra_mr(i)=(mprof(ssdex+1,i)-mprof(ssdex,i))*extra_ssd
                do e = 1, npc
                  extra_eo(i,e) = 
     >              (eofs(ssdex+1,i,e)-eofs(ssdex,i,e))*extra_ssd
                enddo
              endif
              if(extra_ssd.lt.0) then
                extra_mr(i)=(mprof(ssdex,i)-mprof(ssdex-1,i))*extra_ssd
                do e = 1, npc
                  extra_eo(i,e) = 
     >              (eofs(ssdex,i,e)-eofs(ssdex-1,i,e))*extra_ssd
                enddo
              endif
              if(extra_ssd.eq.0) then
                extra_mr(i)   = 0.0
                extra_eo(i,:) = 0.0
              endif
            enddo
            mrmp(:) = mrmp(:) + extra_mr(:)
            peofs(:,:) = peofs(:,:) + extra_eo(:,:)
          endif
          if(oe_sst .gt.303.5) then ! special for extending high SST bin
            do i = 1, nz  
              extra_mr(i) = (mprof(ssdex,i)-mprof(ssdex-1,i))*extra_ssd
              do e = 1, npc
                extra_eo(i,e) = 
     >           (eofs(ssdex,i,e)-eofs(ssdex-1,i,e))*extra_ssd
              enddo
            enddo
            mrmp(:) = mrmp(:) + extra_mr(:)
            peofs(:,:) = peofs(:,:) + extra_eo(:,:)
          endif

      if(useanalysis) then
        xa(4) = pwindmag !10 METER SFC WIND IN M/S
        mrmp(:) = apecmr(:)
      endif

        call opt_est(nz,noretrieval,lapse,
     >                oe_tbs,xa,sa,xr,xmax,xmin,
     >                flag,Amatrix,chisq,sigma,niter)

        ! posterior error on LWP is stddev(log10(LWP)). convert back to
        !  units of mm!
        if(flag.eq.1) then
          ! LWP posterior error approximation
          sigma(5) = (10**(xr(5)+sigma(5))-10**(xr(5)-sigma(5)))/2.0 
          sigma(5) = sigma(5)*1000.0 ! convert into g/m^2

          ! TPW posterior error:
          tpwerr_lyr(:,:) = 0.0
          do i = 1, nz
            do c = 1, npc
              dmr(i) = sigma(c)*peofs(i,c)
              tpwerr_lyr(i,c) = 1.0/(10.0*9.81) * dmr(i)*dp(i)
            enddo
          enddo
          if(npc.eq.3) tpw_posterr = sqrt( (sum(tpwerr_lyr(:,1)))**2 +
     >     (sum(tpwerr_lyr(:,2)))**2 + (sum(tpwerr_lyr(:,3)))**2)
          if(npc.eq.2) tpw_posterr = sqrt( (sum(tpwerr_lyr(:,1)))**2 +
     >     (sum(tpwerr_lyr(:,2)))**2 )
        endif

        chisq = chisq / nch ! normalize ChiSq by nchannels

       ! reminder of variable order: EOF coeff 1-3, WIND, LWP, IWP, SST.
       ! output order of vars: TPW, LWP, WIND, IWP, ChiSq, SST, EOF coeff 1-3
               if (flag .eq. 1 .and. chisq .lt. chisq_out_thresh) then
                oe_output(oepix,oelin,1) = final_tpw ! TPW
                oe_output(oepix,oelin,2) = 1000*10**xr(5) ! CLWP [g/m^2]
                oe_output(oepix,oelin,3) = xr(4) ! WIND
                oe_output(oepix,oelin,4) = miss_flt ! IWP
                if(includeicecloud) oe_output(oepix,oelin,4)=1000*xr(6) ! [g/m^2]
                oe_output(oepix,oelin,5) = chisq ! CHI SQUARED
                oe_output(oepix,oelin,6) = xr(7) ! SST
                oe_output(oepix,oelin,7:6+npc) = xr(1:npc) ! EOF coeffs
                poster(oepix,oelin,1:nvar)  = sigma(:)  ! posterior errors
                poster(oepix,oelin,nvar+1) = tpw_posterr ! separate!
                run_count = run_count+1
                tot_iter = tot_iter + niter ! sum # iterations for successful retrievals
                save_iter(oepix,oelin) = niter
                !write(*,*),' =Out:',oe_output(oepix,oelin,:),oepix,oelin
                last_oeout(:) = xr(:) !save for fg in next pixel
                loeop = oepix !save last converged pix#
                loeos = oelin !save last converged scan#
               endif
               ! dont write out values for high chi squared values:
               if (flag .eq. 1 .and. chisq .ge. chisq_out_thresh) then
                oe_output(oepix,oelin,:) = miss_flt 
                poster(oepix,oelin,:) = sigma(:)
                save_iter(oepix,oelin) = niter
                noret_count = noret_count+1
                !write(*,*),'=ChiSqr too big!'
                cycle
               endif
               if (flag .eq. 0 .or. final_tpw.gt.75) then 
                oe_output(oepix,oelin,:) = miss_flt
                poster(oepix,oelin,:)  = miss_flt  ! posterior errors
                noret_count = noret_count+1
                save_iter(oepix,oelin) = niter
                !print*,'=Max iterations or TPW>75',oepix,oelin
                cycle
               endif
               ! for successful retrievals only, save Tbs
               Tb_diff(oepix,oelin,:) = oe_tbs(:)-test_out(:) 
               Tb_sim(oepix,oelin,:)  = test_out(:)
               !print*,'tb diff: ',Tb_diff(oepix,oelin,:)
              enddo
            enddo
            avg_iter = real(tot_iter) / real(run_count) ! iterations per retrieval
            !write(*,*)'Pixels soln found, soln not found: '
            !write(*,*) run_count,noret_count!,notrun_count
            !print*,'avg iter: ',avg_iter

          call cleanup_crtm

        return
      end subroutine retrieval_nr
!--------------------------------------------------- 


      subroutine noscatter(lapse,X,TBOUT)

!     ROUTINES CALLED:
!     run_crtm:  calls crtm RT model

      INTEGER     I,J,K,M,n,z
      REAL        TEMPAV(nz)
      REAL        BTEMP,WIND,LWP,LR
      REAL        icthick, icbase, IWP, ciwc(nz), clwc(nz)
      REAL        X(nvar),TBOUT(nch)
      real :: teebs(nch), lapse, tpw_lyr(nz)
      real :: mixrat(nz)
      integer :: cldbindex, cldtindex, icbindex,ictindex
      logical :: ice
      real :: eco(3) ! EOF coefficients
      real :: rh(nz), reh, mixr

      eco(:) = 0.0 ! initialize, then define eco(1:npc) below
       ! reminder of variable order: EOFC 1-3, WIND, LWP, IWP, SST.
      WIND      = X(4)
      LWP       = 10**(X(5)) !LWP WAS IN LOG FORM; INVERT IT
      IWP       = X(6)
      BTEMP     = X(7) ! ///JUST SST////
      LR        = lapse !T lapse rate, defined above
      eco(1:npc)= X(1:npc)

        cldbindex = nz-5 ! 16-5 = 11 -- 850-800mb. so 850-750mb cld lyr.
        cldtindex = nz-6 ! P-levels! look at def file for corresponding heights!
        icbindex  = 3 ! 400-300mb ice cld lyr
        ictindex  = 3 
        ciwc(:) = 0.0
        clwc(:) = 0.0
        ciwc(ictindex:icbindex)   = iwp/(icbindex-ictindex+1)
        clwc(cldtindex:cldbindex) = lwp/(cldbindex-cldtindex+1)

      if(useanalysis) then
        tempav(:) = apect(:) ! calculated above!
      else
        do i = 1, nz
          tempav(i) = btemp - lr*(lz(i)-dz(i)/2.0)/1000.0
        enddo
      endif

      if(includeicecloud) then
        ciwc(:) = ciwc(:) 
      else
        ciwc = 0.0
      endif

! START NEW WV PROFILE CODE HERE:
      ! the methodology: 
        ! 1. read in mean WV Mixing Ratio profile for given SST bin.
        ! 2. upon iterating, allow first X EOFs to vary in
        !    magnitude to modify MR profile
        ! 3. output final profile in MR space, plus TPW value
      do i = 1, nz
        mixrat(i) = mrmp(i) +
     >            eco(1) * peofs(i,1) +  ! changed to peofs (BCMZ)
     >            eco(2) * peofs(i,2) +  
     >            eco(3) * peofs(i,3) 
        if(mixrat(i).lt.0) mixrat(i)=0.0 ! cant have <zero mass!
      enddo

        save_tprof(oepix,oelin,:)  = tempav(:)
        save_mrprof(oepix,oelin,:) = mixrat(:)
        save_clwc(oepix,oelin,:)   = clwc(:)
        save_ciwc(oepix,oelin,:)   = ciwc(:)


      final_tpw=sum(1.0/(1000.0*9.81) *
     >               mixrat(:)*(dp(:))*100.0)
      !print*,'TPWnew: ',final_tpw
      ! tpw_layer = 1/(rho_water*g) * MR * dP
      ! UNITS:  mm = (m^3/kg)(s^2/m)(g/kg)(kg/m*s^2)  -- need P in Pa!

!----- pass lyr temperature, wv mixing ratio,
!       LWP/IWP per layer, SST, wind speed

      call run_crtm(tempav,mixrat,clwc,ciwc,
     >              btemp,wind,teebs)

        TBOUT(:) = teebs(:)
        !print*,'tbout: ',tbout(:)
        
      RETURN
      END subroutine noscatter     
!-----------------------------------------------------------

      SUBROUTINE opt_est(nz,noretrieval,lp,
     >                    y,xa2,sa2,x2,xmax,xmin,
     >                    flag,Amatrix2,chisq,sigma2,niter)

      implicit none

      integer nz   !Number of layers (nz) in atmosphere.  
      integer flag !flag indicates whether convergence was reached 
      real :: lp  ! lapse rate

      ! Variables used in optimal estimation framework .
      integer :: n,niter ! defined above too!!
      integer    a,b,c,d,i,count,index,checkpol,end_flag
      real noretrieval
      real x(nretvar),x2(nvar),xnext(nretvar),xprime2(nvar)
      real xa(nretvar),xa2(nvar),xmax(nretvar),xmin(nretvar)
      real sa(nretvar,nretvar),sa2(nvar,nvar),sa_i(nretvar,nretvar)
      real K(nch,nretvar),K_t(nretvar,nch)
      real sx(nretvar,nretvar) ! posteriori error matrix
      real sx_i(nretvar,nretvar)
      real sigma2(nvar)
      real AP(nretvar,nretvar),Amatrix2(nvar)
      real F(nch),Fout(nch),Foutprime(nch) ! changed!
      real Fprime(nch),Fdbprime(nch),y(nch),dF(nch),dx(nretvar)
      real sum1(nretvar,1)         !xa-x
      real sum2(nch,1)           !y-F
      real sum3(nretvar,1)         !prod4+prod3
      real prod1(nretvar,nch)  !K_t*sy_i
      real prod2(nretvar,nretvar)!K_t*sy_i*K
      real prod3(nretvar)        !sa_i*(sum1)
      real prod4(nretvar)        !prod1*sum2
      real prod5(nretvar)        !sx*sum3
      real xdiff(nretvar,1)        !xnext-x
      real xdiff_t(1,nretvar)    !transposed
      real prod6(1,nretvar)      !xdiff_t*sx_i
      real prod7(1,1)            !prod6*xdiff ! CHANGED!!
      real sum4(nch,1)           !F-y
      real sum4_t(1,nch)       !(F-y) transpose
      real sum5(nretvar,1)         !x-xa
      real sum5_t(1,nretvar)     !(x-xa) transpose
      real prod8(1,nch)        !sum4_t*sy_i
      real prod9(1,1)              !prod8*sum2
      real prod10(1,nretvar)     !sum5_t*sa_i
      real prod11(1,1)           !prod10*sum5 ! CHANGED!!!
      real chisqtot(1,1)         !chi squared (prod9 + prod11) ! CHANGED!!
      real chisq                 !cost function (apriori + measurement fit)
      real sum6(nretvar,nretvar) !(IM-AP)
      real prod12(nretvar,nch) !Sx*K_t
      real prod13(nretvar,nch) !Sx*K_t*Sy_i
      real prod14(nretvar)       !prod13*sum2 - contribution from obs
      real prod15(nretvar)       !sum6*sum1 - contribution from apriori
      real prod16(nretvar)       !prod14+prod16 = xnext-x

      real Fold(nch)           !old fwd modeled Tbs
      real diffF(nch,1),diffF_t(1,nch)  !F-F_old, transpose
      real Sdy(nch,nch),Sdy_i(nch,nch) ! (Rodgers Eq5.27) 
      real Ksa(nch,nretvar), Ksakt(nch,nch) ! K*Sa, K*Sa*K_t
      real Sdysum(nch,nch), Sdysum_i(nch,nch) ! (K*Sa*K_t+Sy), *_i
      real almost(nch,nch)! Sy*(Ksakt+Sy)_i (most of 5.27)
      real almost3(1,nch) ! diffF_t*Sdy_t
      real*8 disq(1,1)   !di squared (Eq5.33), used for conv test

      integer :: good(nretvar)
      real :: psdiff ! pix/scan distance from last convergence reached 

      !real :: savejake(nch,nretvar) ! experimental!
      !integer :: moditer


!-----Declare Apriori errors 
      sa(:,:) = 0.0
      b = 1
      do a=1,nvar ! less complex if not using covariances...
        if(sa2(a,a).ne.noretrieval) then 
          good(b) = a
          b=b+1
        endif
      enddo
      do a=1,nretvar
       do c=1,nretvar
         sa(a,c)=sa2(good(a),good(c))
       enddo
      enddo
        
!-----GENERATE FIRST GUESS FORWARD MODEL
! - - instead of using a priori as first guess, use last (good) retrieval
      psdiff = sqrt(real(loeop-oepix)**2 + real(loeos-oelin)**2)
      if(psdiff.le.3.0) then
        b=1
        do a = 1, nvar
          x2(a) = xa2(a) ! should be a priori
          if (sa2(a,a).NE.noretrieval) then
            x(b)  = last_oeout(a) ! use last 'good' OE pixel output
            xa(b) = xa2(a)
            x2(a) = last_oeout(a)
            b=b+1
          endif
        enddo
      else
        b=1
        do a = 1, nvar
          x2(a) = xa2(a)
          if (sa2(a,a).NE.noretrieval) then
            x(b)  = xa2(a)  ! simply use a priori values
            xa(b) = xa2(a)
            b=b+1
          endif
        enddo
      endif 

!-----NEWTONIAN ITERATION-----
      end_flag=0
      do n = 1, nn
        !write(*,*)' n iter: ',n
        !write(*,*)'X: ',x2(:)
        call noscatter(lp,x2,Fout)
        !write(*,*)' Fout:',Fout,n
        F(:) = Fout(:)

!-----CALCULATE JACOBIAN K-----
       !moditer = mod(n,2)
       !print*,'n mod 2: ',moditer
       ! uncomment to speed up code slightly (some downsides!)
       !if(moditer.eq.1) then
        do a = 1, nvar
          xprime2(a) = x2(a)              !Initialize xprime
        enddo
        a = 1
        do d = 1, nvar
         if (sa2(d,d).NE.noretrieval) then
          xprime2(d) = x2(d) + 0.02*x2(d)!//2% perturb to get a slope// 
          if(x2(d).eq.0) xprime2(d) = x2(d) + 0.02*xa2(d)
          dx(a) = xprime2(d) - x2(d)
          ! Forward model call -- EVERY TIME! (major resource hog!)
          call noscatter(lp,xprime2,Foutprime)
          Fprime(:) = Foutprime(:) 
          dF(:) = Fprime(:) - F(:)
          K(:,a) = dF(:)/dx(a) ! Jacobian calculation!
          a = a + 1
          xprime2(d) = x2(d) ! reset for next iteration
         endif
        enddo
        !savejake(:,:) = K(:,:)
       !endif
       !K(:,:) = savejake(:,:)
!          print*,'xpr2:',xprime2(1),xprime2(3),
!     >       xprime2(2),xprime2(4)

c-----EVALUATE NEWTONIAN MATRICES-----   
        !write(*,*)'1 args',K
        K_t(:,:) = transpose(K(:,:))
        !write(*,*)'1 args',K_t
        !write(*,*)'2 args',sa
        call inverse(sa,sa_i,nretvar)
        !write(*,*)'2 args',sa_i
        !write(*,*)'3 args',K_t,sy_i,
        prod1 = matmul(K_t,sy_i)
        !write(*,*)'5 args',prod1
        prod2 = matmul(prod1,K)
        !write(*,*)'6 args',prod2
        sx_i = sa_i + prod2
        !write(*,*)'7 args',sx_i
        call inverse(sx_i,sx,nretvar)
        !write(*,*)'8 args',sx

c-----PERFORM NEWTONIAN STEP-----
        sum1(:,1) = xa(:) - x(:)
        sum2(:,1) = y(:) - F(:)
        prod3 = matmul(sa_i,sum1(:,1))
        prod4 = matmul(prod1,sum2(:,1))
        sum3(:,1) = prod4(:) + prod3(:)
        prod5 = matmul(sx(:,:),sum3(:,1))
        xnext = x + prod5

c-----ERROR DIAGNOSTICS-----
        AP = matmul(sx,prod2)
        sum4(:,1) = F(:) - y(:)
        sum4_t = transpose(sum4)
        sum5(:,1) = x(:) - xa(:)
        sum5_t = transpose(sum5)
        prod8(:,:) = matmul(sum4_t(:,:),sy_i(:,:))
        prod9(:,:) = matmul(prod8(:,:),sum4(:,:))
        prod10(:,:) = matmul(sum5_t(:,:),sa_i(:,:))
        prod11(:,:) = matmul(prod10(:,:),sum5(:,:))
        chisqtot(1,1) = prod9(1,1) + prod11(1,1)
        !Check limits on xnext
        do a=1,nretvar
          if (xnext(a) .gt. xmax(a)) then
            xnext(a) = xmax(a)
          else if (xnext(a) .lt. xmin(a)) then
            xnext(a) = xmin(a)
          endif
        end do
        
        ! Evaluate closeness of Xnext and X
        xdiff(:,1) = xnext(:) - x(:)
        xdiff_t(:,:) = transpose(xdiff(:,:))
        prod6(1,:) = matmul(xdiff_t(1,:),sx_i(:,:))
        prod7(:,:) = matmul(prod6(:,:),xdiff(:,:))

        b = 1
        do a = 1, nvar
         if (sa2(a,a).NE.noretrieval) then
          x(b)  = xnext(b)
          x2(a) = xnext(b)
          b = b + 1
         endif
        end do

! ***As final retrieved state is found, TBs need to be
! ***recomputed.
        call noscatter(lp,x2,Fout)
        F(:) = Fout(:)
        !print*,'Fout: ',Fout

        if (n .gt. 1) then
          diffF(:,1) = F(:) - Fold(:)
          !print*,'Fout: ',Fout
          !print*,'Fold: ',Fold
          !print*,'diffF: ',diffF
          diffF_t = transpose(diffF)

          Ksa = matmul(K,sa)
          Ksakt = matmul(Ksa,K_t)
          Sdysum = Ksakt + Sy
          call inverse(Sdysum,Sdysum_i,nch)
          almost = matmul(sy,Sdysum_i)
          Sdy = matmul(almost,sy)
          call inverse(Sdy,Sdy_i,nch)
          almost3(1,:) = matmul(diffF_t(1,:),Sdy_i(:,:)) ! try this????
          disq = matmul(dble(almost3),dble(diffF))

        ! first is old convergence criteria, 5.33 is new.
          !if (prod7(1,1) .lt. nretvar/10.) then ! Eq 5.29 in Rodgers 2000 
          if (disq(1,1) .lt. dble(nch/conv_factor)) then ! Eq 5.33 in Rodgers 2000 
             end_flag=1
          endif
!        ! secondary convergence crit. for 10+ iter, low chisq
!          if(disq(1,1) .ge. nch/conv_factor .and.
!     >        n.ge.10 .and. chisq.lt.chisq_out_thresh) then 
!             end_flag=1
!             !print*,'secondary conv criteria met'
!          endif
          if(disq(1,1) .lt. 0.0) then ! should be impossible- need real*8?
            !print*,'convergence error -- negative value not allowed!'
            end_flag=0
          endif
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (end_flag .eq. 1) then
	  b = 1
	  do a = 1,nvar 
           if (sa2(a,a).NE.noretrieval) then	  
            Amatrix2(a)=AP(b,b) 
	    sigma2(a)  =sqrt(sx(b,b)) ! -- posterior error standard deviations
	    b = b + 1
	   else
	    Amatrix2(a) = 0.0
	    sigma2(a)   = miss_flt 
	   endif
          enddo
          flag=1
          niter = n
          return
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! now outside of flag check, so chisq for no convergence too
        chisq = chisqtot(1,1) - prod11(1,1) 

        Fold(:) = Fout(:)

      end do  !!!/// end of do loop starting way up above

!      write(*,*) 'Solution not found in ',nn,' steps'
!      No solution found upon reaching max iterations?  Set the following to
!      -99.99.
      do a = 1,nvar
       Amatrix2(a) = -99.99
       sigma2(a)   = -99.99
      enddo
      chisq      = -99.99
      flag       = 0
      niter      = n
      return

      end subroutine opt_est

!----------------------------------------------------------------------
      subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
        implicit none 
        integer n
        real a(n,n), c(n,n)
        real*8 aa(n,n),cc(n,n)
        real*8 L(n,n), U(n,n), b(n), d(n), x(n)
        real*8 coeff
        integer i, j, k

        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0
        U=0.0
        b=0.0

        aa(:,:) = dble(a(:,:)) ! do double precision!
        ! step 1: forward elimination
        do k=1, n-1
           do i=k+1,n
              coeff=a(i,k)/a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                 a(i,j) = a(i,j)-coeff*a(k,j)
              end do
           end do
        end do

        ! Step 2: prepare L and U matrices 
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
          L(i,i) = 1.0
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
          do i=1,j
            U(i,j) = a(i,j)
          end do
        end do

        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
          b(k)=1.0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
            cc(i,k) = x(i)
          end do
          b(k)=0.0
        end do
        c(:,:) = real(cc(:,:))
      end subroutine inverse

      end module optest
