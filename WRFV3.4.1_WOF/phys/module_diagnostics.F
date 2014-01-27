!WRF:MEDIATION_LAYER:PHYSICS
!

MODULE module_diagnostics
CONTAINS
   SUBROUTINE diagnostic_output_calc(                                 &
                      ids,ide, jds,jde, kds,kde,                      &
                      ims,ime, jms,jme, kms,kme,                      &
                      ips,ipe, jps,jpe, kps,kpe,                      & ! patch  dims
                      i_start,i_end,j_start,j_end,kts,kte,num_tiles   &
                     ,dpsdt,dmudt                                     &
                     ,p8w,pk1m,mu_2,mu_2m                             &
                     ,u,v                                             &
                     ,raincv,rainncv,rainc,rainnc                     &
                     ,i_rainc,i_rainnc                                &
                     ,hfx,sfcevp,lh                                   &
                     ,ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC               & ! Optional
                     ,ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC               & ! Optional
                     ,ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC               & ! Optional
                     ,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC               & ! Optional
                     ,I_ACSWUPT,I_ACSWUPTC,I_ACSWDNT,I_ACSWDNTC       & ! Optional
                     ,I_ACSWUPB,I_ACSWUPBC,I_ACSWDNB,I_ACSWDNBC       & ! Optional
                     ,I_ACLWUPT,I_ACLWUPTC,I_ACLWDNT,I_ACLWDNTC       & ! Optional
                     ,I_ACLWUPB,I_ACLWUPBC,I_ACLWDNB,I_ACLWDNBC       & ! Optional
                     ,dt,xtime,sbw,t2                                 &
                     ,diag_print                                      &
                     ,bucket_mm, bucket_J                             &
                     ,prec_acc_c, prec_acc_nc, snow_acc_nc            &
                     ,snowncv, prec_acc_dt, curr_secs                 &
                                                                      )
!----------------------------------------------------------------------

  USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval

   IMPLICIT NONE
!======================================================================
! Definitions
!-----------
!-- DIAG_PRINT    print control: 0 - no diagnostics; 1 - dmudt only; 2 - all
!-- DT            time step (second)
!-- XTIME         forecast time
!-- SBW           specified boundary width - used later
!
!-- P8W           3D pressure array at full eta levels
!-- MU            dry column hydrostatic pressure
!-- RAINC         cumulus scheme precipitation since hour 0
!-- RAINCV        cumulus scheme precipitation in one time step (mm)
!-- RAINNC        explicit scheme precipitation since hour 0
!-- RAINNCV       explicit scheme precipitation in one time step (mm)
!-- SNOWNCV       explicit scheme snow in one time step (mm)
!-- HFX           surface sensible heat flux
!-- LH            surface latent heat flux
!-- SFCEVP        total surface evaporation
!-- U             u component of wind - to be used later to compute k.e.
!-- V             v component of wind - to be used later to compute k.e.
!-- PREC_ACC_C    accumulated convective precip over accumulation time prec_acc_dt
!-- PREC_ACC_NC   accumulated explicit precip over accumulation time prec_acc_dt
!-- SNOW_ACC_NC   accumulated explicit snow precip over accumulation time prec_acc_dt
!-- PREC_ACC_DT   precip accumulation time, default is 60 min
!-- CURR_SECS     model time in seconds
!
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- ips           start index for i in patch
!-- ipe           end index for i in patch
!-- jps           start index for j in patch
!-- jpe           end index for j in patch
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!
!======================================================================

   INTEGER,      INTENT(IN   )    ::                             &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      ips,ipe, jps,jpe, kps,kpe, &
                                                        kts,kte, &
                                                      num_tiles

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                  &
     &           i_start,i_end,j_start,j_end

   INTEGER,      INTENT(IN   )    ::   diag_print
   REAL,      INTENT(IN   )    ::   bucket_mm, bucket_J

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         INTENT(IN ) ::                                       u  &
                                                    ,         v  &
                                                    ,       p8w

   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN) ::           &
                                                           MU_2  &
                                                    ,   RAINNCV  &
                                                    ,    RAINCV  &
                                                    ,   SNOWNCV  &
                                                    ,       HFX  &
                                                    ,        LH  &
                                                    ,    SFCEVP  &  
                                                    ,        T2     

   REAL, DIMENSION( ims:ime , jms:jme ),                         &
          INTENT(INOUT) ::                                DPSDT  &
                                                    ,     DMUDT  &
                                                    ,    RAINNC  &
                                                    ,     RAINC  &
                                                    ,     MU_2M  &
                                                    ,      PK1M
 
   REAL,  INTENT(IN   ) :: DT, XTIME
   INTEGER,  INTENT(IN   ) :: SBW
   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) ::     &
                                                       I_RAINC,  &
                                                       I_RAINNC
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                      ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC,          &
                      ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,          &
                      ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC,          &
                      ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC
   INTEGER, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                      I_ACSWUPT,I_ACSWUPTC,I_ACSWDNT,I_ACSWDNTC,  &
                      I_ACSWUPB,I_ACSWUPBC,I_ACSWDNB,I_ACSWDNBC,  &
                      I_ACLWUPT,I_ACLWUPTC,I_ACLWDNT,I_ACLWDNTC,  &
                      I_ACLWUPB,I_ACLWUPBC,I_ACLWDNB,I_ACLWDNBC

   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                      PREC_ACC_C, PREC_ACC_NC, SNOW_ACC_NC

   REAL, OPTIONAL, INTENT(IN)::  PREC_ACC_DT, CURR_SECS

   INTEGER :: i,j,k,its,ite,jts,jte,ij
   INTEGER :: idp,jdp,irc,jrc,irnc,jrnc,isnh,jsnh
   INTEGER :: prfreq

   REAL              :: no_points
   REAL              :: dpsdt_sum, dmudt_sum, dardt_sum, drcdt_sum, drndt_sum
   REAL              :: hfx_sum, lh_sum, sfcevp_sum, rainc_sum, rainnc_sum, raint_sum
   REAL              :: dmumax, raincmax, rainncmax, snowhmax
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   CHARACTER*256     :: outstring
   CHARACTER*6       :: grid_str

!-----------------------------------------------------------------
! Handle accumulations with buckets to prevent round-off truncation in long runs
! This is done every 360 minutes assuming time step fits exactly into 360 minutes
   IF(bucket_mm .gt. 0. .AND. MOD(NINT(XTIME),360) .EQ. 0)THEN
! SET START AND END POINTS FOR TILES
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )

   DO ij = 1 , num_tiles

      IF (xtime .eq. 0.0)THEN
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
          i_rainnc(i,j) = 0
          i_rainc(i,j) = 0
        ENDDO      
        ENDDO
      ENDIF
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
        IF(rainnc(i,j) .gt. bucket_mm)THEN
          rainnc(i,j) = rainnc(i,j) - bucket_mm
          i_rainnc(i,j) =  i_rainnc(i,j) + 1
        ENDIF
        IF(rainc(i,j) .gt. bucket_mm)THEN
          rainc(i,j) = rainc(i,j) - bucket_mm
          i_rainc(i,j) =  i_rainc(i,j) + 1
        ENDIF
      ENDDO      
      ENDDO

      IF (xtime .eq. 0.0 .and. bucket_J .gt. 0.0 .and. PRESENT(ACSWUPT))THEN
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
          i_acswupt(i,j) = 0
          i_acswuptc(i,j) = 0
          i_acswdnt(i,j) = 0
          i_acswdntc(i,j) = 0
          i_acswupb(i,j) = 0
          i_acswupbc(i,j) = 0
          i_acswdnb(i,j) = 0
          i_acswdnbc(i,j) = 0
        ENDDO      
        ENDDO
      ENDIF
      IF (xtime .eq. 0.0  .and. bucket_J .gt. 0.0 .and. PRESENT(ACLWUPT))THEN
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
          i_aclwupt(i,j) = 0
          i_aclwuptc(i,j) = 0
          i_aclwdnt(i,j) = 0
          i_aclwdntc(i,j) = 0
          i_aclwupb(i,j) = 0
          i_aclwupbc(i,j) = 0
          i_aclwdnb(i,j) = 0
          i_aclwdnbc(i,j) = 0
        ENDDO      
        ENDDO
      ENDIF
      IF (PRESENT(ACSWUPT) .and. bucket_J .gt. 0.0)THEN
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
        IF(acswupt(i,j) .gt. bucket_J)THEN
          acswupt(i,j) = acswupt(i,j) - bucket_J
          i_acswupt(i,j) =  i_acswupt(i,j) + 1
        ENDIF
        IF(acswuptc(i,j) .gt. bucket_J)THEN
          acswuptc(i,j) = acswuptc(i,j) - bucket_J
          i_acswuptc(i,j) =  i_acswuptc(i,j) + 1
        ENDIF
        IF(acswdnt(i,j) .gt. bucket_J)THEN
          acswdnt(i,j) = acswdnt(i,j) - bucket_J
          i_acswdnt(i,j) =  i_acswdnt(i,j) + 1
        ENDIF
        IF(acswdntc(i,j) .gt. bucket_J)THEN
          acswdntc(i,j) = acswdntc(i,j) - bucket_J
          i_acswdntc(i,j) =  i_acswdntc(i,j) + 1
        ENDIF
        IF(acswupb(i,j) .gt. bucket_J)THEN
          acswupb(i,j) = acswupb(i,j) - bucket_J
          i_acswupb(i,j) =  i_acswupb(i,j) + 1
        ENDIF
        IF(acswupbc(i,j) .gt. bucket_J)THEN
          acswupbc(i,j) = acswupbc(i,j) - bucket_J
          i_acswupbc(i,j) =  i_acswupbc(i,j) + 1
        ENDIF
        IF(acswdnb(i,j) .gt. bucket_J)THEN
          acswdnb(i,j) = acswdnb(i,j) - bucket_J
          i_acswdnb(i,j) =  i_acswdnb(i,j) + 1
        ENDIF
        IF(acswdnbc(i,j) .gt. bucket_J)THEN
          acswdnbc(i,j) = acswdnbc(i,j) - bucket_J
          i_acswdnbc(i,j) =  i_acswdnbc(i,j) + 1
        ENDIF
      ENDDO      
      ENDDO
      ENDIF
      IF (PRESENT(ACLWUPT) .and. bucket_J .gt. 0.0)THEN
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
        IF(aclwupt(i,j) .gt. bucket_J)THEN
          aclwupt(i,j) = aclwupt(i,j) - bucket_J
          i_aclwupt(i,j) =  i_aclwupt(i,j) + 1
        ENDIF
        IF(aclwuptc(i,j) .gt. bucket_J)THEN
          aclwuptc(i,j) = aclwuptc(i,j) - bucket_J
          i_aclwuptc(i,j) =  i_aclwuptc(i,j) + 1
        ENDIF
        IF(aclwdnt(i,j) .gt. bucket_J)THEN
          aclwdnt(i,j) = aclwdnt(i,j) - bucket_J
          i_aclwdnt(i,j) =  i_aclwdnt(i,j) + 1
        ENDIF
        IF(aclwdntc(i,j) .gt. bucket_J)THEN
          aclwdntc(i,j) = aclwdntc(i,j) - bucket_J
          i_aclwdntc(i,j) =  i_aclwdntc(i,j) + 1
        ENDIF
        IF(aclwupb(i,j) .gt. bucket_J)THEN
          aclwupb(i,j) = aclwupb(i,j) - bucket_J
          i_aclwupb(i,j) =  i_aclwupb(i,j) + 1
        ENDIF
        IF(aclwupbc(i,j) .gt. bucket_J)THEN
          aclwupbc(i,j) = aclwupbc(i,j) - bucket_J
          i_aclwupbc(i,j) =  i_aclwupbc(i,j) + 1
        ENDIF
        IF(aclwdnb(i,j) .gt. bucket_J)THEN
          aclwdnb(i,j) = aclwdnb(i,j) - bucket_J
          i_aclwdnb(i,j) =  i_aclwdnb(i,j) + 1
        ENDIF
        IF(aclwdnbc(i,j) .gt. bucket_J)THEN
          aclwdnbc(i,j) = aclwdnbc(i,j) - bucket_J
          i_aclwdnbc(i,j) =  i_aclwdnbc(i,j) + 1
        ENDIF
      ENDDO      
      ENDDO
      ENDIF
   ENDDO
!  !$OMP END PARALLEL DO
   ENDIF

! Compute precipitation accumulation in a given time window: prec_acc_dt
   IF (prec_acc_dt .gt. 0.) THEN

!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )

   DO ij = 1 , num_tiles

      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         IF (mod(curr_secs, 60.* prec_acc_dt) == 0.) THEN
            prec_acc_c(i,j)  = 0.
            prec_acc_nc(i,j) = 0.
            snow_acc_nc(i,j)  = 0.
         ENDIF
         prec_acc_c(i,j)  = prec_acc_c(i,j)  +  RAINCV(i,j)
         prec_acc_nc(i,j) = prec_acc_nc(i,j) + RAINNCV(i,j)
         prec_acc_c(i,j)  = MAX (prec_acc_c(i,j), 0.0)
         prec_acc_nc(i,j) = MAX (prec_acc_nc(i,j), 0.0)
         snow_acc_nc(i,j)   = snow_acc_nc(i,j) + SNOWNCV(I,J)
! add convective precip to snow bucket if t2 < 273.15
         IF ( t2(i,j) .lt. 273.15 ) THEN
         snow_acc_nc(i,j)   = snow_acc_nc(i,j) +  RAINCV(i,j)
         snow_acc_nc(i,j)   = MAX (snow_acc_nc(i,j), 0.0)
         ENDIF
      ENDDO     
      ENDDO     

   ENDDO     

!  !$OMP END PARALLEL DO
   ENDIF

   if (diag_print .eq. 0 ) return

   IF ( xtime .ne. 0. ) THEN

    if(diag_print.eq.1) then
       prfreq = dt
!      prfreq = max(2,int(dt/60.))   ! in min
    else
       prfreq=10                   ! in min
    endif
   
    IF (MOD(nint(dt),prfreq) == 0) THEN

! COMPUTE THE NUMBER OF MASS GRID POINTS
   no_points = float((ide-ids)*(jde-jds))

! SET START AND END POINTS FOR TILES
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )

   dmumax = 0.
   DO ij = 1 , num_tiles

!     print *, i_start(ij),i_end(ij),j_start(ij),j_end(ij)
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         dpsdt(i,j)=(p8w(i,kms,j)-pk1m(i,j))/dt
         dmudt(i,j)=(mu_2(i,j)-mu_2m(i,j))/dt
         if(abs(dmudt(i,j)*dt).gt.dmumax)then
           dmumax=abs(dmudt(i,j)*dt)
           idp=i
           jdp=j
         endif
      ENDDO      
      ENDDO

   ENDDO
!  !$OMP END PARALLEL DO

! convert DMUMAX from (PA) to (bars) per time step
   dmumax = dmumax*1.e-5
! compute global MAX
   CALL wrf_dm_maxval ( dmumax,  idp, jdp )

!  print *, 'p8w(30,1,30),pk1m(30,30) : ', p8w(30,1,30),pk1m(30,30)
!  print *, 'mu_2(30,30),mu_2m(30,30) : ', mu_2(30,30),mu_2m(30,30)
   dpsdt_sum = 0.
   dmudt_sum = 0.

   DO j = jps, min(jpe,jde-1)
     DO i = ips, min(ipe,ide-1)
       dpsdt_sum = dpsdt_sum + abs(dpsdt(i,j))
       dmudt_sum = dmudt_sum + abs(dmudt(i,j))
     ENDDO
   ENDDO

! compute global sum
   dpsdt_sum = wrf_dm_sum_real ( dpsdt_sum )
   dmudt_sum = wrf_dm_sum_real ( dmudt_sum )

!  print *, 'dpsdt, dmudt : ', dpsdt_sum, dmudt_sum

   IF ( diag_print .eq. 2 ) THEN
   dardt_sum = 0.
   drcdt_sum = 0.
   drndt_sum = 0.
   rainc_sum = 0.
   raint_sum = 0.
   rainnc_sum = 0.
   sfcevp_sum = 0.
   hfx_sum = 0.
   lh_sum = 0.
   raincmax = 0.
   rainncmax = 0.

   DO j = jps, min(jpe,jde-1)
     DO i = ips, min(ipe,ide-1)
       drcdt_sum = drcdt_sum + abs(raincv(i,j))
       drndt_sum = drndt_sum + abs(rainncv(i,j))
       dardt_sum = dardt_sum + abs(raincv(i,j)) + abs(rainncv(i,j))
       rainc_sum = rainc_sum + abs(rainc(i,j))
! MAX for accumulated conv precip
       IF(rainc(i,j).gt.raincmax)then
          raincmax=rainc(i,j)
          irc=i
          jrc=j
       ENDIF
       rainnc_sum = rainnc_sum + abs(rainnc(i,j))
! MAX for accumulated resolved precip
       IF(rainnc(i,j).gt.rainncmax)then
          rainncmax=rainnc(i,j)
          irnc=i
          jrnc=j
       ENDIF
       raint_sum = raint_sum + abs(rainc(i,j)) + abs(rainnc(i,j))
       sfcevp_sum = sfcevp_sum + abs(sfcevp(i,j))
       hfx_sum = hfx_sum + abs(hfx(i,j))
       lh_sum = lh_sum + abs(lh(i,j))
     ENDDO
   ENDDO

! compute global MAX
   CALL wrf_dm_maxval ( raincmax, irc, jrc )
   CALL wrf_dm_maxval ( rainncmax, irnc, jrnc )

! compute global sum
   drcdt_sum = wrf_dm_sum_real ( drcdt_sum )
   drndt_sum = wrf_dm_sum_real ( drndt_sum )
   dardt_sum = wrf_dm_sum_real ( dardt_sum )
   rainc_sum = wrf_dm_sum_real ( rainc_sum )
   rainnc_sum = wrf_dm_sum_real ( rainnc_sum )
   raint_sum = wrf_dm_sum_real ( raint_sum )
   sfcevp_sum = wrf_dm_sum_real ( sfcevp_sum )
   hfx_sum = wrf_dm_sum_real ( hfx_sum )
   lh_sum = wrf_dm_sum_real ( lh_sum )

   ENDIF

! print out the average values

   CALL get_current_grid_name( grid_str )

#ifdef DM_PARALLEL
   IF ( wrf_dm_on_monitor() ) THEN
#endif
     WRITE(outstring,*) grid_str,'Domain average of dpsdt, dmudt (mb/3h): ', xtime, &
           dpsdt_sum/no_points*108., &
           dmudt_sum/no_points*108.
     CALL wrf_message ( TRIM(outstring) )

     WRITE(outstring,*) grid_str,'Max mu change time step: ', idp,jdp,dmumax
     CALL wrf_message ( TRIM(outstring) )

     IF ( diag_print .eq. 2) THEN
     WRITE(outstring,*) grid_str,'Domain average of dardt, drcdt, drndt (mm/sec): ', xtime, &
           dardt_sum/dt/no_points, &
           drcdt_sum/dt/no_points, &
           drndt_sum/dt/no_points
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Domain average of rt_sum, rc_sum, rnc_sum (mm): ', xtime, &
           raint_sum/no_points, &
           rainc_sum/no_points, &
           rainnc_sum/no_points
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Max Accum Resolved Precip,   I,J  (mm): '               ,&
           rainncmax,irnc,jrnc
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Max Accum Convective Precip,   I,J  (mm): '             ,&
           raincmax,irc,jrc
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Domain average of sfcevp, hfx, lh: ', xtime, &
           sfcevp_sum/no_points, &
           hfx_sum/no_points, &
           lh_sum/no_points
     CALL wrf_message ( TRIM(outstring) )
     ENDIF
#ifdef DM_PARALLEL
   ENDIF
#endif

    ENDIF        ! print frequency
   ENDIF

! save values at this time step
   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij,i,j )
   DO ij = 1 , num_tiles

      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         pk1m(i,j)=p8w(i,kms,j)
         mu_2m(i,j)=mu_2(i,j)
      ENDDO
      ENDDO

      IF ( xtime .lt. 0.0001 ) THEN
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         dpsdt(i,j)=0.
         dmudt(i,j)=0.
      ENDDO
      ENDDO
      ENDIF

   ENDDO
   !$OMP END PARALLEL DO

   END SUBROUTINE diagnostic_output_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE clwrf_output_calc(                                      &
                      ids,ide, jds,jde, kds,kde,                      &
                      ims,ime, jms,jme, kms,kme,                      &
                      ips,ipe, jps,jpe, kps,kpe,                      & ! patch  dims
                      i_start,i_end,j_start,j_end,kts,kte,num_tiles   &
                     ,dpsdt,dmudt                                     &
                     ,p8w,pk1m,mu_2,mu_2m                             &
                     ,u,v                                             &
                     ,is_restart                                      & ! CLWRF
                     ,clwrfH,t2,q2,u10,v10, skintemp                  & ! CLWRF
                     ,t2clmin,t2clmax,tt2clmin,tt2clmax               & ! CLWRF
                     ,t2clmean,t2clstd                                & ! CLWRF
                     ,q2clmin,q2clmax,tq2clmin,tq2clmax               & ! CLWRF
                     ,q2clmean,q2clstd                                & ! CLWRF
                     ,u10clmax,v10clmax,spduv10clmax,tspduv10clmax    & ! CLWRF
                     ,u10clmean,v10clmean,spduv10clmean               & ! CLWRF
                     ,u10clstd,v10clstd,spduv10clstd                  & ! CLWRF
                     ,raincclmax,rainncclmax,traincclmax,trainncclmax & ! CLWRF
                     ,raincclmean,rainncclmean,raincclstd,rainncclstd & ! CLWRF
                     ,skintempclmin,skintempclmax                     & ! CLWRF
                     ,tskintempclmin,tskintempclmax                   & ! CLWRF
                     ,skintempclmean,skintempclstd                    & ! CLWRF
                     ,raincv,rainncv,rainc,rainnc                     &
                     ,i_rainc,i_rainnc                                &
                     ,hfx,sfcevp,lh                                   &
                     ,dt,xtime,sbw                                    &
                     ,diag_print                                      &
                     ,bucket_mm, bucket_J                             &
                                                                      )
!----------------------------------------------------------------------

  USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
  USE module_configure 

   IMPLICIT NONE
!======================================================================
! Definitions
!-----------
!-- DT            time step (second)
!-- XTIME         forecast time
!-- SBW           specified boundary width - used later
!
!-- P8W           3D pressure array at full eta levels
!-- MU            dry column hydrostatic pressure
!-- RAINC         cumulus scheme precipitation since hour 0
!-- RAINCV        cumulus scheme precipitation in one time step (mm)
!-- RAINNC        explicit scheme precipitation since hour 0
!-- RAINNCV       explicit scheme precipitation in one time step (mm)
!-- HFX           surface sensible heat flux
!-- LH            surface latent heat flux
!-- SFCEVP        total surface evaporation
!-- U             u component of wind - to be used later to compute k.e.
!-- V             v component of wind - to be used later to compute k.e.
!
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- ips           start index for i in patch
!-- ipe           end index for i in patch
!-- jps           start index for j in patch
!-- jpe           end index for j in patch
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!
! CLWRF-UC May.09 definitions
!-----------
! is_restart: whether if simulation is a restart
! clwrfH: Interval (hour) of accumulation for computations 
! [var]cl[min/max]: [minimum/maximum] of variable [var] during interval
! t[var]cl[min/max]: Time (minutes) of [minimum/maximum] of variable 
!    [var] during interval 
! [var]clmean: mean of variable [var] during interval
! [var]clstd: standard dev. of variable [var] during interval
!    Variables are written on aux_hist_out7 (established
!    in Registry)
!
!======================================================================

   INTEGER,      INTENT(IN   )                     ::            &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      ips,ipe, jps,jpe, kps,kpe, &
                                                        kts,kte, &
                                                      num_tiles

   INTEGER, DIMENSION(num_tiles), INTENT(IN)       :: i_start,   &
                                      i_end,j_start,j_end

   INTEGER,      INTENT(IN   )                     :: diag_print
   REAL,      INTENT(IN   )                        :: bucket_mm, &
                                      bucket_J

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
                                       INTENT(IN ) :: u,v,p8w    

   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: MU_2,     & 
                                      RAINNCV, RAINCV, HFX,      &
                                      SFCEVP, LH, SKINTEMP 

   REAL, DIMENSION( ims:ime , jms:jme ),                         &
                                     INTENT(INOUT) :: DPSDT,     &
                                     DMUDT, RAINNC, RAINC,       &
                                     MU_2M, PK1M
 
   REAL,  INTENT(IN   )                            :: DT, XTIME
   INTEGER,  INTENT(IN   )                         :: SBW
   INTEGER, DIMENSION( ims:ime , jms:jme ),                      & 
                                     INTENT(INOUT) :: I_RAINC,   &
                                     I_RAINNC

! LOCAL  VAR

   INTEGER                                   :: i,j,k,its,ite,jts,jte,ij
   INTEGER                                   :: idp,jdp,irc,jrc,irnc,jrnc,isnh,jsnh
   INTEGER                                   :: prfreq

   REAL                                      :: dpsdt_sum, dmudt_sum, dardt_sum,  & 
                          drcdt_sum, drndt_sum
   REAL                                      :: hfx_sum, lh_sum, sfcevp_sum,      &
                          rainc_sum, rainnc_sum, raint_sum
   REAL                                      :: dmumax, raincmax, rainncmax,      &
                          snowhmax
   REAL                                      :: xtimep
   LOGICAL, EXTERNAL                         :: wrf_dm_on_monitor
   CHARACTER*256                             :: outstring
   CHARACTER*6                               :: grid_str

!!-------------------
!! CLWRF-UC Nov.09

   CHARACTER (LEN=80)                        :: timestr

   REAL, DIMENSION( ims:ime , jms:jme ),                                          & 
                          INTENT(IN)         :: t2, q2, u10, v10 
   REAL, DIMENSION( ims:ime , jms:jme ),                                          &
                          INTENT(OUT)        :: t2clmin, t2clmax, tt2clmin,       &
                          tt2clmax, t2clmean, t2clstd,                            & 
                          q2clmin, q2clmax, tq2clmin, tq2clmax, q2clmean, q2clstd,&
                          u10clmax, v10clmax, spduv10clmax, tspduv10clmax,        &
                          u10clmean, v10clmean, spduv10clmean,                    &
                          u10clstd, v10clstd, spduv10clstd, skintempclmin,        &
                          skintempclmax, tskintempclmin, tskintempclmax,          &
                          skintempclmean, skintempclstd
   REAL, DIMENSION( ims:ime , jms:jme ),                                          &
                          INTENT(OUT)        :: raincclmax, rainncclmax,          &
                          traincclmax, trainncclmax, raincclmean, rainncclmean,   & 
                          raincclstd, rainncclstd 
   REAL, PARAMETER                           :: minimum0= 1000000.,               &
                          maximum0= -1000000. 
   REAL                                      :: value
   INTEGER, INTENT(IN)                       :: clwrfH
   CHARACTER (LEN=1024)                      :: message
   REAL, SAVE                                :: nsteps
   LOGICAL                                   :: is_restart

!-----------------------------------------------------------------
! Compute minutes from reference times clwrfH

! Initialize [var] values
! SET START AND END POINTS FOR TILES
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )

!  IF ( MOD(NINT(XTIME), clwrfH) == 0 ) THEN
  IF (( MOD(NINT(XTIME*60./dt),NINT(clwrfH*60./dt)) == 0) .AND. (.NOT.is_restart)) THEN
    DO ij = 1 , num_tiles
      IF  ( wrf_dm_on_monitor() ) THEN
        WRITE(message, *)'CLWRFdiag - T2; tile: ',ij,' T2clmin:',           & 
          t2clmin(i_start(ij)+(i_end(ij)-i_start(ij))/2,                    &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' T2clmax:',               &
          t2clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                    &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TT2clmin:',              &
          tt2clmin(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TT2clmax:',              &
          tt2clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' T2clmean:',              &
          t2clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' T2clstd:',               &
          t2clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2,                    &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
        WRITE(message, *)'CLWRFdiag - Q2; tile: ',ij,' Q2clmin:',           &
          q2clmin(i_start(ij)+(i_end(ij)-i_start(ij))/2,                    &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' Q2clmax:',               &
          q2clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                    &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TQ2clmin:',              &
          tq2clmin(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TQ2clmax:',              &
          tq2clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' Q2clmean:',              &
          q2clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' Q2clstd:',               &
          q2clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2,                    &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
        WRITE(message, *)'CLWRFdiag - WINDSPEED; tile: ',ij,' U10clmax:',   &
          u10clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' V10clmax:',              &
          v10clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SPDUV10clmax:',          &
          spduv10clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,               &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TSPDUV10clmax:',         &
          tspduv10clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,              &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' U10clmean:',             &
          u10clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2,                  &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' V10clmean:',             &
          v10clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2,                  &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SPDUV10clmean:',         &
          spduv10clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2,              &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' U10clstd:',              &
          u10clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' V10clstd:',              &
          v10clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2,                   &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SPDUV10clstd:',          &
          spduv10clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2,               &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
        WRITE(message, *)'CLWRFdiag - RAIN; tile: ',ij,' RAINCclmax:',      &
          raincclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                 &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINNCclmax:',           &
          rainncclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TRAINCclmax:',           &
          traincclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,                &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TRAINNCclmax:',          &
          trainncclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,               &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINCclmean:',           &
          raincclmean(i_start(ij)+(i_end(ij)-i_start(ij))/2,                &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINNCclmean:',          &
          rainncclmean(i_start(ij)+(i_end(ij)-i_start(ij))/2,               &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINCclstd:',            &
          raincclstd(i_start(ij)+(i_end(ij)-i_start(ij))/2,                 &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINNCclstd:',           &
          rainncclstd(i_start(ij)+(i_end(ij)-i_start(ij))/2,                &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
        WRITE(message,*)'CLWRFdiag - SKINTEMP; tile: ',ij,' SKINTEMPclmin:',&
          skintempclmin(i_start(ij)+(i_end(ij)-i_start(ij))/2,              &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SKINTEMPclmax:',         &
          skintempclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,              &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TSKINTEMPclmin:',        &
          tskintempclmin(i_start(ij)+(i_end(ij)-i_start(ij))/2,             &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TSKINTEMPclmax:',        &
          tskintempclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2,             &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SKINTEMPclmean:',        &
          skintempclmean(i_start(ij)+(i_end(ij)-i_start(ij))/2,             &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SKINTEMPclstd:',         &
          skintempclstd(i_start(ij)+(i_end(ij)-i_start(ij))/2,              &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
      ENDIF
      DO j = j_start(ij), j_end(ij)
        DO i = i_start(ij), i_end(ij)
          t2clmin(i,j)=t2(i,j)
          t2clmax(i,j)=t2(i,j)
          t2clmean(i,j)=t2(i,j)
          t2clstd(i,j)=t2(i,j)*t2(i,j)
          q2clmin(i,j)=q2(i,j)
          q2clmax(i,j)=q2(i,j)
          q2clmean(i,j)=q2(i,j)
          q2clstd(i,j)=q2(i,j)*q2(i,j)
          spduv10clmax(i,j)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
          u10clmean(i,j)=u10(i,j)
          v10clmean(i,j)=v10(i,j)
          spduv10clmean(i,j)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
          u10clstd(i,j)=u10(i,j)*u10(i,j)
          v10clstd(i,j)=v10(i,j)*v10(i,j)
          spduv10clstd(i,j)=u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j)
          raincclmax(i,j)=raincv(i,j)/dt
          rainncclmax(i,j)=rainncv(i,j)/dt
          raincclmean(i,j)=raincv(i,j)/dt
          rainncclmean(i,j)=rainncv(i,j)/dt
          raincclstd(i,j)=(raincv(i,j)/dt)*(raincv(i,j)/dt)
          rainncclstd(i,j)=(rainncv(i,j)/dt)*(rainncv(i,j)/dt)
          skintempclmin(i,j)=skintemp(i,j)
          skintempclmax(i,j)=skintemp(i,j)
          skintempclmean(i,j)=skintemp(i,j)
          skintempclstd(i,j)=skintemp(i,j)*skintemp(i,j)
!          nsteps=0.
        ENDDO
      ENDDO
  ENDDO
    nsteps=clwrfH*60./dt
  ELSE
    xtimep = xtime + dt/60.   ! value at end of timestep for time info
!    nsteps=nsteps+1.
    nsteps=clwrfH*60./dt
!      DO j = j_start(ij), j_end(ij)
!        DO i = i_start(ij), i_end(ij)
!     DO j = jps, jpe
!       DO i = ips, ipe
! Temperature
          CALL varstatistics(t2,xtimep,ime-ims+1,jme-jms+1,t2clmin,t2clmax,   &
            tt2clmin,tt2clmax,t2clmean,t2clstd)
! Water vapor mixing ratio
          CALL varstatistics(q2,xtimep,ime-ims+1,jme-jms+1,q2clmin,q2clmax,   &
            tq2clmin,tq2clmax,q2clmean,q2clstd)
! Wind speed
          CALL varstatisticsWIND(u10,v10,xtimep,ime-ims+1,jme-jms+1,u10clmax, &
            v10clmax,spduv10clmax,tspduv10clmax,u10clmean,v10clmean,         &
            spduv10clmean,u10clstd,v10clstd,spduv10clstd)
! Precipitation flux
          CALL varstatisticsMAX(raincv/dt,xtimep,ime-ims+1,jme-jms+1,         &
            raincclmax,traincclmax,raincclmean,raincclstd) 
          CALL varstatisticsMAX(rainncv/dt,xtimep,ime-ims+1,jme-jms+1,        &
            rainncclmax,trainncclmax,rainncclmean,rainncclstd)
! Skin Temperature 
          CALL varstatistics(skintemp,xtimep,ime-ims+1,jme-jms+1,skintempclmin,&
            skintempclmax, tskintempclmin,tskintempclmax,skintempclmean,      &
            skintempclstd)

!          IF (MOD(NINT(XTIME),clwrfH) == 0) THEN
!          IF (MOD(NINT(XTIME+dt/60.),clwrfH) == 0) THEN
           IF ((MOD(NINT((XTIME+dt/60.)*60./dt),NINT(clwrfH*60./dt)) == 0)) THEN
             IF  ( wrf_dm_on_monitor() ) PRINT *,'nsteps=',nsteps,' xtime:',  &
               xtime,' clwrfH:',clwrfH 
               t2clmean=t2clmean/nsteps
               t2clstd=SQRT(t2clstd/nsteps-t2clmean**2.)
               q2clmean=q2clmean/nsteps
               q2clstd=SQRT(q2clstd/nsteps-q2clmean**2.)
               u10clmean=u10clmean/nsteps
               v10clmean=v10clmean/nsteps
               spduv10clmean=spduv10clmean/nsteps
               u10clstd=SQRT(u10clstd/nsteps-u10clmean**2.)
               v10clstd=SQRT(v10clstd/nsteps-v10clmean**2.)
               spduv10clstd=SQRT(spduv10clstd/nsteps-                        &
                 spduv10clmean**2)
               raincclmean=raincclmean/nsteps
               rainncclmean=rainncclmean/nsteps
               raincclstd=SQRT(raincclstd/nsteps-raincclmean**2.)
               rainncclstd=SQRT(rainncclstd/nsteps-rainncclmean**2.)
               skintempclmean=skintempclmean/nsteps
              skintempclstd=SQRT(skintempclstd/nsteps-skintempclmean**2.)
            END IF
!        ENDDO
!      ENDDO
  ENDIF
!  !$OMP END PARALLEL DO

   END SUBROUTINE clwrf_output_calc

! UC.CLWRF Nov.09
SUBROUTINE varstatisticsWIND(varu, varv, tt, dx, dy, varumax, varvmax,       &
  varuvmax, tvaruvmax, varumean, varvmean, varuvmean, varustd, varvstd,     & 
  varuvstd) 
! Subroutine to compute variable statistics for a wind somponents 

IMPLICIT NONE

INTEGER                                                        :: i, j
INTEGER, INTENT(IN)                                            :: dx, dy
REAL, DIMENSION(dx,dy), INTENT(IN)                             :: varu, varv
REAL, INTENT(IN)                                               :: tt
REAL, DIMENSION(dx,dy), INTENT(INOUT)                          :: varumax,   &
  varvmax, varuvmax, tvaruvmax, varumean, varvmean, varuvmean, varustd,      & 
  varvstd, varuvstd
REAL                                                           :: varuv

DO i=1,dx
  DO j=1,dy
    varuv=sqrt(varu(i,j)*varu(i,j)+varv(i,j)*varv(i,j))
      IF (varuv > varuvmax(i,j)) THEN
        varumax(i,j)=varu(i,j)
        varvmax(i,j)=varv(i,j)
        varuvmax(i,j)=varuv
        tvaruvmax(i,j)=tt
      END IF
    varuvmean(i,j)=varuvmean(i,j)+varuv
    varuvstd(i,j)=varuvstd(i,j)+varuv**2
  END DO
END DO
varumean=varumean+varu
varvmean=varvmean+varv
varustd=varustd+varu**2
varvstd=varvstd+varv**2

END SUBROUTINE varstatisticsWIND

SUBROUTINE varstatisticsMAX(var, tt, dx, dy, varmax, tvarmax, varmean,       &
   varstd)
! Subroutine to compute variable statistics for a max only variable values

IMPLICIT NONE

INTEGER                                                        :: i,j
INTEGER, INTENT(IN)                                            :: dx, dy
REAL, DIMENSION(dx,dy), INTENT(IN)                             :: var
REAL, INTENT(IN)                                               :: tt
REAL, DIMENSION(dx,dy), INTENT(INOUT)                          :: varmax,    &
  tvarmax, varmean, varstd

DO i=1,dx
  DO j=1,dy
    IF (var(i,j) > varmax(i,j)) THEN
      varmax(i,j)=var(i,j)
      tvarmax(i,j)=tt
    END IF
  END DO
END DO
varmean=varmean+var
varstd=varstd+var**2

END SUBROUTINE varstatisticsMAX 

SUBROUTINE varstatistics(var, tt, dx, dy, varmin, varmax, tvarmin, tvarmax,  & 
  varmean, varstd) 
! Subroutine to compute variable statistics

IMPLICIT NONE

INTEGER                                                        :: i,j
INTEGER, INTENT(IN)                                            :: dx, dy
REAL, DIMENSION(dx,dy), INTENT(IN)                             :: var
REAL, INTENT(IN)                                               :: tt
REAL, DIMENSION(dx,dy), INTENT(INOUT)                          :: varmin,    &
  varmax, tvarmin, tvarmax, varmean, varstd

DO i=1,dx
  DO j=1,dy
    IF (var(i,j) < varmin(i,j)) THEN
      varmin(i,j)=var(i,j)
      tvarmin(i,j)=tt
    END IF
    IF (var(i,j) > varmax(i,j)) THEN
      varmax(i,j)=var(i,j)
      tvarmax(i,j)=tt
    END IF
  END DO
END DO
varmean=varmean+var
varstd=varstd+var**2

END SUBROUTINE varstatistics

SUBROUTINE pld ( u,v,w,t,qv,zp,zb,pp,pb,p,pw,                    &
                 msfux,msfuy,msfvx,msfvy,msftx,msfty,            &
                 f,e,                                            &
                 use_tot_or_hyd_p,missing,                       &
                 num_press_levels,max_press_levels,press_levels, &
                 p_pl,u_pl,v_pl,t_pl,rh_pl,ght_pl,s_pl,td_pl,    &
                 ids,ide, jds,jde, kds,kde,                      &
                 ims,ime, jms,jme, kms,kme,                      &
                 its,ite, jts,jte, kts,kte                       )

   USE module_model_constants

   IMPLICIT NONE


   !  Input variables

   INTEGER, INTENT(IN   )                                          :: ids,ide, jds,jde, kds,kde, &
                                                                      ims,ime, jms,jme, kms,kme, &
                                                                      its,ite, jts,jte, kts,kte
   REAL   , INTENT(IN   ) , DIMENSION(ims:ime , jms:jme)           :: msfux,msfuy,msfvx,msfvy,msftx,msfty, &
                                                                      f,e
   INTEGER, INTENT(IN   )                                          :: use_tot_or_hyd_p
   REAL   , INTENT(IN   )                                          :: missing
   REAL   , INTENT(IN   ) , DIMENSION(ims:ime , kms:kme , jms:jme) :: u,v,w,t,qv,zp,zb,pp,pb,p,pw
   INTEGER, INTENT(IN   )                                          :: num_press_levels, max_press_levels
   REAL   , INTENT(IN   ) , DIMENSION(max_press_levels)            :: press_levels

   !  Output variables

   REAL   , INTENT(  OUT) ,  DIMENSION(num_press_levels)                     :: p_pl
   REAL   , INTENT(  OUT) ,  DIMENSION(ims:ime , num_press_levels , jms:jme) :: u_pl,v_pl,t_pl,rh_pl,ght_pl,s_pl,td_pl

   !  Local variables

   REAL, PARAMETER :: eps = 0.622, t_kelvin = svpt0 , s1 = 243.5, s2 = svp2 , s3 = svp1*10., s4 = 611.0, s5 = 5418.12

   INTEGER :: i, j, ke, kp, ke_h, ke_f
   REAL    :: pu, pd, pm , &
              tu, td     , &
              su, sd     , &
              uu, ud     , &
              vu, vd     , &
              zu, zd     , &
              qu, qd, qm , &
              eu, ed, em , &
              du, dd
   REAL    :: es, qs

   !  Silly, but transfer the small namelist.input array into the grid structure for output purposes.

   DO kp = 1 , num_press_levels
      p_pl(kp) = press_levels(kp)
   END DO

   !  Initialize pressure level data to un-initialized

   DO j = jts , jte
      DO kp = 1 , num_press_levels
         DO i = its , ite
            u_pl  (i,kp,j) = missing
            v_pl  (i,kp,j) = missing
            t_pl  (i,kp,j) = missing
            rh_pl (i,kp,j) = missing
            ght_pl(i,kp,j) = missing
            s_pl  (i,kp,j) = missing
            td_pl (i,kp,j) = missing
         END DO
      END DO
   END DO

   !  Loop over each i,j location

   j_loop : DO j = jts , MIN(jte,jde-1)
      i_loop : DO i = its , MIN(ite,ide-1)

         !  For each i,j location, loop over the selected pressure levels to find

         ke_h = kts
         ke_f = kts
         kp_loop : DO kp = 1 , num_press_levels

            !  For this particular i,j and pressure level, find the eta levels that surround this point
            !  on half-levels.

            ke_loop_half : DO ke = ke_h , kte-2

               IF      ( use_tot_or_hyd_p .EQ. 1 ) THEN     !  total pressure
                  pu = pp(i,ke+1,j)+pb(i,ke+1,j)
                  pd = pp(i,ke  ,j)+pb(i,ke  ,j)
               ELSE IF ( use_tot_or_hyd_p .EQ. 2 ) THEN     !  hydrostatic pressure
                  pu = p(i,ke+1,j)
                  pd = p(i,ke  ,j)
               END IF
               pm = p_pl(kp)

               IF ( ( pd .GE. pm ) .AND. &
                    ( pu .LT. pm ) ) THEN

                  !  Found trapping pressure: up, middle, down.  We are doing first order interpolation.  
                  !  Now we just put in a list of diagnostics for this level.

                  !  1. Temperature (K)

                  tu = (t(i,ke+1,j)+t0)*(pu/p1000mb)**rcp
                  td = (t(i,ke  ,j)+t0)*(pd/p1000mb)**rcp
                  t_pl(i,kp,j) = ( tu * (pm-pd) + td * (pu-pm) ) / (pu-pd)

                  !  2. Speed (m s-1)

                  su = 0.5 * SQRT ( ( u(i,ke+1,j)+u(i+1,ke+1,j) )**2 + ( v(i,ke+1,j)+v(i,ke+1,j+1) )**2 ) 
                  sd = 0.5 * SQRT ( ( u(i,ke  ,j)+u(i+1,ke  ,j) )**2 + ( v(i,ke  ,j)+v(i,ke  ,j+1) )**2 ) 
                  s_pl(i,kp,j) = ( su * (pm-pd) + sd * (pu-pm) ) / (pu-pd)

                  !  3. U and V (m s-1)

                  uu = 0.5 *        ( u(i,ke+1,j)+u(i+1,ke+1,j) )
                  ud = 0.5 *        ( u(i,ke  ,j)+u(i+1,ke  ,j) )
                  u_pl(i,kp,j) = ( uu * (pm-pd) + ud * (pu-pm) ) / (pu-pd)

                  vu = 0.5 *                                           ( v(i,ke+1,j)+v(i,ke+1,j+1) )
                  vd = 0.5 *                                           ( v(i,ke  ,j)+v(i,ke  ,j+1) )
                  v_pl(i,kp,j) = ( vu * (pm-pd) + vd * (pu-pm) ) / (pu-pd)

                  !  4. Dewpoint (K) - Use Bolton's approximation

                  qu = MAX(qv(i,ke+1,j),0.)
                  qd = MAX(qv(i,ke  ,j),0.)
                  eu = qu * pu * 0.01 / ( eps + qu )       ! water vapor pressure in mb.
                  ed = qd * pd * 0.01 / ( eps + qd )       ! water vapor pressure in mb.
                  eu = max(eu, 0.001)
                  ed = max(ed, 0.001)

                  du = t_kelvin + ( s1 / ((s2 / log(eu/s3)) - 1.0) )
                  dd = t_kelvin + ( s1 / ((s2 / log(ed/s3)) - 1.0) )
                  td_pl(i,kp,j) = ( du * (pm-pd) + dd * (pu-pm) ) / (pu-pd)

                  !  5. Relative humidity (%)

                  qm = ( qu * (pm-pd) + qd * (pu-pm) ) / (pu-pd)                           ! qvapor at the pressure level.
                  es = s4 * exp(s5 * (1.0 / 273.0 - 1.0 / t_pl(i,kp,j)))
                  qs = eps * es / (pm - es)
                  rh_pl(i,kp,j)   = qm / qs * 100.

                  !em = qm * pm * 0.01 / ( eps + qm )                                       ! water vapor pressure at the level.
                  !es = s3 * exp( s2 * (t_pl(i,kp,j) - t_kelvin)/(t_pl(i,kp,j) - s4) )      ! sat vapor pressure over liquid water in mb.
                  !rh_pl(i,kp,j) = 100. * em * ( pm * 0.01 - es ) / ( es * ( pm * 0.01 - em ) )

                  ke_h = ke
                  EXIT ke_loop_half
               END IF
            END DO ke_loop_half

            ke_loop_full : DO ke = ke_f , kte-1
               IF ( ( pw(i,ke  ,j) .GE. p_pl(kp) ) .AND. &
                    ( pw(i,ke+1,j) .LT. p_pl(kp) ) ) THEN

                  !  Found trapping pressure: up, middle, down.  We are doing first order interpolation.

                  pu = LOG(pw(i,ke+1,j))
                  pm = LOG(p_pl(kp))
                  pd = LOG(pw(i,ke  ,j))

                  !  Now we just put in a list of diagnostics for this level.

                  !  1. Geopotential height (m)

                  zu = ( zp(i,ke+1,j)+zb(i,ke+1,j) ) / g
                  zd = ( zp(i,ke  ,j)+zb(i,ke  ,j) ) / g
                  ght_pl(i,kp,j) = ( zu * (pm-pd) + zd * (pu-pm) ) / (pu-pd)

                  ke_f = ke
                  EXIT ke_loop_full
               END IF
            END DO ke_loop_full

         END DO kp_loop
      END DO i_loop
   END DO j_loop

   

END SUBROUTINE pld


! David Dowell 5/28/11
! This subroutine uses dbzcalc2 to compute radar reflectivity and/or precipitation
! fall speed.  Currently, only the Thompson microphysics scheme is supported.
!
! A variable Thompson_scheme_version (not in the standard WRF release) must
! be defined in module_mp_thompson.F.
   SUBROUTINE radar_diagnostic_calc(                                  &
                      ims,ime, jms,jme, kms,kme,                      &
                      i_start,i_end,j_start,j_end,num_tiles           &
                     ,mp_physics, qc_curr, qnr_curr, t                &
                     ,qr_curr, qs_curr, qg_curr, rho_curr             &
                     ,refl_10cm, composite_refl_10cm, vt_dbz_wt       &
                                                                      )
!----------------------------------------------------------------------

  USE module_state_description, ONLY : THOMPSON
  USE module_mp_thompson, ONLY : Thompson_scheme_version

   IMPLICIT NONE
!======================================================================
! Definitions
!-----------
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- num_tiles     number of tiles
!
!======================================================================

   INTEGER,      INTENT(IN   )    ::                             &
                                      ims,ime, jms,jme, kms,kme, &
                                                      num_tiles

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                  &
     &           i_start,i_end,j_start,j_end

   INTEGER, INTENT(IN) :: mp_physics

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) ::   &
                                        qr_curr,qs_curr,qg_curr  &
                                            ,t,qc_curr,qnr_curr  &
                                                      ,rho_curr

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL, INTENT(INOUT) :: refl_10cm
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), OPTIONAL, INTENT(INOUT) :: vt_dbz_wt
   REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL, INTENT(INOUT) :: composite_refl_10cm

   INTEGER :: i,j,k,ij

   REAL, PARAMETER :: dbzmin=-20.

   REAL, ALLOCATABLE, DIMENSION(:,:,:) :: qra, qnr, qsn, qgr, qcw, tmk, rho, dbz, fallspd
   INTEGER :: mi, mj, mk, k1, k2


   IF ( mp_physics .EQ. THOMPSON ) THEN

     print *, 'microphysics scheme is Thompson'

     k1 = kms
     k2 = kme

     DO ij = 1 , num_tiles

       mi = i_end(ij) - i_start(ij) + 1
       mj = j_end(ij) - j_start(ij) + 1
       mk = k2 - k1 + 1

!       print *, 'ij, num_tiles, mi, mj, mk = ', ij, num_tiles, mi, mj, mk

       allocate (qra(mi+1,mj+1,mk))
       allocate (qnr(mi+1,mj+1,mk))
       allocate (qsn(mi+1,mj+1,mk))
       allocate (qgr(mi+1,mj+1,mk))
       allocate (qcw(mi+1,mj+1,mk))
       allocate (tmk(mi+1,mj+1,mk))
       allocate (rho(mi+1,mj+1,mk))
       allocate (dbz(mi+1,mj+1,mk))
       allocate (fallspd(mi+1,mj+1,mk))

!       print *, 'arrays allocated'

       ! Reorder indices and flip k index for arrays passed to dbzcalc2.
       DO j=j_start(ij),j_end(ij)
         DO k=k1,k2
           DO i=i_start(ij),i_end(ij)
             qra(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1) = qr_curr(i,k,j)
             qnr(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1) = qnr_curr(i,k,j)
             qsn(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1) = qs_curr(i,k,j)
             qgr(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1) = qg_curr(i,k,j)
             qcw(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1) = qc_curr(i,k,j)
             tmk(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1) = t(i,k,j)
             rho(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1) = rho_curr(i,k,j)
           ENDDO
         ENDDO
       ENDDO

!       print *, 'calling dbzcalc2'

       call dbzcalc2(Thompson_scheme_version,qra,qnr,qsn,qgr,qcw,tmk,rho,dbz,fallspd,mi+1,mj+1,mk)

!       print *, 'updating reflectivity and/or fallspeed arrays'

       DO i=i_start(ij),i_end(ij)
         DO j=j_start(ij),j_end(ij)
           IF ( PRESENT (composite_refl_10cm) ) THEN
             composite_refl_10cm(i,j) = dbzmin
           ENDIF
           DO k=k1,k2
             IF ( PRESENT (refl_10cm) ) THEN
               refl_10cm(i,k,j) = dbz(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1)
             ENDIF
             IF ( PRESENT (composite_refl_10cm) ) THEN
               IF ( dbz(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1) .GT. composite_refl_10cm(i,j) ) THEN
                 composite_refl_10cm(i,j) = dbz(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1)
               ENDIF
             ENDIF
             IF ( PRESENT (vt_dbz_wt) ) THEN
               vt_dbz_wt(i,k,j) = fallspd(i-i_start(ij)+1, j-j_start(ij)+1, k2-k+1)
             ENDIF
           ENDDO
         ENDDO
       ENDDO

!       print *, 'deallocating arrays'

       deallocate (qra)
       deallocate (qnr)
       deallocate (qsn)
       deallocate (qgr)
       deallocate (qcw)
       deallocate (tmk)
       deallocate (rho)
       deallocate (dbz)
       deallocate (fallspd)

     ENDDO

   ELSE

     print*, 'mp_physics = ', mp_physics
     CALL wrf_error_fatal('This microphysics scheme is not supported in radar_diagnostic_calc.')

   ENDIF   ! IF ( mp_physics .EQ. THOMPSON )

   END SUBROUTINE radar_diagnostic_calc


!+---+-----------------------------------------------------------------+
! code below provided by Greg Thompson 4/4/11
! conversion from Fortran 77 to Fortran 90 by David Dowell 4/14/11
! addition of Greg Thompson's fall-speed computations by David Dowell 6/1/11
!+---+-----------------------------------------------------------------+
!
!+---+-----------------------------------------------------------------+
!
      subroutine dbzcalc2(Thompson_scheme_version,                      &
                          qra,qnr,qsn,qgr,qcw,tmk,rho,dbz,vt_dbz,miy,mjx,mkzh)
!                         in0r,in0s,in0g,iliqskin)
!
!     This routine computes equivalent reflectivity factor (in dBZ) at
!     each model grid point.  In calculating Ze, RIP uses the same
!     assumptions as found in the new Thompson bulk microphysics
!     scheme of 2006.  The only difficult part is computing for
!     partially melted snow and graupel particles.  For this, the
!     code from Uli Blahak is simplified by assuming Rayleigh
!     scattering contribution only for 10 cm wavelength (NEXRAD) and
!     a simple function of meltwater coating the ice particles.

!+---+-----------------------------------------------------------------+
!.. Here are a few serious warnings.
!.. 1) Those users that modify parameters set at the top of my
!..  mp_thompson WRF routine MUST change to the same in this code or a
!..  mismatch will occur.  Big changes would likely result if users alter
!..  gamma shape parameters (mu_x) or intercept parameters.
!.. 2) This code is not intended for my mp_thompson WRF routine that was
!..  released prior to WRF-V3.1.  See preceeding subroutine for versions
!..  2.2, 2.2.1, 3.0 or 3.0.1.1 (before 2-moment rain).
!.. 3) The treatment of rain is simplest - just an integration of the
!..  explicit size distrib.  Snow & graupel without melting are also
!..  rather simple integrations if one assumes no water coats the
!..  particles when T<0C - which is what I did.  When snow/graupel exist
!..  above T=0C, it gets interesting (next item).
!.. 4) This code makes some rather simple assumptions about water
!..  coating on outside of frozen species (snow/graupel).  Fraction of
!..  meltwater is simply the ratio of mixing ratio below melting level
!..  divided by mixing ratio at level just above highest T>0C.  Also,
!..  immediately 70% of the melted water exists on the ice's surface
!..  and 30% is embedded within ice.  No water is "shed" at all in these
!..  assumptions.  The work-horse subroutines for reflectivity of melted
!..  particles comes from Ulrich Blahak (Univ. Karlsrhule in Germany):
!..  "rayleigh_soak_wetgraupel" (and related functions).  There are a
!..  large number of assumptions with his scheme that I simplified.
!.. 5) The code is quite slow because it does the reflectivity
!..  calculations based on 100 individual size bins of the distributions.
!.. 6) Ulrich has the code for treating Mie scattering as well as
!..  Rayleigh (which is all I'm using) so that other radar wavelengths
!..  could be considered in a future version.  Right now, only 10-cm
!..  (NEXRAD) is used herein.
!+---+-----------------------------------------------------------------+
!
      implicit none

! DCD 5/31/11
      CHARACTER*8 Thompson_scheme_version
      logical allow_wet_graupel
! DCD

      integer miy,mjx,mkzh
!      real in0r,in0s,in0g
!      integer iliqskin
      real qra(miy,mjx,mkzh), qnr(miy,mjx,mkzh), qcw(miy,mjx,mkzh), &
         qsn(miy,mjx,mkzh),qgr(miy,mjx,mkzh),tmk(miy,mjx,mkzh),     &
         rho(miy,mjx,mkzh),dbz(miy,mjx,mkzh)
      real vt_dbz(miy,mjx,mkzh)        ! reflectivity weighted precipitation fall speed (m/s)

      integer maxk
      parameter (maxk=100)
      double precision ilamg(maxk), N0_g(maxk), &
           mvd_r(maxk), ilamr(maxk), N0_r(maxk)
      real rs(maxk), rg(maxk), rr(maxk), nr(maxk), rc(maxk)
      real temp(maxk), smoz(maxk), smob(maxk), smoc(maxk)
      real ze_rain(maxk), ze_snow(maxk), ze_graupel(maxk)
      real rhof(maxk)
      real vts_dbz_wt, vtr_dbz_wt, vtg_dbz_wt

      real sa(10), sb(10), cse(12), csg(12), oams, obms
!      real WGAMMA
      real cre(13), crg(13), ore1, org1, org2, org3, oamr, obmr
      real cge(5), cgg(5), oge1, ogg1, ogg2, ogg3, oamg, obmg
      real ocms, ocmg
      real tc0, smo2, loga_, a_, b_
      real xDs, Mrat, Mrat2, ils1, ils2
      real t1_vts, t2_vts, t3_vts, t4_vts
      real xDg, Ds_m, Dg_m, r_mvd1, r_mvd2, xkrat
      double precision N0_exp, N0_min, lam_exp, lamg, lamr
      real z_e
      integer k_0, kbot, i, j, k, n
      logical melti

      REAL PI
      PARAMETER (PI = 3.1415926536)
      REAL R1
      PARAMETER (R1 = 1.e-12)
      REAL R2
      PARAMETER (R2 = 1.e-8)
      REAL rho_w
      PARAMETER (rho_w = 1000.0)
      REAL rho_i
      PARAMETER (rho_i = 890.0)
      REAL rho_g
      PARAMETER (rho_g = 400.0)
      REAL mu_r
      PARAMETER (mu_r = 0.0)
      REAL mu_g
      PARAMETER (mu_g = 0.0)
      REAL mu_s
      PARAMETER(mu_s = 0.6357)

      REAL gonv_min

      REAL gonv_max

      REAL am_r
      PARAMETER (am_r = PI*rho_w/6.0)
      REAL bm_r
      PARAMETER (bm_r = 3.0)
      REAL am_g
      PARAMETER (am_g = PI*rho_g/6.0)
      REAL bm_g
      PARAMETER (bm_g = 3.0)
      REAL am_s
      PARAMETER(am_s = 0.069)
      REAL bm_s
      PARAMETER(bm_s = 2.0)

      REAL Kap0
      PARAMETER(Kap0 = 490.6)
      REAL Kap1
      PARAMETER(Kap1 = 17.46)
      REAL Lam0
      PARAMETER(Lam0 = 20.78)
      REAL Lam1
      PARAMETER(Lam1 = 3.29)

!..Fallspeed power laws relations:  v = (av*D**bv)*exp(-fv*D)
!.. Rain from Ferrier (1994), ice, snow, and graupel from
!.. Thompson et al (2008). Coefficient fv is zero for graupel/ice.
      REAL, PARAMETER :: av_r = 4854.0
      REAL, PARAMETER :: bv_r = 1.0
      REAL, PARAMETER :: fv_r = 195.0
      REAL, PARAMETER :: av_s = 40.0
      REAL, PARAMETER :: bv_s = 0.55
      REAL, PARAMETER :: fv_s = 125.0
      REAL, PARAMETER :: av_g = 442.0
      REAL, PARAMETER :: bv_g = 0.89
!     REAL, PARAMETER :: av_i = 1847.5
!     REAL, PARAMETER :: bv_i = 1.0

!..Rho_not used in fallspeed relations (rho_not/rho)**.5 adjustment.
      REAL, PARAMETER :: rho_not = 101325.0/(287.05*298.0)

!..Various radar related variables
      INTEGER nbins
      PARAMETER(nbins = 50)
      DOUBLE PRECISION xDx(nbins+1)
      DOUBLE PRECISION Ds(nbins), Dg(nbins), dts(nbins), dtg(nbins)
      DOUBLE PRECISION fmelt_s, fmelt_g, cback, xx, eta, f_d
      REAL oM3, M0, slam1, slam2
      REAL D0s, D0g
      PARAMETER (D0s = 200.E-6, D0g = 250.E-6)
      REAL xslw1, ygra1, zans1

      DOUBLE PRECISION lamda_radar
      PARAMETER(lamda_radar = 0.10)                                     ! in meters
      DOUBLE PRECISION K_w, PI5, lamda4
      COMPLEX*16 m_w_0, m_i_0
      DOUBLE PRECISION simpson(nbins+1)
      DOUBLE PRECISION basis(3)
      DOUBLE PRECISION melt_outside_s, melt_outside_g
      PARAMETER (melt_outside_s = 0.7d0, melt_outside_g = 0.7d0)
                    
      INTEGER slen
      PARAMETER(slen = 20)
      CHARACTER*(slen)  mixingrulestring_s,  matrixstring_s,       &
                        inclusionstring_s,  hoststring_s,          &
                        hostmatrixstring_s, hostinclusionstring_s, &
                        mixingrulestring_g, matrixstring_g,        &
                        inclusionstring_g,  hoststring_g,          &
                        hostmatrixstring_g, hostinclusionstring_g

!      COMPLEX*16 m_complex_water_ray, m_complex_ice_maetzler

!..For snow moments conversions (from Field et al. 2005)
      DATA sa / 5.065339, -0.062659, -3.032362, 0.029469, -0.000285,  &
                0.31255,   0.000204,  0.003199, 0.0,      -0.015952/
      DATA sb / 0.476221, -0.015896,  0.165977, 0.007468, -0.000141,  &
                0.060366,  0.000079,  0.000594, 0.0,      -0.003577/

      if (mkzh .gt. maxk) then
         print*, ' Cannot continue since mkzh is greater than maxk ', &
                 mkzh, maxk
         stop ' ABORT; Increase maxk'
      endif


! DCD 5/31/11

      allow_wet_graupel = .false.

      if ( (Thompson_scheme_version .eq. '3.2') .or. (Thompson_scheme_version .eq. '3.3') ) then

!       (DCD) gonv_min = 1.E4 in the WRF 3.2 and 3.3 Thompson schemes, but the value has been increased here
!             for the reflectivity calculation to give more realistic (lower) reflectivity values in
!             regions of high graupel mixing ratio.
!       gonv_min = 1.E4
        gonv_min = 1.E5

      else if (Thompson_scheme_version .eq. '3.1.9') then

!       Version 3.1.9 is a non-official release but refers here to code provided by Greg Thompson
!       to Glen Romine on 19 October 2010 for STEP experiments.

        gonv_min = 1.E4
        allow_wet_graupel = .true.

      else

        print*, 'Thompson_scheme_version = ', Thompson_scheme_version
        print*, 'This version is not supported in dbzcalc2.'
        stop

      endif

      if ( (Thompson_scheme_version .eq. '3.2') .or. (Thompson_scheme_version .eq. '3.1.9') ) then
        gonv_max = 3.E6
      else if (Thompson_scheme_version .eq. '3.3') then
        gonv_max = 1.E6
      else
        print*, 'Thompson_scheme_version = ', Thompson_scheme_version
        print*, 'This version is not supported in dbzcalc2.'
        stop
      endif
! DCD 5/31/11

!..Precompute gamma values and exponents used later (rain, snow, graupel).
      cre(1) = bm_r + 1.
      cre(2) = mu_r + 1.
      cre(3) = bm_r + mu_r + 1.
      cre(4) = bm_r*2. + mu_r + 1.
      cre(5) = mu_r + bv_r + 1.
      cre(6) = bm_r + mu_r + bv_r + 1.
      cre(7) = bm_r*0.5 + mu_r + bv_r + 1.
      cre(8) = bm_r + mu_r + bv_r + 3.
      cre(9) = mu_r + bv_r + 3.
      cre(10) = mu_r + 2.
      cre(11) = 0.5*(bv_r + 5. + 2.*mu_r)
      cre(12) = bm_r*0.5 + mu_r + 1.
      cre(13) = bm_r*2. + mu_r + bv_r + 1
      do n = 1, 13
         crg(n) = WGAMMA(cre(n))
      enddo
      ore1 = 1./cre(1)
      org1 = 1./crg(1)
      org2 = 1./crg(2)
      org3 = 1./crg(3)
      oamr = 1./am_r
      obmr = 1./bm_r

      cse(1) = bm_s + 1.
      cse(2) = bm_s + 2.
      cse(3) = bm_s*2.
      cse(4) = bm_s + bv_s + 1.
      cse(5) = bm_s*2. + bv_s + 1.
      cse(6) = bm_s*2. + 1.
      cse(7) = bm_s + mu_s + 1.
      cse(8) = bm_s + mu_s + 2.
      cse(9) = bm_s + mu_s + 3.
      cse(10) = bm_s + mu_s + bv_s + 1.
      cse(11) = bm_s*2. + mu_s + bv_s + 1.
      cse(12) = bm_s*2. + mu_s + 1.
      do n = 1, 12
         csg(n) = WGAMMA(cse(n))
      enddo

      oams = 1./am_s
      obms = 1./bm_s
      ocms = oams**obms

      cge(1) = bm_g + 1.
      cge(2) = mu_g + 1.
      cge(3) = bm_g + mu_g + 1.
      cge(4) = bm_g*2. + mu_g + 1.
      cge(5) = bm_g*2. + mu_g + bv_g + 1.
      do n = 1, 5
         cgg(n) = WGAMMA(cge(n))
      enddo
      oge1 = 1./cge(1)
      ogg1 = 1./cgg(1)
      ogg2 = 1./cgg(2)
      ogg3 = 1./cgg(3)
      oamg = 1./am_g
      obmg = 1./bm_g
      ocmg = oamg**obmg

      basis(1) = 1.d0/3.d0
      basis(2) = 4.d0/3.d0
      basis(3) = 1.d0/3.d0

      PI5 = PI*PI*PI*PI*PI
      lamda4 = lamda_radar*lamda_radar*lamda_radar*lamda_radar
      m_w_0 = m_complex_water_ray (lamda_radar, 0.0d0)
      m_i_0 = m_complex_ice_maetzler (lamda_radar, 0.0d0)
      K_w = (ABS( (m_w_0*m_w_0 - 1.0) /(m_w_0*m_w_0 + 2.0) ))**2
      
      do n = 1, nbins+1
         simpson(n) = 0.0d0
      enddo
      do n = 1, nbins-1, 2
         simpson(n) = simpson(n) + basis(1)
         simpson(n+1) = simpson(n+1) + basis(2)
         simpson(n+2) = simpson(n+2) + basis(3)
      enddo

      do n = 1, slen
         mixingrulestring_s(n:n) = char(0)
         matrixstring_s(n:n) = char(0)
         inclusionstring_s(n:n) = char(0)
         hoststring_s(n:n) = char(0)
         hostmatrixstring_s(n:n) = char(0)
         hostinclusionstring_s(n:n) = char(0)
         mixingrulestring_g(n:n) = char(0)
         matrixstring_g(n:n) = char(0)
         inclusionstring_g(n:n) = char(0)
         hoststring_g(n:n) = char(0)
         hostmatrixstring_g(n:n) = char(0)
         hostinclusionstring_g(n:n) = char(0)
      enddo

      mixingrulestring_s = 'maxwellgarnett'
      hoststring_s = 'air'
      matrixstring_s = 'water'
      inclusionstring_s = 'spheroidal'
      hostmatrixstring_s = 'icewater'
      hostinclusionstring_s = 'spheroidal'

      mixingrulestring_g = 'maxwellgarnett'
      hoststring_g = 'air'
      matrixstring_g = 'water'
      inclusionstring_g = 'spheroidal'
      hostmatrixstring_g = 'icewater'
      hostinclusionstring_g = 'spheroidal'

      do k = 1, mkzh
      do j = 1, mjx-1
      do i = 1, miy-1
         dbz(i,j,k) = 0.0
      enddo
      enddo
      enddo

!     print*, ' DEBUG: ', tmk(140,112,1), tmk(140,112,mkzh)

      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1

      do k = 1, mkzh
        rc(k) = MAX(R1, qcw(i,j,k)*rho(i,j,k))
        rr(k) = MAX(R1, qra(i,j,k)*rho(i,j,k))
        nr(k) = MAX(1., qnr(i,j,k)*rho(i,j,k))
        rg(k) = MAX(R1, qgr(i,j,k)*rho(i,j,k))
        rs(k) = MAX(R1, qsn(i,j,k)*rho(i,j,k))
        temp(k) = tmk(i,j,k)
        rhof(k) = sqrt(RHO_NOT/rho(i,j,k))
        ze_rain(k) = 1.e-22
        ze_snow(k) = 1.e-22
        ze_graupel(k) = 1.e-22
      enddo

!..Contribution to reflectivity from graupel.
!..The line below for N0_exp should change according to which version WRF
!.. and should match how the intercept parameter is set within
!.. module_mp_thompson.F

      N0_min = gonv_max
      do k = 1, mkzh

!..      N0_exp = 100.0*rho(i,j,k)/rg(k)

         if (Thompson_scheme_version .eq. '3.2') then

           N0_exp = (gonv_max-gonv_min)*0.5D0      &
                  * tanh((0.15E-3-rg(k))/0.15E-3)  &
                  + (gonv_max+gonv_min)*0.5D0

         else if (Thompson_scheme_version .eq. '3.1.9') then

           N0_exp = (gonv_max-gonv_min)*0.5D0                &
                  * tanh((0.01E-3-(rc(k)+rr(k)))/0.75E-3)    &
                  + (gonv_max+gonv_min)*0.5D0

         else if (Thompson_scheme_version .eq. '3.3') then

           if (temp(k).lt.273.15 .and. (rc(k)+rr(k)).gt.1.E-5) then
              xslw1 = 5. + alog10(max(1.E-5, min(1.E-2, (rc(k)+rr(k)))))
           else
              xslw1 = 0.
           endif
           ygra1 = 5. + alog10(max(1.E-5, min(1.E-2, rg(k))))
           zans1 = 3.324 + (3./(5.*xslw1*ygra1               &
                 / (5.*xslw1+1.+0.25*ygra1)+1.+0.25*ygra1))
           N0_exp = 10.**(zans1)
           N0_exp = DMAX1(DBLE(gonv_min), DMIN1(N0_exp, DBLE(gonv_max)))
           N0_min = DMIN1(N0_exp, N0_min)
           N0_exp = N0_min

         else

           print*, 'Thompson_scheme_version = ', Thompson_scheme_version
           print*, 'This version is not supported in dbzcalc2.'
           stop

         endif

         lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
         lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
      enddo
      do k = 1, mkzh
         if (rg(k).gt. R1)                                       &
               ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                                    * (am_g/rho_i)*(am_g/rho_i)  &
                                    * N0_g(k)*cgg(4)*ilamg(k)**cge(4)
      enddo

!..Contribution to reflectivity from snow.
      do k = 1, mkzh
         if (.not. (rs(k).gt. R1) ) goto 992
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams 
         
!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2 = smob(k)
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s         &
               + sa(4)*tc0*bm_s + sa(5)*tc0*tc0            &
               + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s      &
               + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0   &
               + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_  
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s            &
               + sb(4)*tc0*bm_s + sb(5)*tc0*tc0            &
               + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s      &
               + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0   &
               + sb(10)*bm_s*bm_s*bm_s
            smo2 = (smob(k)/a_)**(1./b_)
         endif

!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1)             &
               + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0             &
               + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1)   &
               + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0  &
               + sa(10)*cse(1)*cse(1)*cse(1) 
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1)  &
              + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1)               &
              + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1)    &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a_ * smo2**b_

!..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3)               &
               + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0               &
               + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3)     &
               + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0    &
               + sa(10)*cse(3)*cse(3)*cse(3) 
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3)   &
              + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3)                &
              + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3)     &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
         smoz(k) = a_ * smo2**b_

         ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)             &
                                 * (am_s/rho_i)*(am_s/rho_i)*smoz(k)
 992  continue
      enddo

!..Contribution to reflectivity from rain.
      do k = 1, mkzh
         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
         ilamr(k) = 1./lamr
         N0_r(k) = nr(k)*org2*lamr**cre(2)
      enddo

      do k = 1, mkzh
         if (rr(k).gt. R1)                                        &
                ze_rain(k) = N0_r(k)*crg(4)*ilamr(k)**cre(4)
      enddo


!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+ 

      k_0 = mkzh
      melti = .false.
      do k = 2, mkzh
         if ( (temp(k).gt. 273.15) .and. (rr(k).gt. 0.0001e-3)    &
                   .and. ((rs(k-1)+rg(k-1)).gt. 0.01e-3) ) then
            k_0 = MIN(k-1, k_0)
            melti=.true.
            goto 135
         endif
      enddo
 135  continue

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+
      
      if ((.not. melti) .or. k_0.ge.mkzh) goto 991

!..Create bins of snow (from min diameter up to 2 cm).
      xDx(1) = D0s*1.0d0
      xDx(nbins+1) = 0.02d0
      do n = 2, nbins
         xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbins)               &
                  *DLOG(xDx(nbins+1)/xDx(1)) +DLOG(xDx(1)))
      enddo
      do n = 1, nbins
         Ds(n) = DSQRT(xDx(n)*xDx(n+1))
         dts(n) = xDx(n+1) - xDx(n)
      enddo

!..Create bins of graupel (from min diameter up to 5 cm).
      xDx(1) = D0g*1.0d0
      xDx(nbins+1) = 0.05d0
      do n = 2, nbins
         xDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbins)              &
                  *DLOG(xDx(nbins+1)/xDx(1)) +DLOG(xDx(1)))
      enddo
      do n = 1, nbins
         Dg(n) = DSQRT(xDx(n)*xDx(n+1))
         dtg(n) = xDx(n+1) - xDx(n)
      enddo

      do k = k_0+1, mkzh
      
!..Reflectivity contributed by melting snow
         fmelt_s = DMIN1(1.0d0-rs(k)/rs(k_0), 1.0d0)
         if (fmelt_s.gt.0.01d0 .and. fmelt_s.lt.0.99d0 .and. &
                        rs(k).gt.R1) then
          eta = 0.d0
          oM3 = 1./smoc(k)
          M0 = (smob(k)*oM3)
          Mrat = smob(k)*M0*M0*M0
          slam1 = M0 * Lam0
          slam2 = M0 * Lam1
          do n = 1, nbins
             xx = am_s * Ds(n)**bm_s
             call rayleigh_soak_wetgraupel(xx, DBLE(ocms), DBLE(obms),  &
                   fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar,  &
                   CBACK, mixingrulestring_s, matrixstring_s,           &
                   inclusionstring_s, hoststring_s,                     &
                   hostmatrixstring_s, hostinclusionstring_s,           &
                   PI, rho_w, rho_i, PI5, lamda4)
             f_d = Mrat*(Kap0*DEXP(-slam1*Ds(n))                        &
                   + Kap1*(M0*Ds(n))**mu_s * DEXP(-slam2*Ds(n)))
             eta = eta + f_d * CBACK * simpson(n) * dts(n)
          enddo
          ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
         endif


!..Reflectivity contributed by melting graupel

! DCD 6/1/11  Wet versus dry graupel is a significant source of uncertainty
!              in reflectivity calculations.  Graupel with a wet surface
!              has much higher reflectivity than graupel with a dry surface.
!              The allow_wet_graupel parameter controls whether the possibility
!              of wet graupel is considered in the reflectivity calculation.

         if (allow_wet_graupel) then

           fmelt_g = DMIN1(1.0d0-rg(k)/rg(k_0), 1.0d0)
           if (fmelt_g.gt.0.01d0 .and. fmelt_g.lt.0.99d0 .and.            &
                          rg(k).gt.R1) then
            eta = 0.d0
            lamg = 1./ilamg(k)
            do n = 1, nbins
               xx = am_g * Dg(n)**bm_g
               call rayleigh_soak_wetgraupel(xx, DBLE(ocmg), DBLE(obmg),  &
                     fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar,  &
                     CBACK, mixingrulestring_g, matrixstring_g,           &
                     inclusionstring_g, hoststring_g,                     &
                     hostmatrixstring_g, hostinclusionstring_g,           &
                     PI, rho_w, rho_i, PI5, lamda4)
               f_d = N0_g(k)*Dg(n)**mu_g * DEXP(-lamg*Dg(n))
               eta = eta + f_d * CBACK * simpson(n) * dtg(n)
            enddo
            ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
           endif

         endif

      enddo

 991  continue

!..Scale answer to typical dBZ.

!     if (i.eq.140 .and. j.eq.112) print*, ' DEBUG: test column'
      do k = 1, mkzh
         dbz(i,j,k) = MAX(dbz(i,j,k),                                     &
                 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.e18))
!       if (i.eq.140 .and. j.eq.112) then
!        print*, 10.*log10(ze_rain(k)*1.e18),     &
!                10.*log10(ze_snow(k)*1.e18),     &
!                10.*log10(ze_graupel(k)*1.e18),  &
!                dbz(i,j,k)
!       endif
      enddo

!..Reflectivity-weighted terminal velocity (snow, rain, graupel, mix).

     do k = 1, mkzh
       vt_dbz(i,j,k) = 1.E-3
       if (rs(k).gt.R2) then
         xDs = smoc(k) / smob(k)
         Mrat = 1./xDs
         ils1 = 1./(Mrat*Lam0 + fv_s)
         ils2 = 1./(Mrat*Lam1 + fv_s)
!mass wgt           t1_vts = Kap0*csg(4)*ils1**cse(4)
!mass wgt           t2_vts = Kap1*Mrat**mu_s*csg(10)*ils2**cse(10)
         t1_vts = Kap0*csg(5)*ils1**cse(5)
         t2_vts = Kap1*Mrat**mu_s*csg(11)*ils2**cse(11)
         ils1 = 1./(Mrat*Lam0)
         ils2 = 1./(Mrat*Lam1)
         t3_vts = Kap0*csg(6)*ils1**cse(6)
         t4_vts = Kap1*Mrat**mu_s*csg(12)*ils2**cse(12)
!mass wgt           t3_vts = Kap0*csg(1)*ils1**cse(1)
!mass wgt           t4_vts = Kap1*Mrat**mu_s*csg(7)*ils2**cse(7)
         vts_dbz_wt = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
       else
         vts_dbz_wt = 1.E-3
       endif

       if (rr(k).gt.R2) then
         lamr = 1./ilamr(k)
         vtr_dbz_wt = rhof(k)*av_r*crg(13)*(lamr+fv_r)**(-cre(13))      &
                    / (crg(4)*lamr**(-cre(4)))
       else
         vtr_dbz_wt = 1.E-3
       endif

       if (rg(k).gt.R2) then
         lamg = 1./ilamg(k)
         vtg_dbz_wt = rhof(k)*av_g*cgg(5)*lamg**(-cge(5))               &
                    / (cgg(4)*lamg**(-cge(4)))
       else
         vtg_dbz_wt = 1.E-3
       endif

       vt_dbz(i,j,k) = (vts_dbz_wt*ze_snow(k) + vtr_dbz_wt*ze_rain(k)      &
                        + vtg_dbz_wt*ze_graupel(k))                        &
                        / (ze_rain(k)+ze_snow(k)+ze_graupel(k))

      enddo

!
 1000 continue

      return
      end subroutine dbzcalc2



!+---+-----------------------------------------------------------------+

      REAL FUNCTION GAMMLN(XX)
!     --- RETURNS THE VALUE LN(GAMMA(XX)) FOR XX > 0.
      IMPLICIT NONE
      REAL XX
      DOUBLE PRECISION STP
      PARAMETER (STP = 2.5066282746310005D0)
      DOUBLE PRECISION SER,TMP,X,Y
      INTEGER J
      DOUBLE PRECISION COF(6)
      DATA COF /76.18009172947146D0, -86.50532032941677D0,          &
                       24.01409824083091D0, -1.231739572450155D0,   &
                      .1208650973866179D-2, -.5395239384953D-5/

      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO 11 J=1,6
        Y=Y+1.D0
        SER=SER+COF(J)/Y
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER/X)
      END FUNCTION GAMMLN
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+
      REAL FUNCTION WGAMMA(y)

      IMPLICIT NONE
      REAL y
!      REAL GAMMLN

      WGAMMA = EXP(GAMMLN(y))

      END FUNCTION WGAMMA
!+---+-----------------------------------------------------------------+

!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION m_complex_water_ray(lambda,T)

!      Complex refractive Index of Water as function of Temperature T
!      [deg C] and radar wavelength lambda [m]; valid for
!      lambda in [0.001,1.0] m; T in [-10.0,30.0] deg C
!      after Ray (1972)

      IMPLICIT NONE
      DOUBLE PRECISION T,lambda
      DOUBLE PRECISION epsinf,epss,epsr,epsi
      DOUBLE PRECISION alpha,lambdas,sigma,nenner
      COMPLEX*16 i
      REAL PI
      PARAMETER (PI = 3.1415926536)

      i = (0d0,1d0)

      epsinf  = 5.27137d0 + 0.02164740d0 * T - 0.00131198d0 * T*T
      epss    = 78.54d+0 * (1.0 - 4.579d-3 * (T - 25.0)        &
              + 1.190d-5 * (T - 25.0)*(T - 25.0)               &
              - 2.800d-8 * (T - 25.0)*(T - 25.0)*(T - 25.0))
      alpha   = -16.8129d0/(T+273.16) + 0.0609265d0
      lambdas = 0.00033836d0 * exp(2513.98d0/(T+273.16)) * 1e-2

      nenner = 1.d0+2.d0*(lambdas/lambda)**(1d0-alpha)*sin(alpha*PI*0.5)  &
             + (lambdas/lambda)**(2d0-2d0*alpha)
      epsr = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)     &
           * sin(alpha*PI*0.5)+1d0)) / nenner
      epsi = ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)              &
           * cos(alpha*PI*0.5)+0d0)) / nenner                             &
           + lambda*1.25664/1.88496
      
      m_complex_water_ray = SQRT(CMPLX(epsr,-epsi))
      
      return
      END FUNCTION m_complex_water_ray

!+---+-----------------------------------------------------------------+
      
      COMPLEX*16 FUNCTION m_complex_ice_maetzler(lambda,T)
      
!      complex refractive index of ice as function of Temperature T
!      [deg C] and radar wavelength lambda [m]; valid for
!      lambda in [0.0001,30] m; T in [-250.0,0.0] C
!      Original comment from the Matlab-routine of Prof. Maetzler:
!      Function for calculating the relative permittivity of pure ice in
!      the microwave region, according to C. Maetzler, "Microwave
!      properties of ice and snow", in B. Schmitt et al. (eds.) Solar
!      System Ices, Astrophys. and Space Sci. Library, Vol. 227, Kluwer
!      Academic Publishers, Dordrecht, pp. 241-257 (1998). Input:
!      TK = temperature (K), range 20 to 273.15
!      f = frequency in GHz, range 0.01 to 3000
         
      IMPLICIT NONE
      DOUBLE PRECISION T,lambda
      DOUBLE PRECISION f,c,TK,B1,B2,b,deltabeta,betam,beta,theta,alfa

      c = 2.99d8
      TK = T + 273.16
      f = c / lambda * 1d-9

      B1 = 0.0207
      B2 = 1.16d-11
      b = 335.0d0
      deltabeta = EXP(-10.02 + 0.0364*(TK-273.16))
      betam = (B1/TK) * ( EXP(b/TK) / ((EXP(b/TK)-1)**2) ) + B2*f*f
      beta = betam + deltabeta
      theta = 300. / TK - 1.
      alfa = (0.00504d0 + 0.0062d0*theta) * EXP(-22.1d0*theta)
      m_complex_ice_maetzler = 3.1884 + 9.1e-4*(TK-273.16)
      m_complex_ice_maetzler = m_complex_ice_maetzler            &
                             + CMPLX(0.0d0, (alfa/f + beta*f)) 
      m_complex_ice_maetzler = SQRT(CONJG(m_complex_ice_maetzler))
      
      return
      END FUNCTION m_complex_ice_maetzler

!+---+-----------------------------------------------------------------+

      subroutine rayleigh_soak_wetgraupel (x_g, a_geo, b_geo, fmelt,  &
                     meltratio_outside, m_w, m_i, lambda, C_back,     &
                     mixingrule,matrix,inclusion,                     &
                     host,hostmatrix,hostinclusion,                   &
                     PI, rho_w, rho_i, PI5, lamda4)

      IMPLICIT NONE

      DOUBLE PRECISION x_g, a_geo, b_geo, fmelt, lambda,  &
                                     meltratio_outside
      DOUBLE PRECISION C_back
      COMPLEX*16 m_w, m_i
      CHARACTER*(*)  mixingrule, matrix, inclusion,  &
                                     host, hostmatrix, hostinclusion
      REAL PI, rho_w, rho_i
      DOUBLE PRECISION PI5, lamda4

      COMPLEX*16 m_core, m_air
      DOUBLE PRECISION D_large, D_g, rhog, x_w, xw_a, fm, fmgrenz,  &
                       volg, vg, volair, volice, volwater,          &
                       meltratio_outside_grenz, mra
      INTEGER error
!      COMPLEX*16 get_m_mix_nested

!     refractive index of air:
      m_air = (1.0d0,0.0d0)

!     Limiting the degree of melting --- for safety: 
      fm = DMAX1(DMIN1(fmelt, 1.0d0), 0.0d0)
!     Limiting the ratio of (melting on outside)/(melting on inside):
      mra = DMAX1(DMIN1(meltratio_outside, 1.0d0), 0.0d0)

!    ! The relative portion of meltwater melting at outside should increase
!    ! from the given input value (between 0 and 1)
!    ! to 1 as the degree of melting approaches 1,
!    ! so that the melting particle "converges" to a water drop.
!    ! Simplest assumption is linear:
      mra = mra + (1.0d0-mra)*fm

      x_w = x_g * fm

      D_g = a_geo * x_g**b_geo

      if (D_g .ge. 1d-12) then

       vg = PI/6. * D_g**3
       rhog = DMAX1(DMIN1(x_g / vg, DBLE(rho_i)), 10.0d0)
       vg = x_g / rhog

       meltratio_outside_grenz = 1.0d0 - rhog / rho_w

       if (mra .le. meltratio_outside_grenz) then
!       !..In this case, it cannot happen that, during melting, all the
!       !.. air inclusions within the ice particle get filled with
!       !.. meltwater. This only happens at the end of all melting.
        volg = vg * (1.0d0 - mra * fm)
 
       else
!       !..In this case, at some melting degree fm, all the air
!       !.. inclusions get filled with meltwater.
        fmgrenz=(rho_i-rhog)/(mra*rho_i-rhog+rho_i*rhog/rho_w)

        if (fm .le. fmgrenz) then
!        !.. not all air pockets are filled:
         volg = (1.0 - mra * fm) * vg
        else
!        !..all air pockets are filled with meltwater, now the
!        !.. entire ice sceleton melts homogeneously:
         volg = (x_g - x_w) / rho_i + x_w / rho_w
        endif

       endif

       D_large  = (6.0 / PI * volg) ** (1./3.)
       volice = (x_g - x_w) / (volg * rho_i)
       volwater = x_w / (rho_w * volg)
       volair = 1.0 - volice - volwater
      
!      !..complex index of refraction for the ice-air-water mixture
!      !.. of the particle:
       m_core = get_m_mix_nested (m_air, m_i, m_w, volair, volice,       &
                         volwater, mixingrule, host, matrix, inclusion,  &
                         hostmatrix, hostinclusion, error)
       if (error .ne. 0) then
        C_back = 0.0d0
        return
       endif

!      !..Rayleigh-backscattering coefficient of melting particle: 
       C_back = (ABS((m_core**2-1.0d0)/(m_core**2+2.0d0)))**2  &
                * PI5 * D_large**6 / lamda4

      else
       C_back = 0.0d0
      endif

      return
      END subroutine rayleigh_soak_wetgraupel

!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION get_m_mix_nested (m_a, m_i, m_w, volair,     &
                     volice, volwater, mixingrule, host, matrix,       &
                     inclusion, hostmatrix, hostinclusion, cumulerror)

      IMPLICIT NONE

      DOUBLE PRECISION volice, volair, volwater
      COMPLEX*16 m_a, m_i, m_w
      CHARACTER*(*) mixingrule, host, matrix,              &
                    inclusion, hostmatrix, hostinclusion
      INTEGER cumulerror

      DOUBLE PRECISION vol1, vol2
      COMPLEX*16 mtmp
      INTEGER error
!      COMPLEX*16 get_m_mix

!     !..Folded: ( (m1 + m2) + m3), where m1,m2,m3 could each be
!     !.. air, ice, or water

      cumulerror = 0
      get_m_mix_nested = CMPLX(1.0d0,0.0d0)

      if (host .eq. 'air') then

       if (matrix .eq. 'air') then
!       write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
!       CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volice / MAX(volice+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, 0.0d0, vol1, vol2,       &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error
          
        if (hostmatrix .eq. 'air') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,          &
                         volair, (1.0d0-volair), 0.0d0, mixingrule, &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'icewater') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,          &
                         volair, (1.0d0-volair), 0.0d0, mixingrule, &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error
        else
!        write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',    &
!                          hostmatrix
!        CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'ice') then

       if (matrix .eq. 'ice') then
!       write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
!       CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volair+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, vol1, 0.0d0, vol2,         &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'ice') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,          &
                         (1.0d0-volice), volice, 0.0d0, mixingrule, &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airwater') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,          &
                         (1.0d0-volice), volice, 0.0d0, mixingrule, &
                         'air', hostinclusion, error)
         cumulerror = cumulerror + error          
        else
!        write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',    &
!                          hostmatrix
!        CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'water') then

       if (matrix .eq. 'water') then
!       write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
!       CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volice+volair,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, vol1, vol2, 0.0d0,        &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'water') then
         get_m_mix_nested = get_m_mix (2.0d0*m_a, mtmp, m_w,             &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule,  &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airice') then
         get_m_mix_nested = get_m_mix (2.0d0*m_a, mtmp, m_w,             &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule,  &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error          
        else
!        write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',         &
!                          hostmatrix
!        CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'none') then

       get_m_mix_nested = get_m_mix (m_a, m_i, m_w,                      &
                       volair, volice, volwater, mixingrule,             &
                       matrix, inclusion, error)
       cumulerror = cumulerror + error
        
      else
!      write(mp_debug,*) 'GET_M_MIX_NESTED: unknown matrix: ', host
!      CALL wrf_debug(150, mp_debug)
       cumulerror = cumulerror + 1
      endif

      IF (cumulerror .ne. 0) THEN
!      write(mp_debug,*) 'GET_M_MIX_NESTED: error encountered'
!      CALL wrf_debug(150, mp_debug)
       get_m_mix_nested = CMPLX(1.0d0,0.0d0)    
      endif

      return
      END FUNCTION get_m_mix_nested

!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION get_m_mix (m_a, m_i, m_w, volair, volice,    &
                     volwater, mixingrule, matrix, inclusion, error)

      IMPLICIT NONE

      DOUBLE PRECISION volice, volair, volwater
      COMPLEX*16 m_a, m_i, m_w
      CHARACTER*(*) mixingrule, matrix, inclusion
      INTEGER error
!      COMPLEX*16 m_complex_maxwellgarnett

      error = 0
      get_m_mix = CMPLX(1.0d0,0.0d0)

      if (mixingrule .eq. 'maxwellgarnett') then
       if (matrix .eq. 'ice') then
        get_m_mix = m_complex_maxwellgarnett(volice, volair, volwater,   &
                           m_i, m_a, m_w, inclusion, error)
       elseif (matrix .eq. 'water') then
        get_m_mix = m_complex_maxwellgarnett(volwater, volair, volice,   &
                           m_w, m_a, m_i, inclusion, error)
       elseif (matrix .eq. 'air') then
        get_m_mix = m_complex_maxwellgarnett(volair, volwater, volice,   &
                           m_a, m_w, m_i, inclusion, error)
       else
!       write(mp_debug,*) 'GET_M_MIX: unknown matrix: ', matrix
!       CALL wrf_debug(150, mp_debug)
        error = 1
       endif

      else
!      write(mp_debug,*) 'GET_M_MIX: unknown mixingrule: ', mixingrule
!      CALL wrf_debug(150, mp_debug)
       error = 2
      endif

      if (error .ne. 0) then
!      write(mp_debug,*) 'GET_M_MIX: error encountered'
!      CALL wrf_debug(150, mp_debug)
      endif

      return
      END FUNCTION get_m_mix

!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION m_complex_maxwellgarnett(vol1, vol2, vol3,    &
                     m1, m2, m3, inclusion, error)

      IMPLICIT NONE

      COMPLEX*16 m1, m2, m3
      DOUBLE PRECISION vol1, vol2, vol3
      CHARACTER*(*) inclusion

      COMPLEX*16 beta2, beta3, m1t, m2t, m3t
      INTEGER error

      error = 0

      if (DABS(vol1+vol2+vol3-1.0d0) .gt. 1d-6) then
!      write(mp_debug,*) 'M_COMPLEX_MAXWELLGARNETT: sum of the ',   &
!             'partial volume fractions is not 1...ERROR'
!      CALL wrf_debug(150, mp_debug)
       m_complex_maxwellgarnett=CMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m1t = m1**2
      m2t = m2**2
      m3t = m3**2

      if (inclusion .eq. 'spherical') then
       beta2 = 3.0d0*m1t/(m2t+2.0d0*m1t)
       beta3 = 3.0d0*m1t/(m3t+2.0d0*m1t)
      elseif (inclusion .eq. 'spheroidal') then
       beta2 = 2.0d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1.0d0)
       beta3 = 2.0d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1.0d0)
      else
!      write(mp_debug,*) 'M_COMPLEX_MAXWELLGARNETT: ',         &
!                        'unknown inclusion: ', inclusion
!      CALL wrf_debug(150, mp_debug)
       m_complex_maxwellgarnett=DCMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m_complex_maxwellgarnett =                                          &
       SQRT(((1.0d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) /   &
       (1.0d0-vol2-vol3+vol2*beta2+vol3*beta3))

      return
      END FUNCTION m_complex_maxwellgarnett

!+---+-----------------------------------------------------------------+
! end of code provided by Greg Thompson 4/4/11
!+---+-----------------------------------------------------------------+



END MODULE module_diagnostics
