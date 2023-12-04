; modified to allow for single month rather than season
;
; modified for use with CESM2-LE
; LLandrum
; Jan 2023
;
; calculates and plots latitude-summed AICE stds and trends for observations 
; (SSMI) and simulations and plots them as a function of longitude.
; This works for pop grid only in SH where lat/lon are regular
;
; LLandrum
; April 2014
;
; modified May, 2014 to allow for inclusion of LE control run in plots
;
; LLandrum
; May 2014
;
; modified May, 2015 for CESM High Resolution run
;
;----------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "/home/landrum/ncl/landrum_util.ncl"
load "/glade/u/home/landrum/ncl/landrum_util.ncl"

begin
; ==============================================================
; User defined parameters that specify region of globe to consider
; ==============================================================
  plot_CESM1_LE = False				; True to plot CESM1 data as well

  CESM2_case = "b.e21.BHISTcmip6.f09_g17.LE2-"
  CESM2_case_short = "CESM2-LE"
  case_numbers2 = new(50,integer)
  case_numbers2(:9) = ispan(1001,1181,20)
  case_numbers2(10:19) = 1231
  case_numbers2(20:29) = 1251
  case_numbers2(30:39) = 1281
  case_numbers2(40:49) = 1301

  case_numbers = new(50,integer)
  case_numbers(:9) = ispan(1,10,1) 
  case_numbers(10:19) = case_numbers(:9)
  case_numbers(20:29) = case_numbers(:9)
  case_numbers(30:39) = case_numbers(:9)
  case_numbers(40:49) = case_numbers(:9)

  scnumbers = sprinti("%0.3i",case_numbers)
;  case_names = case + scnumbers
  scnumbers2 = sprinti("%0.4i",case_numbers2)
  CESM2_case_names = CESM2_case + scnumbers2+"."+scnumbers
  nCESM2 = dimsizes(CESM2_case_names)

  CESM2dir = "/glade/scratch/landrum/CESM2-LE/aice/"
  CESM2_yr1 = 1850
  CESM2_yr2 = 2014

  syr1 =  sprinti("%0.4i", CESM2_yr1)   ; left zero fill
  syr2 =  sprinti("%0.4i", CESM2_yr2)   ; left zero fill

  fCESM2 = CESM2dir+"aice_sh_long_"+CESM2_case_names+"."+syr1+"01-"+syr2+"12.nc"

  detrend = True
  shade_STD = True				; True to shade std onto mean plot, false to show STD on separate plot

  plot_yr1 = 1979				; yrs over which to calculate means and trends and then to plot
  plot_yr2 = 2014
  nyrs = plot_yr2 - plot_yr1+ 1
  nmnths = nyrs*12
  nyrs_trnds = nyrs

;  seasons = (/0,3,6,9/)
  seasons = "JAS"
;  seasons = (/"ANN","DJF","MAM","JJA","SON"/)
;  seasons = (/"JFM","AMJ","JAS","OND"/)
  nseasons = dimsizes(seasons)

  ssmi_yr1 = 1979
  ssmi_yr2 = 2020
  obs_type = "CDR"

  grid = "gx1v5"
  ssmi_long_file = "/glade/p/cesm/pcwg/landrum/SSMI/SH/ifrac_sh_long_SSMI_CDR."+ssmi_yr1+"01-"+ssmi_yr2+"12.nc"                              ; this script assumes ssmi ifrac are on same grid as model
  cesm_grid_case = "b.e11.B1850C5CN.f09_g16.005"
  fgrid = "/glade/p/cesm/pcwg/landrum/"+cesm_grid_case+"/"+cesm_grid_case+".pop.h.SST.190001-199912.nc"			; needs to be ocn file so no missing TAREA. Here I am taking cesmL
	
  month_names = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
														; as grid no longer in omwg directory
  if (detrend)
     psout = get_script_prefix_name()+".dtrnd"
  else
     psout = get_script_prefix_name()
  end if 
  if (shade_STD) then
     psout = psout+".ShadeSTD"
  end if 
  psout = psout+"."+plot_yr1+"-"+plot_yr2
  psout_season = new(nseasons,typeof(psout))
;------------------------------------------------------------------------------------------------------
; CESM1_LE data
;------------------------------------------------------------------------------------------------------
  LE_20th_case = "b.e11.B20TRC5CNBDRD.f09_g16."
  LE_RCP_case = "b.e11.BRCP85C5CNBDRD.f09_g16."
  LENScase_numbers = ispan(1,40,1)
  LENScase_numbers(35:) = ispan(101,105,1)
  ncases = dimsizes(LENScase_numbers)

  scnumbers = sprinti("%0.3i",case_numbers)

  LEdatadir = "/glade/p/cesm/pcwg/landrum/LE/ice/analysis/monthly/"

  CESM1_yr1 = 1920
  CESM1_yr2 = 2005
  CESM12_yr1 = CESM1_yr2 + 1
  CESM12_yr2 = 2080

  LEcase1 = LEdatadir+"aice_sh_long_"+LE_20th_case+scnumbers+"."+CESM1_yr1+"01-"+CESM1_yr2+"12.nc"
  LEcase2 = LEdatadir+"aice_sh_long_"+LE_RCP_case+scnumbers+"."+CESM12_yr1+"01-"+CESM12_yr2+"12.nc"
; ==============================================================

  year = ispan(plot_yr1,plot_yr2,1)

;-------------------------------------------------------------------------------------------
; grid info for normalizing. obs and cesm are all on this grid
;------------------------------------------------------------------------------------------

  f = addfile(fgrid,"r")
  tarea = doubletofloat(f->TAREA(:75,:))
  tarea = tarea/(100.*100.)					; cm2 to m2
  tlat = doubletofloat(f->TLAT(:75,:))
  tlon = doubletofloat(f->TLONG(:75,:))

  n10 = max(ind(tlat(:,0).lt.-10.))
  tlon1d = dim_avg_n(tlon(:n10,:),0)
  nlon = dimsizes(tlon1d)
  if (tlon1d(nlon-1).lt.tlon1d(0)) then			; not monotonically increasing
     dlon = tlon1d
     dlon(1:) = dlon(1:) - dlon(:nlon-2)
     ilon = min(ind(dlon.lt.0.))
     if (tlon1d(ilon).gt.0.) then
        tlon1d(0:ilon-1) = tlon1d(0:ilon-1) - 360.
     else
        tlon1d(ilon:) = tlon1d(ilon:) + 360.
     end if
     delete(dlon)
  end if

  dlon_for = tlon1d
  dlon_for(:nlon-2) = tlon1d(1:) - tlon1d(:nlon-2)
  dlon_for(nlon-1) = tlon1d(0) - tlon1d(nlon-1)
  dlon_for = where(dlon_for.lt.0.,dlon_for+360.,dlon_for)

  dlon_back = tlon1d
  dlon_back(1:nlon-1) = tlon1d(1:) - tlon1d(:nlon-2)
  dlon_back(0) = tlon1d(0) - tlon1d(nlon-1)
  dlon_back = where(dlon_back.lt.0.,dlon_back+360.,dlon_back)

  dcesm = dlon_for
  dcesm = 0.5*(dlon_for+dlon_back)                    ; use for normalizing data
;-------------------------------------------------------------------------------------------
; CESM2 LE 
;------------------------------------------------------------------------------------------
print("Loading CESM2-LE data")
  CESM2_aice_long = new((/nCESM2,nyrs,nseasons,nlon/),float)
  CESM2_aice_long!0 = "case"
  CESM2_aice_long!1 = "year"
  CESM2_aice_long!2 = "season"
  CESM2_aice_long!3 = "longitude"
  CESM2_aice_long&longitude = tlon1d

  CESM2_aice_long_trnd = new((/nCESM2,nseasons,nlon/),float)
  CESM2_aice_long_trnd!0 = "case"
  CESM2_aice_long_trnd!1 = "season"
  CESM2_aice_long_trnd!2 = "longitude"
  CESM2_aice_long_trnd&longitude = tlon1d

  tmpaice = new((/nyrs*12,nlon/),float)
  i1 = (plot_yr1 - CESM2_yr1)*12
  i2 = (plot_yr2 - CESM2_yr1 + 1)*12 - 1

  dcesm2d = conform(tmpaice,dcesm,1)

  do i=0,nCESM2-1
     f = addfile (fCESM2(i), "r")
     tmpaice = f->aice_sh_long(i1:i2,:)
     tmpaice = tmpaice/dcesm2d

     do is=0,nseasons-1
        if (typeof(seasons).eq."integer") then
           CESM2_aice_long(i,:,is,:) = (/ tmpaice(seasons(is)::12,:) /)
        else
           season = seasons(is)
           if (season.ne."ANN") then
              tmp3d = new((/nmnths,2,nlon/),typeof(tmpaice))
              tmp3d(:,0,:) = tmpaice
              tmp3d(:,1,:) = tmpaice
              tmp3d!1 = ""  ; assign no name but prevents error message
              tmp_seasonal = month_to_season(tmp3d,season)
              CESM2_aice_long(i,:,is,:) = (/ tmp_seasonal(:,0,:) /)

              delete(tmp3d)
              delete(tmp_seasonal)
           else
              CESM2_aice_long(i,:,is,:) = (/ month_to_annual(tmpaice,1) /)
           end if
        end if

        aiceDtrnd = dtrend_msg_n(year,CESM2_aice_long(i,:,is,:),True,True,0)
        CESM2_aice_long_trnd(i,is,:) = aiceDtrnd@slope
     end do
     delete(aiceDtrnd)
  end do
  delete(tmpaice)

  mn_CESM2_aice_long = dim_avg_n_Wrap(CESM2_aice_long,1)
  ens_avg_CESM2_aice_long = dim_avg_n_Wrap(mn_CESM2_aice_long,0)
  ens_avg_CESM2_aice_long_trnd = dim_avg_n_Wrap(CESM2_aice_long_trnd,0)

  max_CESM2_aice_long = dim_max_n_Wrap(mn_CESM2_aice_long,0)				; largest mean value 
  min_CESM2_aice_long = dim_min_n_Wrap(mn_CESM2_aice_long,0)				; smallest mean value

  if (detrend)
     tmp = dtrend_msg_n(year,CESM2_aice_long,True,False,1)
     std_CESM2_aice_long = dim_stddev_n_Wrap(tmp,1)
     delete(tmp)
     std_CESM2_aice_long!0 = "case"
     std_CESM2_aice_long!1 = "season"
     std_CESM2_aice_long!2 = "longitude"
     std_CESM2_aice_long&longitude = tlon1d
  else
     std_CESM2_aice_long = dim_stddev_n_Wrap(CESM2_aice_long,1)
  end if
  ens_avg_std_CESM2_aice_long = dim_avg_n_Wrap(std_CESM2_aice_long,0)
  max_std_CESM2_aice_long = dim_max_n_Wrap(std_CESM2_aice_long,0)                            ; largest mean value 

  mn_plus_std = mn_CESM2_aice_long
  mn_plus_std = mn_CESM2_aice_long + std_CESM2_aice_long
  mn_minus_std = mn_CESM2_aice_long
  mn_minus_std = mn_CESM2_aice_long - std_CESM2_aice_long

  max_CESM2_aice_long = dim_max_n_Wrap(mn_plus_std,0)                            ; largest mean value
  min_CESM2_aice_long = dim_min_n_Wrap(mn_minus_std,0)                            ; smallest mean value

  delete(mn_plus_std)
  delete(mn_minus_std)
  delete(CESM2_aice_long)

;-------------------------------------------------------------------------------------------
; now SSMI observations
;------------------------------------------------------------------------------------------
print("loading SSMI observations")
  ifrac_long = new((/nyrs_trnds,nseasons,nlon/),float)
  ifrac_long!0 = "year"
  ifrac_long!1 = "season"
  ifrac_long!2 = "longitude"
  ifrac_long&longitude = tlon1d

  ifrac_long_trnd = new((/nseasons,nlon/),float)
  ifrac_long_trnd!0 = "season"
  ifrac_long_trnd!1 = "longitude"
  ifrac_long_trnd&longitude = tlon1d

  i1 = (plot_yr1 - ssmi_yr1)*12
  i2 = (plot_yr2 - ssmi_yr1 + 1)*12 - 1
  if (obs_type.eq."CDR") then						; ifrac as a function of longitude already exists
     f = addfile(ssmi_long_file,"r")
     tmp_ifrac_long_CDR = f->ifrac_sh_long(i1:i2,:)
     tmp3d = new((/nyrs_trnds*12,nlon,1/),typeof (tmp_ifrac_long_CDR))       ; need 3d for month_to_season
     tmp3d(:,:,0) = tmp_ifrac_long_CDR
     do is=0,nseasons-1
        if (typeof(seasons).eq."integer") then
           tmp_ifrac_long = tmp_ifrac_long_CDR(seasons(is)::12,:)
        else
           season = seasons(is)
           if (season.ne."ANN") then
              tmp_season = month_to_season(tmp3d,season)
              tmp_ifrac_long = tmp_season(:,:,0)
              delete(tmp_season)
           else
              tmp_ifrac_long = month_to_annual(tmp_ifrac_long_CDR,1)
           end if
        end if

        if (is.eq.0) then
            dobs2d = conform(tmp_ifrac_long,dcesm,1)
        end if
        tmp_ifrac_long = tmp_ifrac_long/dobs2d
        ifrac_long(:,is,:) = (/ tmp_ifrac_long /)                                  ; aice in m2 as fnctn of longitude

        ifracDtrnd = dtrend_msg_n(year,tmp_ifrac_long,True,True,0)
        ifrac_long_trnd(is,:) = ifracDtrnd@slope

        delete(tmp_ifrac_long)
        delete(ifracDtrnd)
     end do
     delete(tmp_ifrac_long_CDR)
     delete(tmp3d)
  else
     f = addfile(ssmi_file,"r")
     if (obs_type.eq."Bootstrap") then
        tmp_ifrac = f->Bootstrap(i1:i2,:75,:)                                 ; only read in 1st 76 tlat for SH
     else
        tmp_ifrac = f->Nasa_Team(i1:i2,:75,:)                                 ; only read in 1st 76 tlat for SH
     end if
     tmp_ifrac@_FillValue = 9.96921e+36                                           ; SSMI uses 254 for land
     tmp_ifrac = where(tmp_ifrac.ge.200, tmp_ifrac@_FillValue,tmp_ifrac)
;     tmp_ifrac = where(tmp_ifrac.eq.-999, tmp_ifrac@_FillValue,tmp_ifrac)
     tmp_ifrac = where(tmp_ifrac.lt.0, tmp_ifrac@_FillValue,tmp_ifrac)
     tmp_ifrac = 0.01*tmp_ifrac					; from % to concentration

     tarea3d = conform(tmp_ifrac(:nyrs-1,:,:),tarea,(/1,2/))

     do is=0,nseasons-1
        if (typeof(seasons).eq."integer") then
           tmp_seasonal = tmp_ifrac(seasons(is)::12,:)
        else
           season = seasons(is)
           if (season.ne."ANN") then
              tmp_seasonal = month_to_season(tmp_ifrac,season)
           else
              tmp_seasonal = month_to_annual(tmp_ifrac,1)
           end if
        end if

        tmp_ifrac_m2 = tmp_seasonal*tarea3d*1.e-10                                  ; convert ice area from % to m^2 for averaging over latitude
        tmp_ifrac_long = dim_sum_n(tmp_ifrac_m2,1)
        if (is.eq.0) then
            dobs2d = conform(tmp_ifrac_long,dcesm,1)
        end if
        tmp_ifrac_long = tmp_ifrac_long/dobs2d
        ifrac_long(:,is,:) = (/ tmp_ifrac_long /)                                  ; aice in m2 as fnctn of longitude

        ifracDtrnd = dtrend_msg_n(year,tmp_ifrac_long,True,True,0)
        ifrac_long_trnd(is,:) = ifracDtrnd@slope

        delete(tmp_ifrac_long)
        delete(tmp_seasonal)
        delete(ifracDtrnd)
        delete(tmp_ifrac_m2)
     end do
     delete(tmp_ifrac)
  end if

  mn_ifrac_long = dim_avg_n_Wrap(ifrac_long,0)
  if (detrend)
     tmp = dtrend_msg_n(year,ifrac_long,True,False,0)
     std_ifrac_long = dim_stddev_n_Wrap(tmp,0)
     delete(tmp)
     std_ifrac_long!0 = "season"
     std_ifrac_long!1 = "longitude"
     std_ifrac_long&longitude = tlon1d
  else
     std_ifrac_long = dim_stddev_n_Wrap(ifrac_long,0)
  end if

  delete(ifrac_long)
;-------------------------------------------------------------------------------------------
; CESM1 LE
;------------------------------------------------------------------------------------------

  if (plot_CESM1_LE) then
     print("Loading CESM1-LE data")
     LENS_aice_long = new((/ncases,nyrs,nseasons,nlon/),float)
     LENS_aice_long!0 = "case"
     LENS_aice_long!1 = "year"
     LENS_aice_long!2 = "season"
     LENS_aice_long!3 = "longitude"
     LENS_aice_long&longitude = tlon1d

     LENS_aice_long_trnd = new((/ncases,nseasons,nlon/),float)
     LENS_aice_long_trnd!0 = "case"
     LENS_aice_long_trnd!1 = "season"
     LENS_aice_long_trnd!2 = "longitude"
     LENS_aice_long_trnd&longitude = tlon1d

     tmpaice = new((/nyrs*12,nlon/),float)
     i1 = (plot_yr1 - case1_yr1)*12
     i2 = (plot_yr2 - case2_yr1 + 1)*12 - 1
     nt1 = (case1_yr2 - plot_yr1 + 1)*12
     nt2 = (plot_yr2 - case1_yr2)*12

     do i=0,ncases-1
        f1      = addfile (LEcase1(i), "r")
        f2      = addfile (LEcase2(i), "r")
        tmpaice1 = f1->aice_sh_long(i1:,:)
        tmpaice2 = f2->aice_sh_long(:i2,:)

        tmpaice(:nt1-1,:) = tmpaice1
        tmpaice(nt1:,:) = tmpaice2
        tmpaice = tmpaice/dcesm2d

        do is=0,nseasons-1
           if (typeof(seasons).eq."integer") then
              LENS_aice_long(i,:,is,:) = tmpaice(seasons(is)::12,:)
           else
              season = seasons(is)
              if (season.ne."ANN") then
                 tmp3d = new((/nmnths,2,nlon/),typeof(tmpaice))
                 tmp3d(:,0,:) = tmpaice
                 tmp3d(:,1,:) = tmpaice

                 tmp_seasonal = month_to_season(tmp3d,season)
                 LENS_aice_long(i,:,is,:) = (/ tmp_seasonal(:,0,:) /)

                 delete(tmp3d)
                 delete(tmp_seasonal)
              else
                 LENS_aice_long(i,:,is,:) = (/ month_to_annual(tmpaice,1) /)
              end if
           end if
           aiceDtrnd = dtrend_msg_n(year,LENS_aice_long(i,:,is,:),True,True,0)
           LENS_aice_long_trnd(i,is,:) = aiceDtrnd@slope
        end do
     end do
     delete(tmpaice1)
     delete(tmpaice2)

     delete(tmpaice)
     delete(aiceDtrnd)

     mn_LENS_aice_long = dim_avg_n_Wrap(LENS_aice_long,1)
     ens_avg_LENS_aice_long = dim_avg_n_Wrap(mn_LENS_aice_long,0)
     ens_avg_LENS_aice_long_trnd = dim_avg_n_Wrap(LENS_aice_long_trnd,0)

     if (detrend)
        tmp = dtrend_msg_n(year,LENS_aice_long,True,False,1)
        std_LENS_aice_long = dim_stddev_n_Wrap(tmp,1)
        delete(tmp)
        std_LENS_aice_long!0 = "case"
        std_LENS_aice_long!1 = "season"
        std_LENS_aice_long!2 = "longitude"
        std_LENS_aice_long&longitude = tlon1d
     else
        std_LENS_aice_long = dim_stddev_n_Wrap(LENS_aice_long,1)
     end if
     ens_avg_std_LENS_aice_long = dim_avg_n_Wrap(std_LENS_aice_long,0)

     delete(LENS_aice_long)

  end if

;-------------------------------------------------------------------------------------------
; PLOTTING
;------------------------------------------------------------------------------------------
print("Making plots")
; pivot about 0 deg east for plotting
;  decor_pivot = lonPivot(decorrel_lon,1.)

;-------pivot about 0 deg east for plotting if longitude does not start ~0
  if (min(tlon1d).lt.0.) then
     mn_CESM2_aice_long = lonPivot(mn_CESM2_aice_long,1.)
     std_CESM2_aice_long = lonPivot(std_CESM2_aice_long,1.)
     CESM2_aice_long_trnd = lonPivot(CESM2_aice_long_trnd,1.)

     ens_avg_CESM2_aice_long = lonPivot(ens_avg_CESM2_aice_long,1.)
     ens_avg_std_CESM2_aice_long = lonPivot(ens_avg_std_CESM2_aice_long,1.)
     ens_avg_CESM2_aice_long_trnd = lonPivot(ens_avg_CESM2_aice_long_trnd,1.)
     max_CESM2_aice_long = lonPivot(max_CESM2_aice_long,1)                            ; largest mean value
     min_CESM2_aice_long = lonPivot(min_CESM2_aice_long,1)                            ; smallest mean value
     max_std_CESM2_aice_long = lonPivot(max_std_CESM2_aice_long,1)                            ; largest mean value

     mn_ifrac_long = lonPivot(mn_ifrac_long,1.)
     std_ifrac_long = lonPivot(std_ifrac_long,1.)
     ifrac_long_trnd = lonPivot(ifrac_long_trnd,1.)

     if (plot_CESM1_LE) then
        mn_LENS_aice_long = lonPivot(mn_LENS_aice_long,1.)
        ens_avg_LENS_aice_long = lonPivot(ens_avg_LENS_aice_long,1.)
        std_LENS_aice_long = lonPivot(std_LENS_aice_long,1.)
        ens_avg_std_LENS_aice_long = lonPivot(ens_avg_std_LENS_aice_long,1.)
     end if
  end if
;************************************************
; resource list for contours
;************************************************
  res                      = True                 ; plot mods desired
  res@gsnDraw		   = False		
  res@gsnFrame		   = False		

  res@trGridType = "TriangularMesh"    ; allows missing coordinates

  res@gsnStringFontHeightF = 0.015

  res@tiYAxisString= "Mean Sea Ice Area (10~S~4~N~ km~S~2~N~)"
  res@tiXAxisString= "Longitude"

  res@gsnMaximize           = False         ; enlarge plot 
;  res@gsnLeftStringOrthogonalPosF = 0.03

  res@gsnCenterString = ""

   pres = True
   pres@gsnFrame = False
   pres@txFontHeightF = 0.01
   pres@gsnMaximize     = True                 ; fill up the page
   pres@gsnPanelRowSpec  = True         ; Specify plots per row

   xSSMI = (/ mn_ifrac_long&longitude /)
   xp = new(2*nlon,float)
   xp(:nlon-1) = xSSMI
   xp(nlon:) = xSSMI(::-1)
   yp = new(2*nlon,float)

   res@tmXBMode          = "Explicit"              ; explicit labels
   res@tmXBValues        = (/60,120,160,210,240,300/)  ; location of labels 
   res@tmXBLabels        = (/"60E","120E","160E","150W","120W","60W"/); labels themselves
   res@tmLabelAutoStride = True                    ; nice stride on labels
   res@trXMaxF = 360.

   res@vpHeightF= 0.4                    ; change aspect ratio of plot
   res@vpWidthF = 0.8                  

   res@xyLineThicknessF = 3.0
   model_colors = (/"gray","black","black","red"/)

   res_ssmi = True
   res_ssmi = res
   res_CESM2 = True
   res_CESM2 = res
 
   if (plot_CESM1_LE) then
      data = new((/ncases+1,nlon/),float)
      data!0 = "case"
      data!1 = "longitude"
      data&longitude = tlon1d
      LENS_colors = new(ncases+1,typeof(model_colors))
      LENS_colors(:ncases-1) = "gray"
      LENS_colors(ncases) = "black"

      opacities = new(ncases+1,float)
      opacities(:ncases-1) = 0.5
      opacities(ncases) = 1.

      res_LENS = True
      res_LENS = res
      LENS_colors = new(ncases+1,typeof(model_colors))
      LENS_colors(:ncases-1) = "gray"
      LENS_colors(ncases) = "black"
      res_LENS@xyLineColors = LENS_colors
      res_LENS@xyMonoDashPattern = True
      res_LENS@xyDashPattern = 0
      res_LENS@xyLineOpacities = opacities
   end if

   data2 = new((/nCESM2+1,nlon/),float)
   data2!0 = "case"
   data2!1 = "longitude"
   data2&longitude = tlon1d
   CESM2_colors = new(nCESM2+1,typeof(model_colors))
   CESM2_colors(:nCESM2-1) = "gray"
   CESM2_colors(nCESM2) = "black"

   opacities2 = new(nCESM2+1,float)
   opacities2(:nCESM2-1) = 0.5
   opacities2(nCESM2) = 1.

   CESM2_colors = new(nCESM2+1,typeof(model_colors))
   CESM2_colors(:nCESM2-1) = "gray"
   CESM2_colors(nCESM2) = "black"
   res_CESM2@xyLineColors = CESM2_colors
   res_CESM2@xyMonoDashPattern = True
   res_CESM2@xyDashPattern = 0
   res_CESM2@xyLineOpacities = opacities2

   res@xyLineColor = model_colors(0)
   res_ssmi@xyLineColor = model_colors(3)

  lres = True
  lres@gsLineDashPattern = 0

; legend resources
  lgres                    = True
  lgres@lgMonoLineColor    = False
  lgres@lgLineThicknessF   = 2
  lgres@lgMonoDashIndex = True
  lgres@lgDashIndex        = 0
  lgres@lgLabelFontHeightF = .02            ; legend label font thickness
  lgres@vpWidthF           = 0.2           ; width of legend (NDC)
  lgres@vpHeightF          = 0.04           ; height of legend (NDC)
  lgres@lgPerimOn = False

   resM = True
   resM = True
   resM@gsnDraw = False
   resM@gsnFrame = False
   resM@mpMinLatF = -80.
   resM@mpMaxLatF = -40.
   resM@mpCenterLonF =  180.
   resM@tmXBMode          = "Explicit"              ; explicit labels
   resM@tmXBValues        = (/60,120,160,210,240,300/)  ; location of labels
   resM@tmXBLabels        = (/"60E","120E","160E","150W","120W","60W"/); labels themselves
   res@tmLabelAutoStride = True                    ; nice stride on labels

   resM@vpWidthF = 0.8
   resM@vpHeightF = 0.15
;   resM@vpYF = 0.4
;   resM@vpXF = 0.125
   resM@tmXBLabelFontHeightF = 0.015
   resM@tmYROn = False
   resM@tmYLOn = False

   attachres1 = True
   attachres2 = True
   attachres2@gsnAttachPlotsXAxis = True ; attaches along x-axis 

   gsres                   = True                        ; poly res
   gsres@gsFillColor       = "LightGrey"                 ; color chosen
   gsres@gsFillOpacityF = 0.2

   do is=0,nseasons-1
      plot = new(3,graphic)
      plot_ssmi = new(3,graphic)			; on different grid than HRes so will have to overlay
      plot_CESM2 = new(3,graphic)                        ; on different grid than HRes so will have to overlay
      plot_LENS = new(3,graphic)
      poly = new(5,graphic)
      newplot = new(3,graphic)

      if (typeof(seasons).eq."integer") then
         season = month_names(seasons(is))
      else
         season = seasons(is)
      end if
      psout_season(is) = psout+"."+season
;      wks  = gsn_open_wks("png",psout_season(is))             ; open a ps file
      wks  = gsn_open_wks("pdf",psout_season(is))             ; open a ps file
      gsn_define_colormap(wks,"BkBlAqGrYeOrReViWh200")
      ltgry    = NhlNewColor(wks,0.7,0.7,0.7)      ; add gray to map
      dkgry    = NhlNewColor(wks,0.4,0.4,0.4)
      blk    = NhlNewColor(wks,0.,0.,0.)      ; add gray to map


      ;;; Plot the ice area
      plot_ssmi(0) = gsn_csm_xy(wks,xSSMI,mn_ifrac_long(is,:),res_ssmi)

      ;data2(:nCESM2-1,:) = mn_CESM2_aice_long(:,is,:)
      data2(nCESM2,:) = (/ ens_avg_CESM2_aice_long(is,:) /)
      plot_CESM2(0) = gsn_csm_xy(wks,xSSMI,data2,res_CESM2)

      if (shade_STD) then
         gsres@gsFillColor       = model_colors(3)
         yp(:nlon-1) = mn_ifrac_long(is,:)+std_ifrac_long(is,:)
         yp(nlon:) = mn_ifrac_long(is,::-1)-std_ifrac_long(is,::-1)
         poly(0) = gsn_add_polygon (wks,plot_ssmi(0),xp,yp,gsres)

         gsres@gsFillColor       = model_colors(2)
         yp(:nlon-1) = max_CESM2_aice_long(is,:)
         yp(nlon:) = min_CESM2_aice_long(is,::-1)
;         yp(:nlon-1) = max_CESM2_aice_long(is,:)+max_std_CESM2_aice_long(is,:)
;         yp(nlon:) = min_CESM2_aice_long(is,::-1)-max_std_CESM2_aice_long(is,::-1)
         poly(2) = gsn_add_polygon (wks,plot_CESM2(0),xp,yp,gsres)

      end if

      overlay(plot_ssmi(0),plot_CESM2(0))

      res@trYMaxF = 1.8

      res@gsnLeftString = "Sea Ice Area STD"
      res@tiYAxisString= "(10~S~10~N~ m~S~2~N~) per ~S~o~N~longitude"
      plot_ssmi(1) = gsn_csm_xy(wks,xSSMI,std_ifrac_long(is,:),res_ssmi)
      data2(:nCESM2-1,:) = std_CESM2_aice_long(:,is,:)
      data2(nCESM2,:) = (/ ens_avg_std_CESM2_aice_long(is,:) /)
      plot_CESM2(1) = gsn_csm_xy(wks,xSSMI,data2,res_CESM2)

      overlay(plot_ssmi(1),plot_CESM2(1))
      res@trYMaxF = 0.06


      ;;; Plot the Sea ice trends
      res_ssmi@tiYAxisString = "Sea Ice Area Trends (10~S~2~N~ km~S~2~N~ year~S~-1~N~)" 
      res_ssmi@trYMaxF = 3
      res_ssmi@trYMinF = -7
      res_CESM2@trYMaxF = 3
      res_CESM2@trYMinF = -7

      lres@xyLineColor = dkgry
      ytemp = ifrac_long_trnd(is,:)*100
      plot_ssmi(2) = gsn_csm_xy(wks,xSSMI,ytemp,res_ssmi)
      delete(ytemp)
      data2(:nCESM2-1,:) = CESM2_aice_long_trnd(:,is,:)
      ;data2(nCESM2,:) = (/ ens_avg_CESM2_aice_long_trnd(is,:) /) ; don't load ensemble mean
      ytemp = data2*100
      plot_CESM2(2) = gsn_csm_xy(wks,xSSMI,ytemp,res_CESM2)
 
      dum = gsn_add_polyline(wks,plot_ssmi(2),(/0.,360./),(/0.,0./),lres)
      overlay(plot_ssmi(2),plot_CESM2(2))
      if (plot_CESM1_LE) then
         data(:ncases-1,:) = mn_LENS_aice_long(:,is,:)
         data(ncases,:) = (/ ens_avg_LENS_aice_long(is,:) /)
         plot_LENS(0) = gsn_csm_xy(wks,xSSMI,data,res_LENS)

         data(:ncases-1,:) = std_LENS_aice_long(:,is,:)
         data(ncases,:) = (/ ens_avg_std_LENS_aice_long(is,:) /)
         plot_LENS(1) = gsn_csm_xy(wks,xSSMI,data,res_LENS)

         data(:ncases-1,:) = LENS_aice_long_trnd(:,is,:)
         data(ncases,:) = (/ ens_avg_LENS_aice_long_trnd(is,:) /)
         plot_LENS(2) = gsn_csm_xy(wks,xSSMI,data,res_LENS)

         overlay(plot_ssmi(0),plot_LENS(0))
         overlay(plot_ssmi(1),plot_LENS(1))
         overlay(plot_ssmi(2),plot_LENS(2))
      end if

      if (detrend)
         ;pres@txString = season + " lat. avg Ice Area, STD in detrended Ice Area and Ice Area Trends~C~ ~Z85~ "
      else
         ;pres@txString = season + " lat. avg Ice Area, STD in Ice Area and Ice Area Trends~C~ ~Z85~ "
      end if

      ;pres@txString = plot_yr1+"-"+plot_yr2+" "+season + " Ice Area and Ice Area Trends~C~ ~Z85~ "

      map1 = gsn_csm_map_ce(wks,resM)
      map2 = gsn_csm_map_ce(wks,resM)
      map3 = gsn_csm_map_ce(wks,resM)
      newplot(0) = gsn_attach_plots(plot_ssmi(0),map1,attachres1,attachres2)
      newplot(1) = gsn_attach_plots(plot_ssmi(1),map2,attachres1,attachres2)
      newplot(2) = gsn_attach_plots(plot_ssmi(2),map3,attachres1,attachres2)
      if (shade_STD) then
         gsn_panel(wks,(/plot_ssmi(0),plot_ssmi(2)/),(/1,2/),pres)
      else
         gsn_panel(wks,plot_ssmi,(/1,3/),pres)
      end if

      xpleg = 0.3
      ypleg = 0.75
;      if (shade_STD) then
;         ypleg = 0.32
;      end if
      if (plot_CESM1_LE) then
         lgres@lgLineColors       = model_colors(1:)
         lgres@lgItemCount        = 3
         lgres@lgItemOrder = (/2,1,0/)
         labels = (/CESM2_case_short,"CESM1-LE","SSMI"+" "+obs_type/)
         gsn_legend_ndc(wks,3,labels,0.4,ypleg,lgres)
      else
         ;lgres@lgLineColors       = model_colors(2:)
         ;lgres@lgItemCount        = 2
         ;lgres@lgItemOrder = (/1,0/)
         ;lgres@vpWidthF    = 0.15
         ;labels = (/CESM2_case_short,"SSMI"+" "+obs_type/)
         ;gsn_legend_ndc(wks,2,labels,xpleg,ypleg,lgres)
      end if

      frame(wks)
      delete(plot)
      delete(newplot)

  end do
end
