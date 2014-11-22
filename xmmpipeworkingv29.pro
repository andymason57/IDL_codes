pro xmmpipeworkingv29, other_args
	;compile_opt strictarr

!EXCEPT=2

;two parameters passed args[0] = path to observation, arg[1] = observationalID


; -------------------------------------------------------------

; Procedure to extract lightcurves from XMM data

; INPUTS:
;        args[0] - path to observation
;        args[1] - observationID
;        File enrng.dat with energy ranges in eV required for
;        extraction in a list
;        i.e. 200 10000
;             200 1000
;             1000 10000
;        Will extract lightcurves for all sources in the energy ranges
;        0.2-10 keV, 0.1-1.0 keV and 1.0-10 keV 
; --- Qaulity control feature - observations with obs_class=5 are removed in Postgres from list of observations (~265 obs in total)


; -------------------------------------------------------------
; SAS referenced in IDL setup
; Read in the file with the link to the directory
args = command_line_args()
obsID_path = args[0]
obs = obsID_path
obsID = args[1]

hard_coded_dir = '/mnt/4tbdata/limflcorr_test'



; -------------------------------------------------------------
; Read in the file with the energy ranges in it
enfl = hard_coded_dir + '/enrng.dat'
nen = file_lines(enfl)
ens = fltarr(nen, 2)
en = ''
  
openr, snow, enfl, /get_lun
 for i = 0, nen-1 do begin
      readf, snow, en
      bts = strsplit(en, /extract)
      ens[i,0] = bts[0]
      ens[i,1] = bts[1]
 endfor

; -------------------------------------------------------------
; Set the paths
odfpth = obsID_path +'/odf'
pth = obsID_path + '/processed/test/'
ppspth = obsID_path + '/pps/'


; tell SAS the paths
setod = 'export SAS_ODF='+odfpth+' && '
setoc = 'export SAS_ODF='+odfpth+' && export SAS_CCF='+odfpth+ $
        '/ccf.cif && '

tst = file_search(odfpth+'/*SUM.SAS')
if tst[0] ne '' then spawn, 'rm '+odfpth+'/*SUM.SAS'

; -------------------------------------------------------------
; run cifbuild
 ;removed call to cifbuild - use pre-existing ccf.cif file     
     
; run odfingest
spawn, setoc+'echo $SAS_ODF'
spawn, setoc+'odfingest outdir='+odfpth+'/'

; ---------------------------------------------------------------
; Begin lightcurve selection

; Find the eventlist
; - adapted to only deal with pn eventlists
evli = file_search(ppspth+'/*PIEVLI*.fit')

; If there is no eventlist then run epchain
if evli[0] eq '' then begin
   spawn, setoc+'cd '+pth+' && epchain odf='+odfpth
   evli = file_search(pth+'/*PIEVLI*.FIT')
endif 

; --------------------------------------------------------------
; Apply the barycentric corrections If you need to
spawn, setoc+' barycen table='+evli+':EVENTS' 

; --------------------------------------------------------------
; Make a test timeseries to extract the start/stop times
fbken = '(PI in [10000:15000])'
tst = pth+'lctst.fits'
;print, "test path: "
;print, tst
fbk = pth+'FLRBKLC.fits'

if file_search(tst) eq '' then spawn, setoc+'evselect table='+evli+' withrateset=true expression="(PATTERN <= 4)&&(FLAG == 0)&&('+fbken+')" maketimecolumn=true timebinsize=10. rateset='+tst

; --------------------------------------------------------------
; Extract the minimum timebin 
dat = mrdfits(tst, 1, hdr, columns=['TIME'])
tzero = dat[0].time
tmax = dat[n_elements(dat.time)-1].time
dat = 0

print, 'passed extract min timebin'


dat = mrdfits(tst, 0, hdr)
frmtm = sxpar(hdr, 'FRMTIME')
;!!!! Also check that it is a clean observation - check catalog - otherwise skip
if frmtm eq 73 then frmtm = 73.4
if frmtm eq 199 then frmtm = 199.1
if frmtm eq 50 then frmtm = 47.7
if frmtm eq 0 then frmtm = 2600

strdt = strtrim(string(fix(frmtm)),2)
dt = frmtm/1000.
frmtm = strtrim(string(dt),2)
dat = 0

print, 'passed frametime'

;----------------------------------
; used to find unique DETID for each detection the location of the XMM catalog fits file is hardcoded
catalog_dat = mrdfits(hard_coded_dir + '/3XMM_DR4cat_v1.0.fits', 1, hdr, $
          columns=['DETID', 'SRC_NUM', 'OBS_ID', 'SUM_FLAG']) 

print, 'passed find detID' 


; -------------------------------------------------------------------------
; Find the source and background region files
; The number of Objects in an observation can be given by the number
; of region files found
rgfls = file_search(ppspth+'*PN*.ASC')    
nsrcs = n_elements(rgfls)
print, 'number of elements: ', nsrcs 

; LOOP over the number of sources
   for k = 0, nsrcs-1 do begin
          
;          print, rgfls[k] 
          
         ; Find the source number for this object
         dtps = strpos(rgfls[k], '.')
         ext = strmid(rgfls[k], dtps-4, 4)
         ; new lucy code find DETID
         nwnum = 0L
        
;        print, ext
        
         reads, ext, nwnum, format='(Z)'
         
         print, 'source number: ', nwnum
         
         
       
         ; Find detectionID for source from source number found and consult with main XMM fits catalog file. 
         detID = catalog_dat[where(catalog_dat.src_num eq nwnum and catalog_dat.obs_id eq obsID)].detid
         
           print, '###detID: ', detID
        
        
;--------------------------------------------------------------------------------------
; get detID into standard format of 6 digits. detID can be 4 or 5 digits 
        convert_detid = strtrim(string(detID),2)
        length_convert_detid = strlen(convert_detid)
      
        case length_convert_detid of
          3: convert_detid = '000' + convert_detid
          4: convert_detid = '00' + convert_detid
          5: convert_detid = '0' + convert_detid
          else: convert_detid = convert_detid 
        endcase
        

        
;--------------------------------------------------------------------------------------
        
          ; Quality control - test detection does not have sum_flag = 4 if so from definition - 'located in an area where spuriuos detection may occur and possibly spurious' - increment k loop (~5850 detections)
         sum_flag = catalog_dat[where (catalog_dat.detid eq detID)].sum_flag
         if sum_flag eq 4 then continue

         ; Find the regions from the region file
         dat = strarr(file_lines(rgfls[k]))
         
;         print, "REGIONFILE............."
;         print, rgfls[k]
      
         openr, snow, rgfls[k], /get_lun
         readf, snow, dat
         free_lun, snow
        
         prbs1 = strpos(dat, 'circle')
         prbs2 = strpos(dat, 'green') ; Source extraction region
         prbs3 = strpos(dat, 'white') ; Background extraction region
         prbs4 = strpos(dat, 'red')   ; Regions to be excluded from the background
        
        print, 'STRING POS RED'
        print, prbs4

         ; sloc is line number where 'green' appears - i.e. source region
         sloc = where (prbs1 ne -1 and prbs2 ne -1)
         print, "SLOC"
         print, sloc
         ; bloc is line number where 'white' appears - i.e. background region
         bloc = where (prbs1 ne -1 and prbs3 ne -1)
         print, "BLOC"
         print, bloc
         ; beloc is line number where 'red' appears - i.e regions to be excluded
         beloc =  where (prbs1 ne -1 and prbs4 ne -1)
         PRINT,  "beloc"
         PRINT, beloc


         if strpos(dat[sloc], '#') eq -1 then src = dat[sloc] else $
            bsrc = strmid(dat[sloc],0,strpos(dat[sloc],'#')-1)
            ; bsrc graps string detector;circle(x,y, radius) for source region 
            print, "bsrc : ", bsrc
         if strpos(dat[bloc], '#') eq -1 then bkg = dat[bloc] else $
            bbkg = strmid(dat[bloc],0,strpos(dat[bloc],'#')-1)
            ;bbkg graps string detector;circle(x,y, radius) for background region
            print, "bbkg : ", bbkg
         bebkg = strarr(n_elements(beloc))
         print, "BEBKG no_elements : ", n_elements(beloc)
       
          ; I have put in - if no regions to be excluded then just print out else gather excluded regions
          if n_elements(beloc) eq 1 and beloc[0] eq -1 then print, "no regions to exclude  " else $
            for l = 0, n_elements(beloc)-1 do $
              if strpos(dat[beloc[l]], '#') eq -1 then ebkg = dat[beloc[l]] else $
                bebkg[l] = strmid(dat[beloc[l]],0,strpos(dat[beloc[l]],'#')-1)
          
         
            
         ;splits string removes 'detector' 
         src = strsplit(bsrc, ';', /extract)
         print, "src : ", src
         src = src[1]
         bkg = strsplit(bbkg, ';', /extract)
         print, "bkg : ", bkg
         bkg = bkg[1]
         ; instantiates string array to hold all occurances of excluded 'red' background regions 
         ebkg = strarr(n_elements(bebkg))
         
         ;initialise boolean ebkg_empty_flag as false - flag used to tell whether to include excluded regions in built up string for source/bgkd spectrum
         ebkg_empty_flag = 1
         
         ; if no excluded regions exist set flag to true and don't build string
         if n_elements(beloc) eq 1 and beloc[0] eq -1 then ebkg_empty_flag = 0 else $
         for l = 0, n_elements(bebkg)-1 do begin
              tebkg = strsplit(bebkg[l], '-', /extract)
              ebkg[l] = '!((X,Y) in '+tebkg[1]+')'
         endfor
      
        
; -------------------------------------------------------------------------
         ; Loop over the energy ranges required
         cnvens1 = strtrim(string(ens[*,0], format='(F5.2)'),2)
         cnvens2 = strtrim(string(ens[*,1], format='(F5.2)'),2)
         strens = cnvens1+' - '+cnvens2+' keV'

         imgnm = strarr(nen) & fltev = strarr(nen)
         ltcv = strarr(nen) & bkcv = strarr(nen)
         srspc = strarr(nen) & bkspc = strarr(nen)
         emid = fltarr(nen) & strlw = strarr(nen)
         strup = strarr(nen)

      ; LOOP over each energy range
      for j = 0, nen-1 do begin

          ; Image creation 
          emid[j] = ((ens[j,1]-ens[j,0])/2.) + ens[j,0]

          strlw[j] = strtrim(string(fix(ens[j,0])),2)
          strup[j] = strtrim(string(fix(ens[j,1])),2)
          
          ; Set up the filtering expressions
          expr='(PI in ['+strlw[j]+':'+strup[j]+ $
               '])&&(PATTERN <= 12)&&(FLAG == 0)' 
          neexpr = '(PATTERN <= 12)&&(FLAG == 0)'

          timflt = '(TIME in ['+strtrim(string(tzero),2)+':'$
                   +strtrim(string(tmax),2)+'])'
          srexpr = expr+'&&((X,Y) in '+src+')&&'+timflt
          bkexpr = expr+'&&((X,Y) in '+bkg+')&&'+timflt
          srneexpr = neexpr+'&&((X,Y) in '+src+')&&'+timflt
          
          ; split up bkneexpr string concat so that if ebkg_empty_flag is false omit excluded regions, if true include them 
          bkneexpr = neexpr+'&&((X,Y) in '+bkg+ ')'  
          if ebkg_empty_flag eq 1 then rest_bkneexpr = '&&' + ebkg + '&&' +timflt else $
            rest_bkneexpr = '&&' +timflt
          bkneexpr = bkneexpr + rest_bkneexpr
          print, "bkneexpr : ", bkneexpr 
          
          imgnm[j] = pth+'IMAGE'+strlw[j]+'-'+strup[j]+'.img'
          fltev[j] = pth+'FLTLI'+strlw[j]+'-'+strup[j]+'.ev'

          if file_search(imgnm[j]) eq '' then spawn, setoc+'evselect table="'+evli+':EVENTS" imageset='+imgnm[j]+' xcolumn="X" ycolumn="Y" imagebinning="binSize" ximagebinsize=80 yimagebinsize=80 squarepixels="true" expression="'+expr+'" imagedatatype="Int32" withimagedatatype="true" writedss="true" updateexposure="true" keepfilteroutput="Y" withfilteredset="Y" filteredset='+fltev[j] 
          ; -------------------------------------------------------------------------

          ; Extract Lightcurves
          ltcv[j] = pth+obsID+'_'+convert_detid+'_LC_'+strlw[j]+'-'+strup[j]+'_'+strdt+'ms_'+ext+'.fit'
          bkcv[j] = pth+obsID+'_'+convert_detid+'_BKLC_'+strlw[j]+'-'+strup[j]+'_'+strdt+'ms_'+ext+'.fit'

          ; Source Lightcurve
          if file_search(ltcv[j]) eq '' then spawn, setoc+'evselect table='+evli+' withrateset=true expression="'+srexpr+'" maketimecolumn=true timebinsize='+frmtm+' makeratecolumn=true rateset='+ltcv[j]
       
          ; Background Lightcurve
          if file_search(bkcv[j]) eq '' then spawn, setoc+'evselect table='+evli+' withrateset=true expression="'+bkexpr+'" maketimecolumn=true timebinsize='+frmtm+' makeratecolumn=true rateset='+bkcv[j]
       
          ; Flare Background Lightcurve
          fbk = file_search(ppspth+'*PN*FBKTSR*.fit')

; --------------------------------------------------------------------------

          ; EXTRACT Spectra for the backscal correction - added xxxx so easy to remove from compressed results
          srspc[j] = pth+'xxxxSRSPEC_'+strlw[j]+'-'+strup[j]+'_'+ext +  '_' + strtrim(string(detID),2) +'.fit'
          bkspc[j] = pth+'xxxxBKSPEC_'+strlw[j]+'-'+strup[j]+'_'+ext +  '_' + strtrim(string(detID),2) +'.fit'
                    
          if file_search(srspc[j]) eq '' then spawn, setoc+'evselect table='+evli+' withspecranges=true specchannelmin=0 specchannelmax=20479 spectralbinsize=5 energycolumn=PI expression="'+srneexpr+'" spectrumset='+srspc[j]
          spawn, setoc+' arfgen spectrumset='+ srspc[j]+' setbackscale=yes'
     
          if file_search(bkspc[j]) eq '' then spawn, setoc+'evselect table='+evli+' withspecranges=true specchannelmin=0 specchannelmax=20479 spectralbinsize=5 energycolumn=PI expression="'+bkneexpr+'" spectrumset='+bkspc[j] 
          spawn, setoc+' arfgen spectrumset='+bkspc[j]+' setbackscale=yes' + ' badpixlocation=' + evli
        
     
;------------------------------------------------------------------------------
          ; Extract Backscal values
          dat = mrdfits(srspc[j], 1, hdr)
          dat = mrdfits(bkspc[j], 1, bhdr)
        
          sclsr = sxpar(hdr, 'BACKSCAL')
          sclbk = sxpar(bhdr, 'BACKSCAL')
        
          scfc = float(sclsr)/float(sclbk)
          print, "source backscal is : ", sclsr
          print, "background backscal is : ", sclbk
          print, "THE BACKSCAL is : ", scfc
	
; --------------------------------------------------------------------------

          ; Process each of the XMM files 
          dat = mrdfits(ltcv[j],1,hdr)
          time = dat.time
          s_rt = dat.rate       ;/dt
          s_err = dat.error
          
          dat = mrdfits(bkcv[j],1,bhdr)
          b_tm = dat.time
          b_rt = dat.rate       ;/dt
          b_err = dat.error

; --------------------------------------------------------------------------
          ; Correct for the scaling of the background.
          b_rt = b_rt*scfc
          b_err = b_err*scfc

          rate = s_rt-b_rt
          err = sqrt((s_err^2.)+(b_err^2.))
; --------------------------------------------------------------------------

; --------------------------------------------------------------------------
          ; Correct for small gaps in the source lightcurve ---- NOT DONE - PROBLEMS
          corrt = smlgp(s_rt, time, dt, psdt, altm=altm, cgti=cgti)
          time = altm
          gti = cgti
; ------------------------------------------------------------------------------          
          
; -------------------------------------------------------------------
          ; Run the flare correction routine 
          gti = -1
          ; Loop over the flare background light curves
          for l = 0, n_elements(fbk)-1 do begin
             dat = mrdfits(fbk[l],1,fhdr)
             f_tm = dat.time
             f_rt = dat.rate
             
             print, "dat.rate : ", f_rt

	           ;FLCUTTR is the optimised flare cut threshold 
             ctlim = sxpar(fhdr, 'FLCUTTHR')
             
             print, "time res: ", dt
             
             ; ****** flare gap size set here **************
             gpsz = 300.
             igti = limflrcorr(f_rt, f_tm, 10., segs=segs, gpsz=gpsz, $
                               flrlim=ctlim)
             if gti[0,0] eq -1 then gti = igti else gti = [[gti],[igti]]
          endfor
          
;--------------------------------------------------------------------------

          ; GTI correct the lightcurve - creates GTI's where flare rate less than set value
          vm = make_array(n_elements(time), /integer, value=-1)
          
          for i = 0, n_elements(gti[*,0])-1 do begin
             mask = where(time ge gti[i,0] and time le gti[i,1])
             vm[mask] = 0
          endfor

          mask = where (vm eq 0)
          ; **** added when don't correct for small gaps 
          ;rate = back_sub_source_rate
          gti_rate = rate[mask]
          gti_time = time[mask]
          ;added by me to get errors for GTI LCs 
          gti_err = err[mask]

;--------------------------------------------------------------------------------------------------------
;         FORMAT of FINAL PRODUCT : obsid(10digits)_detid(6disgits)+_'corrected_source'_+lowerenergy-upperenergy+'kev' 
;
          ; write out FINAL PRODUCTS - the source background subtracted and GTI corrected lightcurve 
          num_data_rows = size(rate)
          num_data_pts = num_data_rows[1]

          ; create binary table with source background subtracted rate, time and err data
          fxhmake,header,/initialize,/extend,/date

          ; read header fields from uncorrected source LC extension [0] 
          theader = HEADFITS(ltcv[j],exten=0)
          
          obs_id = fxpar(theader, 'OBS_ID')
          exp_id = fxpar(theader, 'EXP_ID')
          object = fxpar(theader, 'OBJECT')
          RA = fxpar(theader, 'RA_OBJ')
          DEC = fxpar(theader, 'DEC_OBJ')
         
          ; add keywords so FTOOLS can read 
          fxaddpar, header, 'TIMESYS','TDB  '
          fxaddpar, header, 'TIMEUNIT','s '
          fxaddpar, header, 'obs_id', obs_id
          fxaddpar, header, 'exp_id', exp_id
          fxaddpar, header, 'object', object
          fxaddpar, header, 'RA', RA
          fxaddpar, header, 'DEC', DEC
          
          ; write out final corrected source extention[0] source fits header  
          fxwrite, hard_coded_dir + '/' + obsID + '/processed/test/' + obsID + '_' + convert_detid + '_corrected_source' + '_' + strlw[j]+'-'+strup[j]+ 'keV' +  '.fit',header
          
; ----------------------------         
          ;Create fits extension [1] to hold source corrected data
          fxbhmake,header,num_data_pts,'RATE','source binary table extension'


          ;read start and stop times for LC from fits info. 
          theader = HEADFITS(ltcv[j],exten=1)
          split_tstart = STRSPLIT(theader[35], '= ', /EXTRACT)
          split_tstop = STRSPLIT(theader[36], '= ', /EXTRACT)
          tstart = split_tstart[1]
          tstop = split_tstop[1]
 
         
 
          ; define mjd here - standard time in all XMM fits files?  
          mjd = 50814.00000000
          
          ; add keywords so FTOOLS can read LC fits file        
          fxaddpar, header, 'TIMEDEL', frmtm
          fxaddpar, header, 'TSTART', tstart
          fxaddpar, header, 'TSTOP', tstop
          fxaddpar, header, 'TUNIT1', 'count/s '
          fxaddpar, header, 'TUNIT2', 'count/s '
          fxaddpar, header, 'TUNIT3', 's    '
          fxaddpar, header, 'HDUCLASS', 'OGIP'
          fxaddpar, header, 'HDUCLAS1', 'GTI'
          fxaddpar, header, 'HDUCLAS2', 'STANDARD'
          fxaddpar, header, 'TIMEUNIT', 's'
          fxaddpar, header, 'TIMESYS', 'TT'
          fxaddpar, header, 'MJDREF', mjd
          fxaddpar, header, 'TIMEREF', 'LOCAL  '
          fxaddpar, header, 'TASSIGN', 'SATELLITE'
          ; add new keywords so can read BACKSCAL values
          fxaddpar, header, 'bkgBACKSCAL', sclbk
          fxaddpar, header, 'srcBACKSCAL', sclsr
          fxaddpar, header, 'calcBACKSCAL', scfc
          
          ; create columns for fits extension [1] data          
          fxbaddcol,1,header,rate[0],'rate'             ;Use first element in each array
          fxbaddcol,2,header,err[0],'err'
          fxbaddcol,3,header,time[0],'time'             ;to determine column properties
          
      
          fxbcreate,unit,hard_coded_dir + '/' + obsID + '/processed/test/' +  obsID + '_' + convert_detid + '_corrected_source' + '_' + strlw[j]+'-'+strup[j]+ 'keV' +  '.fit',header
          
          fxbwritm,unit,['rate','err','time'], rate, err, time
          fxbfinish,unit                 ;Close the file

;--------------------------------------------------------------------------------------------
          
;--------------------------------------------------------------------------------------------
          ; Save GTI to text file
          openw, snow, pth+'/GTI'+ext+'.txt', /get_lun
          trgti = transpose(gti)
          for counter_trg=0, n_elements(trgti[0,*])-1 do begin
             printf, snow, trgti[*,counter_trg], format = '(f0.3, 1X, f0.3)'
          endfor
          free_lun, snow

; write out GTI to fits 
          num_dimensions_gti = size((trgti), /N_DIMENSIONS)
         
         ; find start and stop times of GTIs - either 1 row or multirow
         ;one row
          num_row_index = 0 
          if num_dimensions_gti eq 1 then begin
              tstart = double(trgti[0,0])
              tstop = double(trgti[1,0])
              totaltime = tstop - tstart
              num_data_pts = 1 
              start = trgti[0,0]
              stop = trgti[1,0]
          endif
              
          ; more than one row GTI
          if num_dimensions_gti gt 1 then begin
            totaltime = 0
            num_data_pts = n_elements(trgti[0,*])
            start = make_array(1,num_data_pts, /double)
            stop = make_array(1,num_data_pts, /double)
            for counter_trg=0, n_elements(trgti[0,*])-1 do begin
              tstart = double(trgti[0,counter_trg])
              start[0,counter_trg] = tstart
              tstop = double(trgti[1, counter_trg])
              stop[0,counter_trg] = tstop
              totaltime = totaltime + (tstop - tstart)     
            endfor  
        
          endif
                 
          ; write out GTI       
          fxhmake,header,/initialize,/extend
        
          fxwrite,hard_coded_dir + '/' + obsID + '/processed/test/' + 'GTI_FITS.fit',header
          fxbhmake,header,num_data_pts,'STDGTI','GTI binary table extension'
     
          timezero = 0
          mjd = 50814.00000000
          
          fxaddpar, header, 'TUNIT1', 's '
          fxaddpar, header, 'TUNIT2', 's '
          fxaddpar, header, 'HDUCLASS', 'OGIP'
          fxaddpar, header, 'HDUCLAS1', 'GTI'
          fxaddpar, header, 'HDUCLAS2', 'STANDARD'
          fxaddpar, header, 'ONTIME', totaltime
          fxaddpar, header, 'TSTART', tstart
          fxaddpar, header, 'TSTOP', tstop
          fxaddpar, header, 'TIMEUNIT', 's'
          fxaddpar, header, 'TIMESYS', 'TT'
          fxaddpar, header, 'MJDREF', mjd
          fxaddpar, header, 'TIMEZERO', '0 '
          fxaddpar, header, 'EXTNAME', 'STDGTI '
          fxaddpar, header, 'TIMEREF', 'LOCAL  '
          fxaddpar, header, 'TASSIGN', 'SATELLITE'
          fxaddpar, header, 'TIMEZERO', timezero
          
          fxbaddcol,1,header,start[0],'start'             ;Use first element in each array 
          fxbaddcol,2,header,stop[0],'stop' 
          
                                                         
          fxbcreate,unit,hard_coded_dir + '/' + obsID + '/processed/test/' + 'GTI_FITS.fit',header
          fxbwritm,unit,['start','stop'], start, stop
          fxbfinish,unit                 ;Close the file

; ---------------------------------------------
      endfor
   endfor
   
end
  
