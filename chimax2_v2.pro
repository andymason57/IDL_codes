pro chimax2_v2, other_args

; e.g. of use : chimax2,t,f,err,1000,20,100.0,300.0,per,chi2,fv
;
;   Computes Chi^2 vs. period for the null hypothesis that
;   there is no variability over the trial period. So Max(chi2)
;   givess the most likely period. The fvalue is an analysis of
;   variance version of the same method.
;
;       I/O variables:
;
;       time (inout) -- time vector
;       flux (input) -- flux vector
;       err (input) --  error for flux vector
;       nper (input) -- number of periods to search for
;       nbin (input) -- number of phase bins  (try 10-20)
;       minper (input) -- minimum period to search for
;       maxper (inut) -- maximum period ----,,-----
;       per (output) --  period vector 
;       chisqu (output) --  chi^2 vector, maximum gives best period
;       fvalue (output) -- F-stat vector, maximum gives best period


args = command_line_args()
FILE_PATH = args[0] 
GTI_FILE_PATH = args[1] 
OUTPUT_FILE_PATH = args[2]
obsID = args[3]
nper = args[4]
nbin = args[5]
minper = args[6]
maxper = args[7]
detid = args[8]

; added by andy m. assign x and y and y error
time_dat = mrdfits(FILE_PATH, 1, hdr, columns=['TIME'])
flux_dat = mrdfits(FILE_PATH, 1, hdr, columns=['RATE'])
error_dat = mrdfits(FILE_PATH, 1, hdr, columns=['ERR'])
frmtm = sxpar(hdr, 'TIMEDEL')

print, "frametime passed ", frmtm

temp_X = time_dat.time
temp_Y = flux_dat.rate
temp_Z = error_dat.err

;sort out GTIs 
num_lines = file_lines(GTI_FILE_PATH)
en = ''
openr, snow, GTI_FILE_PATH, /get_lun
  for i = 0, num_lines-1 do begin
      readf, snow, en
      bts = strsplit(en, /extract)
      temp = where(temp_X gt bts[0] and temp_X lt bts[1])
      if i eq 0 then (cumulative_X = temp) else (cumulative_X = [cumulative_X,temp])
      if i eq 0 then (cumulative_Y = temp) else (cumulative_Y = [cumulative_Y,temp])
      if i eq 0 then (cumulative_Z = temp) else (cumulative_Z = [cumulative_Z,temp])
  endfor

gti_time = temp_x(cumulative_x)
gti_flux = temp_y(cumulative_y)
gti_err = temp_z(cumulative_z)

; assign time, flux, error variables
time = gti_time
flux = gti_flux
err = gti_err


nlimit=fix(1)
p=dblarr(nper)
chisqu=dblarr(nper)
pd=dblarr(nper)
aov=dblarr(nper)
fvalue=dblarr(nper)
me=dblarr(nper)
bin=dblarr(nbin)
pfrac=dblarr(nbin)
ebin=dblarr(nbin)
nf=dblarr(nbin)
vbin=dblarr(nbin)
phbin=dblarr(nbin)
minf=1.0/(maxper)
maxf=1.0/(minper)
df=(maxf-minf)/(1.*nper)
old=0.0
tzero=systime(1)
i=long(0)
while i lt nper-1 do begin
  bflag=intarr(nbin)
  bflag(*)=1
  phbin(*)=-1.0
  bin(*)=-1.0
  ebin(*)=-1.0  
  p(i)=1.0/(minf+i*df)
  phase=(time-min(time))/p(i)
  phase=phase-fix(phase)
  for j=0,nbin-1 do begin
    ii=where((phase ge j/(1.0*nbin))and(phase lt (j+1)/(1.*nbin)))
    if((n_elements(ii) ge nlimit)and(ii(0) ne -1)) then begin
      fset=flux(ii)
      eset=err(ii)
      nf(j)=n_elements(ii) 
      phbin(j)=(j+0.5)/(1.*nbin)
      bin(j)=mean(fset)
      ebin(j)=sqrt(total(eset^2))/float(nf(j))
      bflag(j)=0 
      vbin(j)=total((fset-bin(j))^2)
    endif else begin
      bflag(j)=1       
    endelse
  endfor
  jj=where(bflag eq 0)
  nbf=n_elements(jj)
  mean=total(bin(jj)*1.0/ebin(jj)^2)/total(1.0/ebin(jj)^2)
  chisqu(i)=total(((bin(jj)-mean)/ebin(jj))^2)/float(nbf-1)
  pd(i)=mean(ebin(jj))/stdev(bin(jj))
  ssb=total(nf(jj)*(bin(jj)-mean)^2)
  ssw=total(vbin(jj))
  fvalue(i)=(ssb/(nbf-1.0))/(ssw/(total(nf(jj))-nbf)*1.0)
  if(1.*i/(1.*nper)*100.-old ge 1.) then begin 
  tnow=systime(1)
  tpass=(tnow-tzero)/(1.*i)
  tleft=(nper-i)*tpass
  old=1.*i/(1.*nper)*100.
  endif
;
  i=i+1
endwhile


;calculate max_chisqu
max_chi = max(chisqu)
max_fv = max(fvalue)

;print to file results
resultfile = OUTPUT_FILE_PATH + obsID  + "_" + detid + '_epoch_results.dat'
cd, OUTPUT_FILE_PATH
openw,1, resultfile
printf, 1, max_chi
printf, 1, max_fv

;code to write array results to file
fname_per = obsID + "_" + detid + '_per.dat'
fname_chisqu = obsID + "_" + detid + '_chisq.dat'
openw, p_data, fname_per, /get_lun
openw, chi_data, fname_chisqu, /get_lun
printf, p_data, p, format = '(f7.3, 1X, f7.3, 1X, f7.3)' 
printf, chi_data, chisqu, format = '(f7.3, 1X, f7.3, 1X, f7.3)' 

print, "reached end of writes" 

free_lun, p_data
free_lun, chi_data



end
