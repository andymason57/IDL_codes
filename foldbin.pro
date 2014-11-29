pro foldbin,t,f,e,per,tmin,nbin,phbin,bin,berr
;
;!p.multi=[0,1,1,0,0]
bin=fltarr(nbin)
bflag=intarr(nbin)
bflag(*)=0
np=intarr(nbin)
phbin=fltarr(nbin)
berr=fltarr(nbin)
dbin=1.0/(1.0*nbin)
phase=(double(t)-double(tmin))/per & phase=phase-long(phase)
for i=0,nbin-1 do begin
  phbin(i)=(i*dbin+(i+1)*dbin)/2.0
  ii=where((phase ge i*dbin)and(phase lt (i+1)*dbin))
  np(i)=n_elements(ii)
  if((n_elements(ii) eq 1)and(ii(0) eq -1)) then begin
    bin(i)=0.0
    berr(i)=0.0
    bflag(i)=1
  endif else begin
    arr=f(ii)
    bin(i)=mean(arr)
    if (n_elements(ii) gt 1) then begin
      berr(i)=stdev(arr)/sqrt(n_elements(ii))
    endif
    if (n_elements(ii) eq 1) then begin 
      berr(i)=e(ii)
    endif
  endelse
endfor
ii=where(bflag eq 0)
if(ii(0) ne -1) then begin
 np=np(ii)
 phbin=phbin(ii)
 bin=bin(ii)
 berr=berr(ii)
endif
;window,0,title='FOLDED LIGHTCURVE'
;!mtitle='Period: '+strtrim(per,2)
;ploterr,phbin,bin,berr,psym=10
;for i=0,n_elements(ii)-1 do begin
; xyouts,phbin(i),bin(i)+0.03,strtrim(np(i),2)
;endfor
return
end
