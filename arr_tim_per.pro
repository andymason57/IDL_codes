function uni_cdf,dd
return,dd
end 
;
;
pro arr_tim_per,times,pmin,pmax,nper,per,nbin,kdist,kprob,chi2,chipr
;
; NOTE!!! chi2 and chipr not used at the moment!
;
;  Times: Arrival times
;  pmin : Minimum period to be searched
;  pmax : Maximum ....
;  nper : Number of periods 
;  per  : the searched periods (on exit)
;  kdist: The K-S max. distances (on exit)
;  kprob: The K-S probabilities (on exit)
;
t=times-min(times)
nd=n_elements(t)
per=dblarr(nper) & kdist=per & kprob=per & chi2=per & chipr=per
chicell=lonarr(nbin)
fmax=1.0/pmin
fmin=1.0/pmax
df=(fmax-fmin)/double(nper)
;
;csamp=randomu(ss,10000)
;
dpha=1.0/float(nbin)
defn=n_elements(times)/float(nbin)
for i=0l,nper-1 do begin
  chicell(*)=0l
  p=1.0/(fmin+i*df)
  per(i)=p
  phase=t/p
  phase=phase-long(phase)
  ksone,phase,'uni_cdf',d,prob
  kdist(i)=d
  kprob(i)=prob
;  for j=0,nbin-1 do begin
;    ib=where((phase ge j*dpha)and(phase lt (j+1)*dpha))
;    if ib(0) ne -1 then chicell(j)=n_elements(ib)
;  endfor
;  chi2(i)=total((chicell-defn)^2/defn)
;  chipr(i)=lowgamma(nbin/2,chi2(i)/2.0)/gamma(nbin/2)
;  
;  
;  print,i,i/float(nper)*100.0,' %'
endfor
return
end
;
  
