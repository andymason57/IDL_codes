;
; Individual files must be stored in 'filepath' (below)
; After that just run this code. If you want to make a ps-plot
; of the results, just define 
;
; IDL>set_plot,'ps'  and
; IDL>device,file='plot.ps'
;
; run the code and do 
;
; IDL>device,/close 
;
;
filepath='/mnt/4tbdata/xmm_reduce_results/old_results_ntt_FASPER/'
!p.multi=[0,1,3,0,0]
;
; Start by creating a list of files and count them
;
spawn,'ls '+filepath+'*tphot.sav > tobslist'
spawn,'wc -l tobslist > tobsn'
nin=''
openr,1,'tobsn'
readf,1,nin
close,1
nfiles=long(strmid(nin,0,8))
print,nfiles,' detections'
;
; Loop over observations
;
totcount=0l
obsdir=''
openr,1,'tobslist'
openw,4,'tcandidates'
;
fil=''
npmax=5000000l
logfap_arr=dblarr(npmax)
nout_arr=logfap_arr
jmax_arr=lonarr(npmax)
cratearr=lonarr(npmax)
temp=0.0d
inum=0l
;
pmin=1.0
pmax=10000.0
nper=1000
nbin=10
while not eof(1) do begin
    readf,1,fil 
;
;   Restore the arrival time data amd analyze it
;
    restore,fil
    arr_tim_per,tbin,pmin,pmax,nper,per,nbin,kdist,kprob,chi2,chipr 
;
;   Plot the periodogram, the probability and the phase histogram of
;   arrival times if FAP lt 0.001 
;
    if(min(kprob) lt 0.00000001) then begin 
      ii=where(kprob eq min(kprob))
      per0=per(ii(0))
      per0s=per0(0)
      Phase=(tbin-min(tbin))/per0s & phase=phase-long(phase) 
      plot_oi,per,kdist
      plot_oo,per,kprob
      plot,histogram(phase,min=0.0,max=1.0,binsize=0.1),psym=10    
      wait,0.1
    endif
    print,fil
;   
    inum=inum+1
;
endwhile
close,/all
end  

