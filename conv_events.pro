;
; Individual files must be stored in 'filepath' (below)
; foldbin.pro must be in the same directory as this programme
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
;         /mnt/4tbdata/xmm_reduce_results/old_results_ntt_FASPER/
filepath='/mnt/4tbdata/xmm_reduce_results/old_results_ntt_FASPER/'
!p.multi=[0,1,2,0,0]
;
; Start by creating a list of files and count them
;
spawn,'ls '+filepath+'*freq.dat > obslist'
spawn,'wc -l obslist > obsn'
nin=''
openr,1,'obsn'
readf,1,nin
close,1
nfiles=long(strmid(nin,0,8))
print,nfiles,' detections'
;
; Loop over observations
;
totcount=0l
obsdir=''
openr,1,'obslist'
openw,4,'candidates'
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
while not eof(1) do begin
    readf,1,fil 
;   
;   N O T E !!!
;   if you change the filepath, you need to change the 'filbase'
;   length below accordingly!!! (not set to first 58 chars of 'fil')
;
    filbase=strmid(fil,0,72)+'_1000_'
    print, filbase
    powerfil=filbase+'power.dat'
    print, powerfil
    freqfil=filbase+'freq.dat'
    resultfil=filbase+'restresults.dat'
    tbinfil=filbase+'tbin.dat'
    print, tbinfil
    fbinfil=filbase+'fbin.dat'
    print,'file: ',filbase,' ',float(inum)/float(nfiles)*100.0,' % done'
    freq=dblarr(npmax)
    power=freq
    tbin=dblarr(npmax*2)
    fbin=tbin
;
    openr,2,tbinfil 
    openr,3,fbinfil 
    ic=0l
    while not eof(2) do begin
      readf,2,temp
      tbin(ic)=temp
      readf,3,temp
      fbin(ic)=temp
;      print,ic,temp
      ic=ic+1
    endwhile
    tbin=tbin(0:ic)
    fbin=fbin(0:ic)
    close,2
    close,3
;
    ii=where(fbin gt 0.0)
    nii=n_elements(ii)
    if(nii gt 200) then begin
      tbin=tbin(ii)
      save,tbin,file=filbase+'tphot.sav'
    endif
    print,inum,nii
;   
  inum=inum+1
;
endwhile
close,/all
end  

