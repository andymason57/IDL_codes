pro plot_FASPER_results, other_args

args = command_line_args()
filepath = args[0] 




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
;           /mnt/4tbdata/one_results_FASPER/
;filepath='/mnt/4tbdata/six_results_FASPER/'
!p.multi=[0,1,2,0,0]
;
; set up printing
set_plot,'ps' 
device,file='plot.ps'


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
backscal_count = 0l
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

; code so no need to count filename length 
    energy_pos = strpos(fil, '_1000_')
    filbase = strmid(fil,0,energy_pos) + '_1000_'
    
    ;find detID 
;    det_pos = strpos(fil,'_')
    print, "energy_pos: ", energy_pos
    det_string = strmid(filbase, energy_pos-6, 6)
    print, "det string: ", det_string
    
    
;    print, "new file base: ", new_file_base
;    filbase=strmid(fil,0,49)+'_1000_'
    print, "filbase", filbase
    powerfil=filbase+'power.dat'
    freqfil=filbase+'freq.dat'
    resultfil=filbase+'restresults.dat'
    tbinfil=filbase+'tbin.dat'
    fbinfil=filbase+'fbin.dat'
    print,'file: ',filbase,' ',float(inum)/float(nfiles)*100.0,' % done'
    freq=dblarr(npmax)
    power=freq
    tbin=dblarr(npmax*2)
    fbin=tbin
;
    openr,2,resultfil
    readf,2,nout
    jmax=0.0
    fap=1.0
;    if(nout gt 11.0) then begin
      readf,2,jmax
      readf,2,fap
;    endif
    close,2
;
;   Limits the rest to cases where log(fap) < -3.0 - changed 2.0
;
    if(alog10(fap) lt -2) then begin
;
    openr,2,freqfil
    ic=0l
    temp=0d
    while not eof(2) do begin
      readf,2,temp
      freq(ic)=temp
      ic=ic+1
    endwhile
    freq=freq(0:ic)
    close,2
;
    openr,2,powerfil
    ic=0l
    while not eof(2) do begin
      readf,2,temp
      power(ic)=temp
      ic=ic+1
    endwhile
    close,2
    power=power(0:ic)
; 
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
    logfap_arr(inum)=alog10(fap)
    nout_arr(inum)=nout
    jmax_arr(inum)=jmax
    cratearr(inum)=mean(tbin)
    totcount=totcount+1
; 
    print,totcount,nout,jmax,fap,alog10(fap),median(power),min(fbin)
    if(min(fbin) lt -100.) then begin
    print,'BACKSCALE wrong!, min(fbin) = ',min(fbin)
    backscal_count = backscal_count + 1
   endif 
;   Select detections with log(fap) < -3.0 and correct backscale
;   (min(fbin) > -100.0) for plotting
;  
    if(logfap_arr(inum) lt -2.0 and finite(logfap_arr(inum)) and (min(fbin) gt -100.)) then begin
      printf,4,strmid(fil,0,energy_pos),1.0/freq(jmax),power(jmax),logfap_arr(inum)
      save,freq,power,jmax,fap,tbin,fbin,file=fil+'_LS.sav'
      !mtitle= 'detid: # ' + det_string + '  log power: ' + string(logfap_arr(inum))
      !ytitle='Power'
      !xtitle='Period (sec)'
      plot_oi,1.0/freq,power
      foldbin,tbin,fbin,sqrt(fbin),1.0/freq(jmax),min(tbin),20,phbin,bin,berr
      !mtitle='Folded over '+strtrim(1.0/freq(jmax),2)+' sec'
      !xtitle='Phase'
      !ytitle='Flux'
      ploterr,[phbin,phbin+1.0],[bin,bin],[berr,berr]
;      wait,1
    endif
;  
  endif
  inum=inum+1
;
endwhile
print, "total BACKSCAL count : ", backscal_count 

; close printing 
device,/close 

close,/all
end  

