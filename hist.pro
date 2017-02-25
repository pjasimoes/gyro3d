FUNCTION hist,y,nbins,rr,binsize,loc

;+
; hist.pro
; purpose: histograms input y into nbins
; call: histogram=hist(y,nbins,indices,binsize,loc)
; input: y(data), nbins
; output: h(histogram), indices, binsize, loc(start locations of bins)
; by PSimoes, 2013, Glasgow
;-

 miny=min(y)
 maxy=max(y)
 binsize = (maxy - miny) / nbins 
 h=lonarr(nbins)
 loc=fltarr(nbins)

 ind=lonarr(nbins+1)
 pos=0l
 FOR i=0,nbins-2 DO BEGIN 
    loc[i]=binsize*i+miny
    ;print,binsize*i+miny,binsize*(i+1)+miny
    pp=where(y GE binsize*i+miny AND y LT binsize*(i+1)+miny,count)
    h[i]=count
    pos=[pos,pp]
    ind[i+1]=ind[i]+count
 ENDFOR 
 i=nbins-1
 loc[i]=binsize*i+miny
 ;;print,binsize*i+miny,binsize*(i+1)+miny
 pp=where(y GE binsize*i+miny AND y LE binsize*(i+1)+miny,count)
 h[i]=count
 pos=[pos,pp]
 ind[i+1]=ind[i]+count

 pos=pos[1:*]

 rr=[ind+n_elements(ind),pos]

return,h

END 
