function testread_am2pp,inpfile,verbose


; d duncan, csu, 1/25/16 for first version of L1R pp for AMSR2
; updated DD 6/9/16 to include TAI93 time and sunglint calculated internally
; updated DD 12/8/16 for very different pp output file format

if n_elements(verbose) EQ 0 then verbose=0 else verbose=1

; Define Header Record ( 536 Bytes Total)
date = create_struct('Year',     0, $
                     'Month',    0, $
                     'Day',      0, $
                     'Hour',     0, $
                     'Minute',   0, $
                     'Second',   0)

hdr=create_struct('Satellite',            bytarr(12), $ ; 12
                  'Sensor',               bytarr(12), $ ; 12  24
                  'PreProcessorVersion',  bytarr(12), $ ; 12  36
                  'RadFile',              bytarr(128), $ ;128 164
                  'PrfDbFile',            bytarr(128), $ ;128 292
                  'CalibrationFile',      bytarr(128), $ ;128 420
                  'Granule',                      0L, $ ;  4 424
		  'NScans',			  0L, $ ;  4 428
		  'NPixels',			  0L, $ ;  4 432
                  'NChannels',                    0L, $ ;  4 436
                  'ChanFreqs',            fltarr(15), $ ; 60 496
                  'Comment',              bytarr(40))   ; 40 536

; Read in File Header

openr,lun,inpfile,/compress,/get_lun
readu,lun,hdr
close,lun
free_lun,lun
npix = hdr.npixels
if verbose then help,hdr,/structure

;stop

; Define Scan Header Record (24 Bytes)

scanhdr=create_struct('ScanDate',               date, $ ; 12
                      'ScLat',                   0.0, $ ; 4
                      'ScLon',                   0.0, $ ; 4
                      'ScAlt',                   0.0,$   ; 4 
                      'TAI93',                   0.d);, $ ; 8
if verbose then help,scanhdr,/structure

; Define Data Record (160 Bytes)

drec=create_struct('Lat',                        0.0, $ ;  4  
                   'Lon',                        0.0, $ ;  4   8
		   'Tb',		  fltarr(15), $ ; 60  68
		   'EIA',		  fltarr(15), $ ; 60 128 
                   'Twb',                        0.0, $ ;  4 132
                   'TCWV',                       0.0, $ ;  4
                   'Tskin',                      0.0, $ ;  4
                   'T2m',                        0.0, $ ;  4
                   'Azimuth',                    0.0, $ ;  4 152
                   'WindMag',                    0.0, $ ;  4 156
                   'SLP',                        0.0, $ ;  4 160
                   'TProf',               intarr(16), $ ; 32 192
                   'WVProf',              intarr(16), $ ; 32 224
                   'QualityFlag',                  0, $ ;  2 226
                   'WindDir',                      0, $ ;  2 228
                   'LOpct',                bytarr(4), $ ;  4 
		   'SunGlintAngle',               0B, $ ;  1
		   'SfcTypeIndex',                0B, $ ;  1
		   'SnowCoverIndex',              0B, $ ;  1
		   'OroliftIndex',                0B)   ;  1 
if verbose then help,drec,/structure

  openr,lun,inpfile,/compress,/get_lun
    scan = create_struct('hdr', scanhdr, $
                         'scan', replicate(drec,npix))

    data = create_struct('hdr',       hdr, $
                         'data', replicate(scan,hdr.nscans))

  readu,lun,data
  close,lun
  free_lun,lun

return,data
end
