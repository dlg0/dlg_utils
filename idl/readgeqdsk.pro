function readGeqdsk, fileName, $
	plot_ = plot_, $
	pAngle = pAngle, $
	half = half, $
	reWrite = reWrite, $
	bFactor = bFactor, $
	bPolFactor = bPolFactor, $
	cMod_lim_adjust = cMod_lim_adjust, $
	nstx_lim_adjust = nstx_lim_adjust, $
	iter_lim_adjust = iter_lim_adjust, $
	bInterpS = bInterpS, $
	fieldLineIn = fieldLineIn, $
	fieldLineOut = fieldLineOut, $
	use_dlg_bField = use_dlg_bField, $
	noToroidalFlux = noToroidalFlux, $
	create_aorsa_nc_input = create_aorsa_nc_input

    @constants

;   Read in data from g-eqdsk file

openr, lun, fileName, /get_lun

case_    = strArr ( 6 )

f1  = '(6a8,3i4)'
f2  = '(5e16.9)'
f3  = '(2i5)'
f4	= '(i5,e16.9,i5)'

readf, lun, format = f1, case_, idum, nW, nH
readf, lun, format = f2, rdim, zdim, rcentr, rleft, zmid
readf, lun, format = f2, rmaxis, zmaxis, simag, sibry, bcentr
readf, lun, format = f2, current, simag, xdum, rmaxis, xdum
readf, lun, format = f2, zmaxis, xdum, sibry, xdum, xdum

fpol    = fltArr ( nW )
pres    = fltArr ( nW )
ffprim  = fltArr ( nW )
pprime  = fltArr ( nW )
psizr   = fltArr ( nW, nH )
qpsi    = fltArr ( nW )

readf, lun, format = f2, fpol
readf, lun, format = f2, pres
readf, lun, format = f2, ffprim 
readf, lun, format = f2, pprime 
readf, lun, format = f2, psizr 
readf, lun, format = f2, qpsi

readf, lun, format = f3, nbbbs, limitr

bbbs    = fltArr ( 2, nbbbs )
lim    = fltArr ( 2, limitr )

readf, lun, format = f2, bbbs
readf, lun, format = f2, lim

rbbbs   = bbbs[0,*]
zbbbs    = bbbs[1,*]

rlim    = lim[0,*]
zlim    = lim[1,*]

pressw  = fltArr ( nW )
pwprim  = fltArr ( nW )
dmion  	= fltArr ( nW )
rhovn  	= fltArr ( nW )

if ~ keyword_set(noToroidalFlux) then begin
readf, lun, format = f4, kvtor, rvtor, nmass
if kvtor gt 0 then begin
	print, 'kvtor > 0 so reading pressq & pwprim'
	readf, lun, format = f2, pressw
	readf, lun, format = f2, pwprim
endif
if nmass gt 0 then begin
	print, 'nmass > 0 so reading dmion'
	readf, lun, format = f2, dmion
endif
readf, lun, format = f2, rhovn ; sqrt(toroidal flux)
endif

if keyword_set ( half ) then begin

	psizr	= psizr / 2.0
	fpol	= fpol / 2.0
	simag	= simag / 2.0
	sibry	= sibry / 2.0

endif

if keyword_set ( bFactor ) then begin

	psizr	= psizr * bFactor
	fpol	= fpol * bFactor
	simag	= simag * bFactor
	sibry	= sibry * bFactor

endif

if keyword_set ( bPolFactor ) then begin

	psizr	= psizr * bPolFactor
	simag	= simag * bPolFactor
	sibry	= sibry * bPolFactor

endif

if keyword_set ( cMod_lim_adjust ) then begin

	limitrNew = 0
	rLimNew = !NULL
	zLimNew = !NULL

	for i=0,limitr-1 do begin

		keep = 1
		if rlim[i] gt 0.75 and zlim[i] gt 0.22 then keep = 0
		if rlim[i] gt 0.75 and zlim[i] lt -0.22 then keep = 0
		if zlim[i] lt -0.45 then keep = 0
		if zlim[i] gt 0.5 then keep = 0

		if keep eq 1 then begin

			rLimNew = [rLimNew,rLim[i]]
			zLimNew = [zLimNew,zLim[i]]
			limitrNew = limitrNew + 1

		endif

	endfor

	limitr = limitrNew
	lim = fltArr(2,limitr)
	lim[0,*]   = rlimNew
	lim[1,*]   = zlimNew
	rLim = lim[0,*]
	zLim = lim[1,*]
	rLim = lim[0,*]
	zLim = lim[1,*]

endif


if keyword_set ( nstx_lim_adjust ) then begin

	rlim = [ $
		0.185100, $
		0.185100, $
		0.279400, $
		0.279400, $
		0.297900, $
		0.571200, $
		0.571200, $
		0.617000, $
		0.717000, $
		1.14330 , $
		1.41920 , $
		1.43580 , $
		1.58510 , $
		1.60910 , $
		1.64740 , $
		1.66130 , $
		1.67640 , $
		1.69080 , $
		1.69700 , $
		1.69570 , $
		1.68430 , $
		1.66410 , $
		1.64810 , $
		1.61180 , $
		1.58510 , $
		1.43580 , $
		1.41920 , $
		1.14330 , $
		0.717000, $
		0.617000, $
		0.571200, $
		0.571200, $
		0.297900, $
		0.279400, $
		0.279400, $
		0.185100, $
		0.185100, $
		0.185100 ]

	zlim = [ $
      0.00000, $
      1.00810, $
      1.17140, $
      1.57800, $
      1.60340, $
      1.60340, $
      1.62800, $
      1.62800, $
      1.62800, $
      1.43000, $
      1.03970, $
      0.99760, $
      0.54500, $
      0.49950, $
      0.30600, $
      0.23550, $
      0.15860, $
      0.08010, $
      0.00000, $
     -0.01770, $
     -0.11230, $
     -0.22100, $
     -0.30260, $
     -0.48600, $
     -0.54500, $
     -0.99760, $
     -1.03970, $
     -1.43000, $
     -1.62800, $
     -1.62800, $
     -1.62800, $
     -1.60340, $
     -1.60340, $
     -1.57800, $
     -1.17140, $
     -1.00810, $
      0.00000, $
      0.00000 ]

	lim = fltArr ( 2, n_elements ( rLim ) )

	lim[0,*]   = rlim
	lim[1,*]   = zlim

	limitr = n_elements ( lim[0,*] )

endif

if keyword_set ( iter_lim_adjust ) then begin

	rlim = [ $
			4.066, $
			4.325, $
			5.000, $
			5.779, $
			6.951, $
			7.578, $
			7.993, $
			8.270, $
			8.408, $
			8.322, $
			7.924, $
			7.318, $
			6.315, $
			5.796, $
			5.606, $
			5.588, $
			5.294, $
			5.260, $
			5.035, $
			4.671, $
			4.516, $
			4.204, $
			4.516, $
			4.516, $
			4.412, $
			4.083, $
			4.083, $
			4.066 ]



	zlim = [ $
			3.581, $
			4.256, $
			4.671, $
			4.498, $
			3.599, $	
			2.976, $
			2.215, $
			1.626, $
			0.588, $
			-0.519, $
			-1.401, $
			-2.318, $
			-3.253, $
			-3.374, $
			-3.702, $
			-4.550, $
			-4.291, $
			-3.979, $
			-3.754, $
			-3.754, $
			-3.927, $
			-3.893, $
			-3.235, $
			-2.976, $
			-2.768, $
			-2.561, $
			0.000, $
			3.581 ]

	;rlim = [ $
	;	8.24609, $
	;	8.24609, $
	;	8.12097, $
	;	4.45905, $
	;	4.43897, $
	;	4.43897, $
	;	4.45884, $
	;	8.12601, $
	;	8.24609, $
	;	8.24609 ]

	;zlim = [ $
    ;	 0.36876, $
    ;	 3.89886, $
    ;	 3.96721, $
    ;	 3.96721, $
    ;	 3.86538, $
    ;	-3.12342, $
    ;	-3.22969, $
    ;	-3.22969, $
    ;	-3.09334, $
    ;	 0.368762 ]
	
	lim = fltArr ( 2, n_elements ( rLim ) )

	lim[0,*]   = rlim
	lim[1,*]   = zlim

	limitr = n_elements ( lim[0,*] )

endif

if keyword_set ( reWrite ) then begin

	openw, lun, fileName+'.dlgMod', /get_lun
	
	printf, lun, format = f1, case_, idum, nW, nH
	printf, lun, format = f2, rdim, zdim, rcentr, rleft, zmid
	printf, lun, format = f2, rmaxis, zmaxis, simag, sibry, bcentr
	printf, lun, format = f2, current, simag, xdum, rmaxis, xdum
	printf, lun, format = f2, zmaxis, xdum, sibry, xdum, xdum
	printf, lun, format = f2, fpol
	printf, lun, format = f2, pres
	printf, lun, format = f2, ffprim 
	printf, lun, format = f2, pprime 
	printf, lun, format = f2, psizr 
	printf, lun, format = f2, qpsi
	printf, lun, format = f3, nbbbs, limitr
	printf, lun, format = f2, bbbs
	printf, lun, format = f2, lim

	if ~ keyword_set(noToroidalFlux) then begin
		printf, lun, format = f4, kvtor, rvtor, nmass
		if kvtor gt 0 then begin
			printf, lun, format = f2, pressw
			printf, lun, format = f2, pwprim
		endif
		if nmass gt 0 then begin
			printf, lun, format = f2, dmion
		endif
		printf, lun, format = f2, rhovn ; sqrt(toroidal flux)
	endif
	close, lun

endif

;   Calculate other desired quantities

rStep   = rdim / ( nW - 1 )
zStep   = zdim / ( nH - 1 )
fStep   = -( simag - sibry ) / ( nW - 1 )

R   = fIndGen ( nW ) * rStep + rleft
z   = fIndGen ( nH ) * zStep + zmid - zdim / 2.0

fluxGrid    = fIndGen ( nW ) * fStep + simag

;	Adjust to an up down symmetric field

;psizr	= ( psizr + reverse ( psizr, 1 ) ) / 2.0

;	Remember psi = -R * A

bR  = -1.0 * dlg_pDeriv ( psizr, 2, zStep ) / rebin ( R, nW, nH )
bz  = dlg_pDeriv ( psizr, 1, rStep ) / rebin ( R, nW, nH )
print, max(bR), min(bR)
print, max(bz), min(bz)
fPol_spline = spl_init ( fluxGrid, fPol )

fPolRZ  = reform ( spl_interp ( fluxGrid, fPol, fPol_spline, psizr[*] ), nW, nH )
fPolRZ2  = reform ( interpol ( fPol, fluxGrid, (psizr[*])<max(fluxGrid) ), nW, nH )
fPolRZ3  = reform ( interpol ( fPol, fluxGrid, psizr[*], /quad ), nW, nH )
fPolRZ4  = reform ( interpol ( fPol, fluxGrid, psizr[*], /spli ), nW, nH )

bPhi    = fPolRZ2 / rebin ( R, nW, nH )
APhi    = -psizr / rebin ( R, nW, nH )

bMag    = sqrt ( bR^2 + bPhi^2 + bz^2 )

buR = bR / bMag
buPhi   = bPhi / bMag
buZ = bZ / bMag

if keyword_set ( create_aorsa_nc_input ) then begin

	nX	= 1024 
	nY	= 1024 

	rSave	= rLeft + fIndGen ( nX ) * ( rDim / ( nX - 1 ) )
	zSave	= min(z) + fIndGen ( nY ) * ( zDim / ( nY - 1 ) )

	rSave2D	= rebin ( rSave, nX, nY )
	zSave2D	= transpose ( rebin ( zSave, nY, nX ) )

	brSave = interpolate ( br, ( rSave2D - rleft ) / rdim * (nw-1), $
   			( zSave2D + zmaxis - min ( z ) ) / zdim * (nh-1) )
	bzSave = interpolate ( bz, ( rSave2D - rleft ) / rdim * (nw-1), $
   			( zSave2D + zmaxis - min ( z ) ) / zdim * (nh-1) )
	bPhiSave = interpolate ( bPhi, ( rSave2D - rleft ) / rdim * (nw-1), $
   			( zSave2D + zmaxis - min ( z ) ) / zdim * (nh-1) )

	outFileName	= 'dlg_bField.nc'
	nc_id	= nCdf_create ( outFileName, /clobber )
	nCdf_control, nc_id, /fill
	
	nR_id	= nCdf_dimDef ( nc_id, 'nR', nX )
	nz_id	= nCdf_dimDef ( nc_id, 'nz', nY )
	scalar_id	= nCdf_dimDef ( nc_id, 'scalar', 1 )
	
	R_id = nCdf_varDef ( nc_id, 'R', [ nR_id ], /float )
	z_id = nCdf_varDef ( nc_id, 'z', [ nz_id ], /float )
	bR_id = nCdf_varDef ( nc_id, 'bR', [nR_id, nz_id], /float )
	bz_id = nCdf_varDef ( nc_id, 'bz', [nR_id, nz_id], /float )
	bPhi_id = nCdf_varDef ( nc_id, 'bPhi', [nR_id, nz_id], /float )

	nCdf_control, nc_id, /enDef
	
	nCdf_varPut, nc_id, R_id, rSave 
	nCdf_varPut, nc_id, z_id, zSave 
	nCdf_varPut, nc_id, bR_id, brSave 
	nCdf_varPut, nc_id, bz_id, bzSave 
	nCdf_varPut, nc_id, bPhi_id, bPhiSave 

	nCdf_close, nc_id

endif 



pPrime_spline   = spl_init ( fluxGrid, pPrime )
fPrime_spline   = spl_init ( fluxGrid, ffprim )

pPrime_Rz   = reform ( spl_interp ( fluxGrid, pPrime, pPrime_spline, psizr[*] ), nW, nH )
fPrime_Rz   = reform ( spl_interp ( fluxGrid, ffprim, fPrime_spline, psizr[*] ), nW, nH )

R2D = rebin ( R, nW, nH )
z2D = transpose ( rebin ( z, nH, nW ) )

JT_Rz   = R2D * pPrime_Rz + fPrime_Rz / R2D / u0

iiAxis = where ( abs ( R2d - rmaxis ) eq min ( abs ( R2d - rmaxis ) ) $
        and abs ( z2d - zmaxis ) eq min ( abs ( z2d - zmaxis ) ) )

JT_axis = JT_Rz[iiAxis]

;   Find points inside the boundary

mask   = intArr ( size ( psizr, /dim ) )

for i = 0, nW - 1 do begin
    for j = 0, nH - 1 do begin 

        q1_  = n_elements ( where ( ( R[i] - rbbbs gt 0 ) and ( z[j] - zbbbs gt 0 ), q1 ) )
        q2_  = n_elements ( where ( ( R[i] - rbbbs gt 0 ) and ( z[j] - zbbbs le 0 ), q2 ) )
        q3_  = n_elements ( where ( ( R[i] - rbbbs le 0 ) and ( z[j] - zbbbs gt 0 ), q3 ) )
        q4_  = n_elements ( where ( ( R[i] - rbbbs le 0 ) and ( z[j] - zbbbs le 0 ), q4 ) )

        if ( q1 gt 0 ) and ( q2 gt 0 ) and ( q3 gt 0 ) and ( q4 gt 0 ) then $
            mask[i,j]  = 1
               
    endfor
endfor
 
iiInside    = where ( mask gt 0 )
iiOutside   = where ( mask eq 0 )

    bInterpS    = { bR : bR, $
                    rleft : rleft, $
                    rdim : rdim, $
                    nW : nW, $
                    z : z, $
                    zdim : zdim, $
                    nH : nH, $
                    bPhi : bPhi, $
                    bz : bz }   

    ; Trace a field line from coords fieldLine = [r,t,z] 
	; --------------------------------------------------
 
	if keyword_set ( fieldLineIn ) then begin

		if keyword_set ( use_dlg_bField ) then begin

			R_orig	= R
			z_orig	= z
			bR_orig	= bR
			bPhi_orig	= bPhi
			bz_orig	= bz
			bInterpS_orig	= bInterpS

			cdfId = ncdf_open ( 'dlg_bField.nc', /noWrite ) 

				ncdf_varget, cdfId, 'R', R 
				ncdf_varget, cdfId, 'z', z 
				ncdf_varget, cdfId, 'bR', bR 
				ncdf_varget, cdfId, 'bPhi', bPhi 
				ncdf_varget, cdfId, 'bz', bz 

			ncdf_close, cdfId

    		bInterpS    = { bR : bR, $
        	            rleft : r[0], $
        	            rdim : r[-1]-r[0], $
        	            nW : n_elements(r), $
        	            z : z, $
        	            zdim : z[-1]-z[0], $
        	            nH : n_elements(z), $
        	            bPhi : bPhi, $
        	            bz : bz }   
		endif

		rStart	= fieldLineIn[0]
		tStart	= fieldLineIn[1]	
		zStart	= fieldLineIn[2]

        ; RK4 
		; ---

        rPos    = rStart
		tPos	= tStart
        zPos    = zStart
    
        rArray  = rStart
		tArray	= tStart
        zArray  = zStart

        stepCnt = 750
        dPhi    = -2 * !pi / 200.0

        for s=0,stepCnt-1 do begin

            bHere   = interpB ( bInterpS, rPos, zPos )

            K1_R  = dPhi * bHere.bR / bHere.bMag
            K1_z    = dPhi * bHere.bz / bHere.bMag
 
            bHere   = interpB ( bInterpS, rPos + K1_R / 2.0, zPos + K1_z / 2.0 )

            K2_R    = dPhi * bHere.bR / bHere.bMag
            K2_z    = dPhi * bHere.bz / bHere.bMag 
    
            bHere   = interpB ( bInterpS, rPos + K2_R / 2.0, zPos + K2_z / 2.0 )

            K3_R    = dPhi * bHere.bR / bHere.bMag
            K3_z    = dPhi * bHere.bz / bHere.bMag

            bHere   = interpB ( bInterpS, rPos + K3_R, zPos + K3_z )

            K4_R    = dPhi * bHere.bR / bHere.bMag
            K4_z    = dPhi * bHere.bz / bHere.bMag

            rPos    = rPos + ( K1_R + 2.0 * K2_R + 2.0 * K3_R + K4_R ) / 6.0
            zPos    = zPos + ( K1_z + 2.0 * K2_z + 2.0 * K3_z + K4_z ) / 6.0
			tPos	= tPos + dPhi
 
            rArray  = [ rArray, rPos ]
            tArray  = [ tArray, tPos ]
            zArray  = [ zArray, zPos ]
        
		endfor

		xArray	= rArray * cos ( tArray )
		yArray	= rArray * sin ( tArray )

		fieldLineOut	= [[rArray],[tArray],[zArray]]

		if keyword_set ( use_dlg_bField ) then begin

			R		= R_orig
			z		= z_orig
			bR		= bR_orig
			bPhi	= bPhi_orig
			bz		= bz_orig
			bInterpS= bInterpS_orig

		endif	

	endif

if keyword_set ( pAngle ) then begin
    
    ;   Calculate a poloidal angle coordinate which labels
    ;   the position along a a flux surface. This will be on 
    ;   the rz grid so we can create grad Chi.
    
    ;   Trace the poloidal field lines at each r (minor) coord.
 
    rMinor   = R - rmaxis
    rMinorRight = max ( rbbbs ) - rmaxis
    iiPositiveRMinor    = where ( rMinor gt 0 and rMinor le rMinorRight, iirMinorCnt )
    
    dl  = 0.001
   
    window, 0, xSize = 800, ySize = 800
    device, decomposed = 0
    loadct, 39, /silent
    !p.background = 255
    plot, [0,0], [0,0], /noData, $
        xRange = [0.0,3.0], yRange = [-1.4,1.3], xStyle = 9, yStyle = 9, color = 0
   
    rChi_all    = fltArr ( 12, iirMinorCnt )
    zChi_all    = fltArr ( 12, iirMinorCnt )

    rPos_all    = 0.0
    zPos_all    = 0.0
    chi_all     = 0.0
    gradChi_R_all   = 0.0
    gradChi_Z_all   = 0.0
    lengthP = 0.0
    lengthP_fluxGrid_R  = rMinor[iiPositiveRMinor] + rmaxis
    lengthP_fluxGrid_z  = lengthP_fluxGrid_R * 0.0 + zmaxis
    lengthP_fluxGrid    = interpolate ( psizr, ( lengthP_fluxGrid_R - rleft ) / rdim * (nW-1.0), $
        ( lengthP_fluxGrid_z - min ( z ) ) / zdim * (nH-1.0), cubic = -0.5 )


    for i = 0, iirMinorCnt - 1 do begin
    
        rStart   = rMinor[iiPositiveRMinor[i]] + rmaxis 
        zStart  = zmaxis
       
        ;   Try RK4 

        rPos    = rStart
        zPos    = zStart
    
        rArray  = rStart
        zArray  = zStart

        lArray  = 0.0
 
        stepCnt = 0
        thetaOld    = 2.0 * !Pi
        dPhi    = -2 * !pi / 100.0
        keepRunning = 1
        while keepRunning do begin

            bHere   = interpB ( bInterpS, rPos, zPos )

            K1_R  = dPhi * bHere.bR / bHere.bMag
            K1_z    = dPhi * bHere.bz / bHere.bMag
 
            bHere   = interpB ( bInterpS, rPos + K1_R / 2.0, zPos + K1_z / 2.0 )

            K2_R    = dPhi * bHere.bR / bHere.bMag
            K2_z    = dPhi * bHere.bz / bHere.bMag 
    
            bHere   = interpB ( bInterpS, rPos + K2_R / 2.0, zPos + K2_z / 2.0 )

            K3_R    = dPhi * bHere.bR / bHere.bMag
            K3_z    = dPhi * bHere.bz / bHere.bMag

            bHere   = interpB ( bInterpS, rPos + K3_R, zPos + K3_z )

            K4_R    = dPhi * bHere.bR / bHere.bMag
            K4_z    = dPhi * bHere.bz / bHere.bMag

            rPos    = rPos + ( K1_R + 2.0 * K2_R + 2.0 * K3_R + K4_R ) / 6.0
            zPos    = zPos + ( K1_z + 2.0 * K2_z + 2.0 * K3_z + K4_z ) / 6.0
 
            if stepCnt gt 0 then thetaOld    = theta
            theta   = aTan ( zPos - zmaxis, rPos - rmaxis )
            if theta lt 0 then theta = theta + 2.0 * !pi 
            
            rArray  = [ rArray, rPos ]
            zArray  = [ zArray, zPos ]
            
            dlO =  sqrt ( ( rStart - rArray[stepCnt] )^2 + ( zStart - zArray[stepCnt] )^2 )
            dl  = sqrt ( ( rPos - rArray[stepCnt] )^2 + ( zPos - zArray[stepCnt] )^2 )
            lArray  = [ lArray, lArray[stepCnt] + dl ]
 
            ++ stepCnt 
            if stepCnt gt 3 then oPlot, rArray, zArray, color = 0
            if ( theta - thetaOld ) gt 0 and stepCnt gt 10 then keepRunning = 0
        
		endWhile

        chiArray    = lArray / max ( lArray ) * ( theta + 2.0 * !pi )
        bHere   = interpB ( bInterpS, rArray, zArray )
        
        lengthP = [ lengthP, max ( lArray ) * 2.0 * !pi / ( theta + 2.0 * !pi ) ]

        dChi_dl = ( theta + 2.0 * !pi ) / max ( lArray )      
        gradChi_R = dChi_dl * bHere.bR / sqrt ( bHere.bR^2 + bHere.bz^2 )
        gradChi_Z   = dChi_dl * bHere.bz / sqrt ( bHere.bR^2 + bHere.bz^2 )

        rPos_all    = [ rPos_all, rArray ]
        zPos_all    = [ zPos_all, zArray ]
        chi_all = [ chi_all, chiArray ]
        gradChi_R_all = [ gradChi_R_all, gradChi_R ]
        gradChi_Z_all = [ gradChi_Z_all, gradChi_Z ]


        chi = 0 
        rChi    = rArray[0]
        zChi    = zArray[0]
        for j = 1, n_elements ( chiArray ) - 1 do begin
           if abs ( fix ( chiArray[j] * !radeg / 30.0 ) $
					- fix ( chiArray[j-1] * !radeg / 30.0 ) ) gt 0 then begin 

                chi = [ chi, chiArray[j] ]
                rChi    = [ rChi, rArray[j] ]
                zChi    = [ zChi, zArray[j] ]

            endif
        endfor
        
        ;   Interpolate the r,z location of a regular chi spacing

        rChi_   = interpol ( rArray, chiArray, fIndGen ( 12 ) * 30.0 * !dtor )
        zChi_   = interpol ( zArray, chiArray, fIndGen ( 12 ) * 30.0 * !dtor )

        rChi_all[*,i]   = rChi_        
        zChi_all[*,i]   = zChi_

        oPlot, rChi, zChi, psym = 5, color = 14 * 16 - 1
        for j = 1, n_elements ( chi ) - 1 do begin
            plots, [rChi[j],rChi[j]], [zChi[j],zChi[j]], $
				psym = 4, color = 4 * 16 - 1, symSize = abs ( j * 30.0 - chi[j] * !radeg ) * 2.0
        endfor
    endfor

rPos_all    = rPos_all[1:*]
zPos_all    = zPos_all[1:*]
chi_all    = chi_all[1:*]

triangulate, rPos_all, zPos_all, triangles, b

rz_chi  = triGrid ( rPos_all, zPos_all, chi_all, triangles, $
    [ rStep, zStep ], [ min ( R ), min ( z ), max ( R ), max ( z ) ] ); , exptrapolate = b )

gradChi_R_all   = gradChi_R_all[1:*]
gradChi_Z_all   = gradChi_Z_all[1:*]
lengthP = lengthP[1:*]

lengthP_spline = spl_init ( fluxGrid, fPol )
lengthPolRZ  = reform ( spl_interp ( lengthP_fluxGrid, lengthP, lengthP_spline, psizr[*] ), nW, nH )

for i = 0, 11 do begin

    oPlot, rChi_all[i,*], zChi_all[i,*], color = 0

endfor
!p.background = 0

endif
if keyword_set ( plot_ ) then begin

    device, decomposed = 0
    !p.multi = [ 0, 2, 2 ]
    
    contour, psizr, R, z
    oplot, rbbbs, zbbbs
    oplot, rlim, zlim
    
    contour, bPhi * mask, R, z, levels = fIndGen ( 21 ) / 5 - 2
    
    veloVect, bR * mask , bz * mask , R, z
    
    !p.multi = 0

endif

;   Create data structure

if keyword_set ( pAngle ) then begin

	eqdsk   = { case_ : case_, $
	            nW : nW, $
	            nH : nH, $
	            rDim : rdim, $
	            zDim : zdim, $
	            rLeft : rleft, $
	            zMid : zmid, $
	            rMAxis : rmaxis, $
	            zMAxis : zmaxis, $
	            siMag : siMag, $
	            siBry : sibry, $
	            rCentr : rcentr, $
	            bCentr : bcentr, $
	            current : current, $
	            fPol : fpol, $
	            pres : pres, $
	            ffPrim : ffprim, $
	            pPrime : pprime, $
	            psizr : psizr, $
	            qPsi : qpsi, $
	            nbbbs : nbbbs, $
	            limitr : limitr, $
	            rbbbs : rbbbs, $
	            zbbbs : zbbbs, $
	            rLim : rlim, $
	            zLim : zlim, $
	            bR : bR, $
	            bz : bz, $
	            bPhi : bPhi, $
	            bMag : bMag, $
	            mask : mask, $
	            iiInside : iiInside, $
	            iiOutside : iiOutside, $
	            fPolRZ : fPolRZ, $
	            chiRZ : rz_chi, $
	            lengthPolRZ : lengthPolRZ, $
	            rStep : rStep, $
	            zStep : zStep, $
	            fStep : fStep, $
	            R : R, $
	            z : z, $
				R2D : rebin ( R, nW, nH ), $
				z2D : transpose ( rebin ( z, nH, nW ) ), $ 
	            fluxGrid : fluxGrid, $
	            lengthP : lengthP, $
	            APhi : APhi, $
	            buR : buR, $
	            buPhi : buPhi, $
	            buZ : buZ }

endif else begin

	eqdsk   = { case_ : case_, $
	            nW : nW, $
	            nH : nH, $
				idum : idum, $
				xdum : xdum, $
	            rDim : rdim, $
	            zDim : zdim, $
	            rLeft : rleft, $
	            zMid : zmid, $
	            rMAxis : rmaxis, $
	            zMAxis : zmaxis, $
	            siMag : siMag, $
	            siBry : sibry, $
	            rCentr : rcentr, $
	            bCentr : bcentr, $
	            current : current, $
	            fPol : fpol, $
	            pres : pres, $
	            ffPrim : ffprim, $
	            pPrime : pprime, $
	            psizr : psizr, $
				rhovn : rhovn, $ ; This is sqrt(toroidal flux)
	            qPsi : qpsi, $
	            nbbbs : nbbbs, $
	            limitr : limitr, $
	            rbbbs : rbbbs, $
	            zbbbs : zbbbs, $
	            rLim : rlim, $
	            zLim : zlim, $
	            bR : bR, $
	            bz : bz, $
	            bPhi : bPhi, $
	            bMag : bMag, $
	            mask : mask, $
	            iiInside : iiInside, $
	            iiOutside : iiOutside, $
	            fPolRZ : fPolRZ, $
	            rStep : rStep, $
	            zStep : zStep, $
	            fStep : fStep, $
	            R : R, $
	            z : z, $
				R2D : rebin ( R, nW, nH ), $
				z2D : transpose ( rebin ( z, nH, nW ) ), $ 
	            fluxGrid : fluxGrid, $
	            APhi : APhi, $
	            buR : buR, $
	            buPhi : buPhi, $
	            buZ : buZ, $
                Jt_axis : Jt_axis }

endelse

	
return, eqdsk

end
