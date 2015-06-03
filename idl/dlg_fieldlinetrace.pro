function dlg_fieldlinetrace, b, StartPoint_CYL, $
    direction = direction, $
    ds = _ds, $
    nS = nS, $
	FieldLine_XYZ = fieldLine_XYZ, $
	FieldLine_CYL = fieldLine_CYL, $
    B_FieldLine_XYZ = B_FieldLine_XYZ, $
    B_FieldLine_CYL = B_FieldLine_CYL, $
	perp = perp, $
	use_dlg_bField = use_dlg_bField
	
	;r0 = g.rmaxis
	;z0 = g.zmaxis

; Trace a field line from coords fieldLine = [r,t,z] 
; --------------------------------------------------

	rStart	= StartPoint_CYL[0]
	tStart	= StartPoint_CYL[1]	
	zStart	= StartPoint_CYL[2]

    ; RK4 
	; ---

    c_CYL = [rStart,tStart,zStart]
    c_XYZ = Coords_CYL_to_XYZ(c_CYL)    
    b_XYZ = bHere_XYZ (b, c_XYZ)
    b_CYL = vector_XYZ_to_CYL (c_XYZ,b_XYZ)  

    PolAngle = 0
    TorAngle = 0

    c_XYZ_array = c_XYZ
    c_CYL_array = c_CYL

    b_XYZ_array = b_XYZ
    b_CYL_array = b_CYL

    if keyword_set(nS) then stepCnt = nS else stepCnt = 10000
    if keyword_set(_ds) then _dS = _ds else _dS = 0.01
    if keyword_set(direction) then _TraceSign = direction else _TraceSign = 1

    dS = _TraceSign * _dS

    for s=0,stepCnt-1 do begin

        b_XYZ  = bHere_XYZ ( b, c_XYZ, bMag=bMagTrace, Perp=perp )

        K1  = dS * b_XYZ / bMagTrace 

        b_XYZ  = bHere_XYZ ( b, c_XYZ + K1 / 2.0, bMag=bMagTrace, Perp=perp )

        K2  = dS * b_XYZ / bMagTrace 

        b_XYZ  = bHere_XYZ ( b, c_XYZ + K2 / 2.0, bMag=bMagTrace, Perp=perp )

        K3  = dS * b_XYZ / bMagTrace 

        b_XYZ  = bHere_XYZ ( b, c_XYZ + K3, bMag=bMagTrace, Perp=perp )

        K4  = dS * b_XYZ / bMagTrace 

        c_XYZ   = c_XYZ + ( K1 + 2 * K2 + 2 * K3 + K4 ) / 6.0
        c_CYL   = Coords_XYZ_to_CYL(c_XYZ) 
       
        b_XYZ   = bHere_XYZ ( b, c_XYZ, bMag=bMagTrace, Perp=perp )
        b_CYL   = vector_XYZ_to_CYL (c_XYZ,b_XYZ)

        c_XYZ_array  = [ [c_XYZ_array],[c_XYZ] ]
        c_CYL_array  = [ [c_CYL_array],[c_CYL] ]

        ;ThisTorAngle = c_CYL_array[1,s+1]*!radeg
        ;PrevTorAngle = c_CYL_array[1,s]*!radeg
        ;TorStep = ThisTorAngle-PrevTorAngle
        ;if(TorStep gt 20) then TorStep = TorStep-360
        ;if(TorStep lt -20) then TorStep = TorStep+360
        ;TorAngle = TorAngle + TorStep

        ;ThisPolAngle = atan((c_CYL_array[2,s+1]-z0),(c_CYL_array[0,s+1]-r0))*!radeg
        ;PrevPolAngle = atan((c_CYL_array[2,s]-z0),(c_CYL_array[0,s]-r0))*!radeg
        ;PolStep = ThisPolAngle-PrevPolAngle
        ;if(PolStep gt 20) then PolStep = PolStep-360
        ;if(PolStep lt -20) then PolStep = PolStep+360
        ;PolAngle = PolAngle + PolStep 

        b_XYZ_array  = [ [b_XYZ_array],[b_XYZ] ]
        b_CYL_array  = [ [b_CYL_array],[b_CYL] ]

	endfor

    ;SafetyFactor = TorAngle / PolAngle

	FieldLine_XYZ= c_XYZ_array
    B_FieldLine_XYZ = b_XYZ_array
	FieldLine_CYL= c_CYL_array
    B_FieldLine_CYL = b_CYL_array

	return, FieldLine_CYL

end
