@xyz_cyl

function interpB,  b, _r, _z

    bRHere  = interpolate ( b.br, ( _r - b.r[0] ) / b.rSize * (b.nR-1.0), $
        ( _z - b.z[0] ) / b.zSize * (b.nZ-1.0), cubic = -0.5 )
    bTHere  = interpolate ( b.bt, ( _r - b.r[0] ) / b.rsize * (b.nR-1.0), $
        ( _z - b.z[0] ) / b.zSize * (b.nZ-1.0), cubic = -0.5 )
    bzHere  = interpolate ( b.bz, ( _r - b.r[0] ) / b.rSize * (b.nR-1.0), $
        ( _z - b.z[0] ) / b.zSize * (b.nZ-1.0), cubic = -0.5 )
    
    bMag    = sqrt ( bRHere^2 + bTHere^2 + bzHere^2 )

    bOut    = { br : bRHere, $
                bt : bTHere, $
                bz : bzHere, $
                bMag : bMag }

    return, bOut

end 

function bHere_CYL, bInterpS, c_CYL, bMag=bMag

    r = c_CYL[0]
    t = c_CYL[1] 
    z = c_CYL[2]

    bHere_CYL   = interpB ( bInterpS, r, z )
    b_CYL = [bHere_CYL.bR,bHere_CYL.bPhi,bHere_CYL.bZ]

    bMag = sqrt(b_CYL[0]^2+b_CYL[1]^2+b_CYL[2]^2)

    return, b_CYL

end


function bHere_XYZ, b, c_XYZ, bMag=bMag, perp=perp

    c_CYL = Coords_XYZ_to_CYL(c_XYZ)

    x = c_XYZ[0]
    y = c_XYZ[1]
    z = c_XYZ[2]
 
    r = c_CYL[0]
    t = c_CYL[1] 

    bHere_CYL   = interpB ( b, r, z )
    b_CYL = [bHere_CYL.br,bHere_CYL.bt,bHere_CYL.bz]

    b_XYZ = vector_CYL_to_XYZ(c_CYL,b_CYL)

    bMag = sqrt(b_XYZ[0]^2+b_XYZ[1]^2+b_XYZ[2]^2)

	if keyword_set(perp) then begin
		; Add capability to return a vector perp to b
		zU = [0,0,1]
		b_XYZ = [zU[1]*b_XYZ[2]-zU[2]*b_XYZ[1],$
					 -(zU[0]*b_XYZ[2]-zU[2]*b_XYZ[0]),$
					 zU[0]*b_XYZ[1]-zU[1]*b_XYZ[0]]
		bMag = sqrt(b_XYZ[0]^2+b_XYZ[1]^2+b_XYZ[2]^2)
	endif

    return, b_XYZ

end


