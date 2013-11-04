@xyz_cyl

function interpB,  bStruct, rPos, zPos

    bRHere  = interpolate ( bStruct.bR, ( rPos - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( zPos - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )
    bPhiHere  = interpolate ( bStruct.bPhi, ( rPos - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( zPos - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )
    bzHere  = interpolate ( bStruct.bz, ( rPos - bStruct.rleft ) / bStruct.rdim * (bStruct.nW-1.0), $
        ( zPos - min ( bStruct.z ) ) / bStruct.zdim * (bStruct.nH-1.0), cubic = -0.5 )
    
    bMag    = sqrt ( bRHere^2 + bPhiHere^2 + bzHere^2 )

    bOut    = { bR : bRHere, $
                bPhi :bPhiHere, $
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


function bHere_XYZ, bInterpS, c_XYZ, bMag=bMag

    c_CYL = Coords_XYZ_to_CYL(c_XYZ)

    x = c_XYZ[0]
    y = c_XYZ[1]
    z = c_XYZ[2]
 
    r = c_CYL[0]
    t = c_CYL[1] 

    bHere_CYL   = interpB ( bInterpS, r, z )
    b_CYL = [bHere_CYL.bR,bHere_CYL.bPhi,bHere_CYL.bZ]

    b_XYZ = vector_CYL_to_XYZ(c_CYL,b_CYL)

    bMag = sqrt(b_XYZ[0]^2+b_XYZ[1]^2+b_XYZ[2]^2)

    return, b_XYZ

end


