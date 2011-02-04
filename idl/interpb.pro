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


