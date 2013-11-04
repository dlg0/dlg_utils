
function Coords_CYL_to_XYZ, c_CYL 

    r = c_CYL[0]
    t = c_CYL[1]
    z = c_CYL[2]

    x = r * cos(t)
    y = r * sin(t)
    
    return, [x,y,z]

end

function Coords_XYZ_to_CYL, c_XYZ

    x = c_XYZ[0]
    y = c_XYZ[1]
    z = c_XYZ[2]

    r = sqrt(x^2+y^2)
    t = atan(y,x)

    return, [r,t,z]

end


function Get_CYL_to_XYZ_RotMat, c_CYL 

    r = c_CYL[0]
    t = c_CYL[1]
    z = c_CYL[2]

    return, [ $
        [cos(t),-sin(t),0],$
        [sin(t),cos(t),0],$
        [0,0,1]]        

end 

function vector_XYZ_to_CYL, c_XYZ, b_XYZ

    c_CYL = Coords_XYZ_to_CYL(c_XYZ)
    Rot2XYZ = Get_CYL_to_XYZ_RotMat(c_CYL)
    Rot2CYL = transpose(Rot2XYZ)

    return, transpose(Rot2CYL ## transpose(b_XYZ))

end

function vector_CYL_to_XYZ, c_CYL, b_CYL

    c_XYZ = Coords_CYL_to_XYZ(c_CYL)
    Rot2XYZ = Get_CYL_to_XYZ_RotMat(c_CYL)

    return, transpose(Rot2XYZ ## transpose(b_CYL))

end

pro Coords_test

    
    c_CYL = [1.5,65*!dtor,1]
    c_XYZ = Coords_CYL_to_XYZ(c_CYL)

    ; Compare with anlaytic result from Chen. pg 31.
    v_CYL = [3*cos(c_CYL[1]),-2*c_CYL[0],5]
    v_XYZ = [$
            3*c_XYZ[0]^2/(c_XYZ[0]^2+c_XYZ[1]^2)+2*c_XYZ[1],$
            3*c_XYZ[0]*c_XYZ[1]/(c_XYZ[0]^2+c_XYZ[1]^2)-2*c_XYZ[0],$
            5]

    v_XYZ_test = Vector_CYL_to_XYZ(c_CYL,v_CYL)

    v_CYL_test = Vector_XYZ_to_CYL(c_XYZ,v_XYZ)

    print, v_XYZ, v_XYZ_test
    print, v_CYL, v_CYL_test

end


