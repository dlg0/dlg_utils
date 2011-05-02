
pro oversample_boundary,  r, z, newR, newz

	newR	= r[0]
	newZ	= z[0]

	for i=0,n_elements(r)-2 do begin

		m	= ( z[i+1]-z[i] ) $
				/ ( r[i+1] - r[i] )
		b	= z[i] - m * r[i]

		d	= sqrt ( ( r[i+1] - r[i] )^2 $
				+ ( z[i+1] - z[i] )^2 )

		nExtra	= 40
		dStep	= (r[i+1] - r[i]) / nExtra

		for j = 0, nExtra - 1 do begin

			if dStep ne 0 then begin
				newR	= [ newR, r[i] + dStep*j ]
				newZ	= [ newZ, m * (r[i] + dStep*j) + b ]
			endif 

		endfor

	endfor

	newR	= [r[*], newR]
	newz	= [z[*], newz]	
	
end 


