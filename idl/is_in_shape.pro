function is_in_shape, r, z, rBry, zBry

	mask	= intArr ( size(r,/dim) )
	maskA	= mask

	rBry2D	= rebin ( rBry[*], n_elements (rBry[*]), n_elements ( r[*] ) )
	zBry2D	= rebin ( zBry[*], n_elements (zBry[*]), n_elements ( z[*] ) )

	r2D	= transpose ( rebin ( r[*], n_elements ( r[*] ), n_elements (rBry[*]) ) )
	z2D	= transpose ( rebin ( z[*], n_elements ( z[*] ), n_elements (zBry[*]) ) )

	q1_  = where ( ( r2D - rBry2D gt 0 ) and ( z2D - zBry2D gt 0 ), q1 )
	q2_  = where ( ( r2D - rBry2D gt 0 ) and ( z2D - zBry2D le 0 ), q2 )
	q3_  = where ( ( r2D - rBry2D le 0 ) and ( z2D - zBry2D gt 0 ), q3 )
	q4_  = where ( ( r2D - rBry2D le 0 ) and ( z2D - zBry2D le 0 ), q4 )

	q12d	= intArr ( size ( r2D, /dim ) )
	q22d	= intArr ( size ( r2D, /dim ) )
	q32d	= intArr ( size ( r2D, /dim ) )
	q42d	= intArr ( size ( r2D, /dim ) )

	if q1 gt 0 then q12d[q1_]	= 1
	if q2 gt 0 then q22d[q2_]	= 1
	if q3 gt 0 then q32d[q3_]	= 1
	if q4 gt 0 then q42d[q4_]	= 1

	q11d	= total ( q12d, 1 )
	q21d	= total ( q22d, 1 )
	q31d	= total ( q32d, 1 )
	q41d	= total ( q42d, 1 )

	iiIn	= where ( q11d gt 0 and q21d gt 0 and q31d gt 0 and q41d gt 0, iiInCnt )
	if iiInCnt gt 0 then $
			maskA[iiIn]	= 1

	;for i = 0, n_elements(r[*])-1 do begin

	;    	q1_  = where ( ( r[i] - rBry gt 0 ) and ( z[i] - zBry gt 0 ), q1 )
	;    	q2_  = where ( ( r[i] - rBry gt 0 ) and ( z[i] - zBry le 0 ), q2 )
	;    	q3_  = where ( ( r[i] - rBry le 0 ) and ( z[i] - zBry gt 0 ), q3 )
	;    	q4_  = where ( ( r[i] - rBry le 0 ) and ( z[i] - zBry le 0 ), q4 )

	;    	if ( q1 gt 0 ) and ( q2 gt 0 ) and ( q3 gt 0 ) and ( q4 gt 0 ) then begin

	;			mask[i] = 1

	;		endif
	;endfor

	return, reform ( maskA, size(r, /dim ) )

end 


