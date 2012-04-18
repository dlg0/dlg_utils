pro dlg_rootPlotter, x, y, imag=imag
		;x = runData.r
		;y = kRPlot
		for g=0,3 do begin
			plotY = y[0,g]
			lastY = y[0,g]		
			for i=0,n_elements(x)-2 do begin
				distance = sqrt((y[i+1,*]-lastY)^2)
				iiKeep = where(distance eq distance)
				distance = distance[iiKeep]
				nextyi = where(distance eq min(distance),iyCnt)
				nextyi = iiKeep[nextyi]
				lastY = y[i+1,nextyi[0]]
				plotY = [plotY,lastY]
				y[i+1,nextyi[0]] = !values.f_nan
			endfor
			if keyword_set(imag) then begin
				p=plot(x,plotY,thick=2.0,transparency=50,color='r',/over)
			endif else begin
				p=plot(x,plotY,thick=2.0,color='b',/over)
			endelse
		endfor
end


