; Function to read a netcdf file and return
; an IDL hash of key/value pairs with keys 
; equal to the variable names.
; Requires IDL 8.6.1

function dlg_read_netcdf, file

ncdf_list, file, vname=vars, /var, /quiet

nc = dictionary()

nVars = n_elements(vars)

for v=0,nVars-1 do begin

    varStr = vars[v] 
    ncdf_get, file, varStr, var, /quiet 
    nc[vars[v]] = var[varStr,'value']

endfor

; Combine _re and _im vars into single complex variable

for v=0,nVars-1 do begin

    varStr = vars[v] 

    if stRegex(varStr,'\_im$',/bool) then begin

        varName = (strSplit(varStr, '\_im$',/regex,/extract))[0]

        nc[varName] = complex( nc[varName+'_re'], nc[varName+'_im'] )

    endif

endfor

return, nc

end
