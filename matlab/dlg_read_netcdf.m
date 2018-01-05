function [nc] = dlg_read_netcdf(file)

nc = containers.Map();

ncI = ncinfo(file);

[mVars,nVars] = size(ncI.Variables);

for v=1:nVars
    
    thisVarStr = ncI.Variables(v).Name;
    nc(thisVarStr) = ncread(file,thisVarStr);
    
end

% Combine foo_re and foo_im vars as foo=complex(foo_re,foo_im)

for v=1:nVars
   
    varStr = ncI.Variables(v).Name;
    
    reg = regexp(varStr,'_im$');
    
    if sum(size(reg))>0 
        
        cmplxStr = varStr(1:reg-1);
        
        nc(cmplxStr) = complex(nc(strcat(cmplxStr,'_re')),nc(strcat(cmplxStr,'_im')));
        
    end
    
end

end