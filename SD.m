% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program 
% under the terms of the GNU Lesser General Public 
% License, either version 3 of the License, or any 
% later version.

% Steepest descent

function [x,i] = SD(A,b,x0,e,m)
i = 0;
x = x0;
r = m;
recal = round(m/100);

while i<m && norm(r)/norm(b)>e
    r = b-A*x;
    delta = r'*r;
    Ar = A*r;
    r_A_norm = r'*Ar;
    
    if r_A_norm == 0
        return   
    end
    
    if mod(i,recal) == 0
        r = b-A*x;
    else
        r = r - alpha*Ar;
    end
    
    alpha = delta/r_A_norm;
    x = x + alpha*r;
    i = i+1;
end
end


