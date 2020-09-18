% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program 
% under the terms of the GNU Lesser General Public 
% License, either version 3 of the License, or any 
% later version.

% Conjugate Gradient 

% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program 
% under the terms of the GNU Lesser General Public 
% License, either version 3 of the License, or any 
% later version.

% Conjugate Gradient 

function [x,i] = CG(A,b,x0,e,m)
r = b - A*x0;               
% n+2*n^2
d = r;
x = x0;
i = 0;
recal = round(m/100);

while i<m && norm(r)/norm(b)>e
    q = r'*r;                
    %2n
    Ad = A*d;
    d_A_norm = d'*Ad; 
    % 4n^2
    if d_A_norm ==0    
        return
    end
    alpha = q/d_A_norm;             
    %1
    x = x + alpha*d;         
    %2n
    if mod(i,recal) == 0
        r = b-A*x;
    else
        r = r - alpha*Ad;
    end            
    %n+2n^2 or 2n
    d = r + (r'*r/q)*d;      
    %4n+1
    i = i + 1;
    %1
end
end

% each iteration: 4*n^2+11*n+4