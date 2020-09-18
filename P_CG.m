% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program 
% under the terms of the GNU Lesser General Public 
% License, either version 3 of the License, or any 
% later version.

% Preconditioned Conjugate Gradient

function [x,i] = P_CG(A,b,x0,e,m)
r = b - A*x0;

% preconditioning
M_1 = 1./diag(A);
% M_1 = ones(size(r));
d = M_1 .* r;
s =  d;
x = x0;
i = 0;

while i<m && norm(r)/norm(b)>e
    q = r' * s;
    % 2n
    d_A_norm = d' * A * d;
    % 4n^2
    if d_A_norm == 0
        return
    end
    alpha = q / d_A_norm;
    % 1
    x = x + alpha*d;
    % 2n
    r = b - A*x;
    % n+2n^2
    s = M_1.*r;
    % n
    d = s + (r'*s/q)*d;
    % 4n+1
    i = i+1;
    % 1
end
end

