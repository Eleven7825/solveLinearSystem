% This is a linear equation solver using PLU fractorization, be sure that
% LU function is in the same folder

% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of the
% GNU Lesser General Public License, either version 3 of the License, or 
% any later version.

% Syntax x = PLU_solve(A,b)
% For linear equation Ax = b, take A and b as parameters, output the
% solution x. 

function x = PLU_solve(A,b)
sz = size(A);
sz = sz(1,1);

%initialize PLU
[P,L,U] = LU(A);
Pb = P*b;
Ux = zeros(sz,1);
x = zeros(sz,1);

%Solve the first equation
for i = 1 : sz
    Ux(i,1) = (Pb(i,1)-L(i,:)*Ux)/L(i,i);
end

%Solve the second equation
for j = 0 : sz-1
    x(sz-j,1) = (Ux(sz-j,1)-U(sz-j,:)*x)/U(sz-j,sz-j);
end
end