% PLU factorization

% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of the
% GNU Lesser General Public License, either version 3 of the License, or 
% any later version.

% Syntax: [P,L,U] = PLU(A) 
% Take a square matrix as input, gives a permutation matrix P, a lower
% triangular matrix L, and a upper triangular matrix U which satisfy PA =
% LU.

% This is an older version of LU, which require more flops.

function [P,L,U] = PLU(A)
format short
%size of the matrix
global sz
sz = size(A);
sz = sz(1,1);

%initial PLU
P = eye(sz);
L = P;
U = A;

%proccess of PLU
for i = 1 : sz
    [P1,L1,U1] = exchange(P,L,U,i);
    [P2,L2,U2] = replace(P1,L1,U1,i);
    P = P2;
    L = L2;
    U = U2;
end
end

%row exchange function
function [P,L,U]=exchange(P,L,U,i)
global sz
I = i;
Max = abs(U(i,i));
for j = i:sz
    if abs(U(j,i))>Max
        I = j;
        Max = abs(U(I,i));
    end
end 

% [M,I] = max(abs(U(i:sz,i)));
Q = eye (sz);
Q([i I],:)=Q([I i],:);
P = Q*P;
L=Q*L*Q;
U=Q*U;
end

%row replacement function
function [P,L,U] = replace(P,L,U,i)
global sz
E = eye(sz);
F = E;
u = U(i,i);
for j = i+1 : sz
    E(j,i) = -U(j,i)/u;
    F(j,i) = U(j,i)/u;
    
end
L = L*F;
U = E*U;
end
