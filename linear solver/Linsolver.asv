% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of 
% the GNU Lesser General Public License, either version 3 of the License, 
% or any later version.

% linear solver

function x = Linsolver(A,b)
sizeA = size(A);
sizeb = size(b);
n = sizeA(2);
n2 = sizeA(2);
if sizeb(1) ~= n
    error('Size of A must match size of b!')
end

% for n<=32, use PLU_solve function
if n<= 32
    x = PLU_solve(A,b);
    return
end

% determine whether symmetric or not, if not, use least square solution
if n==n2 && issymmetric(A)
    B = A;
    y = b;
else
    B = A'*A;
    y = A'*b;
end

density = nnz(B)/n^2;
x0 = zeros(n,1);
e = 0.01;
m = 5000;
recal = 50;

if (density < 16/n) && all(diag(B) ~= 0)
    x = P_CG(B,y,x0,e,m,recal);
    return
else
    x = CG(B,y,x0,e,m,recal);
    return
end
end

function [x,i] = P_CG(A,b,x0,e,m,recal)
r = b - A*x0;
% preconditioning
M_1 = 1./diag(A);
d = M_1 .* r;
s =  d;
x = x0;
i = 0;

while i<m && norm(r)/norm(b)>e
    q = r' * s;
    Ad = A*d;
    d_A_norm = d' * Ad;
    if d_A_norm == 0
        return
    end
    alpha = q / d_A_norm;
    x = x + alpha * d;
    if mod(i,recal) == 0
       r = b - A * x;
    else
       r = r - alpha * Ad;
    end
    s = M_1 .* r;
    d = s + (r' * s / q) * d;
    i = i + 1;
end
end

function [x,i] = CG(A,b,x0,e,m,recal)
r = b - A*x0;
d = r;
x = x0;
i = 0;

while i<m && norm(r)/norm(b)>e
    q = r'*r;                
    Ad = A*d;
    d_A_norm = d'*Ad; 
    if d_A_norm ==0    
        return
    end
    alpha = q / d_A_norm;             
    x = x + alpha * d;         
    if mod(i,recal) == 0
       r = b - A * x;
    else
       r = r - alpha * Ad;
    end            
    d = r + (r' * r / q) * d;      
    i = i + 1;
end
end

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