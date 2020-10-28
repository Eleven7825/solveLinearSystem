function [eig_value,eig_vector,i] = IPM(A,mu,x0,m,e)
% initialize the method 
x = x0;
B = A - mu*eye(size(A));
i = 0;
ERR = 1+e;
Xp = 0;
for j = 1:length(x)
    if abs(x(j)) > Xp
        Xp = x(j);
    end
end
x = x/Xp;

while i < m && ERR > e
    y = Linsolver(B,x);
    for j = 1 : length(y)
        Yp = 0;
        if abs(y(j)) > Yp
            Yp = y(j);
        end
    end
    ERR = 0;
    errvec = x - (y/Yp);
    for j = 1:length(errvec)
        if abs(errvec(j)) > ERR
            ERR = abs(errvec(j));
        end
    end
    x = y/Yp;
    i = i+1;
end

if i < m
    eig_value = 1/Yp+mu;
    eig_vector = x;
    return
else
    disp('max m exceeded!')
    return
end
end
