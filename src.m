%Hadeel Mabrouk - Mariam Abulela - Sara Ahmed - Nadine Sakr
%Matlab Project

%CODE
clc
clear
close all

disp("Please, enter the values of the matrix A");
A = input("A = ");
disp("Please, define the desired inner product by inserting its positive coefficients as follows [c1 c2 c3 ... cn]");
C = input("C = ");

[mA,nA] = size(A);
[mC,nC] = size(C);

if (rank(A)~=nA)
    disp("The column-vectors of A are not linearly independent!");
elseif (mC~=1 || nC~=mA || any(C(1,:)<=0))
    disp("The inserted coefficient matrix is invalid!")
else
    [Q,R] = QR(A,C);
    disp("Q = ");
    disp(Q);
    disp("R = ");
    disp(R);
    disp("QR =");
    disp(Q*R);
    disp("A = QR");
end


function [Q, R] = QR(X,C)
    [m,n] = size(X);
    R = zeros(n,n);
    Q = zeros(m,n);
    
    for i = 1:n %this is to calculate the orthonormal matrix Q
        Q(:,i) = X(:,i);
        for j=1:i-1
            Q(:,i) = Q(:,i)-(inner(X(:,i),Q(:,j),C)/inner(Q(:,j),Q(:,j),C))*Q(:,j);
        end
        Q(:,i)=Q(:,i)./sqrt(inner(Q(:,i),Q(:,i),C));
    end
    
    for i = 1:n %this is to calculate the upper triangular matrix R
        for j = i:n
            R(i,j) = inner(X(:,j),Q(:,i),C);
        end
    end
end

function y = inner(a,b,d)
    c=0;      % initialize the variable c
    n= length(a);    % get the length of the vector a
    for k=1:n    % start the loop
      c=c+a(k)*b(k)*d(k);  % update c by the k-th product in inner product
    end      % end loop
    y = c;      % print value of c = inner product
end

