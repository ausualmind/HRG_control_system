clc
clear
A=[1 5 1 0 0;
2 1 0 1 0;
1 1 0 0 1];
b=[40,20,12];
c=[-3,-5,0,0,0];
Simplex(A,b,c,[3,4,5],10);
function Simplex(A,b,c,bases,iterations)
% A, A matrix in Ax = b, should be a tableau
% b, b vector in Ax = b,
% c, c vector in object function, c'*x
% bases, index vector for column basis vectors
    [m, n] = size(A);
    if m ~= length(bases)
        error("ERROR: # of basis vectors does not match with the input tableau");
    end
    baseB = bases;
    baseD = setdiff(1:n, bases);
    Simplex1(A, eye(m), b, c, baseB, baseD,iterations);
end

function [B, x, baseB, baseD] = Simplex1(A,B,x,c,baseB,baseD,iterations)
% A, tableau
% B, Basis block in a tableau
% x, b vector in Ax = b,
% c, c vector in object function, c'*X
% baseB, baseD; index vector for column basis vectors
    [m, n] = size(A);
    for k = 1:iterations
        cb = c(baseB);
        cd = c(baseD);
        D = A(:,baseD);
        r = cd-((cb*B)*D)';
        q =SelectQ(r);
        if q == 0 % Optimal solution is found
%             println("Minimum cost = ",cb'*x);
            break;
        end
        xq=B*D(:,q);
        p=SelectP(x, xq);
        if p == 0 % Not-bounded
%             println("Not bounded!");
            break;
        end
        B = PivotalTransform(hcat(B,x,xq),p, m+2);
        x = B(:, m+1);
        B = B(:, 1:m);
%         println("B = ",B);
        temp =baseB(p);
        baseB(p)= baseD(q);% Basis exchange
        baseD(q)=temp;
    end
%     println("WARNING: ", iterations ," iterations have been");
end

function b = PivotalTransform(a, p, q)
    [m,n] = size(a);
    b = zeros(m,n);
    b(p,:) = a(p,:) ./ a(p,q);
    for i = 1:p-1
        b(i,:) = a(i,:) - b(p,:) .* a(i,q) ;%vector .* scalar
    end
    for i = p+1:m
        b(i,:) = a(i,:) - b(p,:) .* a(i,q); % vector .* scalar
    end
end

function q = SelectQ(r) %进基的选择
    n=length(r);
    q=0;
    qv=0;
    for j =1:n
        if r(j)<qv
            q =j;
            qv=r(j);
        end
    end
end

function p = SelectP(x0, xq) % 离基的选择
    m = length(x0);
    if m ~= length(xq)
        error("ERROR: vectors x_q and x_0 does not have the same length.")
    end
    p = 0;
    pv= Inf;
    for i = 1:m
        if xq(i)> 0
            temp =x0(i)/xq(i);
            if temp < pv
                p = i;
                pv= temp;
            end
        end
    end
end
% . . .
% Modified Simplex Algorithm for Linear Programming given an initial tableau
% Input: tableau matrix (a) for the initial basic feasible solution
% basis index vector, where column basis vectors are located
% Output：Minimum solution the objective function
% Note: extra input and output arguments are for the 2-phase algorithm
% ......