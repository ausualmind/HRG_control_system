clc
clear
f=@ffun;g=@gfun;
[xn, fn, gn] = steepest_descent(f,g,[0,0]);
function [output] = ffun(x,y)
%FFUN 此处显示有关此函数的摘要
%   此处显示详细说明
output=100*(y-x^2)^2+(1-x)^2;
end
function [output] = gfun(x,y)
%GFUN 此处显示有关此函数的摘要
%   此处显示详细说明
output=[2*(-1+x+200*x^3-200*x*y),200*(y-x^2)];
end
function [xn, fn, gn] = steepest_descent(f, g, x0)
    ex=0.01; % precision for step size
    ef=0.01;
    eg=0.01;
    maxIterations=10000;
    debug=false;
    xk = x0;
    fk = f(xk(1),xk(2));
    gk = g(xk(1),xk(2));
    d = -gk;
    for i = 1:maxIterations
        % iteration
        % linear search
        [alpha, delta, xn, fn, gn]=search_alpha(f, g, xk, fk, gk, d);
        % convergence？
        if (norm(6)<=ex)&&(abs(fn-fk)<=ef)&&(norm(gn)<=eg)
            fprintf("Convergence is reached after %f iterations.\n",i);
            return ;
        end
        if debug
            fprintf('i=%f,alpha=%f,xk=%f,%f,d=%f,%f,delta=%f,%f\n',i,alpha,xn,d,delta);
        end
        xk = xn;
        fk = fn;
        gk = gn;
        d = -gk;
    end
    fprintf("WARN: %f iterations have been exceeded!\n", maxIterations);
end
%Inexact method, Wolfe condtions
% phi(alpha) <=phi(0) + epsilon alpha phi'(0)
% phi'(0）= d'*g
% phi'(alpha) >= eta phi'(0)
function [alpha, delta , xn, fn, gn] = search_alpha(f, g, xk, fk, gk, d)
    alpha0=1; epsilon=0.1; tau=0.5; eta=0.5; zeta=2.0;
    alpha = alpha0;
    phi0= d'*gk;
    delta = alpha .* d;
    xn = xk + delta;
    fn = f(xn(1),xn(2));
    gn = g(xn(1),xn(2));
    % Armijo condition
    while fn > fk+ epsilon*alpha*phi0
        alpha= tau*alpha;
        delta = alpha .* d;
        xn = xk + delta;
        fn = f(xn(1),xn(2));
        gn = g(xn(1),xn(2));
    end
    % Wolfe condition
    while d'*gn < eta*phi0
        alpha = zeta*alpha;
        delta = alpha .* d;
        xn = xk + delta;
        fn = f(xn(1),xn(2));
        gn = g(xn(1),xn(2));
    end
    return ;
end
