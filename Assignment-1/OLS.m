function [alpha, beta, uhat, yhat, s] = OLS(u,y)
    N = length(u);
    ybar = mean(y);
    ubar = mean(u);
    syy = var(y,1);
    suu = var(u,1);
    syu = 1/N*sum((y-ybar).*(u-ubar));
    uhat = u;
    alpha = syu/suu;
    beta = ybar - alpha*ubar;
    yhat = alpha*uhat + beta;
    s = struct('alpha',alpha,'beta',beta,'ybar',ybar,'ubar',ubar,'syy',syy,'syu',syu,'suu',suu);
end