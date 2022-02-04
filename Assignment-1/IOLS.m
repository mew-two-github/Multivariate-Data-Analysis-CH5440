function [alpha, beta, uhat, yhat, s] = IOLS(u,y)
    N = length(u);
    ybar = mean(y);
    ubar = mean(u);
    syy = var(y,1);
    suu = var(u,1);
    syu = 1/N*sum((y-ybar).*(u-ubar));
    yhat = y;
    alpha = syy/syu;
    beta = ybar - alpha*ubar;
    uhat = (yhat - beta)/alpha;
    s = struct('alpha',alpha,'beta',beta,'ybar',ybar,'ubar',ubar,'syy',syy,'syu',syu,'suu',suu);
end