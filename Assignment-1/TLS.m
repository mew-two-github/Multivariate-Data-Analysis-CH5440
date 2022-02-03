function [alpha, beta, uhat, yhat, s] = TLS(u,y)
    N = length(u);
    ybar = mean(y);
    ubar = mean(u);
    syy = var(y,1);
    suu = var(u,1);
    syu = 1/N*sum((y-ybar).*(u-ubar));
    alpha = (syy - suu + sqrt((syy-suu)^2 + 4*syu^2))/2/syu;
    beta = ybar - alpha*ubar;
    uhat = (alpha*(y-beta)+u)/(alpha^2+1);
    yhat = alpha*uhat + beta;
    s = struct('alpha',alpha,'beta',beta,'ybar',ybar,'ubar',ubar,'syy',syy,'syu',syu,'suu',suu);
end