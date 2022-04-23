function [vecs,scores] = MLPCA(X,sigma,p)
    % Solving objective - eqn 14
    Z = X./sigma;
    [~,~,v] = svd(Z);
    vecs = v(:,1:p);
    % Projection - eqn 12
    scores = zeros(size(Z,1),p);
    for i = 1:size(Z,1)
        covZ = diag(sigma(i,:).^2);
        scores(i,:) = Z(i,:)*inv(covZ)*vecs*inv(vecs'*inv(covZ)*vecs);
    end
end