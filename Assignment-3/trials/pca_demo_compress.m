load iris
%  Arrange the data in a matrix of samples x variables
iris_data = [iris(:,1:4);iris(:,5:8);iris(:,9:12)];
%  Plot attribute 1 vs atribute 2 to see whether the samples form disntinct
%  clusters
plot(iris_data(1:50,1),iris_data(1:50,2),'b+')
hold on
plot(iris_data(51:100,1),iris_data(51:100,2),'g*')
plot(iris_data(101:150,1),iris_data(101:150,2),'rs')
keyboard
%
%  Exercise: Use subplot to try other combinations of variables ie generate the
%  biplots
%
%  Use svd to compress the data
[u s v] = svd(iris_data);
% Compute eigenvalues values
lambda = diag(s).^2;
% Variance explained by firts two PCS
sum(lambda(1:2))/sum(lambda)
%  Choose first two PCs since it explains more than 95% variance.  Obtain the scores matrix 
scores = iris_data*v(:,1:2);
% Examines the scores plot to check whether the samples corresponding to
% different samples are clearly separated
hold off
plot(scores(1:50,1),scores(1:50,2),'b+')
hold on
plot(scores(51:100,1),scores(51:100,2),'g*')
plot(scores(101:150,1),scores(101:150,2),'rs')
