X = [1 1; 1 2; 1 3; 1 4; 1 5];

k = 2;  T = 500;
X = [ones(T,1), (1:1:T)'];
b = [ 1; 2];
b = (1:1:k)';
% X = [ones(T,1) rand(T,k-1)*5];
% k = 5;  T = 200;
% b = magic(k);
b = b(:,1);
xb = X*b;
m = 1000;
tstats = zeros(m,k);
preject = 0;

for i = 1:m
%   Draw the vector of errors  
    e = randn(T,1)*sqrt(0.2);
%   Create y 
    y = xb + e;
%   Compute the ols estimate of beta 
    bols = inv(X'*X)*X'*y;
%   The following three lines do the same thing  
%   They compute the estimate of the variance 
%     s2ols = (y - x*bols)'(y - x*bols)/(T-k);
%     s2ols = (y'*y - y'*x*inv(X'*X)*X'*y)/(T-k);
    s2ols = (y'*y - y'*X*bols)/(T-k);
    bcovols = s2ols*inv(X'*X);
%   The dot '.' in the next line tells MATLAB to do 
%   element by element division of the vectors
%   This creates the t-statistics for both coefficients
%   The null hypotheis is that beta equals the true value
    tstatsi=(bols-b)./sqrt(diag(bcovols));
    tstats(i,:)=tstatsi';
    preject = preject + (abs(tstatsi(2)) > 3.182446)/m;
end


%   Plot a histogram of all draws of the t-stats
% hist(tstats(:,1),50);
hist(tstats(:,2),50);
preject

%   Put all the data into one matrix
% yx = [y X];
%   Save the data in this matrix to a file called lab1.dat
% save -ascii lab1.dat yx
%   Load the data in the file lab1.dat into MATLAB
% load lab1.dat -ascii
