clear all

load lab2.dat -ascii
t = lab2(:,1);
u = lab2(:,2);
r = lab2(:,3);
p = lab2(:,4);
T = size(p,1);
ii = ones(T,1);
['Q.2'] 
plot(t,[u r p]);
legend('u', 'r', 'p');
[T,m]=size(p);

%  Questions 3 & 4
y = p;
X = [ii u];
% The next two lines do the same thing
b = X\y;
b = inv(X'*X)*X'*y;
e =  (y-X*b);
s2 = e'*e/(T-2);
tstats = b./sqrt(diag(s2*inv(X'*X)));
['Q.3'] 
['coeffs  t-stats'] 
[b tstats]
% The next line does the same thing as the two following lines
['Q.4'] 
r = corr(e(2:T),e(1:T-1));
[r  (sqrt(T)*r)]
[e(2:T-1)'*e(1:T-2)/(e'*e)  (sqrt(T)*r)]

%  Questions 5 to 7
% ARDL one lag each
y = p(2:T);
X = [ones(T-1,1) p(1:T-1) u(2:T) u(1:T-1)];
[T1, k1] = size(X);
b = X\y;
b = inv(X'*X)*X'*y;
e =  (y-X*b);
s2 = e'*e/(T1-k1);
tstats = b./sqrt(diag(s2*inv(X'*X)));
['Q.5'] 
['coeffs  t-stats'] 
[b tstats]
r = corr(e(2:T1-1),e(1:T1-2));
['r  t-stats for r']
[r  (sqrt(T1)*r);
e(2:T1-1)'*e(1:T1-2)/(e'*e)  (sqrt(T1)*r)]

% ARDL 2 lags each
% T = size(t,1); 
y = p(3:T);
X = [ones(T-2,1) p(2:T-1) p(1:T-2) u(3:T) u(2:T-1) u(1:T-2)];
[T2, k] = size(X);

['Q.6'] 
b = X\y;
b = inv(X'*X)*X'*y;
e =  (y-X*b);
s2 = e'*e/(T2-k);
tstats = b./sqrt(diag(s2*inv(X'*X)));

['coeffs  t-stats']
[b tstats]
r = corr(e(2:T2-1),e(1:T2-2));

['r    t-stats for r']
[r   (sqrt(T2)*r);
e(2:T2-1)'*e(1:T2-2)/(e'*e)  (sqrt(T2)*r)]

% We prefer the ARDL(2,2) based upon the evidence of
% no autocorrelation in the residual

['Q.7'] 
% Re-estimate the ARDL(2,2) in rearranged form
y = p(3:T) -  p(2:T-1);
X = [ones(T-2,1) p(2:T-1) (p(1:T-2) - p(2:T-1)) u(3:T) u(2:T-1) u(1:T-2)];
[T2, k] = size(X);
b = inv(X'*X)*X'*y;
e =  (y-X*b);
s2 = e'*e/(T2-k);
tstats = b./sqrt(diag(s2*inv(X'*X)));

['coeffs  t-stats']
[b tstats]

['Q.8'] 
[ 'The relevant t-stat is the second one (on the lagged level of prices).']
['This is not significant indicating a vertical LR Phillips curve']


['Q.9'] 
% Re-estimate the Restricted ARDL(2,2) in rearranged form
y = p(3:T) -  p(2:T-1);
X = [ones(T-2,1) (p(1:T-2) - p(2:T-1)) u(3:T) u(2:T-1) u(1:T-2)];
[T2, k] = size(X);
b = inv(X'*X)*X'*y;
e =  (y-X*b);
s2 = e'*e/(T2-k);
tstats = b./sqrt(diag(s2*inv(X'*X)));
['coeffs  t-stats']
[b tstats]
