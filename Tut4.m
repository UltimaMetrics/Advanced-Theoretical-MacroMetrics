% Simulate a sample of data
% T = 50;
% btrue = 3;
% htrue = .5;
% X = ones(T,1);
% y = X*btrue + sqrt(1/htrue)*randn(T,1);

% % Load data from the second lab session load lab2.dat -ascii t =
load lab2.dat -ascii
t = lab2(:,1);
u = lab2(:,2);
r = lab2(:,3);
p = lab2(:,4);
% lab2(:,1); 
% ym = lab2(:,2:4); 
% Recall that y containts unemployment,
%  ... interest rates and inflatoin 
T = size(t,1); 

% Estimate a simple Phillips curve 
% y = p; 
% T = size(y,1);
% X = [ones(T,1) u];

% Estimate a Phillips curve using an ARDL(2,2) 
y = p(3:T) - p(2:T-1); 
% X = [ones(T-2,1) p(2:T-1) (p(2:T-1)-p(1:T-2)) u(3:T) u(2:T-1)];
X = [ones(T-2,1) (p(2:T-1)-p(1:T-2)) u(3:T) u(2:T-1)];

[T, k] = size(X);
XX = X'*X;
Xy = X'*y;
yy = y'*y;

iVprior = eye(k)*.10;
nup = 3;
nups2 = 3;

M = 10000;

keepp = zeros(M,k + 1);
NAIRU = zeros(M, 1);
b = 1000*ones(k,1);
% Use the Gibbs sampler
for i = 1:M
    e = y - X*b;
    znu = randn(T + nup - 2,1);
    h = znu'*znu/(nups2 + e'*e);
    
    Vpost = inv(XX*h + iVprior);
    bpost = Vpost*h*Xy;
    b = bpost + chol(Vpost)*randn(k,1);

    keepp(i,:) = [b' h];
    NAIRU(i,:) = (b(3)+b(4));
end

% discard the burn-in sample
NAIRU = NAIRU(M/2:M,:);
keepp = keepp(M/2:M,:);
plot(NAIRU);
title('NAIRU');
% plot(keepp(:,1:k));
% histogram(keepp(:,2),50);
% title('beta - coefficients');

% plot(keepp(:,2));
% histogram(keepp(:,2),50);
% title('h - precision');
