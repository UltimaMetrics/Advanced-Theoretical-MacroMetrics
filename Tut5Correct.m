% % Load data from the second lab session
load lab2.dat -ascii
t = lab2(:,1);
ym = lab2(:,2:4);

% Recall that y containts unemployment, interest rates and inflatoin
T = size(ym,1);

% Estimate a simple Phillips curve
y = ym(2:T,3) - ym(1:T-1,3);
X = [ones(T-1,1) ym(2:T,3)];
T = size(y,1);

y_1 = y(1:T-1);
y = y(2:T);
T = size(y,1);

X_1 = X(1:T,:);
X = X(2:T+1,:);

[T, k] = size(X);

iVprior = eye(k)*100;
nuh = 3;
eh = 1;

M = 10000;
keepp = zeros(M,k + 2);
bigrho = 0;
burnin = 2;

% Initialize the parameters before the Gibbs sampler
b = 1000*ones(k,1);
rho = -50;
h = 1000;

% Use the Gibbs sampler
for i = 1:M
    
    e = y - X*b;
    e_1=[0;e(1:T-1)];
    nu = e - rho*e_1;
    
    znu = randn(T + nuh - 2,1);
    h = znu'*znu/(nuh/eh + nu'*nu);
    
    Vr = inv(e_1'*e_1);
    rho = Vr*e_1'*e + chol(Vr/h)'*randn(1,1);
    if (i > burnin)*(rho > 0.65)
        bigrho = bigrho +1/(M-burnin);
    end
%     while (rho < 0.55)+(rho > 0.9)
%         rho = Vr*e_1'*e + chol(Vr/h)'*randn(1,1);
%     end
    
    yt = y - rho*y_1;
    Xt = X - rho*X_1;
    
    Vbpost = inv(Xt'*Xt*h + iVprior);
    bpost = Vbpost*h*Xt'*yt;
    b = bpost + chol(Vbpost)*randn(k,1);

    keepp(i,:) = [b' rho h];
end

% discard the burn-in sample
keepp = keepp(burnin+1:M,:);

% % Plot and histogram of each coefficient
plot(keepp(:,1:k));
histogram(keepp(:,1),50);
histogram(keepp(:,2),50);

% % Plot and histogram of rho
% plot(keepp(:,k+1));
% histogram(keepp(:,k+1),50);

% % Plot and histogram of the precision
% plot(keepp(:,4));
% histogram(keepp(:,4),50);
% title('h - precision');

% % Plot and histogram of nairu
% nairu = -keepp(:,1)./keepp(:,2);
% histogram(nairu,50);
% title('beta - coefficients');
% plot(nairu);
