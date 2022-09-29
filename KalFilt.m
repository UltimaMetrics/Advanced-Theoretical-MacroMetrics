function [bhat,Vtt] = KalFilt(y1, X, Ht,   Qt,   m,p,t,cRt,cTt)

% run the Kalman filter and then return what -- mean of w 
% proc (2) = KalFilt(y1,X,Ht,Qt,m,p,t,capRt,capTt);
% y1 - (tpx1) matrix of dependent variables
% X - (tpxm)) matrix of explanatory variables
% Ht - (pxp) matrix of all measurement equation variances
% Qt - (mxm) matrix of all state equation variances
% p - number of dependent variables in the matrix y1 (i.e., the number 
%     of columns of y1
% m - number of explanatory variables. (i.e., the number of columns of X)
% t - number of observations 
% cRt - m x r matrix in state equation next to errors
% cTt - m x m matrix in state equation next to beta

% kalman filter code 
Kkeep=zeros(t*m,p);
Lkeep=zeros(t*m,m);

Fkeep=zeros(t*p,p);
Pkeep=zeros(t*m,m);
b=zeros(m,t+1);

mny = y1(1,:)';
b(:,1) =   2*chol(Qt(1:m,:))*randn(m,1);
v=zeros(t,p);
Pt=zeros(m,m);

for i=1:t;
    ht = Ht((i-1)*p+1:i*p,:);
    qt = Qt((i-1)*m+1:i*m,:);
    xt = X((i-1)*p+1:i*p,:);
%     cRt = capRt((i-1)*m+1:i*m,:);
%     [size(v(i,:)) size(y1(i,:)') size(xt*b(:,i))]
    v(i,:) = y1(i,:)' - xt*b(:,i);
    Ft= xt*Pt*xt' + ht;
    Ftinv=inv(Ft);
%   Save the Ft for use in the Kalman smoother
    Fkeep((i-1)*p+1:i*p,:)=Ftinv;
    Kt = cTt*Pt*xt'*Ftinv;
    Kkeep((i-1)*m+1:i*m,:)=Kt;
    Ltt = cTt*eye(m) - Kt*xt;
%   Save the Lt for use in the Kalman smoother
    Lkeep((i-1)*m+1:i*m,:)=Ltt;

    b(:,i+1) = cTt*b(:,i) + Kt*v(i,:);
    Pt = cTt*Pt*Ltt' + cRt*qt*cRt';
    Pkeep((i-1)*m+1:i*m,:) = Pt;
end

% backward recursion to evaluate rt and, thus, whatt
rt=zeros(m,t+1);
Nt =0;
bhat = zeros(m,t);
Vtt = zeros(m,m*t);
for i = t:-1:1
%     cRt = capRt[(i-1)*m+1:i*m,:];
    xt = X((i-1)*p+1:i*p,:);
    Ft = Fkeep((i-1)*p+1:i*p,:)';
    vt = v(i,:);
    Ltt = Lkeep((i-1)*m+1:i*m,:);
    Pt = Pkeep((i-1)*m+1:i*m,:);

    rt(:,i) = xt'*inv(Ft)*v(i,:) + Ltt'*rt(:,i+1);
    Nt = xt'*inv(Ft)*xt + Ltt'*Nt*Ltt;
    
    Vtt(:,(i-1)*m+1:i*m,:) = Pt - Pt*Nt*Pt;
    bhat(:,i+1) = b(:,i+1) + Pt*rt(:,i);
    
end
bhat = bhat(:,2:end);
end
