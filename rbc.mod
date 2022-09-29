/*
 This is the file for RBC in IDEC8012
*/


var C, K, h, z;
varexo u;

parameters sigma, chi, alpha, phi, delta, rho, beta, chi;

alpha = 0.36;
rho   = 0.95;
beta  = 0.99;
delta = 0.025;
sigma = 1;
phi   = 1;
chi   = 1;

model;

//(1)
chi*h^phi=C^(-sigma)*(1-alpha)*exp(z)*h^(-alpha)*K(-1)^alpha;

//(2)
C^(-sigma)=beta*C(+1)^(-sigma)*((1-delta)+alpha*exp(z(+1))*h(+1)^(1-alpha)*K^(alpha-1));

//(3) 
C+K-(1-delta)*K(-1)=exp(z)*h^(1-alpha)*K^alpha;

//(4)
z=rho*z(-1)+u;

end;

initval;
C = 0.80359242014163;
h = 0.29175631001732;
K = 11.08360443260358;
z = 0;
u = 0;
end;

shocks;
var u; stderr 0.01;
end;

stoch_simul(order=1);
