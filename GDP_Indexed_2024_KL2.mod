//Dynare code for simulating the extended model w/out lump sum taxes.
//Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

//-----------------------------------------
//1. Variable declaration and calibration
//-----------------------------------------

var c1, c2, k, b, y, l, r_i, r, rk, w, tau, utility, EV, x_TR;
varexo e_a;

parameters alfa, betta, chi, eps, gama, thetta, tau_p, v, n, psi, gbar, bstar, c1star, c2star, utilitystar, taustar, rstar, sig_A;

alfa = 0.3;
betta = 0.9;
gama = 5;
eps = 0.5;
n = 0.4;
chi = 1+n;
gbar = 0.05;
thetta = 0.45;
psi = 1;
tau_p = 0;
v = 1;
sig_A = 0.05;

//----------------------------------
//1. Find steady state init vals
//----------------------------------
Steady_state_KL2_insert

bstar = b_ss; kstar = k_ss; 
rstar = chi;
s_star = s_ss;
taustar = tau_ss;
c1star = c1_ss; c2star = c2_ss;
lstar = l_ss;
utilitystar = utility_ss;
utility2star = utility2_ss; 

//--------------------------------
//2. Model
//--------------------------------

model;

//Output
y = exp(e_a)*(k(-1)/(1+n))^alfa*l^(1-alfa);

//Consumption when young 
c1 = (1-tau-tau_p)*w*l - k - b;

//Consumption when old 
c2 = (1-psi*tau)*rk*k(-1) + r_i*b(-1) + x_TR;

//Pension transfer
x_TR = tau_p*(1+n)*w*l;

//Bond supply
b = bstar;

//Determination of taxes 
tau = ( gbar + r_i*b(-1)/(1+n) - b ) / ( w*l + psi*rk*(k(-1)/(1+n)) );

//Consumption Euler equation (bonds)
1 = betta*r*(y(+1)/STEADY_STATE(y))^v*(c2(+1)/c1)^(eps*thetta-1)*( 1/(1-l) )^( (1-thetta)*eps )*( c2(+1)^thetta / EV^(1/(1-gama)) )^(1-gama-eps);

//Consumption Euler equation (capital) 
1 = betta*(1-psi*tau(+1))*rk(+1)*(c2(+1)/c1)^(eps*thetta-1)*( 1/(1-l) )^( (1-thetta)*eps )*( c2(+1)^thetta / EV^(1/(1-gama)) )^(1-gama-eps);

//Labour supply
thetta*(1-l) = (1-thetta)*c1/( (1-tau-tau_p)*w );

//Real interest rate on bonds 
r_i = r(-1)*(y/STEADY_STATE(y))^v;

//Return on capital 
rk = alfa*y/(k(-1)/(1+n));

//Wage
w = (1-alfa)*y/l;

//Lifetime utility 
utility = (1/(1-gama))*( (c1^thetta*(1-l)^(1-thetta) )^eps + betta*( EV )^(eps/(1-gama)) )^((1-gama)/eps);

//Expectation term
EV =  c2(+1)^(thetta*(1-gama));

end;

//----------------------------------------
//3. Initial values and shock calibration
//----------------------------------------

initval;
c1 = c1star;
c2 = c2star;
k = kstar;
b = bstar;
r = chi;
r_i = r;
l = lstar;
y = (kstar/(1+n))^alfa*lstar^(1-alfa);
w = (1-alfa)*y/l;
tau = taustar;
utility = utilitystar;
EV = c2star^(thetta*(1-gama));
rk = r/(1-psi*tau);
x_TR = tau_p*(1+n)*w*l;
end;

steady;

shocks;
var e_a; stderr sig_A;
end;

stoch_simul(order=2, drop=0, periods=0, irf=0);
 