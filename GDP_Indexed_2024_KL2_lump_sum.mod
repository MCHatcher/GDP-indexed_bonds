//Dynare code for simulating the extended model (Figs 4 in the paper). 
//Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

//-----------------------------------------
//1. Variable declaration and calibration
//-----------------------------------------

var c1, c2, k, b, y, l, r, r_i, rk, w, tau, utility, EV, c1s, c2s, ks, ys, ls, rs, r_is, rks, ws, tau_s, utility_s, EVs, x_TR, x_c1, x_c1s;
varexo e_a;

parameters alfa, betta, chi, chi1, v, eps, gama, tau_p, thetta, n, psi, gbar, bstar, c1star, c2star, utilitystar, taustar, rstar, 
sig_A; 

alfa = 0.3;  
betta = 0.9;  
gama = 5;  
eps = 0.5;  
n = 0.4;  
chi = 1+n;
chi1 = chi;
gbar = 0.05;  
thetta = 0.45; 
psi = 1;   
v = 1;
tau_p = 0;  
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

//Steady_state_LS2_insert

//b1star = b1_ss; k1star = k1_ss; 
//r1star = chi; s1_star = s1_ss;
//tau1star = tau1_ss; l1star = l1_ss;
//c11star = c11_ss; c21star = c21_ss;
//utility1star = utility1_ss;

//--------------------------------
//2. Model
//--------------------------------

model;

//------------
//Main model
//------------

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
utility = (1/(1-gama))*( ( c1^thetta*(1-l)^(1-thetta) )^eps + betta*( EV )^(eps/(1-gama)) )^((1-gama)/eps);

//Expectation term
EV =  c2(+1)^(thetta*(1-gama));

//Composite consumption 
x_c1 = c1^thetta*(1-l)^(1-thetta);

//---------------------------
//Model with lump-sum taxes
//---------------------------

//Output
ys = exp(e_a)*(ks(-1)/(1+n))^alfa*ls^(1-alfa);

//Consumption when young 
c1s = ws*ls - ks - b - tau_s - tau_p*ws*ls;

//Consumption when old
c2s = rks*ks(-1) + r_is*b(-1) - psi*tau_s + (1+n)*tau_p*ws*ls;

//Determination of taxes 
(1+psi/(1+n))*tau_s =  gbar + r_is*b(-1)/(1+n) - b;

//Consumption Euler equation (bonds)
1 = betta*rs*(ys(+1)/STEADY_STATE(ys))^v*(c2s(+1)/c1s)^(eps*thetta-1)*( 1/(1-ls) )^( (1-thetta)*eps )*( c2s(+1)^thetta / EVs^(1/(1-gama)) )^(1-gama-eps);

//Consumption Euler equation (capital) 
1 = betta*rks(+1)*(c2s(+1)/c1s)^(eps*thetta-1)*( 1/(1-ls) )^( (1-thetta)*eps )*( c2s(+1)^thetta / EVs^(1/(1-gama)) )^(1-gama-eps);

//Labour supply
thetta*(1-ls) = (1-thetta)*c1s/ws;

//Real interest rate on bonds 
r_is = rs(-1)*(ys/STEADY_STATE(ys))^v;

//Return on capital 
rks = alfa*ys/(ks(-1)/(1+n));

//Wage
ws = (1-alfa)*ys/ls;

//Lifetime utility 
utility_s = (1/(1-gama))*( (c1s^thetta*(1-ls)^(1-thetta) )^eps + betta*( EVs )^(eps/(1-gama)) )^((1-gama)/eps);

//Expectation term
EVs =  c2s(+1)^(thetta*(1-gama));

//Composite consumption 
x_c1s = c1s^thetta*(1-ls)^(1-thetta);

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
x_c1 = c1^thetta*(1-l)^(1-thetta);
//-------------
c1s = c1star;
c2s = c2star;
ks = kstar;
rs = chi;
r_is = rs;
ls = lstar;
ys = (kstar/(1+n))^alfa*lstar^(1-alfa);
ws = (1-alfa)*ys/ls;
tau_s = taustar*(w*lstar + psi*rk*kstar/(1+n));
utility_s = utilitystar;
EVs = c2s^(thetta*(1-gama));
rks = rs;
x_c1s = c1s^thetta*(1-ls)^(1-thetta);
end;

steady;

shocks;
var e_a; stderr sig_A;
end;

//---------------------------
//Find optimal bond supply
//---------------------------

n_loop = 50;
chi_stack = linspace(0.876*(1+n),1.175*(1+n),n_loop);  //For LS taxes
//n_loop = 50; chi_stack = linspace(0.865*(1+n),0.955*(1+n),n_loop); 

Stack_utility = NaN(n_loop,1); Stack_tau = Stack_utility; Stack_b = Stack_utility; Lambda = Stack_utility;
Stack_c1 = Stack_utility; Stack_c2 = Stack_utility; Stack_r = Stack_utility; 
Stack_bss = Stack_utility; Stack_rss = Stack_utility;  Stack_Uss = Stack_utility; Lambda_ss = Stack_utility;
var_c1 = Stack_utility; var_c2 = Stack_utility; var_r = Stack_utility; var_tau = Stack_utility;  

for j=1:n_loop

    chi = chi_stack(j);
    chi1 = chi; 
    
    Steady_state_KL2_insert

    bstar = b_ss; kstar = k_ss; 
    rstar = chi; s_star = s_ss;
    taustar = tau_ss;
    c1star = c1_ss; c2star = c2_ss;
    lstar = l_ss;
    utilitystar = utility_ss;

    //Steady_state_LS2_insert

    //b1star = b1_ss; k1star = k1_ss; 
    //r1star = chi; s1_star = s1_ss;
    //tau1star = tau1_ss; l1star = l1_ss;
    //c11star = c11_ss; c21star = c21_ss;
    //utility1star = utility1_ss;

    steady;   //steady(tolf=1e-10,tolx=1e-10);
    stoch_simul(order=2, drop=0, periods=0, irf=0, noprint);

    Loop_record_KL2

end

[U_ss_max,Index_ss] = max(Stack_Uss);
[U_stoch_max,Index_stoch] = max(Stack_utility)

b_opt = Stack_bss(Index_ss);
k_opt = Stack_kss(Index_ss);
y_opt = Stack_yss(Index_ss);
r_opt = Stack_rss(Index_ss);
xc1_opt = Stack_comp_ss(Index_ss);

GDP_Indexed_KL2_Plotter
 