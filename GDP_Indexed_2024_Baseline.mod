//Dynare code for simulating the simple model (Figs 1-3 in the paper). 
//Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

//-----------------------------------------
//1. Variable declaration and calibration
//-----------------------------------------

var c1, c2, y, r, r_i, tau, f_tau, b, utility, EV;
varexo e_a;

parameters alfa, betta, chi, eps, gama, n, v, ybar, gbar, bstar, c1star, c2star, utilitystar, taustar, phi, sig_A;

alfa = 0.3;  
betta = 0.9;  
chi = 1;
gama = 5;  
eps = 0.5; 
n = 0.4; 
v = 1;
ybar = 1;
gbar = 0.15;  
phi = 0.5;
sig_A = 0.05;  

//----------------------------------
//1. Find steady state init vals
//----------------------------------
GDP_indexed_steady_state

bstar = b_root;
rstar = chi*(1+n);
taustar = ( gbar + (chi-1)*bstar ) / ybar;
c1star = (1-alfa-taustar)*ybar - phi*taustar^2 - bstar;
c2star = alfa*(1+n)*ybar + chi*(1+n)*bstar;
utilitystar = (1/(1-gama))*(c1star^eps + betta*c2star^eps)^((1-gama)/eps);

//--------------------------------
//2. Model
//--------------------------------

model;

//Output
y = ybar*exp(e_a);

//Consumption when young 
c1 = (1-alfa-tau)*y - phi*tau^2 - b;

//Consumption when old
c2 = alfa*(1+n)*y + r_i*b(-1);

//Bond supply
b = bstar;

//Determination of taxes 
tau = ( gbar + r_i*b(-1)/(1+n) - b ) / y;

//Consumption Euler equation 
1 = betta*r*(y(+1)/ybar)^v*(c1/c2(+1))^(1-eps)*( c2(+1)/( EV^(1/(1-gama)) ) )^(1-gama-eps);

//Real interest rate on bonds 
r_i = r(-1)*(y/ybar)^v;

//Tax burden
f_tau = phi*tau^2;

//Lifetime utility 
utility = 1/(1-gama)*( c1^eps + betta*( EV )^(eps/(1-gama)) )^((1-gama)/eps);

//Expectation term
EV = c2(+1)^(1-gama);

end;

//----------------------------------------
//3. Initial values and shock calibration
//----------------------------------------

initval;
c1 = c1star;
c2 = c2star;
b = bstar;
r = rstar;
r_i = rstar;
y = ybar;
tau = taustar;
utility = utilitystar;
EV = c2^(1-gama);

end;

steady;

shocks;
var e_a; stderr sig_A;
end;

//---------------------------
//Find optimal bond supply
//---------------------------

n_loop = 100;
chi_stack = linspace(0.9,1.045,n_loop);
Stack_utility = NaN(n_loop,1); Stack_tau = Stack_utility; Stack_b = Stack_utility; 

for i=1:n_loop

    //----------------------------------------------
    //Find determinstic SS and use in Dynare solver
    //----------------------------------------------
    chi = chi_stack(i);
    GDP_indexed_steady_state

    bstar = b_root;
    rstar = chi*(1+n);
    taustar = ( gbar + (chi-1)*bstar ) / ybar;
    c1star = (1-alfa-taustar)*ybar - phi*taustar^2 - bstar;
    c2star = alfa*(1+n)*ybar + chi*(1+n)*bstar;
    utilitystar = (1/(1-gama))*(c1star^eps + betta*c2star^eps)^((1-gama)/eps);

    //------------------------
    //Stochastic simulations
    //------------------------
    steady;   //steady(tolf=1e-10,tolx=1e-10);

    stoch_simul(order=2, drop=0, periods=0, irf=0, noprint);
    
    Loop_record  //Record results in each loop

end

[U_ss_max,Index_ss] = max(U_ss);
[U_stoch_max,Index_stoch] = max(Stack_utility);

GDP_Indexed_Baseline_Plotter

//---------------------
//Optimal bond policy
//---------------------
b_ss_opt = b_ss(Index_ss)
r_ss_opt = r_ss(Index_ss)
chi_ss_opt = chi_stack(Index_ss)

b_opt = Stack_b(Index_stoch);
r_opt = Stack_r(Index_stoch);
chi_opt = chi_stack(Index_stoch);

//------------------------
//Welfare at the optimum
//------------------------
Utility_max = Stack_utility(Index_stoch);  







