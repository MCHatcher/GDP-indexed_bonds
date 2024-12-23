//Dynare code for simulating the parameter sensitivity analysis of the simple model. 
//Written by Michael Hatcher (m.c.hatcher@soton.ac.uk) and Juyi Lyu (julielyu2014@gmail.com). Any errors are ours.

close all;  
clc;


//-----------------------------------------
//1. Variable declaration and calibration
//-----------------------------------------



var c1, c2, y, r, r_i, tau, f_tau, b, utility, EV;
varexo e_a;




parameters alfa, betta, chi, eps, gama, n, v, ybar, gbar, bstar, c1star, c2star, utilitystar, taustar, phi, sig_A;


alfa = 0.3;     //alfa = 0.275,0.325; baseline = 0.3  
betta = 0.9;  //betta = 0.825,0.975;  baseline = 0.9  
chi = 1;
gama = 5;     //gama = 4,6; baseline = 5 


eps = 0.5;   //eps = 0.3,0.7 ;  baseline = 0.5
n = 0.4;     //n = 0.3,0.5 ; baseline = 0.4
v = 1;
ybar = 1;
gbar = 0.15;  
phi = 0.5;   //phi = 0.25,0.75 ; baseline = 0.5  
sig_A = 0.05;  //sig_A = 0.025,0.075 ; baseline = 0.05  





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


n_loop = 500;
//chi_stack = linspace(0.885,0.985,n_loop);  //Fig 2
chi_stack = linspace(0.8,0.995,n_loop);   //Sensitivity analysis
Stack_utility = NaN(n_loop,1); 

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



//------------------------------
//        Stochastic simulations
//-----------------------------
	
	
    steady(tolf=8e-10);

    stoch_simul(order=2, drop=0, periods=0, irf=0, noprint);

    Stack_utility_inner(i) = oo_.mean(9);
    Stack_b(i) = oo_.mean(8);

end

[U_max,Index_max] = max(Stack_utility_inner);

chi = chi_stack(Index_max);
chi_max = chi;
b_opt = Stack_b(Index_max);
bstar = b_opt;




//---------------------------
//Outer loops for v
//---------------------------



n_loop_v = 250;
v_stack = linspace(0,1,n_loop_v);
Lambda_v = NaN(n_loop_v,1); Stack_c1 = Lambda_v; Stack_c2 = Lambda_v; Stack_tau = Lambda_v; var_tau = Lambda_v; 

for m=1:n_loop_v

    v = v_stack(m);

stoch_simul(order=2, drop=0, periods=0, irf=0, noprint);

Stack_utility_outer(m) = oo_.mean(9);
Lambda_v(m) = 100*( (Stack_utility_outer(m)/Stack_utility_outer(1))^(1/(1-gama)) - 1);  //Consumption equiv. welfare gain 

Stack_c1(m) = oo_.mean(1); Stack_c2(m) = oo_.mean(2);
var_c1(m) = oo_.var(1,1); var_c2(m) = oo_.var(2,2);
Stack_tau(m) = oo_.mean(6); var_tau(m) = oo_.var(6,6);
Stack_ftau(m) = oo_.mean(7);

Stack_b(m) = oo_.mean(8); 
Stack_r(m) = oo_.mean(4); var_r(m) = oo_.var(4,4);
Stack_ri(m) = oo_.mean(5); var_ri(m) = oo_.var(5,5);

Stack_utility(m) = oo_.mean(9);
U_mean(m) = (1/(1-gama))*( Stack_c1(m)^eps + betta*Stack_c2(m)^eps  )^((1-gama)/eps);
Lambda_mean(m) = 100*( (U_mean(m)/U_mean(1))^(1/(1-gama)) - 1);
Lambda_var(m) = Lambda_v(m) - Lambda_mean(m);

Lambda_tau(m) = 100*(sqrt(Stack_ftau(m)/Stack_ftau(1)) -1);
ftau_mean(m) = 100*(sqrt(Stack_tau(m)^2/Stack_tau(1)^2) -1);
ftau_var(m) = Lambda_tau(m) - ftau_mean(m);

b_ss(m) = oo_.steady_state(8);
r_ss(m) = oo_.steady_state(4);
ri_ss(m) = oo_.steady_state(5);
U_ss(m) = oo_.steady_state(9);
Lambda_ss(m) = 100*( (U_ss(m)/U_ss(1))^(1/(1-gama))  - 1);   //Consumption equiv. welfare gain 

end

[U_max_max,Index_max_max] = max(Lambda_v);
vstar = v_stack(Index_max_max)

[min_vri,Index_min] = min(var_ri);
vmin = v_stack(Index_min)

[min_vf,Index_minf] = min(Lambda_tau);
vminf = v_stack(Index_minf)


set(0, 'defaultTextInterpreter', 'latex'); 

h = figure; // Create a new figure

// Open the previously saved figure  (Uncomment below lines as necessary)

//h = openfig('first_run_plot.fig', 'reuse'); // Reuse the existing figure

//h = openfig('second_run_plot.fig', 'reuse'); // Reuse the existing figure

//figure(h); % Make the figure active


subplot(1, 2, 1);
plot(v_stack, Lambda_tau, 'k', 'LineWidth', 1);
xlabel('Indexation parameter ($v$)', 'Interpreter', 'latex');
ylabel('% tax equiv.', 'Interpreter', 'latex');
title('Tax burden vs indexation', 'Interpreter', 'latex');
yline(0, '--k');
legend('Baseline: \sigma_A = 0.05','Box', 'off'); 


subplot(1, 2, 2);
plot(v_stack, Lambda_v, 'k', 'LineWidth', 1);
xlabel('Indexation parameter ($v$)', 'Interpreter', 'latex');
ylabel('% c.e. welfare gain', 'Interpreter', 'latex');
title('Welfare gain of indexation', 'Interpreter', 'latex');
yline(0, '--k');
legend('Baseline: \sigma_A = 0.05','Box', 'off'); 


//Save the first plot

savefig(h, 'first_run_plot.fig');   //Comment out as necessary

//Save the updated plot
//savefig(h, 'second_run_plot.fig');  //Uncomment as necessary

//savefig(h, 'third_run_plot.fig');   //Uncomment as necessary

