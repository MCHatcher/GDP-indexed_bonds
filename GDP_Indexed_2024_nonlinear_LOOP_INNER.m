%GDP_Indexed_2024_nonlinear_LOOP_INNER
%Nonlinear simulation of the simple model (Algorithm in Supp Appendix).
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

%clear

alfa = 0.3;  
betta = 0.9;  
chi = 1;
gama = 5;  
eps = 0.5;  
n = 0.4;
ybar = 1;
gbar = 0.15;  
phi = 0.5; 
sig_A = 0.05; 
%v = 1;

n_sim = 100;  %No. of sims
T_sim = 210; drop = 10;
T_sim2 = T_sim - drop;

N_guess0 = 850;
y_init = 1;

n_loop = 16;  %No. of b values
chi_stack = linspace(0.9,1.045,n_loop);
r_lower = -0.85; r_upper = 0.85;

n_states = 5;  %No. of states
prob = ones(1,n_states); prob = prob / sum(prob);

%Shocks
sigma = sig_A;
Discretization_short
x1 = e_i;
%Grid of values
states = x1;
y_prime = ybar*exp(states);

Stack_utility = NaN(n_loop,1); Stack_tau = Stack_utility; Stack_ftau = Stack_utility; Stack_b = Stack_utility; 
Stack_c1 = Stack_utility; Stack_c2 = Stack_utility; var_tau = Stack_b; var_ftau = Stack_b; U = Stack_b; 
var_c1 = Stack_utility; var_c2 = Stack_utility; U_mean = var_c1; U_var = var_c1; ftau_mean = var_c1;
Stack_ri = Stack_utility; var_ri = Stack_utility;

Resid_max = U; Index_mini = Resid_max; Index_maxi = Resid_max; Resid_max2 = Resid_max;  

for m=1:n_loop

    chi = chi_stack(m);
    GDP_indexed_steady_state

    bstar = b_root;
    rstar = chi*(1+n);
    Stack_b(m) = bstar;

    r_init = rstar; 

    %Guesses for interest rate
    r_guess_stack0 = r_init + linspace(r_lower,r_upper,N_guess0);
    %r_guess_stack = r_guess_stack(randperm(N_guess));

    %---------------------------
    %Stochastic simulations
    %---------------------------
    Dum = zeros(n_sim,1); Max_resid = Dum; Max_resid2 = Dum; Index_loc = Dum; Index_loc2 = Dum;  
    c1_vec = Dum; c2_vec = Dum; tau_vec = Dum; ftau_vec = Dum; ri_vec = Dum; Utility_vec = Dum; 
    c1_vec2 = Dum; c2_vec2 = Dum; tau_vec2 = Dum; ftau_vec2 = Dum; ri_vec2 = Dum;

    for j=1:n_sim

        rng(500+j)
        y = ybar*exp(randn(T_sim,1)*sig_A);

        GDP_Indexed_2024_nonlinear_SIM

        Max_resid(j) = max(Max_Resid);
        Max_resid2(j) = max(Resid_check);

        Index_loc(j) = max(Index);
        Index_loc2(j) = min(Index);
    
        %---------------------------------
        %Store and stack key variables
        %---------------------------------

        %For means
        c1_vec(j) = sum(c1(drop+1:end));
        c2_vec(j) = sum(c2(drop+1:end));
        tau_vec(j) = sum(tau(drop+1:end));
        ftau_vec(j) = sum(ftau(drop+1:end));
        ri_vec(j) = sum(r_i(drop+1:end));
        Utility_vec(j) = sum(Utility(drop+1:end));

        %For variances
        c1_vec2(j) = sum(c1(drop+1:end).^2);
        c2_vec2(j) = sum(c2(drop+1:end).^2);
        tau_vec2(j) = sum(tau(drop+1:end).^2);
        ftau_vec2(j) = sum(ftau(drop+1:end).^2);
        ri_vec2(j) = sum(r_i(drop+1:end).^2);

    end

    Resid_max(m) = max(Max_resid);
    Resid_max2(m) = max(Max_resid2);

    Index_mini(m) = min(Index_loc2);
    Index_maxi(m) = max(Index_loc);

    Stack_utility(m) = sum(Utility_vec)/(n_sim*T_sim2);
    Stack_c1(m) = sum(c1_vec)/(n_sim*T_sim2);
    Stack_c2(m) = sum(c2_vec)/(n_sim*T_sim2);
    Stack_tau(m) = sum(tau_vec)/(n_sim*T_sim2);
    Stack_ftau(m) = sum(ftau_vec)/(n_sim*T_sim2);
    Stack_ri(m) = sum(ri_vec)/(n_sim*T_sim2);
    
    var_c1(m) = sum(c1_vec2)/(n_sim*T_sim2) - (Stack_c1(m))^2;
    var_c2(m) = sum(c2_vec2)/(n_sim*T_sim2) - (Stack_c2(m))^2;
    var_tau(m) = sum(tau_vec2)/(n_sim*T_sim2) - (Stack_tau(m))^2;
    var_ftau(m) = sum(ftau_vec2)/(n_sim*T_sim2) - (Stack_ftau(m))^2;
    var_ri(m) = sum(ri_vec2)/(n_sim*T_sim2) - (Stack_ri(m))^2;

end






