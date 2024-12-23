%GDP_Indexed_2024_nonlinear_LOOP_scatter: plot of W0 for diff initial conditions. 
%Nonlinear simulation of the simple model (Algorithm in Supp Appendix).
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

clear; clc;

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

tic

n_sim = 10;  %No. of sims
T_sim = 510; 
drop = 10; T_sim2 = T_sim - drop;

%For discounted social welfare
omega = 0.98;  %Social discount factor
omega_vec = omega.^(1:T_sim2+1) / omega;

N_guess0 = 800;
y_init = 1;

n_loop = 15;  %no. of b values
chi_stack = linspace(0.9,1.045,n_loop);
r_lower = -0.8; r_upper = 0.8;

n_states = 5;  %No. of states
prob = ones(1,n_states); prob = prob / sum(prob);
n_shocks = 8;

%Shocks
sigma = sig_A;
Discretization_short
x1 = e_i;
%Grid of values
states = x1;
y_prime = ybar*exp(states);

n_v = 11; %Number of v values
v_stack = linspace(0,1,n_v); Lambda_v = NaN(n_sim,n_v);

for pp = 1:n_v

    v = v_stack(pp);

    pp

    %---------------------------
    %Stochastic simulations
    %---------------------------
    b_opt = zeros(n_sim,1); Utility_max = b_opt; Resid_max = b_opt; Resid_max2 = b_opt; Index_maxi = b_opt; Index_mini = b_opt; 

    for j=1:n_sim

        rng(j)
        y = ybar*exp(randn(T_sim,1)*sig_A);

        U_store = NaN(n_shocks,n_loop); Max_resid = U_store; Max_resid2 = U_store; Index_loc = U_store; Index_loc2 = U_store;
        Utility_vec = NaN(n_loop,1); Stack_b = Utility_vec; 

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
            
            for f=1:n_shocks
                
                rng(500+f);
                y(drop+1:T_sim) = ybar*exp(randn(T_sim2,1)*sig_A);

                GDP_Indexed_2024_nonlinear_SIM

                U_store(f,m) = (1-omega)*(omega_vec*Utility(drop:end) + omega^(-1)*Utility(drop-1));
                
                Max_resid(f,m) = max(Max_Resid);
                Max_resid2(f,m) = max(Resid_check);

                Index_loc(f,m) = max(Index);
                Index_loc2(f,m) = min(Index);

            end

            Utility_vec(m) = mean(U_store(:,m));
            %Utility_vec(m) = mean(Utility(drop+1:end));

        end

        [U_max,Index_max] = max(Utility_vec);
        
        b_opt(j) = Stack_b(Index_max);
        Utility_max(j) = U_max;

        Resid_max(j) = max(Max_resid(:,Index_max));
        Resid_max2(j) = max(Max_resid2(:,Index_max));

        Index_mini(j) = min(Index_loc2(:,Index_max));
        Index_maxi(j) = max(Index_loc(:,Index_max));

    end

if pp==1
    U_maxi = mean(Utility_max);
end

Lambda_v(:,pp) = 100*( (Utility_max/U_maxi).^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain 

end

max(Resid_max)
max(Resid_max2)

figure(3)
hold on, plot(v_stack,Lambda_v,'o', 'MarkerSize',5, 'MarkerEdgeColor','k'), xlabel('Indexation parameter (v)'), ylabel('% c.e. welfare gain' )
%%hold on, plot(v,Lambda_v(:,pp),'o','MarkerFaceColor','k', 'MarkerSize',5), xlabel('Indexation parameter (v)'), ylabel('% c.e. welfare gain' )

toc


