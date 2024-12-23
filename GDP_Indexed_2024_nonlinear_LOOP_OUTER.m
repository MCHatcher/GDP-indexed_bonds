%GDP_Indexed_2024_nonlinear_LOOP_OUTER
%Nonlinear simulation of the simple model (Algorithm in Supp Appendix). 
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

clear; clc;

n_v = 20; %Number of v values

Stack_utility_outer = NaN(n_v,1); Lambda_v = Stack_utility_outer; Stack_b_outer = Stack_utility_outer; chi_max = Stack_utility_outer; b_opt = chi_max;
Stack_c1_outer =  Stack_utility_outer; Stack_c2_outer = Stack_utility_outer; var_c1_outer =  Stack_utility_outer; var_c2_outer = Stack_utility_outer;
Stack_tau_outer = Stack_utility_outer; var_tau_outer = Stack_utility_outer; Stack_ri_outer = Stack_utility_outer; var_ri_outer = Stack_utility_outer;
Index_lower = Stack_utility_outer; Index_upper = Stack_utility_outer; Max_Resid1 = Stack_utility_outer; Max_Resid2 = Stack_utility_outer; 
Stack_ftau_outer = Stack_utility_outer;

tic

v_stack = linspace(0,1,n_v);

    for pp = 1:n_v

        v = v_stack(pp);

        pp

        GDP_Indexed_2024_nonlinear_LOOP_INNER

        [Stack_utility_outer(pp),Ind_max] = max(Stack_utility);

        Stack_b_outer(pp) = Stack_b(Ind_max);
        
        Stack_c1_outer(pp) = Stack_c1(Ind_max);
        Stack_c2_outer(pp) = Stack_c2(Ind_max);
        var_c1_outer(pp) = var_c1(Ind_max);
        var_c2_outer(pp) = var_c2(Ind_max);

        Stack_tau_outer(pp) = Stack_tau(Ind_max);
        var_tau_outer(pp) = var_tau(Ind_max);
        Stack_ftau_outer(pp) = Stack_ftau(Ind_max);
        
        Stack_ri_outer(pp) = Stack_ri(Ind_max);
        var_ri_outer(pp) = var_ri(Ind_max);

        Index_lower(pp) = min(Index_mini);
        Index_upper(pp) = max(Index_maxi);
        
        Max_Resid1(pp) = max(Resid_max);
        Max_Resid2(pp) = max(Resid_max2);

        chi = chi_stack(Ind_max);
        chi_max(pp) = chi;
        b_opt(pp) = Stack_b(Ind_max);

        Lambda_v(pp) = 100*( (Stack_utility_outer(pp)/Stack_utility_outer(1))^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain 

    end

U_mean = (1/(1-gama))*( Stack_c1_outer.^eps + betta*Stack_c2_outer.^eps  ).^((1-gama)/eps);

Lambda_mean = 100*( (U_mean/U_mean(1)).^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain
Lambda_var = Lambda_v - Lambda_mean;  

Lambda_tau = 100*( (Stack_ftau_outer/Stack_ftau_outer(1)).^0.5 - 1);
ftau_mean = 100*( sqrt(Stack_tau_outer.^2/Stack_tau_outer(1)^2) - 1);
ftau_var = Lambda_tau - ftau_mean;

%-----------------------------------
%Compute turning points w.r.t. v
%-----------------------------------

[U_max_max,Index_max_max] = max(Lambda_v);
vstar = v_stack(Index_max_max)

bstar_opt = b_opt(Index_max_max)
chi_opt = chi_max(Index_max_max)

[min_vri,Index_min] = min(var_ri);
vmin = v_stack(Index_min)

[min_vf,Index_minf] = min(Lambda_tau);
vminf = v_stack(Index_minf)

max(Max_Resid)
max(Max_Resid2)

figure(1)
subplot(2,2,1), plot(v_stack,Lambda_tau,'k','LineWidth', 1), xlabel('Indexation parameter (v)'), ylabel('% tax equiv.'), title('Tax burden: increase with indexation')
yline(0,'--k')
subplot(2,2,2), plot(v_stack,ftau_mean,'k','LineWidth', 1), hold on, plot(v_stack,ftau_var,'--k','LineWidth', 1), xlabel('Indexation parameter (v)'), ylabel('% tax equiv.'), title('Tax burden: decomposition')
subplot(2,2,3), plot(v_stack,Lambda_v,'k','LineWidth', 1), xlabel('Indexation parameter (v)'), ylabel('% c.e. welfare gain'), title('Welfare gain of indexation')
yline(0,'--k')
subplot(2,2,4), plot(v_stack,Lambda_mean,'k','LineWidth', 1), hold on, plot(v_stack,Lambda_var,'--k','LineWidth', 1), xlabel('Indexation parameter (v)'), ylabel('% c.e. welfare gain'), title('Welfare gain: decomposition')

figure(2)
subplot(2,3,1), plot(v_stack, Stack_c1_outer/Stack_c1_outer(1),'k','LineWidth', 1), hold on, plot(v_stack, Stack_c2_outer/Stack_c2_outer(1),'--k','LineWidth', 1), xlabel('Indexation parameter (v)'), title('Consumption means: $E[c_{j,t}]$')
subplot(2,3,2), plot(v_stack, var_c1_outer/var_c1_outer(1),'k','LineWidth', 1), hold on, plot(v_stack, var_c2_outer/var_c2_outer(1),'--k','LineWidth', 1), xlabel('Indexation parameter (v)'), title('Consumption variance: $var[c_{j,t}]$')
subplot(2,3,3), plot(v_stack, Stack_tau_outer/Stack_tau_outer(1),'k','LineWidth', 1), xlabel('Indexation parameter (v)'), title('Mean taxes: $E[\tau_t]$')
subplot(2,3,4), plot(v_stack, var_tau_outer/var_tau_outer(1),'k','LineWidth', 1), xlabel('Indexation parameter (v)'), title('Variance of taxes: $var[\tau_t]$')
subplot(2,3,5), plot(v_stack, Stack_ri_outer/Stack_ri_outer(1),'k','LineWidth', 1), xlabel('Indexation parameter (v)'), title('Mean return: $E[r^i_t]$')
subplot(2,3,6), plot(v_stack, var_ri_outer/var_ri_outer(1),'k','LineWidth', 1), xlabel('Indexation parameter (v)'), title('Return risk: $var[r^i_t]$')

toc


