%Steady_state_KL2_insert
%Computes numerical approximation to steady state which is later fed to the
%Dynare steady state solver (Algorithm described in Supp Appendix)

n_guess = 1E7;
tol = 8e-7;
l_min = 0.25; l_max = 0.85;
tau_min = 0.085; tau_max = 0.195; Dum = 0;

for i = 1:n_guess

        lab_guess = l_min + rand*(l_max-l_min);
        tau_guess = tau_min + rand*(tau_max-tau_min);
        tau_k_guess = psi*tau_guess;
        
        w = (1-alfa)*(alfa*(1-tau_k_guess)/chi)^(alfa/(1-alfa));
        TR = tau_p*(1+n)*w*lab_guess;

        s = ( ( betta*chi*(1/(1-lab_guess))^( (1-thetta)*eps ) )^(1/(1-thetta*eps))*(1-tau_guess-tau_p)*w*lab_guess - TR) / ( chi + ( betta*chi*(1/(1-lab_guess))^( (1-thetta)*eps ) )^(1/(1-thetta*eps)) );
        k_tild = (alfa*(1-tau_k_guess)/chi)^(1/(1-alfa))*lab_guess;
        b = s - (1+n)*k_tild; 

        c1 = (1-tau_guess-tau_p)*w*lab_guess - s;
        c2 = chi*s + TR;

        Resid1 = (1-thetta)/thetta*c1 / ((1-tau_guess-tau_p)*w) - (1-lab_guess);
        Resid2 = tau_guess*(w*lab_guess + psi*(chi/(1-psi*tau_guess))*k_tild ) - gbar - (chi/(1+n) -1)*b; 

        Loss = Resid1^2 + Resid2^2; 

        if Loss < tol
            Dum = 1;
        end

        if Dum==1
            break
        end

end

b_ss = b; 
k_ss = (1+n)*k_tild; 
s_ss = s; 
tau_ss = tau_guess;
c1_ss = c1; c2_ss = c2;
l_ss = lab_guess;
utility_ss = (1/(1-gama))*(  (c1_ss^thetta*(1-l_ss)^(1-thetta))^eps + betta*(c2_ss^thetta)^eps )^((1-gama)/eps);
utility2_ss = (1/(1-gama))*(c2_ss^thetta)^(1-gama); 
          