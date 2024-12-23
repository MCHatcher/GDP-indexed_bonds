%GDP_Indexed_2024_nonlinear_SIM: inner loop simulatios.
%Insert for the file GDP_Indexed_2024_nonlinear_LOOP.m
%Nonlinear simulation of the simple model (Algorithm in Supp Appendix).
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). Any errors are my own.

tau = NaN(T_sim,1); c1 = tau; c2 = tau; Index = tau; Max_Resid = tau; r_i  = tau;
r = tau; EV = tau; EV_reverse = tau; Utility = tau; ftau = tau; Resid_check = tau;

for t=1:T_sim

    if t==1

        r_lag = r_init;
        y_lag = y_init;

    else

        r_lag = r(t-1);
        y_lag = y(t-1);

    end

    r_i(t) = r_lag*(y(t)/ybar)^v;
    tau(t) = ( gbar + (r_i(t)/(1+n) - 1)*bstar ) / y(t);
    ftau(t) = phi*tau(t)^2;
    c1(t) = (1-alfa)*y(t) - tau(t)*y(t) - phi*tau(t)^2 - bstar;
    c2(t) = alfa*(1+n)*y(t) + r_i(t)*bstar;

    ret_prime = (y_prime/ybar).^v;
    Resid0 = NaN(N_guess0,1);  

    for k0=1:N_guess0

        Dum = 0;
        r_guess = r_guess_stack0(k0);    
        cprime0 = alfa*(1+n)*y_prime + r_guess*bstar*ret_prime;
        c_power0 = cprime0.^(1-gama);

        E_c_power0 = prob*c_power0;
        E_reverse0 = E_c_power0^(1/(1-gama));

        sdf_adj0 = betta*c1(t)^(1-eps).*cprime0.^(eps-1)*( 1 / E_reverse0 )^(1-gama-eps).*cprime0.^(1-gama-eps).*ret_prime;

        Resid0(k0) = abs(1 - r_guess*prob*sdf_adj0);

        if k0 > 1 && Resid0(k0) > Resid0(k0-1)
            Dum = 1;
            break 
        end  

    end

    [Resid_mini0,Index_min0] = min(Resid0);
    r(t) = r_guess_stack0(Index_min0);

    Max_Resid(t) = Resid_mini0;

    Index(t) = Index_min0;
        
    cprime = alfa*(1+n)*y_prime + r(t)*bstar*ret_prime;
    c_power = cprime.^(1-gama);
    E_c_power = prob*c_power;
    E_reverse = E_c_power^(1/(1-gama));

    SDF_check = betta*c1(t)^(1-eps).*cprime.^(eps-1)*( 1 / E_reverse )^(1-gama-eps).*cprime.^(1-gama-eps).*ret_prime;
    Resid_check(t) = abs(1 - r(t)*prob*SDF_check);

    EV(t) = prob*c_power;

    Utility(t) =  1/(1-gama)*( c1(t)^eps + betta*( EV(t) )^(eps/(1-gama)) )^((1-gama)/eps);

end




