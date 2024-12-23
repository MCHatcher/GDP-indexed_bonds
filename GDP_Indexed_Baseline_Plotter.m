%GDP_Indexed_Baseline_Plotter
%Written by M. Hatcher
%Used in GDP_Indexed_2024_Baseline.mod and its variants

%--------------
%Plot results
%--------------
Lambda = NaN(n_loop,1);   
  for i=1:n_loop
        Lambda(i) = 100*( (Stack_utility(i)/Stack_utility(1))^(1/(1-gama)) - 1);  %Consumption equiv. welfare gain 
        if phi == 0
             Stack_utility0 = -0.024854062;  %Update if parameters change
             Lammba0_ss = Lambda_ss;
             Lambda_ss = 100*( (Stack_utility/Stack_utility0).^(1/(1-gama)) - 1);
        end
  end

figure(1)
hold on, subplot(1,2,1), plot(r_ss, b_ss, '--k','DisplayName', '\Phi = 0.5', 'LineWidth', 1), hold on, 
title('Bond supply vs interest rate'), xlabel('Steady state interest rate'), ylabel('Steady state bond supply')
hold on, subplot(1,2,2), plot(r_ss, Lambda_ss, '--k','LineWidth', 1) 
title('Interest rate vs steady state welfare'), xlabel('Steady state interest rate'), ylabel('% c.e. welfare gain')
xline(1.4,'--k')

figure(2)
hold on, plot(b_ss, Lambda, '--k','LineWidth', 1) 
title('Bond supply vs expected welfare'), xlabel('Fixed bond supply'), ylabel('% c.e. welfare gain')






