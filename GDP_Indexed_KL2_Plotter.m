%GDP_Indexed_KL2_Plotter

%--------------
%Plot results
%--------------
figure(1)
subplot(1,2,1), plot(Stack_b,Stack_Lambda_ss,'k','LineWidth', 1), hold on, plot(Stack_b,Stack_Lambda_ss1,'--k','LineWidth', 1)
title('Welfare'), hold on, xlabel('Steady state bond supply'), ylabel('% cons. equiv.'), hold on,
subplot(1,2,2), plot(Stack_b,Tax_burden_ss,'k','LineWidth', 1)
title('Tax burden'), hold on, xlabel('Steady state bond supply'), ylabel('% cons. equiv.')








