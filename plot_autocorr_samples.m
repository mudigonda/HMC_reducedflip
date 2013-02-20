function [h] = plot_autocorr_samples(Xstandard,Xstandard_persist,...
    Xreduced_flip,X_2_momentum)

ac_std=calc_autocorr(Xstandard);
ac_std_per=calc_autocorr(Xstandard_persist);
ac_red_flip=calc_autocorr(Xreduced_flip);
ac_2mom=calc_autocorr(X_2_momentum);

h=figure(222);
plot(ac_std,'r');
hold on;
plot(ac_std_per,'g');
plot(ac_red_flip,'b');
plot(ac_2mom,'y');
legend('standard','standard, persistent mom','reduced flip','two momentum');
xlabel('Auto correlation windows');
ylabel('Auto correlation values');

end