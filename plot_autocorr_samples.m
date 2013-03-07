function [h] = plot_autocorr_samples(X, names)

h=figure(222);
clf();

for ii in length(X)
    ac{ii} = calc_autocorr(X{ii})
    plot(ac{ii});
    hold on;
end

legend(names);
xlabel('Auto correlation windows');
ylabel('Auto correlation values');

end