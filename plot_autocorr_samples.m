function [h_scaled,h_notscaled,ac] = plot_autocorr_samples(X, names,avg_fevals, max_shift)
colorlist=['r','g','b','k','m','y'];
h_notscaled=figure(222);
clf();

num_steps = 500;
max_shift = min([max_shift, size(X{1},3)]);
target_shifts = ceil(logspace(0,log10(max_shift),num_steps));
target_shifts = target_shifts - min(target_shifts);
target_shifts = unique(target_shifts);

%acn=500; %The number of windows we are looking at, stick to 500 samples as default
for ii=1:length(X)
    ac{ii} = calc_autocorr(X{ii},target_shifts);
    ac{ii} = ac{ii}/ac{ii}(1);
    plot(target_shifts,ac{ii},colorlist(ii));
    hold on;
end

legend(names);
xlabel('Sampling Steps');
ylabel('Autocorrelation');

h_scaled=figure(333);
clf();

for ii=1:length(X)
    plot(target_shifts*avg_fevals{ii},ac{ii},colorlist(ii));
    hold on;
end
legend(names);
xlabel('Gradient Evaluations');
ylabel('Autocorrelation');

end
