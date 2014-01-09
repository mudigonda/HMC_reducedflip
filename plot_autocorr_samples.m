function [h_scaled,h_notscaled,ac] = plot_autocorr_samples(X, names,avg_fevals, max_shift)

fig_shape = [1,1,4,3];

colorlist={'r:','b--','m-.','k-','y-','r:'};
h_notscaled=figure(222);
clf();
set(h_notscaled, 'Units','inches');
set(h_notscaled, 'Position', fig_shape);
set(h_notscaled, 'Color', 'w');

num_steps = 500;
max_shift = min([max_shift, size(X{1},3)]);
target_shifts = ceil(logspace(0,log10(max_shift),num_steps));
target_shifts = target_shifts - min(target_shifts);
target_shifts = unique(target_shifts);

%acn=500; %The number of windows we are looking at, stick to 500 samples as default
for ii=1:length(X)
    ac{ii} = calc_autocorr(X{ii},target_shifts);
    ac{ii} = ac{ii}/ac{ii}(1);
    plot(target_shifts,ac{ii},colorlist{ii});
    hold on;
end

grid on;
legend(names);
xlabel('Sampling Steps');
ylabel('Autocorrelation');

h_scaled=figure(333);
clf();
set(h_scaled, 'Units','inches');
set(h_scaled, 'Position', fig_shape);
set(h_scaled, 'Color', 'w');

min_fevals = inf;

for ii=1:length(X)
    min_fevals = min([max(target_shifts*avg_fevals{ii}), min_fevals]);
    plot(target_shifts*avg_fevals{ii},ac{ii},colorlist{ii});
    hold on;
end
axis([0, min_fevals, -0.1, 1.01]);
grid on;
legend(names);
xlabel('Gradient Evaluations');
ylabel('Autocorrelation');

end
