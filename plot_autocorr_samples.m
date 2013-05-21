function [h_scaled,h_notscaled] = plot_autocorr_samples(X, names,avg_fevals)
colorlist=['r','g','b','k','m','y'];
h_scaled=figure(222);
clf();

for ii=1:length(X)
    ac{ii} = calc_autocorr(X{ii})
%     ac{ii} = ac{ii}*avg_fevals{ii};
    plot(ac{ii},colorlist(ii));
    hold on;
end

legend(names);
xlabel('Auto correlation windows');
ylabel('Auto correlation values');

h_notscaled=figure(333);
clf();

for ii=1:length(X)
    plot((1:length(ac{ii}))*avg_fevals{ii},ac{ii},colorlist(ii));
    hold on;
end
legend(names);
xlabel('Auto correlation windows scaled by avg fevals');
ylabel('Auto correlation values');


end
