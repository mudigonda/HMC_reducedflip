function [h_scaled,h_notscaled,ac] = plot_autocorr_samples(X, names,avg_fevals,Mu)
colorlist=['r','g','b','k','m','y'];
h_notscaled=figure(222);
clf();
acn=500; %The number of windows we are looking at, stick to 500 samples as default
for ii=1:length(X)
    if exist('Mu','var')
        ac{ii} = calc_autocorr_mog(X{ii},acn,Mu)
    else
        ac{ii} = calc_autocorr(X{ii},acn)
    end
    plot(ac{ii},colorlist(ii));
    hold on;
end

legend(names);
xlabel('Auto correlation windows');
ylabel('Auto correlation values');

h_scaled=figure(333);
clf();

for ii=1:length(X)
    plot((1:length(ac{ii}))*avg_fevals{ii},ac{ii},colorlist(ii));
    hold on;
end
legend(names);
xlabel('Auto correlation windows scaled by avg fevals');
ylabel('Auto correlation values');

end
