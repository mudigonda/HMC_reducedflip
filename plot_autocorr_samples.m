%%Plots the auto-correlation of the samples from the various models. It
%%uses a sliding window scale. If the sampler works well for a set of
%%parameters, the autocorrelation reaches zero.
%%Input = Xstandard, Xstandard_persist, Xreduced_flip, X_2mom
%%Output = handle to plot
function h_ac = plot_autocorr_samples(Xstandard,Xstandard_persist,Xreduced_flip,X_2_momentum)

acn=500;

acstandard = calc_autocorr(Xstandard,acn);

acstandard_persist = calc_autocorr(Xstandard_persist,acn);

acreduced_flip = calc_autocorr(Xreduced_flip,acn);

if nargin <4
    %Actually plot
    h_ac=figure(12); clf;
    hold on;
    plot(1:acn,acstandard,'r');
    plot(1:acn,acstandard_persist,'g');
    plot(1:acn,acreduced_flip,'b');    
    title( 'Autocorrelation' );
    xlabel( 'Autcorrel windows' );
    ylabel( 'Correlation' );
    legend( 'Standard', 'Standard, persistent momentum', 'Reduced flip');
    hold off;
else
    %Actually plot   
    ac2mom_reduced = calc_autocorr(X_2_momentum,acn);
    h_ac=figure(12); clf;
    hold on;
    plot(1:acn,acstandard,'r');
    plot(1:acn,acstandard_persist,'g');
    plot(1:acn,acreduced_flip,'b');
    plot(1:acn,ac2mom_reduced,'y');    
    title( 'Autocorrelation' );
    xlabel( 'Autcorrel windows' );
    ylabel( 'Correlation' );
    legend( 'Standard', 'Standard, persistent momentum', 'Reduced flip','Two momentum reduced flip' );
    hold off;
end
end
