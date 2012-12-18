%function make_figures(epsilon,LeapSize)
    %clc;
    %close all;
%     clear;
    opts = [];
    % the energy function and gradient for the circle distribution in the arXiv
    % TODO: I think I should have used a very long narrow Gaussian instead!
    % opts.E = @E_circle;
    opts.E = @E_gauss;
    % opts.dEdX = @dEdX_circle;
    opts.dEdX = @dEdX_gauss;
    % make this 1 for more output
    opts.Debug = 0;
    % step size for HMC
    if 1 %nargin<1
        opts.epsilon = 1.7;
    else
        opts.epsilon = epsilon;
    end
    %opts.epsilon = 0.05;
    %opts.beta = 0.03;

    % number of times to call the sampler
%     Nsamp = 5e5;
%     Nsamp = 1e4;
    Nsamp = 100;
    % number of sampling steps to take in each sampler call
    opts.T = 1;
    opts.BatchSize = 3;
    % number of data dimensions
    opts.DataSize = 2;
    if 1 %nargin<2
        opts.LeapSize = 10;
    else
        opts.LeapSize = LeapSize;
    end
    opts.beta = 0.03;

    % scaling factor for energy function
    % theta = 100;
%     theta = [1,0;0,0.00001];
    theta = [1,0;0,1e-3];

    %Initalize Options
    optsstandard = opts;
    optsstandard.FlipOnReject = 0;
    optsstandard.beta = 1;

    optsstandard_persist = opts;
    optsstandard_persist.FlipOnReject = 0;

    optsreduced_flip = opts;
    optsreduced_flip.FlipOnReject = 1;

    opts_2_momemntum = opts;
    opts_2_momemntum.FlipOnReject = 2;

    
    %Initialize States
    statestandard = [];
    statestandard_persist = [];
    state_2_momentum = [];
    statereduced_flip = [];

    % arrays to keep track of the samples
    Xstandard = zeros(opts.DataSize,Nsamp);
    fevalsstandard = zeros(1,Nsamp);

    Xstandard_persist = zeros(opts.DataSize,Nsamp);
    fevalsstandard_persist = zeros(1,Nsamp);

    Xreduced_flip = zeros(opts.DataSize,Nsamp);
    fevalsreduced_flip = zeros(1,Nsamp);

    X_2_momentum = zeros(opts.DataSize,Nsamp);
    feval_2_momentum = zeros(1,Nsamp);

    ttt = tic();
    % call the sampling algorithm Nsamp times
    for ii = 1:Nsamp  

         % the standard sampler
        [Xs, statestandard] = rf2vHMC( optsstandard, statestandard, theta );
        if ii == 1
            Xstandard = Xs;
        else
            Xstandard = horzcat(Xstandard, Xs);
        end
        fevalsstandard(ii) = statestandard.steps.total;

        % the standard sampler with persistent momentum
        [Xs, statestandard_persist] = rf2vHMC( optsstandard_persist, statestandard_persist, theta );
        if ii == 1
            Xstandard_persist = Xs;
        else
            Xstandard_persist = horzcat(Xstandard_persist, Xs);
        end
        fevalsstandard_persist(ii) = statestandard_persist.steps.total;

        % the reduced flip sampler
        [Xs, statereduced_flip] = rf2vHMC( optsreduced_flip, statereduced_flip, theta );
        if ii == 1
            Xreduced_flip = Xs;
        else
            Xreduced_flip = horzcat(Xreduced_flip, Xs);
        end
        fevalsreduced_flip(ii) = statereduced_flip.steps.total;

        %the two momentum sampler
        %[Xs, state_2_momentum] = rf2vHMC( opts_2_momemntum, state_2_momentum,theta );
        %X_2_momentum(:,ii) = Xs;
        %feval_2_momentum(ii) = state_2_momentum.steps.total;
        

        if mod( ii, 100 ) == 1
            fprintf('%d / %d in %f sec (%f sec remaining)\n', ii, Nsamp, toc(ttt), toc(ttt)*Nsamp/ii - toc(ttt) );
        end
    end
    toc(ttt);

    save
    
    %% and now a bunch of code to display the results

    % DEBUG
    state_2_momentum = statestandard;
    
    fprintf('\n');
    % print out some info
    statestandard
    statestandard.steps
    statestandard_persist
    statestandard_persist.steps
    statereduced_flip
    statereduced_flip.steps
    state_2_momentum
    state_2_momentum.steps

    % make an image out of the samples
    figure(4); clf;
    subplot( 4, 1, 1 );
    imagesc( Xstandard ); colorbar;
    title( 'circle - Xstandard' );
    subplot( 4, 1, 2 );
    imagesc( Xstandard_persist ); colorbar;
    title( 'circle - Xstandard_persist' );
    subplot( 4, 1, 3 );
    imagesc( Xreduced_flip ); colorbar;
    title( 'circle - Xreduced_flip' );
    subplot( 4, 1, 4 );
    imagesc( X_2_momentum ); colorbar;
    title( 'circle - X_2_momentum' );


    % auto-correlation plot
    figure(12); clf;
    mn = mean( [Xreduced_flip, Xstandard], 2 );
    mn(:) = 0;
    acn = ceil(size(Xreduced_flip,2)/2);
    acn = 5000;

    tmp = Xstandard;
    tmp = bsxfun( @plus, tmp, -mn );
    C = tmp*tmp'/size(tmp,2);
    %tmp = sqrtm(inv(C))*tmp;
    autocorr = zeros(acn+1,1);
    for ii = 0:acn
        autocorr(ii+1) = mean(sum(tmp(:,1:end-ii).*tmp(:,1+ii:end)));
    end
    autocorr = autocorr / autocorr(1);
    acstandard = autocorr;

    tmp = Xstandard_persist;
    tmp = bsxfun( @plus, tmp, -mn );
    C = tmp*tmp'/size(tmp,2);
    %tmp = sqrtm(inv(C))*tmp;
    autocorr = zeros(acn+1,1);
    for ii = 0:acn
        autocorr(ii+1) = mean(sum(tmp(:,1:end-ii).*tmp(:,1+ii:end)));
    end
    autocorr = autocorr / autocorr(1);
    acstandard_persist = autocorr;

    tmp = Xreduced_flip;
    mn(:) = 0;
    %mn(:) = 2;
    tmp = bsxfun( @plus, tmp, -mn );
    C = tmp*tmp'/size(tmp,2);
    %tmp = sqrtm(inv(C))*tmp;
    autocorr = zeros(acn+1,1);
    for ii = 0:acn
        autocorr(ii+1) = mean(sum(tmp(:,1:end-ii).*tmp(:,1+ii:end)));
    end
    autocorr = autocorr / autocorr(1);
    acreduced_flip = autocorr;

    tmp=X_2_momentum;
    tmp = bsxfun(@plus,tmp,-mn);
    C = tmp*tmp'/size(tmp,2);
    autocorr = zeros(acn+1,1);
    for ii=0:acn
        autocorr(ii+1) = mean(sum(tmp(:,1:end-ii).*tmp(:,1+ii:end)));
    end
    autocorr = autocorr / autocorr(1);
    ac2mom_reduced = autocorr;

    figure(12); clf;
    subplot( 2, 1, 1 );
    %plot( [0:(acn)]' * [statestandard.funcevals, statestandard_persist.funcevals, statereduced_flip.funcevals, state_2_momentum.funcevals] / size(Xreduced_flip,2), [acstandard, acstandard_persist, acreduced_flip,ac2mom_reduced] );
    plot( [0:(acn)]' * [statestandard.funcevals, statestandard_persist.funcevals, statereduced_flip.funcevals] / size(Xreduced_flip,2), [acstandard, acstandard_persist, acreduced_flip] );
    title( 'Autocorrelation' );
    xlabel( 'Gradient evaluations' );
    ylabel( 'Correlation' );
    %legend( 'Standard', 'Standard, persistent momentum', 'Reduced flip','Two momentum reduced flip' );
    legend( 'Standard', 'Standard, persistent momentum', 'Reduced flip', 'Location', 'best' );

    subplot( 2, 1, 2 );
    %plot( [0:(acn)]' * [statestandard.funcevals, statestandard_persist.funcevals, statereduced_flip.funcevals, state_2_momentum.funcevals] / size(Xreduced_flip,2), [acstandard, acstandard_persist, acreduced_flip,ac2mom_reduced] );
    plot( [0:(acn)]' * [statestandard.funcevals, statestandard_persist.funcevals, statereduced_flip.funcevals] / size(Xreduced_flip,2), [acstandard, acstandard_persist, acreduced_flip] );
    %semilogx( [1:(acn+1)]' * [statestandard.funcevals, statestandard_persist.funcevals, statereduced_flip.funcevals] / size(Xreduced_flip,2), [acstandard, acstandard_persist, acreduced_flip] );
    title( 'Autocorrelation' );
    xlabel( 'Gradient evaluations' );
    ylabel( 'Correlation' );
    %legend( 'Standard', 'Standard, persistent momentum', 'Reduced flip','Two momentum reduced flip' );
    legend( 'Standard', 'Standard, persistent momentum', 'Reduced flip', 'Location', 'best' );

    
    % plot the sample positions
    figure(70); clf;
    fmax = 200;
    gd = (fevalsstandard >= (max(fevalsstandard)-fmax));
    plot( Xstandard(1,gd), Xstandard( 2, gd ), '.g' );
    hold on;
    gd = (fevalsstandard_persist >= (max(fevalsstandard_persist)-fmax));
    plot( Xstandard_persist(1,gd), Xstandard_persist( 2, gd ), '.y' );
    gd = (fevalsreduced_flip >= (max(fevalsreduced_flip)-fmax));
    plot( Xreduced_flip(1,gd), Xreduced_flip( 2, gd ), '.b' );
    gd = (feval_2_momentum >= (max(feval_2_momentum)-fmax));
    plot( X_2_momentum(1,gd),X_2_momentum(2,gd),'.r');
    legend( 'standard', 'standard, persistent momentum', 'reduced_flip','two momentum reduced flip' );
    title( 'samples' );

    % plot the sample positions on top of an image of the probability density

    mnx = min(Xreduced_flip(1,:));
    mxx = max(Xreduced_flip(1,:));
    mny = min(Xreduced_flip(2,:));
    mxy = max(Xreduced_flip(2,:));
    xxc = linspace( mnx, mxx, 256 );
    yyc = linspace( mny, mxy, 256 );
    [ex,ey] = meshgrid(xxc,yyc);
    ee = opts.E( [ex(:)';ey(:)'], theta );
    ee = reshape(ee,size(ex));
    %ee(ee>4) = 4;
    ee = exp(-ee);
    figure(77); clf;
    %imshow(ee, 'DisplayRange', [-1,4], 'Xdata', [-2,2], 'Ydata', [-2,2]); colormap gray; colorbar;
    imagesc(ee,'Xdata', [mnx,mxx], 'Ydata', [mny,mxy]); colormap gray; colorbar;
    hold on;
    fmax = 1000;
    gd = (fevalsstandard >= (max(fevalsstandard)-fmax));
    plot( Xstandard(1,gd), Xstandard( 2, gd ), '.g' );
    hold on;
    gd = (fevalsstandard_persist >= (max(fevalsstandard_persist)-fmax));
    plot( Xstandard_persist(1,gd), Xstandard_persist( 2, gd ), '.y' );
    gd = (fevalsreduced_flip >= (max(fevalsreduced_flip)-fmax));
    plot( Xreduced_flip(1,gd), Xreduced_flip( 2, gd ), '.b' );
    gd = (feval_2_momentum >= (max(feval_2_momentum)-fmax));
    plot( X_2_momentum(1,gd),X_2_momentum(2,gd),'.r');
    legend( 'standard', 'standard, persistent momentum', 'reduced_flip','two momentum reduced flip' );
    title( 'samples' );

    % a histogram to compare the two sampling conditions
    figure(87); clf;
%     hist( [Xstandard(:), Xstandard_persist(:), Xreduced_flip(:), X_2_momentum(:)], 50 ); 
    hist( [Xstandard(:), Xstandard_persist(:), Xreduced_flip(:)], 50 ); 
    legend( 'standard', 'standard, persistent momentum', 'reduced_flip','two momentum reduced flip' );
    
%     
%         figure(237); clf;
%     imshow(ee, 'DisplayRange', [-1,4], 'Xdata', [-2,2], 'Ydata', [-2,2]); colormap gray; colorbar;
%     imagesc(ee','Xdata', [mny,mxy], 'Ydata', [mnx,mxx]); colormap gray; colorbar;

    
%end
