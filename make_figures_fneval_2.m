function make_figures_fneval_2(epsilon,LeapSize,beta)
			clc;
			close all;
			clear;
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
					opts.epsilon = 1.2;
			else
					opts.epsilon = epsilon;
			end
			opts.beta = beta;

			% number of times to call the sampler
			Nsamp = 2000;
			% number of sampling steps to take in each sampler call
			opts.T = 1;
			opts.BatchSize = 1;
			% number of data dimensions
			opts.DataSize = 2;
			if 1 %nargin<2
					opts.LeapSize = 1;
			else
					opts.LeapSize = LeapSize;
			end
			opts.beta = 0.2;
			opts.skew_dims = 2;

			% scaling factor for energy function
			theta = [1,0;0,1e-5];

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
			Xstandard_persist = zeros(opts.DataSize,Nsamp);
			Xreduced_flip = zeros(opts.DataSize,Nsamp);
			X_2_momentum = zeros(opts.DataSize,Nsamp);

			ttt = tic();
			% call the sampling algorithm Nsamp times
			for ii = 1:Nsamp  

					 % the standard sampler
					[Xs, statestandard] = rf2vHMC( optsstandard, statestandard, theta );
					if ii == 1
							Xstandard = Xs;
					else
							Xstandard = cat(3,Xstandard, Xs);
					end
	%         fevalsstandard(ii) = statestandard.steps.total;
					%let's store a two-D array one for fevals and the other for diff in
					%cov(samples,model)
					fevalsstandard(ii,1) = statestandard.funcevals;
					fevalsstandard(ii,2) = calc_samples_err(Xstandard,theta);

	% %         the standard sampler with persistent momentum
					[Xs, statestandard_persist] = rf2vHMC( optsstandard_persist, statestandard_persist, theta );
					if ii == 1
							Xstandard_persist = Xs;
					else
							Xstandard_persist = cat(3,Xstandard_persist, Xs);
					end
					fevalsstandard_persist(ii) = statestandard_persist.steps.total;
					fevalsstandard_persist(ii,1) = statestandard_persist.funcevals;
					fevalsstandard_persist(ii,2) = calc_samples_err(Xstandard_persist,theta);       

	% %         the reduced flip sampler
					[Xs, statereduced_flip] = rf2vHMC( optsreduced_flip, statereduced_flip, theta );
					if ii == 1
							Xreduced_flip = Xs;
					else
							Xreduced_flip = cat(3,Xreduced_flip, Xs);
					end
					fevalsreduced_flip(ii) = statereduced_flip.steps.total;
					fevalsreduced_flip(ii,1) = statereduced_flip.funcevals;
					fevalsreduced_flip(ii,2) = calc_samples_err(Xreduced_flip,theta);

	% %         %the two momentum sampler
					[Xs, state_2_momentum] = rf2vHMC( opts_2_momemntum, state_2_momentum,theta );
					if ii == 1
							X_2_momentum = Xs;
					else
							X_2_momentum = cat(3,X_2_momentum,Xs);
					end
					feval_2_momentum(ii,1) = state_2_momentum.funcevals;
					feval_2_momentum(ii,2) = calc_samples_err(X_2_momentum,theta);

					if mod( ii, 100 ) == 1
							fprintf('%d / %d in %f sec (%f sec remaining)\n', ii, Nsamp, toc(ttt), toc(ttt)*Nsamp/ii - toc(ttt) );
					end
			end
			toc(ttt);

			%% and now a bunch of code to display the results    
	sprintf('\n');
	% print out some info
	statestandard
	statestandard_persist
	statereduced_flip
	state_2_momentum

	h1=plot_autocorr_samples(Xstandard,Xstandard_persist,Xreduced_flip,X_2_momentum);
	h2=plot_fevals(fevalsstandard,fevalsstandard_persist,fevalsreduced_flip,feval_2_momentum);
	% batch_2d_plot(Xstandard);
	savestr = strcat('LeapSize',int2str(LeapSize),'Epsilon',int2str(epsilon),'Beta',int2str(beta));
	savepath = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/',savestr);
	save(savepath);
end
