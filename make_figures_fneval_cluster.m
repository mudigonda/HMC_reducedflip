function make_figures_fneval_2(LeapSize,epsilon,beta)
whos
			clc;
			close all;
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
			if nargin<1
					opts.LeapSize = 1;
			else
					LeapSize=str2num(LeapSize);
					opts.LeapSize = LeapSize
			end
			if nargin<2
					opts.epsilon = 1.2;
			else
					epsilon =str2num(epsilon);
					opts.epsilon = epsilon
			end
			if nargin<3
					opts.beta = .03;
			else
					beta    =str2num(beta);
					opts.beta = beta
			end            

			% number of times to call the sampler
			Nsamp = 5000;
			% number of sampling stpdf to take in each sampler call
% 			opts.T = 1;
			opts.BatchSize = 100;
			% number of data dimensions
			opts.DataSize = 2;
			opts.funcevals = 0

			% scaling factor for energy function
			theta = [1,0;0,1e-3];
			FEVAL_MAX = 20000

			%Initalize Options
			optsstandard = opts;
			optsstandard.FlipOnReject = 0;
			optsstandard.beta = 1;

			optsstandard_persist = opts;
			optsstandard_persist.FlipOnReject = 0;

			optsreduced_flip = opts;
			optsreduced_flip.FlipOnReject = 1;

			opts_2_momentum = opts;
			opts_2_momentum.FlipOnReject = 2;

			
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
				tic()
			if ii > 1
				if statestandard.funcevals < FEVAL_MAX 
						[Xs, statestandard] = rf2vHMC( optsstandard, statestandard, theta );
						Xstandard = cat(3,Xstandard, Xs);
						statestandard.funcevals
				end
			else
						[Xs, statestandard] = rf2vHMC( optsstandard, statestandard, theta );
						Xstandard = Xs;
			end
				toc()
					fevalsstandard(ii,1) = statestandard.funcevals;
					fevalsstandard(ii,2) = calc_samples_err(Xstandard,theta);
	% %         the standard sampler with persistent momentum
				tic()
				if ii > 1
					if statestandard_persist.funcevals < FEVAL_MAX 
						[Xs, statestandard_persist] = rf2vHMC( optsstandard_persist, statestandard_persist, theta );
						Xstandard_persist = cat(3,Xstandard_persist, Xs);
						statestandard_persist.funcevals
				 end
			 else
						[Xs, statestandard_persist] = rf2vHMC( optsstandard_persist, statestandard_persist, theta );
						Xstandard_persist = Xs;
			 end
				toc()
						fevalsstandard_persist(ii,1) = statestandard_persist.funcevals;
						fevalsstandard_persist(ii,2) = calc_samples_err(Xstandard_persist,theta);      

	% %         the reduced flip sampler
				tic()
			if ii > 1
				if statereduced_flip.funcevals < FEVAL_MAX 
					[Xs, statereduced_flip] = rf2vHMC( optsreduced_flip, statereduced_flip, theta );
					Xreduced_flip = cat(3,Xreduced_flip, Xs);
					statereduced_flip.funcevals
				end
			else
					[Xs, statereduced_flip] = rf2vHMC( optsreduced_flip, statereduced_flip, theta );
					Xreduced_flip = Xs;
			end
				toc()
				fevalsreduced_flip(ii,1) = statereduced_flip.funcevals;
				fevalsreduced_flip(ii,2) = calc_samples_err(Xreduced_flip,theta);
	% %         %the two momentum sampler
				tic()
				if ii > 1
				 if state_2_momentum.funcevals < FEVAL_MAX 
					[Xs, state_2_momentum] = rf2vHMC( opts_2_momentum, state_2_momentum,theta );
					X_2_momentum = cat(3,X_2_momentum,Xs);
					state_2_momentum.funcevals
				 end
				else
					[Xs, state_2_momentum] = rf2vHMC( opts_2_momentum, state_2_momentum,theta );
					X_2_momentum = Xs;
				end
				toc()
					feval_2_momentum(ii,1) = state_2_momentum.funcevals;
					feval_2_momentum(ii,2) = calc_samples_err(X_2_momentum,theta);

				%Display + Saving 
					if mod( ii, 100 ) == 1
							fprintf('%d / %d in %f sec (%f sec remaining)\n', ii, Nsamp, toc(ttt), toc(ttt)*Nsamp/ii - toc(ttt) );
							h1=plot_autocorr_samples(Xstandard,Xstandard_persist,Xreduced_flip,X_2_momentum);
							h2=plot_fevals(fevalsstandard,fevalsstandard_persist,fevalsreduced_flip,feval_2_momentum);
							LeapSize
							int2str(epsilon)
							int2str(beta)
							savestr = strcat('LeapSize-',int2str(LeapSize),'epsilon-',int2str(epsilon*10),'Beta-',int2str(beta*100));
							savepath = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/',savestr);
							figpath1 = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/figures/',savestr,'autocor');
							figpath2 = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/figures/',savestr,'fneval');
							saveas(h1,figpath1,'pdf');
							saveas(h2,figpath2,'pdf');
							save(savepath);
					end
			end
			toc(ttt);

			%% and now a bunch of code to display the results    
	sprintf('\n');
	fevalsstandard(end,:)
	fevalsstandard_persist(end,:)
	fevalsreduced_flip(end,:)
	feval_2_momentum(end,:)

	h1=plot_autocorr_samples(Xstandard,Xstandard_persist,Xreduced_flip,X_2_momentum);
	h2=plot_fevals(fevalsstandard,fevalsstandard_persist,fevalsreduced_flip,feval_2_momentum);
	% batch_2d_plot(Xstandard);
	savestr = strcat('LeapSize-',int2str(LeapSize),'epsilon-',int2str(epsilon*10),'Beta-',int2str(beta*100));
	savepath = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/',savestr);
	figpath1 = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/figures/',savestr,'autocor');
	figpath2 = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/figures/',savestr,'fneval');
	saveas(h1,figpath1,'pdf');
	saveas(h2,figpath2,'pdf');
	save(savepath);
end
