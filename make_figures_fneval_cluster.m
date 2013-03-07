function make_figures_fneval_2(LeapSize,epsilon,beta)

whos
clc;
close all;
opts_init = [];
% the energy function and gradient for the circle distribution in the arXiv
% TODO: I think I should have used a very long narrow Gaussian instead!
% opts_init.E = @E_circle;
opts_init.E = @E_gauss;
% opts_init.dEdX = @dEdX_circle;
opts_init.dEdX = @dEdX_gauss;
% make this 1 for more output
opts_init.Debug = 0;
% step size for HMC
if nargin<1
    opts_init.LeapSize = 1;
else
    LeapSize=str2num(LeapSize);
    opts_init.LeapSize = LeapSize
end
if nargin<2
    opts_init.epsilon = 1.2;
else
    epsilon =str2num(epsilon);
    opts_init.epsilon = epsilon
end
if nargin<3
    opts_init.beta = .03;
else
    beta    =str2num(beta);
    opts_init.beta = beta
end            

% number of times to call the sampler
Nsamp = 5000;
% number of sampling stpdf to take in each sampler call
% 			opts_init.T = 1;
opts_init.BatchSize = 100;
% number of data dimensions
opts_init.DataSize = 2;
opts_init.funcevals = 0

% scaling factor for energy function
theta = [1,0;0,1e-3];
FEVAL_MAX = 250000

%Initalize Options
ii = 1
names{ii} = 'standard'
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 0;
opts{ii}.beta = 1;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = []

ii = ii + 1
names{ii} = 'persist'
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 0;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = []

ii = ii + 1
names{ii} = 'reduced flip'
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 1;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = []

ii = ii + 1
names{ii} = 'two momentum'
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 2;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = []



ttt = tic();
    % call the sampling algorithm Nsamp times
    for ii = 1:Nsamp  
        for jj = 1:length(names)
            tic()
                if states{jj}.funcevals < FEVAL_MAX 
                    [Xloc, statesloc] = rf2vHMC( opts{jj}, states{jj}, theta );
                    if ii > 1                                        
                        X{jj} = cat(3,X{jj}, Xloc);
                    else
                        X{jj} = Xloc;
                    end
                    
                    fevals{jj}(ii,1) = statestandard.funcevals;
                    fevals{jj}(ii,2) = calc_samples_err(Xstandard,theta);
                end
            toc()
        end
        
        %Display + Saving 
        if (mod( ii, 100 ) == 1) || (ii == Nsamp)
            fprintf('%d / %d in %f sec (%f sec remaining)\n', ii, Nsamp, toc(ttt), toc(ttt)*Nsamp/ii - toc(ttt) );
            h1=plot_autocorr_samples(X, names);
            h2=plot_fevals(fevals, names);
            savestr = strcat('LeapSize-',int2str(LeapSize),'epsilon-',int2str(epsilon*10),'Beta-',int2str(beta*100));
            disp(savestr);
            savepath = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/',savestr);
            figpath1 = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/figures/',savestr,'autocor');
            figpath2 = strcat('/clusterfs/cortex/scratch/mayur/HMC_reducedflip/figures/',savestr,'fneval');
            saveas(h1,figpath1,'pdf');
            saveas(h2,figpath2,'pdf');
            save(savepath);
        end
    end
ttt = toc(ttt);
