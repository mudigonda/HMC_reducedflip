function make_figures_fneval_local(LeapSize,epsilon,beta)

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

%Model Name
FEVAL_MAX = 5000000
modelname='2dGausSkew10-6'
savestr = strcat('ModelName-',modelname,'-LeapSize-',int2str(opts_init.LeapSize),...
    '-epsilon-',int2str(opts_init.epsilon*10),'-Beta-',int2str(opts_init.beta*100)...
    ,'-fevals-',int2str(FEVAL_MAX));
savepath = strcat('/Users/mudigonda/Data/HMC_reducedflip/2d/',savestr);
figpath1 = strcat('/Users/mudigonda/Data/HMC_reducedflip/2d/figures/',savestr,'autocor');
figpath2 = strcat('/Users/mudigonda/Data/HMC_reducedflip/2d/figures/',savestr,'fneval');
% number of times to call the sampler
Nsamp = 6000;
% number of sampling stpdf to take in each sampler call
% 			opts_init.T = 1;
opts_init.BatchSize = 1000;
% number of data dimensions
opts_init.DataSize = 2;
opts_init.funcevals = 0;

% scaling factor for energy function
theta = [1,0;0,1e-6];


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

if 0
ii = ii + 1
names{ii} = 'two momentum'
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 2;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = []
end

RUN_FLAG=1;
ttt = tic();
ii=1;
    % call the sampling algorithm Nsamp times
    while (ii <=Nsamp && RUN_FLAG == 1)
        for jj = 1:length(names)
            tic()
                if ii == 1 || states{jj}.funcevals < FEVAL_MAX 
                    [Xloc, statesloc] = rf2vHMC( opts{jj}, states{jj}, theta );
                    states{jj} = statesloc
                    if ii > 1                                        
                        X{jj} = cat(3,X{jj}, Xloc);
                    else
                        X{jj} = Xloc;
                    end
                    
                    fevals{jj}(ii,1) = states{jj}.funcevals/opts_init.BatchSize;
                    assert(opts_init.BatchSize == size(Xloc,2));
                    fevals{jj}(ii,2) = calc_samples_err(X{jj},theta);
                else
                    RUN_FLAG = 0;
                    break;
                end
            toc()
        end
        
        %Display + Saving 
        if (mod( ii, 600 ) == 0) || (ii == Nsamp) || RUN_FLAG == 0
            fprintf('%d / %d in %f sec (%f sec remaining)\n', ii, Nsamp, toc(ttt), toc(ttt)*Nsamp/ii - toc(ttt) );
            h1=plot_autocorr_samples(X, names);
						disp('Autocorr plot completed')
            h2=plot_fevals(fevals, names);
						disp('Fevals plot completed')
            disp(savestr)
            saveas(h1,figpath1,'pdf');
            saveas(h2,figpath2,'pdf');
            save(savepath);
        end
        ii = ii + 1;
    end
ttt = toc(ttt);
