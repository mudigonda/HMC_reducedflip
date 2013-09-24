function make_figures_fneval_local(LeapSize,epsilon,beta)

HOME=getenv('HOME');
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
    opts_init.LeapSize = 10;
else
    LeapSize=str2num(LeapSize);
    opts_init.LeapSize = LeapSize
end
if nargin<2
    opts_init.epsilon = 1.3;
else
    epsilon =str2num(epsilon);
    opts_init.epsilon = epsilon
end
if nargin<3
    opts_init.beta = .35;
else
    beta    =str2num(beta);
    opts_init.beta = beta
end            

%Model Name

FEVAL_MAX = 5000000
% modelname='100dCircle100'
modelname='2D'
Nsamp = 10000;
opts_init.BatchSize = 100;
opts_init.DataSize = 10;
opts_init.DataSize = 2;
savestr = strcat('ModelName-',modelname,'-LeapSize-',int2str(opts_init.LeapSize),...
    '-epsilon-',int2str(opts_init.epsilon*10),'-Beta-',int2str(opts_init.beta*100)...
    ,'-fevals-',int2str(FEVAL_MAX),'-Nsamp-',int2str(Nsamp)...
    ,'-BS-',int2str(opts_init.BatchSize),'-DS-',int2str(opts_init.DataSize));

savepath = strcat(HOME,'/Data/HMC_reducedflip/10d-MOG-20/',savestr);
figpath1 = strcat(HOME,'/Data/HMC_reducedflip/10d-MOG-20/figures/',savestr,'autocor');
figpath2 = strcat(HOME,'/Data/HMC_reducedflip/10d-MOG-20/figures/',savestr,'autocor-fevals');
mkdir(savepath);
mkdir(figpath1);
mkdir(figpath2);


opts_init.funcevals = 0;

% % % % scaling factor for energy function
% % % % theta = [1,0;0,1e-6];
% % % % % theta = 100; %%circle
% % % %%Crazy ass MOG with crazy ass everything
% % %  x=4;b=-6;
% % %  for ii=1:opts_init.DataSize
% % %    rng(ii);
% % %    J{ii}=diag(exp(linspace(log(1e-6), log(1), opts_init.DataSize)).*rand(1,opts_init.DataSize));
% % %    Mu{ii} = [zeros(opts_init.DataSize-1,1);b+(ii)*x];
% % % % %   Mu{ii}=randn(opts_init.DataSize,1)*(rand(1));
% % %  end
% % % % % 

% % %%2D MOG with n mixtures
% % % x=4;b=-6;
% %  
% % % for ii=1:MOG
% % %    J{ii}=diag(ones(opts_init.DataSize,1));%Unit covariances
% % %    Mu{ii}=[0;b + (ii)*x];%Means
% % % end

%logalpha = zeros(opts_init.DataSize,1);
%W = eye(opts_init.DataSize);
%theta = [W, logalpha];

theta = diag(exp(linspace(log(1e-5), log(1), opts_init.DataSize)));

opts_init.Xinit = sqrtm(inv(theta))*randn( opts_init.DataSize, opts_init.BatchSize );

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
names{ii} = 'forever forward'
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 3;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = []

ii = ii + 1
names{ii} = 'default + ff'
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 3;
opts{ii}.beta = 1;
%Initialize States
states{ii} = [];
%arrays to keep track of samples
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
                if ii == 1 || states{jj}.funcevals < FEVAL_MAX 
                    [Xloc, statesloc] = rf2vHMC( opts{jj}, states{jj},theta);
                    states{jj} = statesloc;
                    if ii > 1                                        
                        X{jj} = cat(3,X{jj}, Xloc);
                    else
                        X{jj} = Xloc;
                    end
                    
                    fevals{jj}(ii,1) = states{jj}.funcevals;
                    assert(opts_init.BatchSize == size(Xloc,2));
                else
                    RUN_FLAG = 0;
                    break;
                end
        end

        %Display + Saving 
        if (mod( ii, 500 ) == 0) || (ii == Nsamp) || RUN_FLAG == 0
            fprintf('%d / %d in %f sec (%f sec remaining)\n', ii, Nsamp, toc(ttt), toc(ttt)*Nsamp/ii - toc(ttt) );

            for jj = 1:length(names)
                disp(names{jj})
                states{jj}
                states{jj}.steps
                states{jj}.steps.leap'
            end


            %Calculate average fevals by taking total fevals at this point
            %and dividing it by the number of samples we have acquired
            sprintf('calculating average fevals')
            for jj=1:length(names)
                avg_fevals{jj}=fevals{jj}(end,1)/size(X{jj},3);
            end
            if exist('Mu','var')
                [h1,h2]=plot_autocorr_samples(X, names,avg_fevals,Mu);
            else
                [h1,h2]=plot_autocorr_samples(X, names,avg_fevals);
            end
                disp('Autocorr plot completed')
%             h2=plot_fevals(fevals, names);
						disp('Fevals plot completed')
            disp(savestr)
            saveas(h1,figpath1,'pdf');
            saveas(h2,figpath2,'pdf');
            save(savepath);
        end
        ii = ii + 1;
    end
ttt = toc(ttt);
