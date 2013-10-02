HOME=getenv('HOME');
whos
clc;
close all;
opts_init = [];
opts_init.E = @E_gauss;
opts_init.dEdX = @dEdX_gauss;

% make this 1 for more output
opts_init.Debug = 0;
opts_init.LeapSize = 10;
opts_init.epsilon = 1.3;

%Model Name

FEVAL_MAX = 5000000
modelname='2D'
Nsamp = 10000;
opts_init.BatchSize = 10000;
opts_init.DataSize = 2;
opts_init.funcevals = 0;
theta = diag(exp(linspace(log(1e-4), log(1), opts_init.DataSize)));

%opts_init.Xinit = sqrtm(inv(theta))*randn( opts_init.DataSize, opts_init.BatchSize );
opts_init.Xinit = randn( opts_init.DataSize, opts_init.BatchSize );

%Initalize Options
ii = 1;
names{ii} = 'standard';
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 0;
opts{ii}.alpha = 1;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = [];

ii = ii + 1;
names{ii} = 'persist';
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 0;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = [];

ii = ii + 1;
names{ii} = 'reduced flip';
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 1;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = [];

ii = ii + 1;
names{ii} = 'forever forward';
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 3;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = [];

ii = ii + 1;
names{ii} = 'default + ff';
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 3;
opts{ii}.alpha = 1;
%Initialize States
states{ii} = [];
%arrays to keep track of samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = [];

if 0
ii = ii + 1;
names{ii} = 'two momentum';
opts{ii} = opts_init;
opts{ii}.FlipOnReject = 2;
%Initialize States
states{ii} = [];
% arrays to keep track of the samples
X{ii} = zeros(opts{ii}.DataSize,Nsamp);
fevals{ii} = [];
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
        if (mod( ii, 100 ) == 0) || (ii == Nsamp) || RUN_FLAG == 0
            fprintf('%d / %d in %f sec (%f sec remaining)\n', ii, Nsamp, toc(ttt), toc(ttt)*Nsamp/ii - toc(ttt) );

            for jj = 1:length(names)
                disp(names{jj})
                disp(states{jj})
                disp(states{jj}.steps)
                disp(states{jj}.steps.leap')
                disp('most recent covariance')
                xx = X{jj}(:,:,end);
                disp( xx * xx' / size(xx,2))
                disp('total covariance')
                xx = X{jj}(:,:,:);
                xx = reshape(xx, size(xx,1), size(xx,2)*size(xx,3));
                disp( xx * xx' / size(xx,2))
            end
        end
        ii = ii + 1;
    end
ttt = toc(ttt);
