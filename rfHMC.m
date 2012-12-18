% HMC sampler with reduced flip

% X is the current sample
% state holds the complete state of the sampler, and can be passed back in to resume
%   at the same location
% varargin contains arguments that are passed through to E( X, varargin ) and dEdX( X, varargin )
% opts contains options for the sampler

function [X, state] = rfHMC( opts, state, varargin )
    
    % load parameters
    f_E = getField( opts, 'E', 0 );
    f_dEdX = getField( opts, 'dEdX', 0 );
        
    dut = 0.5; % fraction of the momentum to be replaced per unit time
    epsilon = getField( opts, 'epsilon', 0.1 );
    beta = 1 - exp( log( dut ) * epsilon );  % replace a fraction dut of the momentum
    beta = getField( opts, 'beta', beta );
    
    % this controls whether or not to use the reduced flip mode
    % default (0) is reduced flip mode
    flip_on_rej = getField( opts, 'FlipOnReject', 0 );

    
    T = getField( opts, 'T', 1000 );
    szb = getField( opts, 'BatchSize', 1 );
    %szb = 1; % for now, this is always 1...
    szd = getField( opts, 'DataSize', 10 );
    debug = getField( opts, 'Debug', 0 );
    
    % initialize the state variable if not already initialized
    if isempty(state)        
        % the random seed reset will make experiments exactly repeatable
        %reset(RandStream.getDefaultStream);
        state.X = randn( szd, szb ); %rows = dims, cols = no.of.samples
        state.V = randn( szd, szb );
        state.X(:) = 0;
        state.X(1,:) = 1;
        % state.steps provides counters for each kind of transition
        state.steps = [];
        state.steps.leap = 0;
        state.steps.flip = 0;
        state.steps.rej = 0;
        state.steps.total = 0;
        state.funcevals = 0;
        
        % populate the energy and gradient at the initial position
        state.dEdX = f_dEdX( state.X, varargin{:} );
        state.E = f_E( state.X, varargin{:} );
    end


    for t = 1:T % iterate over update steps
        if debug
            fprintf( '.' );
        end
        
        % run a single Langevin dynamics step.
        % TODO: make this a variable number of leapfrog steps
        V0 = state.V;
        X0 = state.X;
        E0 = state.E;
        dEdX0 = state.dEdX;
        Vhalf = V0 - epsilon/2 * dEdX0;
        X1 = X0 + epsilon * Vhalf;
        dEdX1 = f_dEdX( X1, varargin{:} );
        E1 = f_E( X1, varargin{:} );
        V1 = Vhalf - epsilon/2 * dEdX1;
        
        % accept or reject the new state
        delta_E = E1 - E0 + 1/2*sum(V1.^2) - 1/2*sum(V0.^2);
        p_leap = exp( -delta_E );
        p_leap(p_leap>1) = 1;
        % p_leap is the probability of accepting the forward transition
        % compare against a random number to decide whether to accept
        rnd_cmp = rand(1,szb);
        gd = (rnd_cmp < p_leap);
        % update the counter
        state.steps.leap = state.steps.leap + sum(gd);
        % update the current state for the samples where the forward transition
        % was accepted
        state.X(:,gd) = X1(:,gd);
        state.V(:,gd) = V1(:,gd);
        state.E(gd) = E1(gd);
        state.dEdX(:,gd) = dEdX1(:,gd);
        % bd indexes the samples for which the forward transition was rejected
        bd = rnd_cmp > p_leap;
        
        % if the state was rejected, figure out whether or not to flip the momentum
        p_leap_rev = zeros(1,szb);
        if sum(bd) > 0
            % run the rejected states 1 step backwards, for the momentum flip comparison
            V0 = -state.V(:,bd);
            X0 = state.X(:,bd);
            E0 = state.E(:,bd);
            dEdX0 = state.dEdX(:,bd);
            Vhalf = V0 - epsilon/2 * dEdX0;
            X1 = X0 + epsilon * Vhalf;
            dEdX1 = f_dEdX( X1, varargin{:} );
            E1 = f_E( X1, varargin{:} );
            V1 = Vhalf - epsilon/2 * dEdX1;
            delta_E = E1 - E0 + 1/2*sum(V1.^2) - 1/2*sum(V0.^2);
            p_leap_rev(:,bd) = exp( -delta_E );
        end
        p_leap_rev(p_leap_rev>1) = 1;
        % p_flip is the probability of reversing the momentum
        p_flip = p_leap_rev - p_leap;
        p_flip(p_flip < 0) = 0;
        if flip_on_rej == 0
            % if we're not using the reduced momentum flip mode, then always flip
            % the momentum on a rejection
            p_flip(:) = 1 - p_leap;
        end

        % reverse the momentum with probability p_flip
        % note that we use the same random number rnd_cmp
        % used above to decide whether to perform the forward transition
        gd = ((rnd_cmp > p_leap) & (rnd_cmp < (p_leap + p_flip)));
        state.V(:,gd) = -state.V(:,gd);        
        % update the counter
        state.steps.flip = state.steps.flip + sum(gd);

        gd = (rnd_cmp > (p_leap + p_flip));
        % update the counter for the case where the forward transition was rejected but the
        % momentum was not flipped
        state.steps.rej = state.steps.rej + sum(gd);
            
        % slightly randomize the momentum
        N = randn( szd, szb );
        state.V  = real(sqrt(1 - beta)) * state.V + sqrt(beta) * N; % numerical errors if beta == 1
        
        state.steps.total = state.steps.total + 1;    
    end
    
    X = state.X;
end

% to process the fields in our options structure
% this function taken from Mark Schmidt's minFunc
% http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
function [v] = getField(options,opt,default)
options = toUpper(options); % make fields case insensitive
opt = upper(opt);
if isfield(options,opt)
    if ~isempty(getfield(options,opt))
        v = getfield(options,opt); 
    else
        v = default;
    end
else
    v = default;
end
end
function [o] = toUpper(o)
if ~isempty(o)
    fn = fieldnames(o);
    for i = 1:length(fn)
        o = setfield(o,upper(fn{i}),getfield(o,fn{i}));
    end
end
end


