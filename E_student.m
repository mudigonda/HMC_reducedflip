% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)
% copyright 2010 Jascha Sohl-Dickstein

function E = E_student( X, theta ) %x, phi, nu )

        X = [X;ones(1,size(X,2))];

%    E = (nu+1)/2 * sum(log( 1 + ((phi*x).^2)/nu  ));
    
        % product of student-t experts
        ndims = size(X, 1);
        nbatch = size(X, 2);
        nexperts = prod(size(theta)) / (ndims+1);
        W = reshape( theta, [nexperts, ndims+1] );
        alpha = W(:,ndims+1);
        W = W(:, 1:ndims);
        
        ff = W * X;
        ff2 = ff.^2;
        off2 = 1 + ff2;
        lff2 = log( off2 );
        %        alff2 = diag(alpha) * lff2;
        alff2 = bsxfun( @times, lff2, alpha );

        E = sum( alff2, 1 );

        if 0 %nbatch > 500
        [as, ord] = sort(alpha);
        figure(776);
        display_receptive_fields( W );
        title( '\Phi' );
        figure(777);
        display_receptive_fields( W(ord,:) );
        title( '\Phi sorted by alpha' );
        figure(778);
        subplot( 2,1,1 );
        plot( as );
        title( '\alpha sorted' );
        subplot( 2,1,2 );
        plot( sqrt(sum( W(ord,:).^2, 2 )) );
        title( '| \Phi |_2 sorted by alpha' );
        drawnow;
        end