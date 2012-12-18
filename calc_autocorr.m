%Calculate Autocorrelation
%%Calculate autocorr for a bunch of samples
function autocorr = calc_autocorr(Samples,acn)

if nargin < 2
    acn = 500;
end

sz = size(Samples);
autocorr = zeros(sz(2),acn+1);
for ii = 0:acn
    for jj = 1:sz(2)%batch id
        autocorr(jj,ii+1) = mean(sum(Samples(:,jj,1:end-ii).*Samples(:,jj,1+ii:end)));
    end
end

%Normalize by the first window value
autocorr = autocorr./repmat(autocorr(:,1),1,acn+1);
autocorr = mean(autocorr(:,1:(end-1)),1);

end