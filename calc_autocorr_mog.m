%Calculate Autocorrelation
%%Calculate autocorr for a bunch of samples
function autocorr = calc_autocorr_mog(Samples,acn,Mu)

if nargin < 2
    acn = 500;
end

%Calculate mean of all the means from the Model
mean_of_models=zeros(size(Mu{1}));
for ii=1:size(Mu,1)
   mean_of_models=mean_of_models+Mu{ii};
end
mean_of_models=mean(mean_of_models);

%Mean subtract the Samples
Samples=bsxfun(@minus,Samples,mean_of_models);

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