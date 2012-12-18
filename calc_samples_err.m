%%Calculate Error between samples and model
%%Input = samples, model
%%Output = err
function err = calc_samples_err(samples, model)
%To prevent NaNs
model = model + eps;

%Find the size of the samples
sz = size(samples);

err = zeros(sz(2),1);
diag_inv_model = diag(inv(model));

%Transpose for making it nxp
%Looping over each batch
for i=1:sz(2)
    samples_i = squeeze(samples(:,i,:));
    cov_samples = diag(cov(samples_i'));
    
    %Err
    err(i,1) = mean(abs(cov_samples-diag_inv_model)./diag_inv_model);
end
err = mean(err);
end