%%Function to calculate loglikelihood of samples
function lle = calculate_lle(Samples,J,Mu)

%%Loop over all J
logsum=0;
for ii=1:size(J,1)
    mean_subtr_X=bsxfun(@minus,Samples,Mu{ii});
    %We can then reshape all the batches into one matrix and then perform
    %the operation
    reshape_X= reshape(mean_subtr_X,[size(Samples,1),size(Samples,2)*size(Samples,3)]);
    logsum = logsum + log(exp(- reshape_X'*inv(J{ii})*reshape_X));
end

lle= mean(logsum);

end