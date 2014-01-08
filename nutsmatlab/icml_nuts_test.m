%ICML NUTs test
clc; clear; close all;

DataSize=2;
BatchSize=1;
theta = diag(exp(linspace(log(1e-5), log(1), DataSize)));
% temp = mvnrnd(zeros(1, 2), [1, 1.98; 1.98, 4], 500);
Xinit = sqrtm(inv(theta))*randn(DataSize,BatchSize);
% Xinit = Xinit';

M=10000; %Number of Samples I need
% samples = zeros(M,DataSize);
% % nfevals = zeros(M,1);
samples(1,:) = randn(1,DataSize);
% % nfevals(1,:) = 1;
% for ii=2:M
%     [tmp,fevals] = nuts(1.2,@gauss_grad_logp,2,samples(ii-1,:));
%     samples(ii,:)=tmp(end,:); %Pick the second sample because the first one is what we passed in
%     nfevals(ii,:)=fevals+nfevals(ii-1,:);%doing a cummulative sum
%     if mod(ii,100)==0
%         disp(ii)
%     end
% end
[samples,fevals]= nuts(1.2,@gauss_grad_logp,M,samples(1,:));
plot(samples(:, 1), samples(:, 2), 'b.');
title('NUTS on 100D with loglinear covariance');
hold off;
mu = mean(samples)
stddev = sqrt(mean(samples.^2))
correlations = corr(samples)
%Calculate average fevals
avg_fevals = mean(fevals);
sprintf('Average Function evaluations %f',avg_fevals)
acn=500; %acn window
autocorr = zeros(1,acn+1);
for ii = 0:acn   
    autocorr(ii+1) = mean(sum(samples(1:end-ii,:).*samples(1+ii:end,:)));   
end
figure;
plot((1:length(autocorr))*avg_fevals,autocorr,'r');
xlabel('Function Evaluations');
ylabel('Auto correlation');
title('NUTS on 100D with loglinear covariance');
hold on;