function [logp,grad] = gauss_grad_logp(X)
theta = diag(exp(linspace(log(1e-5), log(1), size(X,2))));
% A=inv(theta);
A= theta;
% A = [50.251256, -24.874372; -24.874372, 12.562814];

grad = -X * A;
logp = 0.5 * grad * X';
end