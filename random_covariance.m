%%Generates a n-d random symmetric covariance matrix with a heavy skew on p
%%dimensions
%%Input = n,p (where n is the number of dimensions and p is the the number
%%of dimensions with high skew)
%%Output = random covariance matrix of size n
function rand_cov = random_covariance(n,p)
%Check for nargin
if nargin < 2
    p = 1; %choose a rand dimension for the skew
end

%rand_cov init
rand_cov = ones(n,1);

ind = ceil(n.*rand(p,1));

%Now apply random skew
for i = 1:length(ind)
    rand_cov(ind(i),1) = rand(1)*10e-3;
end

rand_cov = diag(rand_cov);

end