%%Function to plot percent diff of model and sample covariance versus
%%function evaluations
%%Inputs = feval_std,feval_std_persist,feval_redflip,feval_2mom
%%Outputs = handle to plot
function h_feval = plot_fevals(fevals, names)
h_feval = figure(111);

clf();

for ii=1:length(fevals)
    plot(fevals{ii}(2:end,1), fevals{ii}(2:end,2));
    hold on;
end

legend(names);
xlabel('function evaluations');
ylabel('percent error diff');
title('percent error diff of model and sample covariances - diagonoal elements only');
