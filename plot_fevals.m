%%Function to plot percent diff of model and sample covariance versus
%%function evaluations
%%Inputs = feval_std,feval_std_persist,feval_redflip,feval_2mom
%%Outputs = handle to plot
function h_feval = plot_fevals(fevalsstandard,fevalsstandard_persist,fevalsreduced_flip...
    ,feval_2_momentum)
h_feval = figure(111);
if nargin < 4
    plot([fevalsstandard(2:end,1),fevalsstandard_persist(2:end,1),...
    fevalsreduced_flip(2:end,1)],...
    [fevalsstandard(2:end,2),fevalsstandard_persist(2:end,2),...
        fevalsreduced_flip(2:end,2)]);
    legend('standard','standard, persistent mom','reduced flip','two momentum');
    xlabel('function evaluations');
    ylabel('percent error diff');
    title('percent error diff of model and sample covariances - diagonoal elements only');
else
    plot([fevalsstandard(2:end,1),fevalsstandard_persist(2:end,1),...
    fevalsreduced_flip(2:end,1),feval_2_momentum(2:end,1)],...
    [fevalsstandard(2:end,2),fevalsstandard_persist(2:end,2),...
        fevalsreduced_flip(2:end,2),feval_2_momentum(2:end,2)]);
    legend('standard','standard, persistent mom','reduced flip','two momentum');
    xlabel('function evaluations');
    ylabel('percent error diff');
    title('percent error diff of model and sample covariances - diagonoal elements only');
end

end