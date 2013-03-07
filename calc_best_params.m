%This script will take from the list of experiments and calculate the best params and plot that figure separately
%I could have done this easier if I had all the mat files also append and write to a single file the final values
%but I didn't. Such is life.

path='/clusterfs/cortex/scratch/mayur/HMC_reducedflip/2d'

cd(path);
files=dir('*.mat')

min_std=zeros(length(files),1);
min_std_per=zeros(length(files),1);
min_red_flip=zeros(length(files),1);
min_2mom=zeros(length(files),1);


%Loop over them
for i=1:length(files)
	load(files(i).name);
	min_std(i)=fevalsstandard(end,2);
	min_std_per(i)=fevalsstandard_persist(end,2);
	min_red_flip(i)=fevalsreduced_flip(end,2);
	min_2mom(i)=feval_2_momentum(end,2);
end

[min_std_val,min_std_idx]=min(min_std)
[min_per_val,min_per_idx]=min(min_std_per)
[min_redflip_val,min_redflip_idx]=min(min_red_flip)
[min_2_mom_val,min_2_mom_idx]=min(min_2mom)

load(files(min_std_idx).name);
std_plot_fevals = fevalsstandard;
X_plot = Xstandard;
load(files(min_per_idx).name);
std_per_plot_fevals = fevalsstandard_persist;
X_per_plot = Xstandard_persist;
load(files(min_redflip_idx).name);
redflip_plot_fevals = fevalsreduced_flip;
X_red_plot = Xreduced_flip;
load(files(min_2_mom_idx).name);
mom_plot_fevals =feval_2_momentum;
X_2mom_plot = X_2_momentum;

cd('/global/home/users/mayur/HMC_reducedflip')
h_final = plot_fevals(std_plot_fevals,std_per_plot_fevals,redflip_plot_fevals,mom_plot_fevals);
h_ac_final = plot_autocorr_samples(X_plot, X_per_plot, X_red_plot, X_2mom_plot);

saveas(h_final,'/global/home/users/mayur/final.pdf')
saveas(h_ac_final,'/global/home/users/mayur/ac_final.pdf')
