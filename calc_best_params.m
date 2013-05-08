%This script will take from the list of experiments and calculate the best params and plot that figure separately
%I could have done this easier if I had all the mat files also append and write to a single file the final values
%but I didn't. Such is life.

path='/clusterfs/cortex/scratch/mayur/HMC_reducedflip/2d/'

cd(path);
files=dir('*10-6*feval*.mat')

min_std=zeros(length(files),1);
min_std_per=zeros(length(files),1);
min_red_flip=zeros(length(files),1);
min_2mom=zeros(length(files),1);


%Loop over 
for i=1:length(files)
	load(files(i).name);
	min_std(i)=fevals{1}(end,2);
	min_std_per(i)=fevals{2}(end,2);
	min_red_flip(i)=fevals{3}(end,2);
%	min_2mom(i)=fevals{4}(end,2);
end

[min_std_val,min_std_idx]=min(min_std)
[min_per_val,min_per_idx]=min(min_std_per)
[min_redflip_val,min_redflip_idx]=min(min_red_flip)
%[min_2_mom_val,min_2_mom_idx]=min(min_2mom)

load(files(min_std_idx).name);
sprintf('standard model fname %s',files(min_std_idx).name)
std_plot_fevals = fevals{1};
X_plot = X{1};
load(files(min_per_idx).name);
std_per_plot_fevals = fevals{2};
sprintf('standard persistent fname %s',files(min_per_idx).name)
X_per_plot = X{2};
load(files(min_redflip_idx).name);
redflip_plot_fevals = fevals{3};
sprintf('Reduced flip model fname %s',files(min_redflip_idx).name)
X_red_plot = X{3};
%load(files(min_2_mom_idx).name);
%mom_plot_fevals =fevals{4};
%sprintf('Two momentum model fname %s',files(min_2_mom_idx).name)
%X_2mom_plot = X{4};

cd('/global/home/users/mayur/HMC_reducedflip')
%making new cell arrays that reflect the appropriate data
clear X;
clear fevals;
fevals{1}=std_plot_fevals;
fevals{2}=std_per_plot_fevals;
fevals{3}=redflip_plot_fevals;
%fevals{4}=mom_plot_fevals;

X{1}=X_plot;
X{2}=X_per_plot;
X{3}=X_red_plot;
%X{4}=X_2mom_plot;

h_final = plot_fevals(fevals,names)
h_ac_final = plot_autocorr_samples(X,names)

saveas(h_final,'/global/home/users/mayur/final.pdf')
saveas(h_ac_final,'/global/home/users/mayur/ac_final.pdf')
