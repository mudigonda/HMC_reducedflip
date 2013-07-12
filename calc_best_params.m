%This script will take from the list of experiments and calculate the best params and plot that figure separately
%I could have done this easier if I had all the mat files also append and write to a single file the final values
%but I didn't. Such is life.

function calc_best_params(path,projpath)
if nargin <1
	path='/clusterfs/cortex/scratch/mayur/HMC_reducedflip/2d/'
end

if nargin<2
	projpath='/global/home/users/mudigonda/HMC_reducedflip/'
end

HOME=getenv('HOME');

cd(path);
files=dir('*.mat')

min_std=zeros(length(files),1);
min_std_per=zeros(length(files),1);
min_red_flip=zeros(length(files),1);
min_2mom=zeros(length(files),1);

No_of_samples = 100;

%Loop over 
for i=1:length(files)
	sprintf('This is the i %d',i)
	load(files(i).name);
	min_std(i)=ac{1}(1,No_of_samples)*avg_fevals{1};
	min_std_per(i)=ac{2}(1,No_of_samples)*avg_fevals{2};
	min_red_flip(i)=ac{3}(1,No_of_samples)*avg_fevals{3};
%	min_2mom(i)=fevals{4}(end,2);
end

[min_std_val,min_std_idx]=min(min_std)
[min_per_val,min_per_idx]=min(min_std_per)
[min_redflip_val,min_redflip_idx]=min(min_red_flip)
%[min_2_mom_val,min_2_mom_idx]=min(min_2mom)

load(files(min_std_idx).name);
sprintf('standard model fname %s',files(min_std_idx).name)
std_plot_ac = ac{1};
avg_fevals_final{1}= avg_fevals{1};

load(files(min_per_idx).name);
std_per_plot_ac = ac{2};
sprintf('standard persistent fname %s',files(min_per_idx).name)
avg_fevals_final{2}= avg_fevals{2};

load(files(min_redflip_idx).name);
redflip_plot_ac = ac{3};
sprintf('Reduced flip model fname %s',files(min_redflip_idx).name)
avg_fevals_final{3}= avg_fevals{3};

%load(files(min_2_mom_idx).name);
%mom_plot_ac =ac{4};
%sprintf('Two momentum model fname %s',files(min_2_mom_idx).name)
% avg_fevals_final{4}= avg_fevals{4};

cd(projpath);
%making new cell arrays that reflect the appropriate data
clear X;
clear ac;
ac{1}=std_plot_ac;
ac{2}=std_per_plot_ac;
ac{3}=redflip_plot_ac;
%ac{4}=mom_plot_ac;


colorlist=['r','g','b','k','m','y'];
h_scaled=figure(222);
clf();
for ii=1:length(ac)
    plot(ac{ii},colorlist(ii));
    hold on;
end

legend(names);
xlabel('Auto correlation windows');
ylabel('Auto correlation values');

h_notscaled=figure(333);
clf();

for ii=1:length(ac)
    plot((1:length(ac{ii}))*avg_fevals_final{ii},ac{ii},colorlist(ii));
    hold on;
end
legend(names);
xlabel('Auto correlation windows scaled by avg fevals');
ylabel('Auto correlation values');


saveas(h_scaled,strcat(HOME,'/ac_scaled_final.pdf'))
saveas(h_notscaled,strcat(HOME,'/ac_final.pdf'))
