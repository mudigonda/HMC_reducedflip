%This script will summarize the auto correlation for all the ranges specified
%Input = path is the path where we want to summarize
function [err] = ac_summary_stats(path)

%first change the path
cd(path)
pwd
max_shift=7000
%then make a search string based on the other variables
%strsearch=
%search for mat files
%ac_files=dir('*-LeapSize-10-*mat');
ac_files=dir('*-Beta-41-*mat')
ac_files
%for each mat file
for ii_files=1:length(ac_files)
%Load data
	disp('Loading another file')
	ac_files(ii_files).name
	load(ac_files(ii_files).name)
%get ac values
	Data{ii_files}.name=ac_files(ii_files).name
	Data{ii_files}.ac = ac;
	Data{ii_files}.avg_fevals=avg_fevals;
	num_steps = 500;
	max_shift = min([max_shift, size(X{1},3)]);
	target_shifts = ceil(logspace(0,log10(max_shift),num_steps));
	target_shifts = target_shifts - min(target_shifts);
	target_shifts = unique(target_shifts);
	Data{ii_files}.target_shifts= target_shifts;
end
%save('~/Data/ac_summary-LS10','Data','ac_files')
save('~/Data/ac_summary-B41','Data','ac_files')
end
