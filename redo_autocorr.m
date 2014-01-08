function redo_autocorr(path)
%Input - takes all the mat files in the path and recomputes autocorr and saves it back 
MAX_SHIFT=7000
pwd
addpath(pwd)
disp(path)
disp('Changing path')
cd(path)
files=dir('*.mat');
whos
sprintf('Total number of files in this directory %d', length(files))
for filesidx=1:length(files)
	disp(filesidx)
	sprintf('File name is %s',files(filesidx).name)
	load(files(filesidx).name);
	[h1,h2,ac]=plot_autocorr_samples(X,names,avg_fevals,MAX_SHIFT);
	saveas(h1,figpath1,'pdf');
	saveas(h2,figpath2,'pdf');
	disp('saving everything')
	save(files(filesidx).name)
end

end
