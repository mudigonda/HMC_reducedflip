function[summary_stats] = summary_stats(path)

cd(path)

files = dir('*.mat');
summary_stats=zeros(length(files),9);

for filesidx=1:length(files)
   load(files(filesidx).name)
   disp(filesidx)
   for sampleridx=1:length(states)
       summary_stats(filesidx,(sampleridx-1)*3+ 1) = states{sampleridx}.steps.leap;
       summary_stats(filesidx,(sampleridx-1)*3+ 2) = states{sampleridx}.steps.flip;
       summary_stats(filesidx,(sampleridx-1)*3+ 3) = states{sampleridx}.steps.stay;
   end
end