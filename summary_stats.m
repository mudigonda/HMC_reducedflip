function[summary_stats] = summary_stats(path)
disp('changing directory')
cd(path)
disp('calculating stored mat files')
files = dir('*.mat');

%Init

%Summarizing
for filesidx=1:length(files)
   load(files(filesidx).name)
	 STATES=length(states);
	 disp(filesidx)
   for sampleridx=1:length(states)
       summary_stats{sampleridx}.leap(filesidx,:) = states{sampleridx}.steps.leap'; 
			 summary_stats{sampleridx}.flip(filesidx,:) = states{sampleridx}.steps.flip;
       summary_stats{sampleridx}.stay(filesidx,:) = states{sampleridx}.steps.stay;
   end
end
