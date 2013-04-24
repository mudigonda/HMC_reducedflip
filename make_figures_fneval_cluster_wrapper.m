%%This script is to get around the large waiting time and limitation on cluster nodes available with no
%%parallelization I will attempt to run all the jobs in serial

%%Leap Size
LeapSize=[1,5,10,15]
Epsilon=[0.8,0.9,1.0,1.1,1.2,1.3]
Beta=[0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55]
for L=1:length(LeapSize)
	for E=1:length(Epsilon)
		for B=1:length(Beta)
			sprintf('###############################################################')
			sprintf('Leap Size %d, Epsilon %f, Beta %f ',LeapSize(L),Epsilon(E),Beta(B))
			sprintf('###############################################################')
			make_figures_fneval_cluster(int2str(LeapSize(L)),int2str(Epsilon(E)),int2str(Beta(B)));
		end
	end
end
