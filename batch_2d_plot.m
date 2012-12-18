%%%make subplots of each of the batch's data
function h_plt = batch_2d_plot(Samples)
sz=size(Samples);
h_plt = figure(333);
for i = 1:sz(2)
    subplot(1,sz(2),i)
    plot(squeeze(Samples(1,i,:)),squeeze(Samples(2,i,:)),'r.');    
end
end