function out = spike_coactivation(whichPt)

% This calculates spike coactivation networks
do_mmn = 0;
coa_method = 1;
time_thresh = 0.15;


%% File path
locations = comp_nets_files;
main_folder = locations.main_folder;
data_folder = [main_folder,'data/'];
script_folder = [main_folder,'scripts/'];
addpath(genpath(script_folder));
spike_struct_folder = locations.spike_struct_folder;
pt_new_locs_file = [data_folder,'spike_structures/pt.mat'];

%% Load files
pt = load([spike_struct_folder,'long_seq']);
pt = pt.pt;


%% Get patient parameters
locs = pt(whichPt).electrodeData.locs(:,2:4);
nchs = size(locs,1);


%% Get sequence matrix, with bad clusters removed
[seq_matrix,all_times_all,all_spikes] = get_seqs(whichPt);
first_time = min(seq_matrix,[],1);

%% Initialize coA array
window = 1800;
nbins = ceil((first_time(end)-first_time(1))/window);

% Get coactivation matrix
if coa_method == 1
    [coA,times,total_ch_counts] = find_seq_coA(window,seq_matrix,nchs);
else
    [coA,times,total_ch_counts] = ...
    find_coactivation(window,time_thresh,all_times_all,all_spikes,nchs);
end

out.name = pt(whichPt).name;
out.fs = pt(whichPt).fs;
out.coA = coA;
out.times = times;
out.elecs = pt(whichPt).electrodeData;
out.seq_matrix = seq_matrix;
out.total_ch_counts = total_ch_counts;

% optional plot of average matrix
avg_coA = nanmean(coA,1);
avg_adj = flatten_or_expand_adj(avg_coA);
ns = sum(avg_adj,1);

if 0
    figure
    subplot(1,3,1)
    imagesc(avg_adj)

    subplot(1,3,2)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
    hold on
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,ns,'filled')

    subplot(1,3,3)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
    hold on
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,sum(total_ch_counts,1),'filled')

end
    
   
if do_mmn == 1
    m = 3;
    %[W,H,D] = nnmf(coA,m);
    [H,W] = pca(coA);
    sz_times = pt(whichPt).newSzTimes(:,1);
    figure
    for i = 1:m
        plot(times,W(:,i))
        hold on
        for j = 1:length(sz_times)
            plot([sz_times(j) sz_times(j)],get(gca,'ylim'),'k--')
        end
    end
    
end
    



end