function out = spike_coactivation(whichPt)

% This calculates spike coactivation networks

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
seq_matrix = get_seqs(whichPt);
first_time = min(seq_matrix,[],1);

%% Initialize coA array
window = 1800;
nbins = ceil((first_time(end)-first_time(1))/window);
coA = zeros(nbins,nchs*(nchs-1)/2);
times = zeros(nbins,1);
counts = zeros(nbins,1);
total_ch_counts = zeros(nbins,nchs);

% Loop through and find coactivated channels (which I define as
% channels in the same sequence)
for tt = 1:nbins

    % Get appropriate times
    curr_times = [first_time(1) + (tt-1)*window, min(first_time(1) + tt*window,first_time(end))];
    times(tt) = curr_times(2);

    % Get appropriate spikes
    seq_idx = (first_time >= curr_times(1) & first_time <= curr_times(2));
    seqs = seq_matrix(:,seq_idx);

    coA_temp = zeros(nchs,nchs);

    % Loop through sequences in the appropriate times
    for i = 1:size(seqs,2)

        curr_seq = seq_matrix(:,i);
        seq_chs = find(~isnan(curr_seq));
        counts(tt) = counts(tt) + length(seq_chs);
        total_ch_counts(tt,seq_chs) = total_ch_counts(tt,seq_chs) + 1;

        for j = 1:length(seq_chs)
            for k = j+1:length(seq_chs)
                coA_temp(seq_chs(j),seq_chs(k)) = coA_temp(seq_chs(j),seq_chs(k)) + 1;
                coA_temp(seq_chs(k),seq_chs(j)) = coA_temp(seq_chs(k),seq_chs(j)) + 1;
            end
        end
    end

    % Normalize by dividing by sum
    coA_temp = coA_temp/sum(sum(coA_temp));

    % Flatten the adjacency matrix
    adj_flat = zeros(nchs*(nchs-1)/2,1);
    flat_idx = 0;

    for i = 1:nchs
        for j = 1:i-1
            flat_idx = flat_idx + 1;
            adj_flat(flat_idx) = coA_temp(j,i);
        end
    end


    % Add to output
    coA(tt,:) = adj_flat;
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
    
   
    
    



end