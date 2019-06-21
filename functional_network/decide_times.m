function sp_times = decide_times(whichPts)

% This determines the spike times to use for calculating interictal
% functional networks

%% Parameters
n_spikes = 1000;

%% File path
locations = comp_nets_files;
main_folder = locations.main_folder;
data_folder = [main_folder,'data/'];
script_folder = [main_folder,'scripts/'];
addpath(genpath(script_folder));
spike_struct_folder = locations.spike_struct_folder;
out_file = [data_folder,'spike_times/times.mat'];

%% Load pt structure
pt = load([data_folder,'spike_structures/pt.mat']);
pt = pt.pt;

if isempty(whichPts)
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
end


for whichPt = whichPts

    %% Get sequence matrix, with bad clusters removed
    seq_matrix = get_seqs(whichPt);
    first_time = min(seq_matrix,[],1);
    
    if isempty(first_time) == 1
        continue
    end

    %% Pick random sample of sequences
    n_total_seqs = length(first_time);
    n_seqs_to_pick = min([n_spikes,n_total_seqs]);
    seq_idx = randi(n_total_seqs,n_seqs_to_pick,1);
    times = first_time(seq_idx)';

    %% Optional plot
    if 0
        figure
        histogram(sp_times/3600)
    end

    %% Add it to the structure
    sp_times(whichPt).times = times;
    sp_times(whichPt).name = pt(whichPt).name;
end

%% Save the structure
save(out_file,'sp_times');
    
end