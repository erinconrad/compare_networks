function [seq_matrix,all_times_all,all_spikes] = get_seqs(whichPt)

%% File path
locations = spike_network_files;
main_folder = locations.main_folder;
spike_struct_folder = locations.spike_struct_folder;

%% Load files
pt = load([spike_struct_folder,'long_seq']);
pt = pt.pt;
cluster = load([spike_struct_folder,'cluster']);
cluster = cluster.cluster;

seq_matrix = pt(whichPt).seq_matrix;

% Reorder seizure times if out of order
szTimes = pt(whichPt).newSzTimes;
oldSzTimes = szTimes;
szTimes = sort(szTimes,1);
if isequal(oldSzTimes,szTimes) == 0
    fprintf('WARNING!!! %s seizure times out of order\n',pt(whichPt).name);
end

% Combine nearly equal seizure times
newIdx = 2;
newSzTimes = [];
newSzTimes(1,:) = szTimes(1,:);
for j = 2:size(szTimes,1)
    if abs(szTimes(j,1)-szTimes(j-1,1)) < 10 && ...
            abs(szTimes(j,2)-szTimes(j-1,2))
       newIdx = newIdx - 1; 
       newSzTimes(newIdx,1) = min(szTimes(j,1),szTimes(j-1,1));
       newSzTimes(newIdx,2) = max(szTimes(j,2),szTimes(j-1,2));  
    else
       newSzTimes(newIdx,:) = szTimes(j,:);
    end
    newIdx = newIdx + 1;
end

if isequal(newSzTimes,szTimes) == 0
    fprintf('WARNING!!! %s had duplicate seizure times\n',pt(whichPt).name);
end

%% Remove ties
keep = ones(size(seq_matrix,2),1);
for s = 1:size(seq_matrix,2)
   curr_seq = seq_matrix(:,s);
   nonans = curr_seq(~isnan(curr_seq));
   norepeats = unique(nonans);
   if length(norepeats) < 0.5*length(nonans)
       keep(s) = 0;
   end
end
seq_matrix(:,keep==0) = [];
fprintf(['%s had %d sequences (%1.2f of all sequences) deleted'...
'for having >50 percent ties\n%d sequences remain\n'],...
pt(whichPt).name,sum(keep == 0),sum(keep == 0)/length(keep),sum(keep==1));

%% Remove ictal sequences
all_times = seq_matrix(:);
icTimes = find(any(all_times >= (szTimes(:,1)-repmat(60,size(szTimes,1),1))' ...
    & all_times <= szTimes(:,2)',2));
seq_matrix(icTimes) = nan;
fprintf('Removed %d ictal spikes\n',length(icTimes));

%% Get cluster info
all_times_all = cluster(whichPt).all_times_all; % all spike times
all_spikes = cluster(whichPt).all_spikes; % all spike channels
all_locs = cluster(whichPt).all_locs;
k = cluster(whichPt).k; % the number of clusters
idx = cluster(whichPt).idx; % the cluster index for every spike
C = cluster(whichPt).C; % the centroids of the clusters
bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad


%% Compare number of spikes in cluster array and my data
if sum(sum(~isnan(seq_matrix))) ~= length(all_times_all)
    error('Warning, number of spikes do not align\n');
end



%% Find bad spikes
bad_idx = find(ismember(idx,bad_cluster));

% Nx2 array of bad spikes, showing the channel and time
bad_spikes = [all_spikes(bad_idx),all_times_all(bad_idx)];

%% Get all sequences

new_seq_matrix = seq_matrix;
n_removed = 0;

%% Go through sequence matrix and remove bad spikes
for ich = 1:size(seq_matrix,1)
    % loop across electrodes

    % All spike times for this channel
    spikeTimesCh = seq_matrix(ich,:);

    % Get the bad spikes in that channel
    bad_times_for_ch = bad_spikes(bad_spikes(:,1) == ich,2);

    % Make sure I am finding all of them
    Lia = ismember(spikeTimesCh,bad_times_for_ch);
    if sum(Lia) ~= length(bad_times_for_ch)
        error(sprintf('Did not find all bad spikes for channel %d\n',ich));
    end

    %{
    if sum(Lia) > 0
        fprintf('Removed %d spikes for channel %d\n',sum(Lia),ich)
    end
    %}
    n_removed = n_removed + sum(Lia);

    % Make bad spikes nans
    spikeTimesCh(Lia==1) = nan;
    new_seq_matrix(ich,:) = spikeTimesCh;


end

if n_removed~=length(bad_idx)
    error('Incorrect number of bad spikes removed\n');
end
fprintf('Removed %d spikes for being in bad clusters\n',n_removed);

 %% Remove sequences that have fewer than 5 spikes
removeSeq = zeros(size(new_seq_matrix,2),1);
for s = 1:size(new_seq_matrix,2)
    currSeq = new_seq_matrix(:,s);
    currSeq(isnan(currSeq)) = [];
    if length(currSeq) < 5
        removeSeq(s) = 1;
    end
end

fprintf('Removed %d sequences for now being too short\n',sum(removeSeq));
new_seq_matrix(:,removeSeq==1) = [];


seq_matrix = new_seq_matrix;

fprintf('%d sequences remain\n',size(seq_matrix,2));

all_times_all(bad_idx) = [];
all_spikes(bad_idx) = [];

end