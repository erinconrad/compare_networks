function [coA,times,total_ch_counts] = find_seq_coA(window,seq_matrix,nchs)

first_time = min(seq_matrix,[],1);
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

end