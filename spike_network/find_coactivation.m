function [coA,times,counts] = find_coactivation(window,time_thresh,all_times_all,all_spikes,nchs)

nbins = ceil((all_times_all(end)-all_times_all(1))/window);
coA = zeros(nbins,nchs*(nchs-1)/2);
times = zeros(nbins,2);
counts = zeros(nbins,nchs);

% Loop through and find coactivated channels
    for tt = 1:nbins
        
        % Get appropriate times
        curr_times = [all_times_all(1) + (tt-1)*window, min(all_times_all(1) + tt*window,all_times_all(end))];
        times(tt,:) = curr_times;
        
        % Get appropriate spikes
        sp_idx = find(all_times_all >= curr_times(1) & all_times_all <= curr_times(2));
        sp_times = all_times_all(sp_idx);
        sp_chs = all_spikes(sp_idx);
        
        coA_temp = zeros(nchs,nchs);
        
        for i = 1:nchs
            counts(tt,i) = sum(sp_chs==i);
        end
        
        % Loop through spikes
        for i = 1:length(sp_times)
            for j = i+1:length(sp_times)
                
                % if the two spikes are close enough, add to coactivation
                if sp_times(j) - sp_times(i) < time_thresh
                    coA_temp(sp_chs(i),sp_chs(j)) = coA_temp(sp_chs(i),sp_chs(j)) + 1;
                    coA_temp(sp_chs(j),sp_chs(i)) = coA_temp(sp_chs(j),sp_chs(i)) + 1;
                else
                    % if j is too far from i then no subsequent j will be
                    % closer, so need to move to the next i
                    break
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