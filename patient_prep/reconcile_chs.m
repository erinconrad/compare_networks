function to_keep = reconcile_chs(new_elecs,old_elecs)
    %% Remove the rows and columns corresponding to channels that don't exist in the functional network
    
    nchs = length(old_elecs.electrodes);
    
    to_keep = zeros(nchs,1);
    for i = 1:nchs
        old_label = old_elecs.electrodes(i).name;
        for j = 1:length(new_elecs.electrodes)
            if strcmp(old_label,new_elecs.electrodes(j).name) == 1
                to_keep(i) = 1;
            end
        end
    end
    
    if sum(to_keep) ~= length(new_elecs.electrodes)
        error('what\n');
    end

end