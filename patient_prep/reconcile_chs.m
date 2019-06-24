function [to_remove,to_remove_flat] = reconcile_chs(fun_elecs,sp_elecs)
    
    to_keep = zeros(length(sp_elecs),1);
    
    % Loop through sp_elecs and identify which ones are in fun_elecs
    for i = 1:length(sp_elecs)
        for j = 1:length(fun_elecs)
            if strcmp(sp_elecs(i).name,fun_elecs(j).name) == 1
                to_keep(i) = 1;
            end
        end
    end
    
    to_remove = ~to_keep;
    
    % Get the rows of the spike adjacency matrix to remove
    to_remove_flat = find_row_flat(to_remove);
    
    % Compare order of remaining electrodes
    fun_elec_names = cell(length(fun_elecs),1);
    sp_elec_names = cell(length(sp_elecs),1);
    
    for i = 1:length(fun_elecs)
        fun_elec_names{i} = fun_elecs(i).name;
    end
    
    for i = 1:length(sp_elecs)
        sp_elec_names{i} = sp_elecs(i).name;
    end
    
    sp_elec_names(to_remove) = [];
    
    if isequal(fun_elec_names,sp_elec_names) == 0
        error('what\n');
    end
    
    if 0
        table(fun_elec_names,sp_elec_names)
    end

end