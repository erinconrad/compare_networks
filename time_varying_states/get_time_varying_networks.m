function get_time_varying_networks(whichPts)

%% Parameters
which_file = 1;
which_fun_time = 2; % 2 is the non-spike time
which_fun_freq = 2; 
gamma = 1;

%% File path
locations = comp_nets_files;
main_folder = locations.main_folder;
script_folder = [main_folder,'scripts/'];
data_folder = [main_folder,'data/'];
addpath(genpath(script_folder));
addpath(genpath(locations.BCT));
results_folder = [main_folder,'results/'];
adj_folder = [results_folder,'adj/'];

%% Load pt structure
pt = load([data_folder,'spike_structures/pt.mat']);
pt = pt.pt;


%% Decide which patients to do
if isempty(whichPts) == 1
    listing = dir(adj_folder);
    for i = 4:length(listing)
        curr_name = listing(i).name;
        if exist([adj_folder,curr_name,'/adj_',curr_name,'.mat'],'file') ~= 0
            % We have the adjacency matrix for this patient, now find the
            % corresponding number
            
            for j = 1:length(pt)
                temp_name = pt(j).name;
                if strcmp(temp_name,curr_name) == 1
                    whichPts = [whichPts,j];
                    break;
                end
            end
        end
    end
end

% Loop over patients
for whichPt = whichPts
    
    %% Get spike coactivation network
    sp_net = spike_coactivation(whichPt);
    
    %% Get functional network
    name = pt(whichPt).name;
    pt_folder = [adj_folder,name,'/'];
    if which_file == 1
        adj_file = [pt_folder,'adj_',sprintf('%s',name),'.mat'];
    else
        adj_file = [pt_folder,'adj_',sprintf('%s',name),'_nocar.mat'];
    end
    out = load(adj_file);
    fun_net = out.out;
    fun_times = fun_net.times;
    
    % Get the specific functional network time (spike or no) and frequency
    % I want
    fun_adj_flat = fun_net.adj(which_fun_freq).which_adj(which_fun_time).data;
    
    % Remove any columns with complex numbers AND remove them from the
    % times
    [~,c] = find(imag(fun_adj_flat)~=0);
    fun_adj_flat(:,c) = [];
    fun_times(c,:) = [];
   

    if isempty(c) == 0
        fprintf('%s freq %d time %d had %d columns removed due to being imaginary\n.',...
            name,f,t,length(unique(c)));
    end
    
    %% Reconcile electrodes
    fun_elecs = fun_net.elec_data.electrodes;
    sp_elecs = sp_net.elecs.electrodes;
    
    [to_remove,to_remove_flat] = reconcile_chs(fun_elecs,sp_elecs);
    
    % Sanity check
    if size(sp_net.coA,2) ~= length(to_remove_flat)
        error('what\n');
    end
    
    % Remove the appropriate rows from the spike coA 
    sp_coa = sp_net.coA';
    sp_coa(to_remove_flat,:) = [];
    if size(sp_coa,1) ~= size(fun_net.adj(1).which_adj(1).data,1) 
        error('what\n');
    end
    
    % Get locs of spike electrodes, removing those not in functional
    % network
    sp_locs = sp_net.elecs.locs;
    sp_locs(to_remove,:) = [];
    sp_locs = sp_locs(:,2:4);
    
    % Get locs of functional electrodes
    fun_locs = fun_net.elec_data.locs;
    if isequal(sp_locs,fun_locs) == 0
        error('what\n');
    end
    
    locs = sp_locs;
    
    %% Bin the functional network times into the same bins as the spike network times
    
    
    sp_times = sp_net.times;
    n_sp_bins = size(sp_times,1);
    
    % Define edges (the starting index of each spike network bin, and then
    % the end time of the last bin as the closing edge)
    edges = [sp_times(:,1);sp_times(end,2)]; 
    
    % Get the indices of the spike bins tht contain the functional times
    fun_bin_assign = discretize(fun_times,edges);
    
    %% For each spike network bin, average the functional times in that bin
    % It may be inappropriate to average here. THINK ABOUT THIS.
    
    % Initialize rebinned functional networks
    fun_re_bin = zeros(size(fun_adj_flat,1),size(sp_times,1));
    
    for i = 1:n_sp_bins
        
        % Identify the functional times that fall within these bins.
        fun_in_bin = fun_bin_assign == i;
        
        % Get the matrix of flattened functional networks that fall within
        % these bins
        matrix_in_bin = fun_adj_flat(:,fun_in_bin);
        
        % Average these flattened networks
        avg_mat_in_bin = mean(matrix_in_bin,2);
        fun_re_bin(:,i) = avg_mat_in_bin;
        
    end
    
    % Now I have two equal size matrices, aligned appropriately in channels
    % and in time, one representing the functional coherence network and
    % one representing the spike network.
    
    %% Define a single vector of times
    both_times = mean(sp_times,2);
    
    %% Get configuration similarity matrix for both spike and functional networks
    cs_spike = configuration_similarity(sp_coa);
    cs_fun = configuration_similarity(fun_re_bin);
    
    %% Detect communities for both spike and functional networks
    %{
    Need to optimize gamma for each patient and each network
    %}
    M_spike = get_community_affiliations(cs_spike,gamma,0);
    M_fun = get_community_affiliations(cs_fun,gamma,0);
    
    % Take the dot product
    dp = M_spike.*M_fun;
    
    % Plot the community affiliations
    if 1
        figure
        stem(both_times,M_spike)
        hold on
        stem(both_times,M_fun)
        legend
        title(sprintf('Dot product: %d',dp))
    end
    
end


end