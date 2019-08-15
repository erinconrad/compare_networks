function get_time_varying_networks(whichPts)

%% Parameters
which_file = 1;

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
    
    %% Get the times for the functional networks
    
    
    
end


end