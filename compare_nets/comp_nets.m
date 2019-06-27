function comp_nets(whichPts)

%{

%}

which_file = 1;

%% File path
locations = comp_nets_files;
main_folder = locations.main_folder;
script_folder = [main_folder,'scripts/'];
data_folder = [main_folder,'data/'];
addpath(genpath(script_folder));
results_folder = [main_folder,'results/'];
adj_folder = [results_folder,'adj/'];

%% Load pt structure
pt = load([data_folder,'spike_structures/pt.mat']);
pt = pt.pt;


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
    
    
    figure
    set(gcf,'position',[0 100 1300 800])
    [ha,~] = tight_subplot(3,5,[0.03 0.01],[0.04 0.01],[0.03 0.01]);
    
    %% Show average spike coactivation array
    axes(ha(1))
    avg_spike_coa = flatten_or_expand_adj(nanmean(sp_coa,2));
    imagesc(avg_spike_coa);
    xticklabels([])
    yticklabels([])
    ylabel('Spike coactivation network')
    
    %% Show spike node strength of electrodes 
    axes(ha(2))
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    hold on
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,sum(avg_spike_coa,2),'filled');
    xticklabels([])
    yticklabels([])
    zticklabels([])
    
    delete(ha(3)); delete(ha(4)); delete(ha(5));
    
    %% For each functional network, show the node strength and correlation with sp network
    count = 5;
    for t = 1:2
        for f = 1:5
            count = count + 1;
            axes(ha(count))
            
            % Get flattened matrix
            curr_adj = fun_net.adj(f).which_adj(t).data;
            
            % Remove any columns with complex numbers
            [~,c] = find(imag(curr_adj)~=0);
            curr_adj(:,c) = [];
            
            if isempty(c) == 0
                fprintf('%s freq %d time %d had %d columns removed due to being imaginary\n.',...
                    name,f,t,length(unique(c)));
            end
            
            % Expand average to adjacency
            avg_fun_coa = flatten_or_expand_adj(mean(curr_adj,2));
            
            % plot node strength
            scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
            hold on
            scatter3(locs(:,1),locs(:,2),locs(:,3),100,sum(avg_fun_coa,2),'filled');
            
            % Get correlation between full average adj
            [r_full,p_full] = corr(nanmean(sp_coa,2),mean(curr_adj,2),'Type','Spearman');
            
            % Get correlation between node strength
            [r_ns,p_ns] = corr(sum(avg_spike_coa,2),sum(avg_fun_coa,2),'Type','Spearman');
            
            title(sprintf('Full: %1.2f, p = %1.3f, ns: %1.2f, p = %1.3f',...
                r_full,p_full,r_ns,p_ns))
            
            xticklabels([])
            yticklabels([])
            zticklabels([])
            
            if f == 1
                if t == 1
                    zlabel('Functional network at spike')
                else
                    zlabel('Non-spike functional network')
                end  
            end
            
            if t == 2
                xlabel(sprintf('%s',fun_net.adj(f).name))
            end
            
        end
    end
    
    
end


end