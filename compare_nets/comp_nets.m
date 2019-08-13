function comp_nets(whichPts)

%{

%}
do_plots = 0;
which_file = 1;
remove_depth = 0;

%% for table
which_freq = 2;
which_time = 2;

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

all_names = cell(length(whichPts),1);
n_depths = nan(length(whichPts),1);
count = 0;

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
    
    count = count + 1;
    all_names{count} = name;
    
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
    
    %% Get indices of depth electrodes
    is_depth = zeros(length(fun_net.elec_data.electrodes),1);
    for i = 1:length(fun_net.elec_data.electrodes)
        if strcmp(fun_net.elec_data.electrodes(i).type,'D') == 1
            is_depth(i) = 1;
        end
    end
    is_depth = logical(is_depth);
    is_depth_flat = find_row_flat(is_depth);
    if remove_depth == 1
        locs(is_depth,:) = [];
    end
    n_depths(count) = sum(is_depth);
    
    %% Get avg spike coactivation matrix
    if remove_depth == 1
        sp_coa(is_depth_flat,:) = [];
    end
    avg_spike_coa = flatten_or_expand_adj(nanmean(sp_coa,2));
    avg_spike(whichPt).data = avg_spike_coa;
    
    %% Get the 1/distance matrix
    dist_graph(whichPt).A = gen_dist_graph(pt,whichPt);
    dist_lin = flatten_or_expand_adj(dist_graph(whichPt).A);
    
    %% Get spike-distance graph correlation
    [r_dist,p_dist] = corr(dist_lin,nanmean(sp_coa,2),'Type','Spearman');
    avg_spike(whichPt).dist_corr.r = r_dist;
    avg_spike(whichPt).dist_corr.p = p_dist;
    
    %% Get average functional matrices
    for t = 1:2
        for f = 1:5
            
            % Get flattened matrix
            curr_adj = fun_net.adj(f).which_adj(t).data;
            
            % Remove any columns with complex numbers
            [~,c] = find(imag(curr_adj)~=0);
            curr_adj(:,c) = [];
            
            if isempty(c) == 0
                fprintf('%s freq %d time %d had %d columns removed due to being imaginary\n.',...
                    name,f,t,length(unique(c)));
            end
            
            % remove depths
            if remove_depth == 1
                curr_adj(is_depth_flat,:) = [];
            end
            
            % Expand average to adjacency
            avg_fun_coa = flatten_or_expand_adj(mean(curr_adj,2));
            
            % Add it to structure
            avg_fun(whichPt).freq(f).time(t).data = avg_fun_coa;
            
            % Get correlation between full average adj and add to struct
            [r_full,p_full] = corr(nanmean(sp_coa,2),mean(curr_adj,2),'Type','Spearman');
            avg_fun(whichPt).freq(f).time(t).corr_full.r = r_full;
            avg_fun(whichPt).freq(f).time(t).corr_full.p = p_full;
            
            % Get correlation between node strength and add to struct
            [r_ns,p_ns] = corr(nansum(avg_spike_coa,2),nansum(avg_fun_coa,2),'Type','Spearman');
            avg_fun(whichPt).freq(f).time(t).corr_ns.r = r_ns;
            avg_fun(whichPt).freq(f).time(t).corr_ns.p = p_ns;
            
            avg_fun(whichPt).freq(f).name = fun_net.adj(f).name;
            
            
            % Get correlation between distance and functional and add to
            % struct
            
            [r_dist,p_dist] = corr(dist_lin,mean(curr_adj,2),'Type','Spearman');
            avg_fun(whichPt).freq(f).time(t).corr_dist.r = r_dist;
            avg_fun(whichPt).freq(f).time(t).corr_dist.p = p_dist;
            
            % Build a univariate linear regression model
            mdl_uni = fitlm(nanmean(sp_coa,2),nanmean(curr_adj,2));
            avg_fun(whichPt).freq(f).time(t).uni_model.spike.t = mdl_uni.Coefficients.tStat(2);
            avg_fun(whichPt).freq(f).time(t).uni_model.spike.p = mdl_uni.Coefficients.pValue(2);
            
            % Build a linear regression model to see the multivariate
            % contribution of the spike network and distance network to the
            % functional network
            X = [dist_lin,nanmean(sp_coa,2)];
            y = nanmean(curr_adj,2);
            mdl = fitlm(X,y);
            avg_fun(whichPt).freq(f).time(t).multi_model.dist.t = mdl.Coefficients.tStat(2);
            avg_fun(whichPt).freq(f).time(t).multi_model.dist.p = mdl.Coefficients.pValue(2);
            avg_fun(whichPt).freq(f).time(t).multi_model.spike.t = mdl.Coefficients.tStat(3);
            avg_fun(whichPt).freq(f).time(t).multi_model.spike.p = mdl.Coefficients.pValue(3);
            
        end
    end
    
    %% Make plots
    if do_plots == 1
        figure
        set(gcf,'position',[0 100 1300 800])
        [ha,~] = tight_subplot(3,5,[0.03 0.01],[0.04 0.01],[0.03 0.01]);

        %% Show average spike coactivation array
        axes(ha(1))
        imagesc(avg_spike(whichPt).data);
        xticklabels([])
        yticklabels([])
        ylabel('Spike coactivation network')

        %% Show spike node strength of electrodes 
        axes(ha(2))
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        hold on
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,nansum(avg_spike(whichPt).data,2),'filled');
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


                % plot node strength
                scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
                hold on
                scatter3(locs(:,1),locs(:,2),locs(:,3),100,...
                    nansum(avg_fun(whichPt).freq(f).time(t).data,2),'filled');


                title(sprintf('Full: %1.2f, p = %1.3f, ns: %1.2f, p = %1.3f',...
                    avg_fun(whichPt).freq(f).time(t).corr_full.r,...
                    avg_fun(whichPt).freq(f).time(t).corr_full.p,...
                    avg_fun(whichPt).freq(f).time(t).corr_ns.r,...
                    avg_fun(whichPt).freq(f).time(t).corr_ns.p))

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


%% Put together a table with all the correlations
if 0
corr_full = zeros(length(whichPts),2);
corr_ns = zeros(length(whichPts),2);
corr_dist = zeros(length(whichPts),2);
sp_corr_dist = zeros(length(whichPts),2);
count = 0;
for whichPt = whichPts
    count = count + 1;
    freq_name = avg_fun(whichPt).freq(which_freq).name;
    corr_full(count,1) = avg_fun(whichPt).freq(which_freq).time(which_time).corr_full.r;
    corr_full(count,2) = avg_fun(whichPt).freq(which_freq).time(which_time).corr_full.p;
    
    corr_ns(count,1) = avg_fun(whichPt).freq(which_freq).time(which_time).corr_ns.r;
    corr_ns(count,2) = avg_fun(whichPt).freq(which_freq).time(which_time).corr_ns.p;
    
    corr_dist(count,1) = avg_fun(whichPt).freq(which_freq).time(which_time).corr_dist.r;
    corr_dist(count,2) = avg_fun(whichPt).freq(which_freq).time(which_time).corr_dist.p;
    
    sp_corr_dist(count,1) = avg_spike(whichPt).dist_corr.r;
    sp_corr_dist(count,2) = avg_spike(whichPt).dist_corr.p;
    
end

fprintf('For %s frequency and time %d, the correlations are:\n',freq_name,which_time);

T = table(corr_full(:,1),corr_full(:,2),...
    corr_dist(:,1),corr_dist(:,2),sp_corr_dist(:,1),sp_corr_dist(:,2),...
    'VariableNames',...
    {'Full_corr_rho','Full_corr_p','dist_corr_rho',...
    'dist_corr_p','sp_dist_r','sp_dist_p'},'RowNames',all_names);
end

%% Put together a table with the multivariate correlations
for t = 1:2
    for f = 1:5
        tStat = zeros(length(whichPts),3);
        pValue = zeros(length(whichPts),3);
        count = 0;
        for whichPt = whichPts
            count = count + 1;
            freq_name = avg_fun(whichPt).freq(f).name;
            tStat(count,1) = avg_fun(whichPt).freq(f).time(t).multi_model.dist.t;
            pValue(count,1) = avg_fun(whichPt).freq(f).time(t).multi_model.dist.p;

            tStat(count,2) = avg_fun(whichPt).freq(f).time(t).multi_model.spike.t;
            pValue(count,2) = avg_fun(whichPt).freq(f).time(t).multi_model.spike.p;

            tStat(count,3) = avg_fun(whichPt).freq(f).time(t).uni_model.spike.t;
            pValue(count,3) = avg_fun(whichPt).freq(f).time(t).uni_model.spike.p;

        end

        fprintf('For %s frequency and time %d, the multivariate model controlling for distance is:\n',...
            freq_name,which_time);

        T = table(tStat(:,1),pValue(:,1),tStat(:,2),pValue(:,2),tStat(:,3),pValue(:,3),...
            'VariableNames',{'t_dist','p_dist','t_spike','p_spike','t_uni','p_uni'},...
            'RowNames',all_names)

        % Are the pooled t stats different from zeros
        [~,p] = ttest(tStat(:,2));
        if mean(tStat) > 0
            sign_text = 'positive';
        else
            sign_text = 'negative';
        end

        fprintf('For %s and time %d, pooled p value is %1.3f (t is %s).\n\n',...
            freq_name,t,p,sign_text)
    end
end

%% Find highest ns electrode and get distance from nearest SOZ

count = 0;
median_elec_dist = [];
for whichPt = whichPts
    
    count = count + 1;
    
    % electrode locations
    locs = pt(whichPt).new_elecs.locs;
    
    % soz locs
    soz_locs = locs(get_soz_chs(pt,whichPt),:);
    
    % Get distance from each electrode and its closest soz electrode
    closest_soz = zeros(size(locs,1),1);
    for i = 1:length(closest_soz)
        curr_locs = locs(i,:);
        curr_dists = vecnorm(curr_locs-soz_locs,2,2);
        min_dist = min(curr_dists);
        closest_soz(i) = min_dist;
    end
    soz_dist(whichPt).all_elecs = closest_soz;
    soz_dist(whichPt).median_elec = median(closest_soz);
    
    median_elec_dist = [median_elec_dist;median(closest_soz)];
    
    for type = 1:2
        if type == 1
            type_name = 'functional';
            struct_head = avg_fun(whichPt);
            adj = struct_head.freq(which_freq).time(which_time).data;
        elseif type == 2
            type_name = 'spike';
            struct_head = avg_spike(whichPt);
            adj = struct_head.data;
        end
 
        
        % Calculate node strength
        ns = sum(adj,1);
        
        % Calculate eigenvector centrality
        ec = eigenvector_centrality_und(adj);
        
        for m = 1:2
            if m == 1
                meas = ns;
                meas_name = 'ns';
            elseif m == 2
                meas = ec;
                meas_name = 'ec';
            end
            
            % Find the most connected electrode and get its location
            [~,most_con] = max(meas);
            con_loc = locs(most_con,:);
            
            % Get the distance between this electrode and each soz and take
            % the minimum
            dist_con_soz = min(vecnorm(con_loc - soz_locs,2,2));
            
            network(type).name = type_name;
            network(type).measure(m).name = meas_name;
            if isfield(network(type).measure(m),'dist') == 0
                network(type).measure(m).dist = [];
            end
            network(type).measure(m).dist = [network(type).measure(m).dist;dist_con_soz];
            
        end

    end
    
end

%% Paired t test
for i = 1:length(network)
    for j = 1:length(network(i).measure)
        [~,p] = ttest(median_elec_dist,network(i).measure(j).dist);
        fprintf('For %s network, is highest %s elec closer to soz than predicted by chance?:\np = %1.3f\n\n',...
            network(i).name,network(i).measure(j).name,p);
        
    end
end

end