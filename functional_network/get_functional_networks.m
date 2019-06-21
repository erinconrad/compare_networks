function get_functional_networks(whichPt)

%{
To do:
- decide if I want to do pre-whiten step
- notch filter? HP (5 Hz) and LP (115 Hz) filters?

%}

%% Parameters
do_pre_whiten = 0;
do_notch = 0;
do_hp = 0;
do_lp = 0;
window = 1;
non_spike_shift = -10; % 10 seconds before the spike
hp_freq = 1;
lp_freq = 115;


freq_bands = [0 256;... %broadband   
    5 15;... %alpha/theta
    15 25;... %beta
    30 40;... % low gamma
    95 105;... % high gamma
    ]; 


%% File path
locations = comp_nets_files;
main_folder = locations.main_folder;
data_folder = [main_folder,'data/'];
script_folder = [main_folder,'scripts/'];
addpath(genpath(script_folder));
spike_struct_folder = locations.spike_struct_folder;
results_folder = [main_folder,'results/'];
adj_folder = [results_folder,'adj/'];
loginname = 'erinconr';
pwname = locations.pwfile;
mscohere_folder = locations.mscohere_folder;
addpath(genpath(mscohere_folder));


%% Load pt structure
pt = load([data_folder,'spike_structures/pt.mat']);
pt = pt.pt;

%% Load times structure
sp_times = load([data_folder,'spike_times/times.mat']);
sp_times = sp_times.sp_times(whichPt).times;

%% Get pt info
name = pt(whichPt).name;
fs = pt(whichPt).fs;
nchs = length(pt(whichPt).new_elecs.electrodes);
pt_folder = [adj_folder,name,'/'];
if exist(pt_folder,'dir') == 0
    mkdir(pt_folder)
end
out_file = [pt_folder,'adj_',sprintf('%s',name),'.mat'];

%% Get spike times
n_spikes = length(sp_times);

%% Get spike and non-spike windows
sp_windows = [sp_times - repmat(window/2,length(n_spikes),1),...
    sp_times + repmat(window/2,length(n_spikes),1)];

%% Attempt to load existing structure to restart partway through

% Initialize if doesn't exist
if exist(out_file,'file') == 0
    out.name = name;
    out.fs = fs;
    out.elec_data = pt(whichPt).electrodeData;
    out.adj(1).name = 'broadband';
    out.adj(1).band = freq_bands(1,:);
    out.adj(2).name = 'alpha_theta';
    out.adj(2).band = freq_bands(2,:);
    out.adj(3).name = 'beta';
    out.adj(3).band = freq_bands(3,:);
    out.adj(4).name = 'low_gamma';
    out.adj(4).band = freq_bands(4,:);
    out.adj(5).name = 'high_gamma';
    out.adj(5).band = freq_bands(5,:);
    for i = 1:length(out.adj)
        out.adj(i).which_adj(1).data = zeros(nchs*(nchs-1)/2,n_spikes);
        out.adj(i).which_adj(2).data = zeros(nchs*(nchs-1)/2,n_spikes);
        
        out.adj(i).which_adj(1).name = 'spike';
        out.adj(i).which_adj(2).name = 'not_spike';
    end
    
    first_zeros = 1;
    
else
    
    % Load the file
    out = load(out_file);
    out = out.out;
    
    % Find the first column of all zeros. That is the spike to start with.
    sum_columns = sum(out.adj(5).which_adj(2).data,1);
    first_zeros = find(sum_columns == 0,1,'first');
    
end

%% Build filters
notch_filter = designfilt('bandstopiir','FilterOrder',6, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
hp_filter = designfilt('highpassiir','FilterOrder',6, ...
               'PassbandFrequency',hp_freq,'PassbandRipple',0.2, ...
               'SampleRate',fs);
lp_filter = designfilt('lowpassiir','FilterOrder',6, ...
               'PassbandFrequency',lp_freq,'PassbandRipple',0.2, ...
               'SampleRate',fs);

% Loop through the spikes
for s = first_zeros:size(sp_windows,1)
    
    fprintf('Doing spike %d of %d\n',s,size(sp_windows,1));
    
    % Loop through 2 conditions of spike and not spike
    for shift_idx = 1:2
        
        if shift_idx == 1
            shift = 0;
        elseif shift_idx == 2
            shift = non_spike_shift;
        end
        
        % Establish the current time window
        curr_window = sp_windows(s,:) + shift;
        
        % Convert the time window to points
        window_indices = round(curr_window*fs);
        indices = window_indices(1):window_indices(2);
        
        %% Get EEG signal
        ieeg_name = pt(whichPt).ieeg_name;
        data = get_eeg_data(ieeg_name,indices,loginname,pwname);
        values = data.values;
        ch_labels = data.ch_labels(:,1);
        
        %% Get channels to ignore
        ignore = zeros(length(ch_labels),1);
    
        for i = 1:length(ch_labels)
            ch_labels{i} = ieeg_ch_parser(ch_labels{i});
            for j = 1:length(pt(whichPt).ignore.names)
                if strcmp(pt(whichPt).ignore.names(j),ch_labels{i}) == 1
                    ignore(i) = 1;
                end
            end
        end
        
        ignore = logical(ignore);
        
        values(:,ignore) = [];
        ch_labels(ignore) = [];
        
        % Confirm they line up with elec labels
        if isequal(ch_labels,pt(whichPt).new_elecs.names) == 0
            error('what\n');
        end
        
        % optional example plot
        if 0
            ex_ch = 20;
            figure
            plot(linspace(curr_window(1),curr_window(2),size(values,1)),...
                values(:,ex_ch))
            title(sprintf('%s',ch_labels{ex_ch}))
            pause
            close(gcf)
        end
        old_values = values;
        
        %% Perform pre-processing
        
        % Do common average reference
        avg_ref = nanmean(values,2);
        values = values - avg_ref;
        
        % Various filters
        for i = 1:size(values,2)
            if do_notch
                values(:,i) = filtfilt(notch_filter,values(:,i));
            end
            if do_hp
                values(:,i) = filtfilt(hp_filter,values(:,i));
            end
            if do_lp
                values(:,i) = filtfilt(lp_filter,values(:,i));
            end
        end
        
        % Pre-whiten
        if do_pre_whiten == 1
            values = pre_whiten(values);
        end
        
        % optional example to plot
        if 0
            ex_ch = 20;
            figure
            plot(linspace(curr_window(1),curr_window(2),size(old_values,1)),...
                old_values(:,ex_ch))
            hold on
            plot(linspace(curr_window(1),curr_window(2),size(values,1)),...
                values(:,ex_ch))
            legend({'Original','Processed'})
            title(sprintf('%s',ch_labels{ex_ch}))
            pause
            close(gcf)
        end
        
        %% Calculate adjacency matrices

        % Get adjacency matrices for all frequencies
        adj = calc_fun_adj(values,fs,freq_bands);
        
        % optional plot
        if 0
            figure
            imagesc(adj(1).adj)
            pause
            close(gcf)
        end
        
        % Flatten the adjacency matrices for storage
        for f = 1:length(adj)
            curr_adj = adj(f).adj;
            curr_adj = flatten_or_expand_adj(curr_adj);
            
            % Fill up the matrix
            out.adj(f).which_adj(shift_idx).data(:,s) = curr_adj;
        end
            
    end
    
    % Save the structure
    save(out_file,'out');
    
end
    
    
end

