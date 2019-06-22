function visualize_fun_adj(whichPt)

%% Parameters
which_freq = 1; 
which_time = 1;

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

%% Load the adj
name = pt(whichPt).name;
pt_folder = [adj_folder,name,'/'];
adj_file = [pt_folder,'adj_',sprintf('%s',name),'.mat'];
out = load(adj_file);
out = out.out;

adj_freq = out.adj(which_freq).which_adj(which_time).data;

%% Get the average
adj_avg = mean(adj_freq,2);
adj_avg = flatten_or_expand_adj(adj_avg);

%% Get the average node strength
ns = sum(adj_avg,1);

%% Get locs
locs = pt(whichPt).new_elecs.locs;

% Plot stuff
figure
set(gcf,'position',[44 396 1397 402])
subplot(1,3,1)
imagesc(adj_freq)

subplot(1,3,2)
imagesc(adj_avg);

subplot(1,3,3)
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),100,ns,'filled')

end