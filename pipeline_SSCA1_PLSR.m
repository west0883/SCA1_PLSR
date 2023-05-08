% Put all needed paramters in a structure called "parameters", which you
% can then easily feed into your functions. 
% Use correlations, Fisher transformed, mean removed within mice (mean
% removed for at least the cases when you aren't using mice as response
% variables).

clear all; 

% Create the experiment name.
parameters.experiment_name='SCA1 PLSR';

% Output directory name bases
parameters.dir_base='Y:\Sarah\Analysis\Experiments\';
parameters.dir_exper=[parameters.dir_base parameters.experiment_name '\']; 

% Load mice_all, pass into parameters structure
load([parameters.dir_exper '\mice_all.mat']);
parameters.mice_all = mice_all;

% ****Change here if there are specific mice, days, and/or stacks you want to work with**** 
parameters.mice_all = parameters.mice_all;

% Other parameters
parameters.digitNumber = 2;
parameters.yDim = 256;
parameters.xDim = 256;

parameters.conditions = {'withinMouse', 'acrossMice'};
parameters.periods = {'rest', 'walk'}; % for across mice comparisons
parameters.mouse_dates = {parameters.mice_all(:).date}; 
parameters.data_types = {'fluor', 'corrs'};

% Load periods_nametable_PLSR.m, if it exists yet. (Otherwise is created in
% first step).
if isfile([parameters.dir_exper 'PLSR\periods_nametable_forPLSR.mat'])
    load([parameters.dir_exper 'PLSR\periods_nametable_forPLSR.mat']);

end

% Put relevant variables into loop_variables.
parameters.loop_variables.mice_all = parameters.mice_all;
parameters.loop_variables.variable_type = {'response variables', 'correlations'};
parameters.loop_variables.conditions = parameters.conditions; 
parameters.loop_variables.periods = parameters.periods;
parameters.loop_variables.mouse_dates = parameters.mouse_dates;
parameters.loop_variables.data_types = parameters.data_types;

parameters.average_and_std_together = false;

%% WT mouse ICs
load([parameters.dir_exper 'data/mouse_819_CCA_by_behavior.mat'])
figure;
for i = 1:size(IC_files.dom_list, 3)

    subplot(5, 7, i); imagesc(IC_files.dom_list(:, :, i) .* IC_files.atlas_id(i));
    caxis([-30 30]);
end 

sgtitle('WT mouse')

%% SCA1 mouse ICs

load('mouse_007_CCA_by_behavior_07-01-2021.mat');
figure;
for i = 1:size(IC_files.dom_list, 3)

    subplot(5, 7, i); imagesc(IC_files.dom_list(:, :, i) .* IC_files.atlas_id(i));
    caxis([-30 30]);
end 

sgtitle('SCA1 mouse')

%% Overlay WT mouse 
load([parameters.dir_exper 'data/mouse_819_CCA_by_behavior.mat'])
holder = ones(256) .* (size(IC_files.dom_list, 3) + 1) ;
for i = 1:size(IC_files.dom_list, 3)
    indices = find(IC_files.dom_list(:, :, i)~=0);
    holder(indices) = i;
end 

figure; imagesc(holder);
mymap = [ parula(max(holder, [], 'all')); 1 1 1 ];
colormap(mymap);
sgtitle('WT mouse')

%% Overlay SCA1 mouse
load([parameters.dir_exper 'data/mouse_819_CCA_by_behavior.mat'])
holder = ones(256) .* (size(IC_files.dom_list, 3) + 1) ;
for i = 1:size(IC_files.dom_list, 3)
    indices = find(IC_files.dom_list(:, :, i)~=0);
    holder(indices) = i;
end 

figure; imagesc(holder);
mymap = [ parula(max(holder, [], 'all')); 1 1 1 ];
colormap(mymap);
sgtitle('SCA1 mouse')

%% Overlay atlas IDs WT mouse 
load([parameters.dir_exper 'data/mouse_819_CCA_by_behavior.mat'])
holder = ones(256) * (max(IC_files.atlas_id) + 1);
for i = 1:size(IC_files.dom_list, 3)
    indices = find(IC_files.dom_list(:, :, i)~=0);
    holder(indices) = abs(IC_files.atlas_id(i));
end 

figure; imagesc(holder);
mymap = [ parula(max(holder, [], 'all')); 1 1 1 ];
colormap(mymap);
sgtitle('WT mouse')

%% Overlay atlas IDs SCA1 mouse
load([parameters.dir_exper 'data/mouse_007_CCA_by_behavior.mat'])
holder = ones(256) * (max(IC_files.atlas_id) + 1);
for i = 1:size(IC_files.dom_list, 3)
    indices = find(IC_files.dom_list(:, :, i)~=0);
    holder(indices) = abs(IC_files.atlas_id(i));
end 

figure; imagesc(holder);
mymap = [parula(max(holder, [], 'all')); 1 1 1 ];
colormap(mymap);
sgtitle('SCA1 mouse')

%% Concatenate all trials 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                                  'period', {'loop_variables.periods'}, 'period_iterator'};

parameters.concatenate_across_cells = true;
parameters.concatDim = 3;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'data\']};
parameters.loop_list.things_to_load.data.filename= {'mouse_', 'mouse', '_CCA_by_behavior.mat'};
parameters.loop_list.things_to_load.data.variable= {'atlas_trial_', 'period', '_CCA'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
parameters.loop_list.things_to_save.concatenated_data.filename= {'corrs_', 'mouse', '_','period', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'corrs_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'period';

RunAnalysis({@ConcatenateData}, parameters);

%% Within mice, prepare datasets

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};

% Run Fisher transform? (on correlations only)
parameters.data_type = 'correlations';
parameters.run_FisherTransform = true;

% Inputs
% rest
parameters.loop_list.things_to_load.rest_data.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
parameters.loop_list.things_to_load.rest_data.filename= {'corrs_', 'mouse', '_rest.mat'};
parameters.loop_list.things_to_load.rest_data.variable= {'corrs_all'}; 
parameters.loop_list.things_to_load.rest_data.level = 'mouse';
% walk
parameters.loop_list.things_to_load.walk_data.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
parameters.loop_list.things_to_load.walk_data.filename= {'corrs_', 'mouse', '_walk.mat'};
parameters.loop_list.things_to_load.walk_data.variable= {'corrs_all'}; 
parameters.loop_list.things_to_load.walk_data.level = 'mouse';

% Outputs 
parameters.loop_list.things_to_save.dataset.dir = {[parameters.dir_exper 'data\for PLSR\']};
parameters.loop_list.things_to_save.dataset.filename= {'dataset_info_', 'mouse', '.mat'};
parameters.loop_list.things_to_save.dataset.variable= {'dataset'}; 
parameters.loop_list.things_to_save.dataset.level = 'mouse';

RunAnalysis({@DatasetPrep_SCA1Exper}, parameters)

parameters.run_FisherTransform = false;

%% Within mice, run PLSR
% Will look at the outputs from 10 calculated components.

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'}; 

% Parameters for calculating best number of components. If
% "findBestNComponents" = false, just run the ncomponents_max
parameters.findBestNComponents = true;
parameters.ncomponents_max = 10; 
parameters.contiguous_partitions = true; 
parameters.kFolds = 10;
parameters.MonteCarloReps = 10;
parameters.comparison_type = 'categorical';
parameters.stratify = true;

% Do you want permutations?
parameters.permutationGeneration = true;
parameters.n_permutations = 5000;

% Input 
parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'data\for PLSR\']};
parameters.loop_list.things_to_load.dataset.filename= {'dataset_info_', 'mouse', '.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset'}; 
parameters.loop_list.things_to_load.dataset.level = 'mouse';

% Output
% single result
parameters.loop_list.things_to_save.results.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_save.results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_save.results.level = 'mouse';
% random permutations
parameters.loop_list.things_to_save.Covs_randomPermutations.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.Covs_randomPermutations.filename= {'PLSR_randomPermutations.mat'};
parameters.loop_list.things_to_save.Covs_randomPermutations.variable= {'PLSR_randomPermutations'}; 
parameters.loop_list.things_to_save.Covs_randomPermutations.level = 'mouse';

RunAnalysis({@PLSR_forRunAnalysis}, parameters);  

parameters.findBestNComponents = false;
parameters.permutationGeneration = true;

%% Within mice -- check components
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'}; 
               
parameters.max_response_vars = 2;

parameters.plot_weights = true;
parameters.plot_MSEPs = true;
parameters.plot_BICs = true;
parameters.plot_percentVars = true;
parameters.plot_PCTVAR_response = true;

parameters.define_number_of_sources = true;
parameters.isCorrelationMatrix = true;

% Input
parameters.loop_list.things_to_load.results.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_load.results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_load.results.level = 'mouse';

parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'data\for PLSR\']};
parameters.loop_list.things_to_load.dataset.filename= {'dataset_info_', 'mouse', '.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset'}; 
parameters.loop_list.things_to_load.dataset.level = 'mouse';

% Output
parameters.loop_list.things_to_save.fig_weights.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_weights.filename= {'PLSR_weights.fig'};
parameters.loop_list.things_to_save.fig_weights.variable= {'fig_weights'}; 
parameters.loop_list.things_to_save.fig_weights.level = 'mouse';

parameters.loop_list.things_to_save.fig_MSEPs_explanatory.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_MSEPs_explanatory.filename= {'PLSR_MSEPs_explanatory.fig'};
parameters.loop_list.things_to_save.fig_MSEPs_explanatory.variable= {'fig_MSEPs_explanatory'}; 
parameters.loop_list.things_to_save.fig_MSEPs_explanatory.level = 'mouse';

parameters.loop_list.things_to_save.fig_MSEPs_response.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_MSEPs_response.filename= {'PLSR_MSEPs_response.fig'};
parameters.loop_list.things_to_save.fig_MSEPs_response.variable= {'fig_MSEPs_response'}; 
parameters.loop_list.things_to_save.fig_MSEPs_response.level = 'mouse';

parameters.loop_list.things_to_save.fig_BICs_explanatory.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_BICs_explanatory.filename= {'PLSR_BICs_explanatory.fig'};
parameters.loop_list.things_to_save.fig_BICs_explanatory.variable= {'fig_BICs_explanatory'}; 
parameters.loop_list.things_to_save.fig_BICs_explanatory.level = 'mouse';

parameters.loop_list.things_to_save.fig_BICs_response.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_BICs_response.filename= {'PLSR_BICs_response.fig'};
parameters.loop_list.things_to_save.fig_BICs_response.variable= {'fig_BICs_response'}; 
parameters.loop_list.things_to_save.fig_BICs_response.level = 'mouse';

parameters.loop_list.things_to_save.fig_PCTVARs_explanatory.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_PCTVARs_explanatory.filename= {'PLSR_PCTVARs_explanatory.fig'};
parameters.loop_list.things_to_save.fig_PCTVARs_explanatory.variable= {'fig_PCTVARs_explanatory'}; 
parameters.loop_list.things_to_save.fig_PCTVARs_explanatory.level = 'mouse';

parameters.loop_list.things_to_save.fig_PCTVARs_response.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_PCTVARs_response.filename= {'PLSR_PCTVARs_response.fig'};
parameters.loop_list.things_to_save.fig_PCTVARs_response.variable= {'fig_PCTVARs_response'}; 
parameters.loop_list.things_to_save.fig_PCTVARs_response.level = 'mouse';

RunAnalysis({@CheckComponents}, parameters);

close all;

%% Within mice -- plot Covs ("betas")
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'
               };
% Adjust beta values based on zscore sigmas?
parameters.adjust_beta = true;
parameters.caxis = [-1.5 1.5];

% Input 
parameters.loop_list.things_to_load.results.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_load.results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_load.results.level = 'mouse';

% Also load in dataset values for the zscore sigma.
parameters.loop_list.things_to_load.dataset_info.dir = {[parameters.dir_exper 'data\for PLSR\']};
parameters.loop_list.things_to_load.dataset_info.filename= {'dataset_info_', 'mouse', '.mat'};
parameters.loop_list.things_to_load.dataset_info.variable= {'dataset'}; 
parameters.loop_list.things_to_load.dataset_info.level = 'mouse';

% Output
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'results within mice corrs\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig.filename= {'PLSR_Cov.fig'};
parameters.loop_list.things_to_save.fig.variable= {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'mouse';

RunAnalysis({@PlotBetas}, parameters); 

close all;


%% Within mice fluorescence -- Concatenate all trials 

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                                  'period', {'loop_variables.periods'}, 'period_iterator'};

parameters.concatenate_across_cells = true;
parameters.concatDim = 3;

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'data\']};
parameters.loop_list.things_to_load.data.filename= {'mouse_', 'mouse', '_CCA_by_behavior.mat'};
parameters.loop_list.things_to_load.data.variable= {'period', '_DFF'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Outputs
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
parameters.loop_list.things_to_save.concatenated_data.filename= {'fluor_', 'mouse', '_','period', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable= {'fluor_all'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'period';

RunAnalysis({@ConcatenateData}, parameters);

%% Within mice fluorescence -- average by atlas area

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                                  'period', {'loop_variables.periods'}, 'period_iterator'};
parameters.IC_dim = 2; 
% Inputs
% concatenated trials
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
parameters.loop_list.things_to_load.data.filename= {'fluor_', 'mouse', '_','period', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'fluor_all'}; 
parameters.loop_list.things_to_load.data.level = 'period';
% atlas info
parameters.loop_list.things_to_load.mouse_regions_ordered.dir = {[parameters.dir_exper 'data\']};
parameters.loop_list.things_to_load.mouse_regions_ordered.filename= {'mouse_', 'mouse', '_CCA_by_behavior.mat'};
parameters.loop_list.things_to_load.mouse_regions_ordered.variable= {'mouse_regions_ordered'}; 
parameters.loop_list.things_to_load.mouse_regions_ordered.level = 'mouse';
% IC atlas ids
parameters.loop_list.things_to_load.atlas_ids.dir = {[parameters.dir_exper 'data\']};
parameters.loop_list.things_to_load.atlas_ids.filename = {'mouse_', 'mouse', '_CCA_by_behavior.mat'};
parameters.loop_list.things_to_load.atlas_ids.variable = {'IC_files.atlas_id'}; 
parameters.loop_list.things_to_load.atlas_ids.level = 'mouse';

% Outputs
parameters.loop_list.things_to_save.data_averaged.dir = {[parameters.dir_exper 'data\fluoresce atlas averaged\']};
parameters.loop_list.things_to_save.data_averaged.filename= {'fluor_', 'mouse','_', 'period', '.mat'};
parameters.loop_list.things_to_save.data_averaged.variable= {'period', '_DFF'}; 
parameters.loop_list.things_to_save.data_averaged.level = 'period';

RunAnalysis({@FluorescenceByAtlasRegion}, parameters)

%% Within mice fluorescence -- average per window/roll

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                                  'period', {'loop_variables.periods'}, 'period_iterator'};

parameters.averageDim = 1; 

% Inputs
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'data\fluoresce atlas averaged\']};
parameters.loop_list.things_to_load.data.filename= {'fluor_', 'mouse','_', 'period', '.mat'};
parameters.loop_list.things_to_load.data.variable= {'period', '_DFF'}; 
parameters.loop_list.things_to_load.data.level = 'period';

% Outputs
parameters.loop_list.things_to_save.average.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
parameters.loop_list.things_to_save.average.filename= {'fluor_', 'mouse', '_','period', '_average.mat'};
parameters.loop_list.things_to_save.average.variable= {'fluor_all'}; 
parameters.loop_list.things_to_save.average.level = 'period';

RunAnalysis({@AverageData}, parameters);

%% Within mice fluorescence -- prepare datasets

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};

% Run Fisher transform? (on correlations only)
parameters.data_type = 'fluorescence';

% Inputs
% rest
parameters.loop_list.things_to_load.rest_data.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
parameters.loop_list.things_to_load.rest_data.filename= {'fluor_', 'mouse', '_rest_average.mat'};
parameters.loop_list.things_to_load.rest_data.variable= {'fluor_all'}; 
parameters.loop_list.things_to_load.rest_data.level = 'mouse';
% walk
parameters.loop_list.things_to_load.walk_data.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
parameters.loop_list.things_to_load.walk_data.filename= {'fluor_', 'mouse', '_walk_average.mat'};
parameters.loop_list.things_to_load.walk_data.variable= {'fluor_all'}; 
parameters.loop_list.things_to_load.walk_data.level = 'mouse';

% Outputs 
parameters.loop_list.things_to_save.dataset.dir = {[parameters.dir_exper 'data\for PLSR\fluorescence\']};
parameters.loop_list.things_to_save.dataset.filename= {'dataset_info_', 'mouse', '.mat'};
parameters.loop_list.things_to_save.dataset.variable= {'dataset'}; 
parameters.loop_list.things_to_save.dataset.level = 'mouse';

RunAnalysis({@DatasetPrep_SCA1Exper}, parameters)

%% Within mice fluorescence, run PLSR
% Will look at the outputs from 10 calculated components.

% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'}; 

% Parameters for calculating best number of components. If
% "findBestNComponents" = false, just run the ncomponents_max
parameters.findBestNComponents = true;
parameters.ncomponents_max = 10; 
parameters.contiguous_partitions = true; 
parameters.kFolds = 10;
parameters.MonteCarloReps = 10;
parameters.comparison_type = 'categorical';
parameters.stratify = true;

% Do you want permutations?
parameters.permutationGeneration = true;
parameters.n_permutations = 5000;

% Input 
parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'data\for PLSR\fluorescence\']};
parameters.loop_list.things_to_load.dataset.filename= {'dataset_info_', 'mouse', '.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset'}; 
parameters.loop_list.things_to_load.dataset.level = 'mouse';

% Output
% single result
parameters.loop_list.things_to_save.results.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_save.results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_save.results.level = 'mouse';
% random permutations
parameters.loop_list.things_to_save.Covs_randomPermutations.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.Covs_randomPermutations.filename= {'PLSR_randomPermutations.mat'};
parameters.loop_list.things_to_save.Covs_randomPermutations.variable= {'PLSR_randomPermutations'}; 
parameters.loop_list.things_to_save.Covs_randomPermutations.level = 'mouse';

RunAnalysis({@PLSR_forRunAnalysis}, parameters);  

parameters.findBestNComponents = false;
parameters.permutationGeneration = false;

%% Within mice -- check components
% Always clear loop list first. 
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'}; 
               
parameters.max_response_vars = 2;

parameters.plot_weights = true;
parameters.plot_MSEPs = true;
parameters.plot_BICs = true;
parameters.plot_percentVars = true;
parameters.plot_PCTVAR_response = true;

parameters.isCorrelationMatrix = false;
parameters.define_number_of_sources = false;

% Input
parameters.loop_list.things_to_load.results.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_load.results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_load.results.level = 'mouse';

parameters.loop_list.things_to_load.dataset.dir = {[parameters.dir_exper 'data\for PLSR\']};
parameters.loop_list.things_to_load.dataset.filename= {'dataset_info_', 'mouse', '.mat'};
parameters.loop_list.things_to_load.dataset.variable= {'dataset'}; 
parameters.loop_list.things_to_load.dataset.level = 'mouse';

% Output
parameters.loop_list.things_to_save.fig_weights.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_weights.filename= {'PLSR_weights.fig'};
parameters.loop_list.things_to_save.fig_weights.variable= {'fig_weights'}; 
parameters.loop_list.things_to_save.fig_weights.level = 'mouse';

parameters.loop_list.things_to_save.fig_MSEPs_explanatory.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_MSEPs_explanatory.filename= {'PLSR_MSEPs_explanatory.fig'};
parameters.loop_list.things_to_save.fig_MSEPs_explanatory.variable= {'fig_MSEPs_explanatory'}; 
parameters.loop_list.things_to_save.fig_MSEPs_explanatory.level = 'mouse';

parameters.loop_list.things_to_save.fig_MSEPs_response.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_MSEPs_response.filename= {'PLSR_MSEPs_response.fig'};
parameters.loop_list.things_to_save.fig_MSEPs_response.variable= {'fig_MSEPs_response'}; 
parameters.loop_list.things_to_save.fig_MSEPs_response.level = 'mouse';

parameters.loop_list.things_to_save.fig_BICs_explanatory.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_BICs_explanatory.filename= {'PLSR_BICs_explanatory.fig'};
parameters.loop_list.things_to_save.fig_BICs_explanatory.variable= {'fig_BICs_explanatory'}; 
parameters.loop_list.things_to_save.fig_BICs_explanatory.level = 'mouse';

parameters.loop_list.things_to_save.fig_BICs_response.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_BICs_response.filename= {'PLSR_BICs_response.fig'};
parameters.loop_list.things_to_save.fig_BICs_response.variable= {'fig_BICs_response'}; 
parameters.loop_list.things_to_save.fig_BICs_response.level = 'mouse';

parameters.loop_list.things_to_save.fig_PCTVARs_explanatory.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_PCTVARs_explanatory.filename= {'PLSR_PCTVARs_explanatory.fig'};
parameters.loop_list.things_to_save.fig_PCTVARs_explanatory.variable= {'fig_PCTVARs_explanatory'}; 
parameters.loop_list.things_to_save.fig_PCTVARs_explanatory.level = 'mouse';

parameters.loop_list.things_to_save.fig_PCTVARs_response.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig_PCTVARs_response.filename= {'PLSR_PCTVARs_response.fig'};
parameters.loop_list.things_to_save.fig_PCTVARs_response.variable= {'fig_PCTVARs_response'}; 
parameters.loop_list.things_to_save.fig_PCTVARs_response.level = 'mouse';

RunAnalysis({@CheckComponents}, parameters);

close all;

%% Within mice fluorescence -- plot Covs ("betas")
if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {
               'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'
               };
% Adjust beta values based on zscore sigmas?
parameters.adjust_beta = true;
%parameters.caxis = [-10 10];

% Input 
parameters.loop_list.things_to_load.results.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_load.results.filename= {'PLSR_results.mat'};
parameters.loop_list.things_to_load.results.variable= {'PLSR_results'}; 
parameters.loop_list.things_to_load.results.level = 'mouse';

% Also load in dataset values for the zscore sigma.
parameters.loop_list.things_to_load.dataset_info.dir = {[parameters.dir_exper 'data\for PLSR\fluorescence\']};
parameters.loop_list.things_to_load.dataset_info.filename= {'dataset_info_', 'mouse', '.mat'};
parameters.loop_list.things_to_load.dataset_info.variable= {'dataset'}; 
parameters.loop_list.things_to_load.dataset_info.level = 'mouse';

% Output
parameters.loop_list.things_to_save.fig.dir = {[parameters.dir_exper 'results within mice fluorescence\'], 'mouse', '\'};
parameters.loop_list.things_to_save.fig.filename= {'PLSR_Cov.fig'};
parameters.loop_list.things_to_save.fig.variable= {'fig'}; 
parameters.loop_list.things_to_save.fig.level = 'mouse';

RunAnalysis({@PlotBetas}, parameters); 

close all;

%% Across mice -- find the brain regions you can use
load('Y:\Sarah\Analysis\Experiments\SCA1 PLSR\data\mouse_007_CCA_by_behavior.mat', 'mouse_regions_ordered');
m007_regions = mouse_regions_ordered;

load('Y:\Sarah\Analysis\Experiments\SCA1 PLSR\data\mouse_819_CCA_by_behavior.mat', 'mouse_regions_ordered');
m819_regions = mouse_regions_ordered;

[common_regions, m007_indices, m819_indices] = intersect(m007_regions, m819_regions, 'stable');

save([parameters.dir_exper 'common_regions.mat'], 'common_regions');
save([parameters.dir_exper 'm007_region_indices.mat'],'m007_indices'); 
save([parameters.dir_exper 'm819_region_indices.mat'],'m819_indices'); 

clear m007_regions m819_regions common_regions m007_indices m819_indices;

%% Across mice -- organize the brain regions
% fluorescence

for datai = 1:numel(parameters.data_types)
    data_type = parameters.data_types{datai};

    parameters.data_type = data_type;

    if isfield(parameters, 'loop_list')
    parameters = rmfield(parameters,'loop_list');
    end
    
    % Iterators
    parameters.loop_list.iterators = {'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator';
                                      'period', {'loop_variables.periods'}, 'period_iterator'};
    
    % Inputs
    % region indices
    parameters.loop_list.things_to_load.region_indices.dir = {[parameters.dir_exper]};
    parameters.loop_list.things_to_load.region_indices.filename= {'m', 'mouse', '_region_indices.mat'};
    parameters.loop_list.things_to_load.region_indices.variable= {'m', 'mouse', '_indices'}; 
    parameters.loop_list.things_to_load.region_indices.level = 'mouse';
    % data 
    parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'data\concatenated across trials\']};
    if strcmp(data_type, 'fluor')
        parameters.loop_list.things_to_load.data.filename= { [data_type '_'], 'mouse', '_', 'period', '_average.mat'};
    else
        parameters.loop_list.things_to_load.data.filename= { [data_type '_'], 'mouse', '_', 'period', '.mat'};
    end
    parameters.loop_list.things_to_load.data.variable= {[data_type '_all']}; 
    parameters.loop_list.things_to_load.data.level = 'period';

    % Outputs 
    parameters.loop_list.things_to_save.data_new.dir = {[parameters.dir_exper 'data\shared regions\'], 'mouse', '\'};
    parameters.loop_list.things_to_save.data_new.filename= {data_type, '_', 'period', '.mat'};
    parameters.loop_list.things_to_save.data_new.variable = {'data_shared'}; 
    parameters.loop_list.things_to_save.data_new.level = 'period';

    RunAnalysis({@OnlySharedData}, parameters);

end

%% Across mice -- concatenate data within periods, across mice 
% [no, you don't need to do this-- can put it in the dataset prep code]

if isfield(parameters, 'loop_list')
parameters = rmfield(parameters,'loop_list');
end

% Iterators
parameters.loop_list.iterators = {'data_type', {'loop_variables.data_types'}, 'data_type_iterator';
                                  'period', {'loop_variables.periods'}, 'period_iterator';
                                  'mouse', {'loop_variables.mice_all(:).name'}, 'mouse_iterator'};

parameters.data_type_looping = true;

% Inputs 
parameters.loop_list.things_to_load.data.dir = {[parameters.dir_exper 'data\shared regions\'], 'mouse', '\'};
parameters.loop_list.things_to_load.data.filename= {data_type, '_', 'period', '.mat'};
parameters.loop_list.things_to_load.data.variable = {'data_shared'}; 
parameters.loop_list.things_to_load.data.level = 'mouse';

% Outputs 
parameters.loop_list.things_to_save.concatenated_data.dir = {[parameters.dir_exper 'data\shared regions\concatenated'], 'mouse', '\'};
parameters.loop_list.things_to_save.concatenated_data.filename= {data_type, '_', 'period', '.mat'};
parameters.loop_list.things_to_save.concatenated_data.variable = {'data_shared'}; 
parameters.loop_list.things_to_save.concatenated_data.level = 'mouse';

RunAnalysis({@DatasetPrep_SCA1Exper}, parameters)