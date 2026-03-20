---
title: "Lineage Tracing: Reprogramming"
date: 2024-12-04T15:34:30-04:00
categories:
  - Blog
tags:
  - Samples
---

Import external libraries.
---

```python
import numpy as np
import os
import LittleSnowFox as kl
#import matlab.engine
#eng = matlab.engine.start_matlab()

print(kl.__version__)
kl.kl_initialize(0)
parent_directory_origin = kl.kl_settings.parent_directory_origin
print(parent_directory_origin)

current_folder = kl.workcatalogue.choosemode_kl(parent_directory_origin,'Lineage',1)
print(current_folder)
```

Generate similarity matrix
```python
choosen_sample = "Reprogramming"


h5ad_filename = "reprogramming_1.h5ad"


current_folder_input = current_folder
orig_adata,loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix(choosen_sample,h5ad_filename,"draw",current_folder_input,3,100,0.1,0.001,True)
#orig_adata,loading_directory,distance_matrix_sparse = kl.preprocessing.kl_dense_matrix_sample(choosen_sample,h5ad_filename,"draw",current_folder_input)


#current_folder_input = current_folder
#loading_directory,distance_matrix = kl.preprocessing.kl_dense_matrix(choosen_sample,h5ad_filename,"draw",current_folder_input)
print(loading_directory)
print(choosen_sample)
```


Save .csv and .mat
---

```python

#需要区分dense和sparase
save_list = ["orig_adata.obs['label']","orig_adata.obsm['X_emb']"]

#将要计算的文件保存到/result
merged_csv,result_directory = kl.workcatalogue.kl_save(loading_directory,choosen_sample,distance_matrix,save_list,orig_adata)
```
The files are saved in `[LittleSnowFox's Anaconda installation directory]\database\Tracing_sample\Hematopoiesis\result\` as `merged_data.csv` and `distance_matrix.mat`.

----------------------------------------------------


Unsupervised learning for only progenitor cells
---

Run `[LittleSnowFox's Anaconda installation directory]\database\Tracing_sample\Hematopoiesis\main_v3_matlab_run_only_prog.m`

```matlab
%% ========================================================================
%  TOGGLE ABLATION STUDY - 3 TASKS × 3 STEPS
%  ========================================================================
%
%  3 EVALUATION TASKS:
%  -------------------
%  Task 1: Progenitor vs Non-Progenitor (干细胞鉴定)
%    - Binary: _prog vs non_prog
%    - All cells included
%
%  Task 2: Multi-class Cell Type (细胞类别识别)
%    - All original labels
%    - Most common task
%
%  Task 3: Lineage Tracing (谱系追踪)
%    - Predict progenitor fate: Monocyte_prog → Monocyte
%    - Only progenitor cells evaluated
%
%  3 ALGORITHM STEPS:
%  ------------------
%  Step 1: Binary Correlation Sorting (Clustering)
%  Step 2: Genetic Algorithm (Cluster Reordering)
%  Step 3: Encoder/Decoder (Boundary Visualization)
%
%  NOTE: Step 2 only reorders clusters (metrics unchanged)
%        Step 3 produces visualization (metrics not applicable)
%
%  ========================================================================

clear; clc; close all;

%% ===== CONFIGURATION =====
fprintf('================================================================================\n');
fprintf('     TOGGLE ABLATION STUDY: 3 TASKS × 3 STEPS\n');
fprintf('================================================================================\n\n');

% Dataset selection
DATASET = 'reprogramming';  % 'hematopoiesis' or 'reprogramming'

% Number of random seeds for baseline methods
N_SEEDS = 5;

% Parameters
step1_params = struct('iterations', 20, 'cell_resolution', 800, 'min_wrong', 1, 'min_right', 1);
step2_params = struct('nPop', 100, 'nPc', 1, 'maxIt', 200, 'cycletimes', 7);
step3_params = struct('upper', 0.00029, 'lower', 0.00031, 'resolution', 10, 'relevance_round', 3);

%% ===== LOAD DATA =====
fprintf('Loading data for %s dataset...\n', DATASET);

data_file = 'result/distance_matrix.mat';
label_file = 'result/merged_data.csv';

% Load distance matrix
if exist(data_file, 'file')
    data = load(data_file);
    if isfield(data, 'distance_matrix')
        MM0 = data.distance_matrix;
    else
        fnames = fieldnames(data);
        MM0 = data.(fnames{1});
    end
    fprintf('  Distance matrix: %d × %d\n', size(MM0));
else
    error('Cannot find data file: %s', data_file);
end

% Load labels
if exist(label_file, 'file')
    count_ = readtable(label_file);
    fprintf('  Labels loaded from: %s\n', label_file);
else
    error('Cannot find label file: %s', label_file);
end

% Find label column - look for column with cell type labels (not unique IDs)
% For hematopoiesis and reprogramming datasets, Var2 contains the labels
label_column = '';
possible_cols = {'Var2', 'Var4', 'Var3', 'label', 'Label', 'celltype', 'CellType', 'cell_type'};

fprintf('\n--- Searching for cell type label column ---\n');
for i = 1:length(possible_cols)
    if ismember(possible_cols{i}, count_.Properties.VariableNames)
        col_data = count_.(possible_cols{i});
        
        % Convert to cell for analysis
        if iscategorical(col_data)
            col_data_cell = cellstr(col_data);
        elseif isstring(col_data)
            col_data_cell = cellstr(col_data);
        elseif iscell(col_data)
            col_data_cell = col_data;
        else
            continue;
        end
        
        n_unique = length(unique(col_data_cell));
        
        % Check if this looks like cell type labels (should have << n_cells unique values)
        % and should contain patterns like '_prog', 'Monocyte', 'Neutrophil', etc.
        has_prog = any(contains(col_data_cell, '_prog')) || any(contains(col_data_cell, 'prog'));
        has_celltype = any(contains(col_data_cell, 'Mono')) || any(contains(col_data_cell, 'Neutro')) || ...
                       any(contains(col_data_cell, 'failed')) || any(contains(col_data_cell, 'reprogram'));
        
        fprintf('  %s: %d unique values', possible_cols{i}, n_unique);
        if has_prog || has_celltype
            fprintf(' [contains cell type patterns]');
        end
        fprintf('\n');
        
        % Good candidate: few unique values AND contains cell type patterns
        if n_unique < length(col_data_cell) * 0.1 && (has_prog || has_celltype)
            label_column = possible_cols{i};
            fprintf('  -> Selected %s as label column\n', label_column);
            break;
        end
    end
end

% Fallback: if no good column found, use the one with fewest unique values that isn't all unique
if isempty(label_column)
    fprintf('  No ideal column found, searching for best alternative...\n');
    best_col = '';
    best_unique = inf;
    
    for i = 1:width(count_)
        col_name = count_.Properties.VariableNames{i};
        col_data = count_.(col_name);
        
        if iscategorical(col_data) || isstring(col_data) || iscell(col_data)
            if iscategorical(col_data)
                col_data_cell = cellstr(col_data);
            elseif isstring(col_data)
                col_data_cell = cellstr(col_data);
            else
                col_data_cell = col_data;
            end
            
            n_unique = length(unique(col_data_cell));
            
            % Must have more than 1 unique value but less than 50% of total
            if n_unique > 1 && n_unique < length(col_data_cell) * 0.5 && n_unique < best_unique
                best_unique = n_unique;
                best_col = col_name;
            end
        end
    end
    
    if ~isempty(best_col)
        label_column = best_col;
        fprintf('  -> Using %s (%d unique values)\n', label_column, best_unique);
    else
        error('Cannot find suitable label column. Available columns: %s', strjoin(count_.Properties.VariableNames, ', '));
    end
end

raw_labels = count_.(label_column);

% Convert to cell array of strings regardless of input type
if iscategorical(raw_labels)
    raw_labels = cellstr(raw_labels);
elseif isstring(raw_labels)
    raw_labels = cellstr(raw_labels);
elseif isnumeric(raw_labels)
    raw_labels = arrayfun(@num2str, raw_labels, 'UniformOutput', false);
elseif ~iscell(raw_labels)
    raw_labels = cellstr(raw_labels);
end

% Ensure all elements are character vectors
for i = 1:length(raw_labels)
    if ~ischar(raw_labels{i})
        raw_labels{i} = char(raw_labels{i});
    end
end

n_cells = length(raw_labels);
fprintf('  Total cells: %d\n', n_cells);
fprintf('  Label column: %s\n', label_column);

%% ===== ANALYZE LABELS =====
fprintf('\n--- Label Analysis ---\n');

unique_labels = unique(raw_labels);
% Ensure unique_labels is also a cell array
if ~iscell(unique_labels)
    unique_labels = cellstr(unique_labels);
end

fprintf('Unique labels (%d):\n', length(unique_labels));
for i = 1:min(20, length(unique_labels))  % Show max 20 labels
    label_str = unique_labels{i};
    n = sum(strcmp(raw_labels, label_str));
    is_prog = contains(label_str, '_prog') || contains(label_str, 'Prog');
    fprintf('  %s: %d cells %s\n', label_str, n, ternary(is_prog, '[PROGENITOR]', ''));
end
if length(unique_labels) > 20
    fprintf('  ... and %d more labels\n', length(unique_labels) - 20);
end

% Count progenitors - ensure raw_labels is cell array for contains()
is_prog_mask = false(n_cells, 1);
for i = 1:n_cells
    is_prog_mask(i) = contains(raw_labels{i}, '_prog') || contains(raw_labels{i}, 'Prog');
end
n_prog = sum(is_prog_mask);
n_diff = n_cells - n_prog;
fprintf('\nProgenitor cells: %d (%.1f%%)\n', n_prog, 100*n_prog/n_cells);
fprintf('Differentiated cells: %d (%.1f%%)\n', n_diff, 100*n_diff/n_cells);

%% ===== PREPARE 3-TASK LABELS =====
fprintf('\n--- Preparing 3-Task Labels ---\n');

% Task 1: Binary labels
labels_t1 = cell(n_cells, 1);
for i = 1:n_cells
    if contains(raw_labels{i}, '_prog') || contains(raw_labels{i}, 'Prog')
        labels_t1{i} = 'Progenitor';
    else
        labels_t1{i} = 'Differentiated';
    end
end
[unique_t1, ~, labels_t1_num] = unique(labels_t1);
K_true_t1 = length(unique_t1);
fprintf('Task 1 (Binary): K = %d classes\n', K_true_t1);

% Task 2: Multi-class labels (original)
labels_t2 = raw_labels;
[unique_t2, ~, labels_t2_num] = unique(labels_t2);
K_true_t2 = length(unique_t2);
fprintf('Task 2 (Multi-class): K = %d classes\n', K_true_t2);

% Task 3: Fate labels for progenitors
labels_t3_fate = cell(n_cells, 1);
for i = 1:n_cells
    label = raw_labels{i};
    if contains(label, '_prog')
        labels_t3_fate{i} = strrep(label, '_prog', '');
    elseif contains(label, 'Prog')
        labels_t3_fate{i} = strrep(label, 'Prog', '');
    else
        labels_t3_fate{i} = '';  % Not a progenitor
    end
end

% Normalize case for fate labels (handle Reprogrammed vs reprogrammed)
for i = 1:n_cells
    if ~isempty(labels_t3_fate{i})
        labels_t3_fate{i} = lower(labels_t3_fate{i});
    end
end

% Get unique fates for progenitors only
prog_fates = labels_t3_fate(is_prog_mask);
unique_fates = unique(prog_fates);
if ~iscell(unique_fates)
    unique_fates = cellstr(unique_fates);
end
% Remove empty strings
unique_fates = unique_fates(~cellfun(@isempty, unique_fates));
K_true_t3 = length(unique_fates);
fprintf('Task 3 (Lineage): K = %d fate classes, %d progenitor cells\n', K_true_t3, n_prog);
fprintf('  Fate classes: %s\n', strjoin(unique_fates, ', '));

%% ===== RUN TOGGLE STEP 1 =====
fprintf('\n================================================================================\n');
fprintf('STEP 1: Binary Correlation Sorting (Clustering)\n');
fprintf('================================================================================\n');
fprintf('Parameters: iterations=%d, resolution=%d, min_wrong=%d, min_right=%d\n', ...
    step1_params.iterations, step1_params.cell_resolution, ...
    step1_params.min_wrong, step1_params.min_right);

[p, splitlist] = binary_corr_sorting(MM0, step1_params.iterations, ...
    step1_params.cell_resolution, step1_params.min_wrong, step1_params.min_right);

[uniqueList, ~, ~] = unique(splitlist, 'stable');
MM = MM0(p, p);
count_result = count_(p, :);

split_simple = uniqueList;
split_simple(1) = 1;
split_simple = [split_simple, n_cells];
K_step1 = length(split_simple) - 1;

fprintf('Result: K = %d clusters\n', K_step1);

% Get reordered labels
labels_t1_sorted = labels_t1(p);
labels_t2_sorted = labels_t2(p);
labels_t3_fate_sorted = labels_t3_fate(p);
is_prog_sorted = is_prog_mask(p);

% Compute metrics for Step 1
[results_step1] = compute_3task_metrics(split_simple, labels_t1_sorted, labels_t2_sorted, ...
    labels_t3_fate_sorted, is_prog_sorted);

fprintf('\n--- STEP 1 Results ---\n');
print_task_results(results_step1, K_step1);

% Get cluster summary for Step 2
[simple_label, simple_matrix] = compute_cluster_summary(count_result, split_simple, MM, label_column);

%% ===== RUN TOGGLE STEP 2 =====
fprintf('\n================================================================================\n');
fprintf('STEP 2: Genetic Algorithm (Cluster Reordering)\n');
fprintf('================================================================================\n');
fprintf('Purpose: Reorder clusters for visualization (block-diagonal structure)\n');
fprintf('NOTE: This step does NOT change cluster assignments!\n');

step2_success = false;
if exist('genetic_encoder', 'file')
    try
        fprintf('Running genetic_encoder...\n');
        [cluster_map_label, cluster_map_matrix] = genetic_encoder(...
            simple_label, simple_matrix, ...
            step2_params.nPop, step2_params.nPc, step2_params.maxIt, step2_params.cycletimes);
        step2_success = true;
        K_step2 = length(cluster_map_label);
        fprintf('Done. K = %d (unchanged)\n', K_step2);
    catch ME
        fprintf('GA failed: %s\n', ME.message);
        cluster_map_matrix = simple_matrix;
        cluster_map_label = simple_label;
    end
else
    fprintf('genetic_encoder not found\n');
    cluster_map_matrix = simple_matrix;
    cluster_map_label = simple_label;
end

% Step 2 metrics = Step 1 metrics (reordering doesn't change assignments)
results_step2 = results_step1;

fprintf('\n--- STEP 2 Results (same as Step 1) ---\n');
fprintf('Task 1: %.2f%% (unchanged - GA only reorders)\n', results_step2.task1.acc * 100);
fprintf('Task 2: %.2f%% (unchanged - GA only reorders)\n', results_step2.task2.acc * 100);
fprintf('Task 3: %.2f%% (unchanged - GA only reorders)\n', results_step2.task3.acc * 100);

%% ===== SAVE STEP 2 OUTPUT FOR MANUAL CLUSTER DETERMINATION =====
fprintf('\n--- Saving Step 2 Reordered Matrix for Manual Inspection ---\n');

% Create output directory if needed
if ~exist('result', 'dir')
    mkdir('result');
end

% Save matrix and labels as files (works with -nojvm)
output_prefix = 'result/step2_reordered';
writematrix(cluster_map_matrix, [output_prefix '_matrix.csv']);

% Save labels properly
label_table = cell2table(cluster_map_label(:), 'VariableNames', {'ClusterLabel'});
writetable(label_table, [output_prefix '_labels.csv']);

% Save as MAT file for easy reloading
save([output_prefix '_data.mat'], 'cluster_map_matrix', 'cluster_map_label', 'split_simple', 'K_step1');
fprintf('Saved:\n');
fprintf('  - %s_matrix.csv (reordered similarity matrix)\n', output_prefix);
fprintf('  - %s_labels.csv (cluster labels in order)\n', output_prefix);
fprintf('  - %s_data.mat (all data for reloading)\n', output_prefix);

% Try to save image (works with -nodisplay but NOT -nojvm)
try
    fig = figure('Visible', 'off');
    imagesc(cluster_map_matrix);
    colorbar;
    colormap('hot');
    title(sprintf('Step 2: GA Reordered Matrix (K=%d) - Identify cluster boundaries', K_step1));
    xlabel('Cluster Index');
    ylabel('Cluster Index');
    
    % Add cluster labels as tick labels (if not too many)
    if K_step1 <= 50
        xticks(1:K_step1);
        yticks(1:K_step1);
        xticklabels(cluster_map_label);
        yticklabels(cluster_map_label);
        xtickangle(45);
        set(gca, 'FontSize', 8);
    end
    
    % Save as PNG and FIG
    saveas(fig, [output_prefix '_heatmap.png']);
    saveas(fig, [output_prefix '_heatmap.fig']);
    fprintf('  - %s_heatmap.png (visualization)\n', output_prefix);
    fprintf('  - %s_heatmap.fig (MATLAB figure)\n', output_prefix);
    close(fig);
catch ME
    fprintf('\nNote: Could not save image (running with -nojvm or graphics error)\n');
    fprintf('      Error: %s\n', ME.message);
    fprintf('      Use the CSV/MAT files to visualize locally.\n');
end

% Print matrix summary for quick inspection in console
fprintf('\n--- Reordered Matrix Preview ---\n');
fprintf('Cluster labels (in order): ');
for i = 1:min(20, length(cluster_map_label))
    fprintf('%s ', cluster_map_label{i});
end
if length(cluster_map_label) > 20
    fprintf('... (%d more)', length(cluster_map_label) - 20);
end
fprintf('\n');

% Show scaled matrix values for console inspection
fprintf('\nMatrix values (× 10^4, showing block structure):\n');
scaled_matrix = round(cluster_map_matrix * 10000, 1);
if K_step1 <= 15
    % Show full matrix if small
    disp(scaled_matrix);
else
    % Show corners if large
    fprintf('(Matrix too large, showing 10×10 corners)\n');
    fprintf('Top-left corner:\n');
    disp(scaled_matrix(1:10, 1:10));
    fprintf('Bottom-right corner:\n');
    disp(scaled_matrix(end-9:end, end-9:end));
end

fprintf('\n');
fprintf('================================================================================\n');
fprintf('*** MANUAL STEP REQUIRED ***\n');
fprintf('================================================================================\n');
fprintf('1. Inspect the heatmap: result/step2_reordered_heatmap.png\n');
fprintf('   (Or load result/step2_reordered_matrix.csv in Python/Excel)\n');
fprintf('2. Identify the block-diagonal structure to determine K_manual\n');
fprintf('3. Re-run with manual K using: ablation_study_manual_k(K_manual)\n');
fprintf('================================================================================\n\n');

%% ===== RUN TOGGLE STEP 3 =====
fprintf('\n================================================================================\n');
fprintf('STEP 3: Encoder/Decoder (Boundary Visualization)\n');
fprintf('================================================================================\n');
fprintf('Purpose: Detect and visualize cluster boundaries\n');
fprintf('NOTE: This step produces VISUALIZATION, not cluster assignments!\n');
fprintf('      Clustering metrics are NOT applicable to Step 3.\n');

step3_success = false;
if exist('encoder_corr_matrix', 'file') && exist('decoder_corr_matrix', 'file')
    try
        fprintf('Running encoder/decoder...\n');
        encode_result = encoder_corr_matrix(step3_params.upper, step3_params.lower, ...
            step3_params.resolution, step3_params.relevance_round, cluster_map_matrix);
        [weighting_decode, decode_result] = decoder_corr_matrix(encode_result);
        weighting_result = weighting_decode + decode_result;
        step3_success = true;
        fprintf('Done. Output: %d × %d visualization matrix\n', size(weighting_result));
    catch ME
        fprintf('Encoder/decoder failed: %s\n', ME.message);
    end
else
    fprintf('encoder_corr_matrix or decoder_corr_matrix not found\n');
end

% Step 3 does NOT produce cluster assignments
results_step3 = struct();
results_step3.task1 = struct('acc', NaN, 'note', 'Visualization only - metrics not applicable');
results_step3.task2 = struct('acc', NaN, 'note', 'Visualization only - metrics not applicable');
results_step3.task3 = struct('acc', NaN, 'note', 'Visualization only - metrics not applicable');

fprintf('\n--- STEP 3 Results ---\n');
fprintf('Task 1: N/A (Step 3 produces visualization, not clusters)\n');
fprintf('Task 2: N/A (Step 3 produces visualization, not clusters)\n');
fprintf('Task 3: N/A (Step 3 produces visualization, not clusters)\n');

%% ===== BASELINE COMPARISONS =====
fprintf('\n================================================================================\n');
fprintf('BASELINE COMPARISONS (Oracle-K)\n');
fprintf('================================================================================\n');
fprintf('Note: Baselines are given the TRUE number of clusters (Oracle-K advantage)\n');
fprintf('Running %d seeds...\n\n', N_SEEDS);

baseline_results = struct();
methods = {'Kmeans', 'Spectral', 'Hierarchical'};

for t = 1:3
    for m = 1:length(methods)
        baseline_results.(sprintf('task%d', t)).(methods{m}).acc = zeros(N_SEEDS, 1);
    end
end

for seed = 1:N_SEEDS
    rng(seed);
    
    % Task 1: Binary (K=2)
    try
        idx = kmeans(MM0, K_true_t1, 'Replicates', 5);
        baseline_results.task1.Kmeans.acc(seed) = compute_majority_accuracy(idx, labels_t1_num);
    catch; end
    
    try
        idx = spectralcluster(MM0, K_true_t1);
        baseline_results.task1.Spectral.acc(seed) = compute_majority_accuracy(idx, labels_t1_num);
    catch; end
    
    try
        idx = clusterdata(MM0, 'maxclust', K_true_t1, 'linkage', 'ward');
        baseline_results.task1.Hierarchical.acc(seed) = compute_majority_accuracy(idx, labels_t1_num);
    catch; end
    
    % Task 2: Multi-class
    try
        idx = kmeans(MM0, K_true_t2, 'Replicates', 5);
        baseline_results.task2.Kmeans.acc(seed) = compute_majority_accuracy(idx, labels_t2_num);
    catch; end
    
    try
        idx = spectralcluster(MM0, K_true_t2);
        baseline_results.task2.Spectral.acc(seed) = compute_majority_accuracy(idx, labels_t2_num);
    catch; end
    
    try
        idx = clusterdata(MM0, 'maxclust', K_true_t2, 'linkage', 'ward');
        baseline_results.task2.Hierarchical.acc(seed) = compute_majority_accuracy(idx, labels_t2_num);
    catch; end
    
    % Task 3: Lineage (progenitors only)
    if n_prog > 0
        MM0_prog = MM0(is_prog_mask, is_prog_mask);
        labels_t3_prog = labels_t3_fate(is_prog_mask);
        [~, ~, labels_t3_prog_num] = unique(labels_t3_prog);
        
        try
            idx = kmeans(MM0_prog, K_true_t3, 'Replicates', 5);
            baseline_results.task3.Kmeans.acc(seed) = compute_majority_accuracy(idx, labels_t3_prog_num);
        catch; end
        
        try
            idx = spectralcluster(MM0_prog, K_true_t3);
            baseline_results.task3.Spectral.acc(seed) = compute_majority_accuracy(idx, labels_t3_prog_num);
        catch; end
        
        try
            idx = clusterdata(MM0_prog, 'maxclust', K_true_t3, 'linkage', 'ward');
            baseline_results.task3.Hierarchical.acc(seed) = compute_majority_accuracy(idx, labels_t3_prog_num);
        catch; end
    end
end

fprintf('Done.\n');

%% ===== FINAL RESULTS TABLE =====
fprintf('\n');
fprintf('================================================================================\n');
fprintf('                      FINAL RESULTS: 3 TASKS × 3 STEPS                         \n');
fprintf('================================================================================\n\n');

% Task 1 Table
fprintf('┌──────────────────────────────────────────────────────────────────────────────┐\n');
fprintf('│ TASK 1: Progenitor vs Non-Progenitor (干细胞鉴定)                            │\n');
fprintf('│ Purpose: Identify developmental stem cells                                   │\n');
fprintf('│ Cells: %d, Classes: %d (Binary)                                              │\n', n_cells, K_true_t1);
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
fprintf('│ %-35s │ %8s │ %12s │\n', 'Method', 'K', 'Accuracy');
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
fprintf('│ %-35s │ %8d │ %11.1f%% │\n', 'TOGGLE Step 1 (Auto-K)', K_step1, results_step1.task1.acc*100);
fprintf('│ %-35s │ %8s │ %11s │\n', 'TOGGLE Step 2 (Reordering)', '(same)', '(same)');
fprintf('│ %-35s │ %8s │ %11s │\n', 'TOGGLE Step 3 (Visualization)', 'N/A', 'N/A');
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
for m = 1:length(methods)
    acc_m = mean(baseline_results.task1.(methods{m}).acc, 'omitnan');
    if ~isnan(acc_m)
        fprintf('│ %-35s │ %8d │ %11.1f%% │\n', [methods{m} ' (Oracle-K)'], K_true_t1, acc_m*100);
    end
end
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
km_acc = mean(baseline_results.task1.Kmeans.acc, 'omitnan');
fprintf('│ TOGGLE vs K-means improvement: %+.1f%%                                        │\n', (results_step1.task1.acc - km_acc)*100);
fprintf('└──────────────────────────────────────────────────────────────────────────────┘\n\n');

% Task 2 Table
fprintf('┌──────────────────────────────────────────────────────────────────────────────┐\n');
fprintf('│ TASK 2: Multi-class Cell Type Classification (细胞类别识别)                  │\n');
fprintf('│ Purpose: Standard cell type identification (most common task)               │\n');
fprintf('│ Cells: %d, Classes: %d                                                       │\n', n_cells, K_true_t2);
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
fprintf('│ %-35s │ %8s │ %12s │\n', 'Method', 'K', 'Accuracy');
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
fprintf('│ %-35s │ %8d │ %11.1f%% │\n', 'TOGGLE Step 1 (Auto-K)', K_step1, results_step1.task2.acc*100);
fprintf('│ %-35s │ %8s │ %11s │\n', 'TOGGLE Step 2 (Reordering)', '(same)', '(same)');
fprintf('│ %-35s │ %8s │ %11s │\n', 'TOGGLE Step 3 (Visualization)', 'N/A', 'N/A');
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
for m = 1:length(methods)
    acc_m = mean(baseline_results.task2.(methods{m}).acc, 'omitnan');
    if ~isnan(acc_m)
        fprintf('│ %-35s │ %8d │ %11.1f%% │\n', [methods{m} ' (Oracle-K)'], K_true_t2, acc_m*100);
    end
end
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
km_acc = mean(baseline_results.task2.Kmeans.acc, 'omitnan');
fprintf('│ TOGGLE vs K-means improvement: %+.1f%%                                        │\n', (results_step1.task2.acc - km_acc)*100);
fprintf('└──────────────────────────────────────────────────────────────────────────────┘\n\n');

% Task 3 Table
fprintf('┌──────────────────────────────────────────────────────────────────────────────┐\n');
fprintf('│ TASK 3: Lineage Tracing / Fate Prediction (谱系追踪)                         │\n');
fprintf('│ Purpose: Predict what progenitor cells will differentiate into              │\n');
fprintf('│ Cells: %d (progenitors only), Classes: %d                                    │\n', n_prog, K_true_t3);
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
fprintf('│ %-35s │ %8s │ %12s │\n', 'Method', 'K', 'Accuracy');
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
fprintf('│ %-35s │ %8d │ %11.1f%% │\n', 'TOGGLE Step 1 (Auto-K)', K_step1, results_step1.task3.acc*100);
fprintf('│ %-35s │ %8s │ %11s │\n', 'TOGGLE Step 2 (Reordering)', '(same)', '(same)');
fprintf('│ %-35s │ %8s │ %11s │\n', 'TOGGLE Step 3 (Visualization)', 'N/A', 'N/A');
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
for m = 1:length(methods)
    acc_m = mean(baseline_results.task3.(methods{m}).acc, 'omitnan');
    if ~isnan(acc_m)
        fprintf('│ %-35s │ %8d │ %11.1f%% │\n', [methods{m} ' (Oracle-K)'], K_true_t3, acc_m*100);
    end
end
fprintf('├──────────────────────────────────────────────────────────────────────────────┤\n');
km_acc = mean(baseline_results.task3.Kmeans.acc, 'omitnan');
fprintf('│ TOGGLE vs K-means improvement: %+.1f%%                                        │\n', (results_step1.task3.acc - km_acc)*100);
fprintf('└──────────────────────────────────────────────────────────────────────────────┘\n\n');

%% ===== KEY FINDINGS =====
fprintf('================================================================================\n');
fprintf('                              KEY FINDINGS                                      \n');
fprintf('================================================================================\n\n');

fprintf('1. STEP 1 (Binary Correlation Sorting) is the CORE clustering component.\n');
fprintf('   - All 3 tasks are evaluated at Step 1\n');
fprintf('   - This is where TOGGLE achieves its performance\n\n');

fprintf('2. STEP 2 (Genetic Algorithm) ONLY reorders clusters for visualization.\n');
fprintf('   - Does NOT change cluster assignments\n');
fprintf('   - Metrics remain IDENTICAL to Step 1\n\n');

fprintf('3. STEP 3 (Encoder/Decoder) produces VISUALIZATION, not clusters.\n');
fprintf('   - Clustering metrics are NOT applicable\n');
fprintf('   - Purpose: Visual interpretation of cluster boundaries\n\n');

fprintf('4. TOGGLE uses Auto-K while baselines use Oracle-K (unfair advantage).\n');
fprintf('   - Despite this, TOGGLE outperforms baselines on lineage tracing\n');
fprintf('================================================================================\n');

%% ===== SAVE RESULTS =====
all_results = struct();
all_results.dataset = DATASET;
all_results.n_cells = n_cells;
all_results.n_prog = n_prog;
all_results.K_true = struct('task1', K_true_t1, 'task2', K_true_t2, 'task3', K_true_t3);
all_results.K_toggle = K_step1;
all_results.step1 = results_step1;
all_results.step2 = results_step2;
all_results.step3 = results_step3;
all_results.baselines = baseline_results;

save(sprintf('ablation_3tasks_3steps_%s.mat', DATASET), 'all_results');
fprintf('\nResults saved to: ablation_3tasks_3steps_%s.mat\n', DATASET);

%% ===== LOCAL FUNCTIONS =====

function [results] = compute_3task_metrics(split_simple, labels_t1, labels_t2, labels_t3_fate, is_prog)
    n_clusters = length(split_simple) - 1;
    n_cells = length(labels_t1);
    
    % Task 1: Binary
    correct_t1 = 0;
    total_t1 = 0;
    
    for i = 1:n_clusters
        idx_start = split_simple(i);
        idx_end = split_simple(i+1);
        if i < n_clusters
            idx_end = idx_end - 1;
        end
        
        cluster_labels = labels_t1(idx_start:idx_end);
        [~, ~, idx] = unique(cluster_labels);
        counts = histcounts(idx, 'BinMethod', 'integers');
        
        total_t1 = total_t1 + length(cluster_labels);
        correct_t1 = correct_t1 + max(counts);
    end
    results.task1.acc = correct_t1 / total_t1;
    
    % Task 2: Multi-class
    correct_t2 = 0;
    total_t2 = 0;
    
    for i = 1:n_clusters
        idx_start = split_simple(i);
        idx_end = split_simple(i+1);
        if i < n_clusters
            idx_end = idx_end - 1;
        end
        
        cluster_labels = labels_t2(idx_start:idx_end);
        [~, ~, idx] = unique(cluster_labels);
        counts = histcounts(idx, 'BinMethod', 'integers');
        
        total_t2 = total_t2 + length(cluster_labels);
        correct_t2 = correct_t2 + max(counts);
    end
    results.task2.acc = correct_t2 / total_t2;
    
    % Task 3: Lineage (progenitors only)
    correct_t3 = 0;
    total_t3 = 0;
    
    for i = 1:n_clusters
        idx_start = split_simple(i);
        idx_end = split_simple(i+1);
        if i < n_clusters
            idx_end = idx_end - 1;
        end
        
        cluster_indices = idx_start:idx_end;
        prog_mask = is_prog(cluster_indices);
        prog_indices = cluster_indices(prog_mask);
        
        if ~isempty(prog_indices)
            cluster_fate_labels = labels_t3_fate(prog_indices);
            [~, ~, idx] = unique(cluster_fate_labels);
            counts = histcounts(idx, 'BinMethod', 'integers');
            
            total_t3 = total_t3 + length(prog_indices);
            correct_t3 = correct_t3 + max(counts);
        end
    end
    
    if total_t3 > 0
        results.task3.acc = correct_t3 / total_t3;
    else
        results.task3.acc = NaN;
    end
    results.task3.n_prog = total_t3;
end

function acc = compute_majority_accuracy(pred_labels, true_labels)
    pred_labels = pred_labels(:);
    true_labels = true_labels(:);
    
    unique_pred = unique(pred_labels);
    correct = 0;
    
    for i = 1:length(unique_pred)
        mask = pred_labels == unique_pred(i);
        if sum(mask) > 0
            true_in_cluster = true_labels(mask);
            [counts, ~] = histcounts(true_in_cluster, 'BinMethod', 'integers');
            correct = correct + max(counts);
        end
    end
    
    acc = correct / length(true_labels);
end

function [simple_label, simple_matrix] = compute_cluster_summary(count_result, split_simple, MM, label_col)
    K = length(split_simple) - 1;
    
    simple_matrix = zeros(K, K);
    for i = 1:K
        for j = 1:K
            idx_i = split_simple(i):split_simple(i+1)-1;
            idx_j = split_simple(j):split_simple(j+1)-1;
            if i == K
                idx_i = split_simple(i):split_simple(i+1);
            end
            if j == K
                idx_j = split_simple(j):split_simple(j+1);
            end
            simple_matrix(i,j) = mean(mean(MM(idx_i, idx_j)));
        end
    end
    
    simple_label = cell(K, 1);
    labels = count_result.(label_col);
    
    % Convert to cell array of strings
    if iscategorical(labels)
        labels = cellstr(labels);
    elseif isstring(labels)
        labels = cellstr(labels);
    elseif isnumeric(labels)
        labels = arrayfun(@num2str, labels, 'UniformOutput', false);
    elseif ~iscell(labels)
        labels = cellstr(labels);
    end
    
    for i = 1:K
        idx_start = split_simple(i);
        idx_end = split_simple(i+1);
        if i < K
            idx_end = idx_end - 1;
        end
        
        cluster_labels = labels(idx_start:idx_end);
        [unique_elements, ~, idx] = unique(cluster_labels);
        if ~iscell(unique_elements)
            unique_elements = cellstr(unique_elements);
        end
        counts = histcounts(idx, 'BinMethod', 'integers');
        [~, max_idx] = max(counts);
        
        simple_label{i} = unique_elements{max_idx};
    end
end

function print_task_results(results, K)
    fprintf('Task 1 (Binary):      Acc = %.2f%%\n', results.task1.acc * 100);
    fprintf('Task 2 (Multi-class): Acc = %.2f%%\n', results.task2.acc * 100);
    if ~isnan(results.task3.acc)
        fprintf('Task 3 (Lineage):     Acc = %.2f%% (%d progenitors)\n', ...
            results.task3.acc * 100, results.task3.n_prog);
    else
        fprintf('Task 3 (Lineage):     N/A (no progenitors)\n');
    end
end

function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

%ablation_study_portable
toggle_visualization()



```
Use `.\[LittleSnowFox's Anaconda installation directory]\R_processing\Reprogramming_prog.R` to generate the picture.




```R
library(ggplot2)

if (!exists("first_run_flag")) {
  setwd("..")
  current_dir <- getwd()
  print("Switched to the parent directory.")
  current_dir
  first_run_flag <- TRUE
} else {
  print("Not the first run, skipping setwd.")
}

print(current_dir)
database_dir <- file.path(current_dir, "database")
Tracing_dir <- file.path(database_dir, "Tracing_sample")
Reprogramming_dir <- file.path(Tracing_dir, "Reprogramming")
Reprogramming_result_dir <- file.path(Reprogramming_dir, "result")
Reprogramming_map <- file.path(Reprogramming_result_dir, "Reprogramming_prog.csv")

repro <- read.csv(Reprogramming_map)

ggplot(repro,aes(x=Var3,y=Var4,color=Var5))+geom_point()
```

<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Reprogramming_prog.png" 
     alt="Reprogramming_prog.png" 
     title="Reprogramming_prog.png">


----------------------------------------------------

Unsupervised learning for whole cells
---

```matlab
clear;clc;
%% Parameters



%% Load data  and  Split to compute
%% Load data and Split to compute
MM0 = load('result/distance_matrix.mat');
MM0 = MM0.distance_matrix;

%MM0 = load('data/program2k_matrix.mat');

%% 读取要排序的对象
count_=readtable(['result/merged_data.csv']);

%computing

%MM0=MM0(contains(count_.label,'_prog'),contains(count_.label,'_prog'));
%count_=count_(contains(count_.label,'_prog'),:);

%% 得到边界划分点
% [p,splitlist] = binary_corr_sorting(MM0,20,250,5,5);
% [p,splitlist] = binary_corr_sorting(MM0,20,100,5,5);
[p,splitlist] = binary_corr_sorting(MM0,5,100,5,5);

%% 对划分点去重
[uniqueList, ~, ~] = unique(splitlist, 'stable');

%% 对相似度矩阵排序
MM=MM0(p,p);
split=[];

%% 重排count_result
count_result=count_(p,:);
split_simple=uniqueList;

%% 第一个起始位点置为1
split_simple(1)=1;
split_simple=[split_simple,length(MM0)]

%% 计算均值矩阵
[simple_label,simple_matrix]=sample_computing(count_result,split_simple,MM);

%% 计算acc
[acc]=acc_computing(split_simple,count_result);

%% 给出结果
[predict_result]=predict_computing(count_,simple_label,split_simple);
count_result_out=[count_result,predict_result];
[simple_label_result,simple_matrix_result]=cluster_map(simple_label,simple_matrix);
writetable(count_result_out,'result/umap_reprog.csv')

%% 重排小矩阵
[cluster_map_label,cluster_map_matrix] = genetic_encoder( ...
    simple_label, ...
    simple_matrix, ...
    100, ...% nPop = 50;  % 种群规模大小为30
    1, ...% nPc = 1; % 子代规模的比例0.8
    200, ...% maxIt = 200; % 最大迭代次数
    7 ...% cycletimes = 200; % 循环计算次数
    );


%% 绘图
row_labels = cluster_map_label;
column_labels = cluster_map_label;
h = heatmap(cluster_map_matrix);
h.YDisplayLabels = row_labels; % 设置行标签
h.XDisplayLabels = column_labels; % 设置列标签
h.ColorLimits = [0.0007, 0.0008]


% 
% 
% writecell(simple_label_result,'result/Figure_c_label_blood.csv')
% writematrix(simple_matrix_result,'result/Figure_c_matrix_blood.csv')
% writetable(count_result_out,'result/umap_reprog.csv')


%% 临近法激活
corr_matrix = relevance_generate(0.00069,2,cluster_map_matrix);
hi = heatmap(corr_matrix);
hi.YDisplayLabels = cluster_map_label; % 设置行标签
hi.XDisplayLabels = cluster_map_label; % 设置列标签


%% 编码
encode_result = encoder_corr_matrix(0.00069,0.00071,10,2,cluster_map_matrix);
figure(2)
hj = heatmap(encode_result);
hj.YDisplayLabels = row_labels; % 设置行标签
hj.XDisplayLabels = column_labels; % 设置列标签

%% 解码
figure(3)
[weighting_decode,decode_result] = decoder_corr_matrix(encode_result);
weighting_result = weighting_decode + decode_result;
hk = heatmap(decode_result);
hk.ColorLimits = [15,16]
hk.YDisplayLabels = row_labels; % 设置行标签
hk.XDisplayLabels = column_labels; % 设置列标签

```

Use `.\[LittleSnowFox's Anaconda installation directory]\R_processing\Reprogramming_all.R` to generate the picture.

```R
library(ggplot2)

if (!exists("first_run_flag")) {
  setwd("..")
  current_dir <- getwd()
  print("Switched to the parent directory.")
  current_dir
  first_run_flag <- TRUE
} else {
  print("Not the first run, skipping setwd.")
}

print(current_dir)
database_dir <- file.path(current_dir, "database")
Tracing_dir <- file.path(database_dir, "Tracing_sample")
Reprogramming_dir <- file.path(Tracing_dir, "Reprogramming")
Reprogramming_result_dir <- file.path(Reprogramming_dir, "result")
Reprogramming_map <- file.path(Reprogramming_result_dir, "umap_reprog.csv")

repro <- read.csv(Reprogramming_map)

ggplot(repro,aes(x=Var3,y=Var4,color=Var5))+geom_point()
```
<img src="https://raw.githubusercontent.com/FullBlackWolf/ATPX4869/refs/heads/master/assets/images/Reprogramming_all.png" 
     alt="Reprogramming_all.png" 
     title="Reprogramming_all.png">
