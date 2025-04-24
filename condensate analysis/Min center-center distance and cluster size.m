% Set folder paths
folderA = 'F:\TY_Yoon_cowork\20240806_SHAM\TS_csv\DMM_C9_75_month\3779\sox9'; % Folder path for the first dataset
folderB = 'F:\TY_Yoon_cowork\20240806_SHAM\TS_csv\DMM_C9_75_month\3779\polII'; % Folder path for the second dataset

% Set input CSV file names
fileA = 'Cell_2744_1_analysis_cluster_stats.csv'; % First data file
fileB = 'Cell_2745_1_analysis_cluster_stats.csv'; % Second data file

% Combine file paths
pathA = fullfile(folderA, fileA);
pathB = fullfile(folderB, fileB);

output_filename = 'nearest_points_with_size.csv'; % Output file name (unused here)

% Load CSV files
dataA = readmatrix(pathA);
dataB = readmatrix(pathB);

% Extract data from the first dataset (A)
A_x = dataA(:, 2); % x coordinates
A_y = dataA(:, 3); % y coordinates
A_size = dataA(:, 5); % size values

% Extract data from the second dataset (B)
B_x = dataB(:, 2); % x coordinates
B_y = dataB(:, 3); % y coordinates
B_size = dataB(:, 5); % size values

% Create coordinate matrices
A = [A_x, A_y]; % (Nx2)
B = [B_x, B_y]; % (Mx2)

%% Find nearest B cluster for each point in A
distances_A_to_B = pdist2(A, B);
[minDistances_A_to_B, minIdx_A_to_B] = min(distances_A_to_B, [], 2);
closestB_x = B_x(minIdx_A_to_B);
closestB_y = B_y(minIdx_A_to_B);
closestB_size = B_size(minIdx_A_to_B);

% Create result table (A -> B)
result_A_to_B = table(A_x, A_y, A_size, closestB_x, closestB_y, closestB_size, minDistances_A_to_B, ...
    'VariableNames', {'A_x', 'A_y', 'A_size', 'ClosestB_x', 'ClosestB_y', 'ClosestB_size', 'Min_Distance'});

% Save result
output_filename_A_to_B = 'nearest_points_A_to_B.csv';
writetable(result_A_to_B, output_filename_A_to_B);

fprintf('Saved nearest B for each A: %s\n', output_filename_A_to_B);

%% Find nearest A cluster for each point in B
distances_B_to_A = pdist2(B, A);
[minDistances_B_to_A, minIdx_B_to_A] = min(distances_B_to_A, [], 2);
closestA_x = A_x(minIdx_B_to_A);
closestA_y = A_y(minIdx_B_to_A);
closestA_size = A_size(minIdx_B_to_A);

% Create result table (B -> A)
result_B_to_A = table(B_x, B_y, B_size, closestA_x, closestA_y, closestA_size, minDistances_B_to_A, ...
    'VariableNames', {'B_x', 'B_y', 'B_size', 'ClosestA_x', 'ClosestA_y', 'ClosestA_size', 'Min_Distance'});

% Save result
output_filename_B_to_A = 'nearest_points_B_to_A.csv';
writetable(result_B_to_A, output_filename_B_to_A);

fprintf('Saved nearest A for each B: %s\n', output_filename_B_to_A);
