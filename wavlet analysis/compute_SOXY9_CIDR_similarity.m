%% SOX9 C-IDR Mass Spectrometry Wavelet Analysis
% This script analyzes similarity between the SOX9 C-terminal IDR aromatic barcode
% and mass-spec–identified IDR barcodes via continuous wavelet transform (CWT) and
% 2D cross-correlation, focusing on magnitude only.
%
% Workflow:
% 1. Load binary-encoded IDR barcodes from an Excel file.
% 2. Compute CWT for each target barcode and select first 30 scales.
% 3. Compute CWT for the original SOX9 C-IDR barcode (positions 244–313).
% 4. Perform cross-correlation (magnitude only) between original SOX9 C-IDR wavelet
%    and each target wavelet code.
% 5. Export results (magnitude) to an Excel sheet.

warning('off');

%% Load and parse barcode Excel
excelFile = 'significant_IDRs_aromatic_barcodes.xlsx' ;
dataTable = readtable(excelFile, 'Sheet', 'Aromatic');
[targetBarcodes, targetIDs, targetNames, targetDiffs, targetPvals, targetFCs] = binary_codes_combination(dataTable);

%% Compute wavelet codes for each target
analyzedWaveletList = {};
filteredIDs       = {};
filteredNames     = {};
filteredDiffs     = {};
filteredPvals     = {};
filteredFCs       = {};
idx = 0;
for i = 1:numel(targetBarcodes)
    code = targetBarcodes{i};
    padLen = 313 - 244;  % length of region 244–313
    padded = [zeros(1,padLen), code, zeros(1,padLen)];
    wt = cwt(padded);
    if size(wt,1) >= 30
        idx = idx + 1;
        analyzedWaveletList{idx} = wt(1:30, :);
        filteredIDs{idx}   = targetIDs{i};
        filteredNames{idx} = targetNames{i};
        filteredDiffs{idx} = targetDiffs{i};
        filteredPvals{idx} = targetPvals(i);
        filteredFCs{idx}   = targetFCs(i);
    end
end

%% Compute original SOX9 C-IDR wavelet code
SOX9_CIDR_barcode = [0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	1	0	0	1	0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	0	0	0	0	0] ;
origWT = cwt(SOX9_CIDR_barcode);
origWavelet = origWT(1:30, 244:313);

%% Cross-correlation (magnitude only)
total_crosscorr_coeff = zeros(idx,1);
for k = 1:idx
    [~,~,mag] = cross_corr2_comp(origWavelet, analyzedWaveletList{k});
    totalMag(k) = mag;
end

disp('Analysis complete');

%% Write results to Excel
outputFile = 'outputFileName.xlsx';
resultTable = table(filteredNames.', filteredIDs.', filteredDiffs.', filteredPvals.', filteredFCs.', ...
                    totalMag.', ...
                    'VariableNames', {'Gene','UniprotID','Difference','P_value','Log2FC','CrossCorr_Coeff'});
writetable(resultTable, outputFile, 'Sheet', 'SOX9_CIDR');

%% Utility Functions
function [codes, ids, names, diffs, pvals, fcs] = binary_codes_combination(tbl)
    % Extract binary barcodes and metadata from table
    rawCodes = tbl{:,5};
    n = numel(rawCodes);
    codes = cell(1,n); ids = cell(1,n); names = cell(1,n);
    diffs = cell(1,n); pvals = zeros(1,n); fcs = zeros(1,n);
    for j = 1:n
        binStr = rawCodes{j};
        codes{j} = double(binStr) - 48;  % ASCII to 0/1
        names{j} = string(tbl.Name(j));
        ids{j}   = string(tbl.UniprotID(j));
        diffs{j} = string(tbl.Difference(j));
        pvals(j) = tbl{j,7};
        fcs(j)   = tbl{j,8};
    end
end

% cross_corr2_comp: compute maximum similarity between two wavelet matrices
% Inputs:
%   wt1: reference wavelet (e.g., SOX9 C-IDR)
%   wt2: target wavelet code
% Outputs:
%   horiz: horizontal slice of the correlation map (central row)
%   maxIdx: index of maximum magnitude
%   max_crosscorr_coeff: value of maximum magnitude
function [horiz, maxIdx, maxMag] = cross_corr2_comp(wt1, wt2)
    ccMap = xcorr2(wt1, wt2);
    row = size(wt2,1);
    horiz = ccMap(row,:);
    mags = abs(horiz);
    [maxMag, maxIdx] = max(mags);
end
