%% SOX9-eGFP and Cy5-DNA colocalization analysis
% This script analyzes signal intensities between SOX9-eGFP (green channel)
% and Cy5-DNA (red channel) puncta in confocal images. It performs:
%   1. Frequency-domain filtering (FFT-based) for each channel
%   2. Intensity-based thresholding to extract puncta
%   3. Cross-channel masking to calculate signal colocalization
%   4. Summary statistics (sum and mean intensity) for each channel's colocalized signal
%   5. Visualization and export of results as PNG and Excel files
% Created by Tae Gyun Kim

%% Set path and search for TIFF files
addpath(cd)
fold_path = uigetdir;                   % Select folder with images
cd(fold_path);

f_info = dir('*.tif');                  % Find .tif images
fname_list = {f_info.name};            % Store filenames

%% Initialize result storage
sum_green_list = [];
mean_green_list = [];
sum_red_list = [];
mean_red_list = [];
fname_list_selec = {};

%% Process images
warning('off')

for i = 1:length(fname_list) 
    fname_tmp = fname_list{i};
    tmp = imread(fname_tmp);
    red_ch = tmp(:,:,1);             % Cy5-DNA (red)
    green_ch = tmp(:,:,2);           % SOX9-eGFP (green)

    % FFT filtering for each channel
    filt_green = image_FFT(green_ch);
    filt_red = image_FFT(red_ch);

    % Extract green signal within red puncta, and vice versa
    [green_in_red, green_bg] = select_puncta(filt_red, green_ch, 700);
    [red_in_green, red_bg] = select_puncta(filt_green, red_ch, 780);

    % Visualization
    close all
    figure_generator(fname_tmp, filt_green, filt_red, green_in_red, red_in_green);

    % Summary values
    sum_green_list(end+1) = sum(green_in_red(:));
    mean_green_list(end+1) = mean(nonzeros(green_in_red));
    sum_red_list(end+1) = sum(red_in_green(:));
    mean_red_list(end+1) = mean(nonzeros(red_in_green));
    fname_list_selec{end+1} = fname_tmp;
end

%% Save summary table
save_fname = 'summarized_result.xlsx';
temp_table = table(fname_list_selec.', sum_green_list.', sum_red_list.', mean_green_list.', mean_red_list.', ...
    'VariableNames', {'FileName', 'Sum_Green_in_Red', 'Sum_Red_in_Green', 'Mean_Green', 'Mean_Red'});
writetable(temp_table, save_fname);

%% Functions 
function figure_generator(fname, green_img, red_img, green_in_red, red_in_green)
    % Save filtered channels and colocalized signals
    figure(1); set(gcf, 'Position', [80, 300, 1800, 400]);
    ax1 = subplot(1,3,1); imagesc(red_img); title('Cy5-DNA'); colorbar;
    colormap(ax1, [linspace(0,1,256)', zeros(256,1), zeros(256,1)]);

    ax2 = subplot(1,3,2); imagesc(green_img); title('SOX9-eGFP'); colorbar;
    colormap(ax2, [zeros(256,1), linspace(0,1,256)', zeros(256,1)]);

    ax3 = subplot(1,3,3); imagesc(green_in_red); title('eGFP colocalized with Cy5'); colorbar;
    colormap(ax3, [zeros(256,1), linspace(0,1,256)', zeros(256,1)]);
    saveas(gcf, [fname(1:end-4),'_eGFP.png']);

    figure(2); set(gcf, 'Position', [80, 300, 1800, 400]);
    ax1 = subplot(1,3,1); imagesc(green_img); title('SOX9-eGFP'); colorbar;
    colormap(ax1, [zeros(256,1), linspace(0,1,256)', zeros(256,1)]);

    ax2 = subplot(1,3,2); imagesc(red_img); title('Cy5-DNA'); colorbar;
    colormap(ax2, [linspace(0,1,256)', zeros(256,1), zeros(256,1)]);

    ax3 = subplot(1,3,3); imagesc(red_in_green); title('Cy5 colocalized with eGFP'); colorbar;
    colormap(ax3, [linspace(0,1,256)', zeros(256,1), zeros(256,1)]);
    saveas(gcf, [fname(1:end-4),'_cy5.png']);
end

function [signal_in_mask, signal_outside] = select_puncta(mask_source, signal_channel, threshold_factor)
    % Detect bright puncta using thresholded version of mask_source
    level = graythresh(uint8(mask_source));
    int_threshold = level * threshold_factor;
    binary_mask = imbinarize(mask_source, int_threshold);
    signal_in_mask = double(signal_channel) .* double(binary_mask);
    signal_outside = double(signal_channel) .* double(~binary_mask);
end

function filt_img = image_FFT(Img_input)
    % Apply low-pass filter in Fourier domain using circular mask
    Img_fft = fft2(Img_input);
    Img_fft_shifted = fftshift(Img_fft);
    log_fft = log(1 + Img_fft);

    figure(1); subplot(1,3,2);
    imagesc(abs(log_fft)); title('Log-FFT');

    circle_radius = 350;
    windowSize = 750;
    center_pos = size(Img_input,1)/2;

    circ = drawcircle('Center', [center_pos, center_pos], 'Radius', circle_radius);
    circ_mask = createMask(circ);
    kernel = ones(windowSize) / windowSize^2;
    circ_mask = conv2(single(circ_mask), kernel, 'same');

    subplot(1,3,2); imshow(log_fft .* circ_mask, []);
    filt_fft = Img_fft_shifted .* double(circ_mask);
    filt_img = real(ifft2(ifftshift(filt_fft)));
    close
end
