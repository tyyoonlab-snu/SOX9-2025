function [num_cond, size_cond, peri_cond, circ_cond, avg_int_cond, sum_int_cond] = extract_condensate_features(fname, circle_radius, windowSize, threshold_factor)
% confocal_cond_analyzer_V5 - Quantifies condensates from a confocal image using FFT filtering and segmentation
% 
% Inputs:
%   fname - Filename of the TIFF image to be analyzed
%   circle_radius - Radius of the circular filter in Fourier domain
%   windowSize - Smoothing kernel size for Fourier mask
%   threshold_factor - Multiplier for Otsu thresholding
%
% Outputs:
%   num_cond - Number of condensates detected
%   size_cond - Area of each condensate
%   peri_cond - Perimeter of each condensate
%   circ_cond - Circularity of each condensate
%   avg_int_cond - Mean pixel intensity per condensate
%   sum_int_cond - Total intensity per condensate

warning('off')

%% Load the image and select green fluorescence channel
tmp = imread(fname);
Img_green = tmp(:,:,2); % Assuming GFP or similar

%% FFT transformation and mask application
Img_fft = fft2(Img_green);
Img_fft_shifted = fftshift(Img_fft);
log_fft = log(1 + Img_fft_shifted);

figure(1); subplot(1,3,2);
imagesc(abs(log_fft)); title('Log-FFT Image');

% Create circular mask
center_pos = size(tmp,1)/2;
circ = drawcircle('Center', [center_pos, center_pos], 'Radius', circle_radius);
circ_mask = createMask(circ);

% Apply smoothing kernel to circular mask
kernel = ones(windowSize) / windowSize^2;
circ_mask = conv2(single(circ_mask), kernel, 'same');

subplot(1,3,2);
imshow(log_fft .* circ_mask, []); title('Masked Log-FFT');

% Apply filtered FFT and inverse transform
Img_fft_filtered = Img_fft_shifted .* double(circ_mask);
Img_fft_ishift = ifftshift(Img_fft_filtered);
Img_filtered = real(ifft2(Img_fft_ishift));

subplot(1,3,1);
imshow(Img_green, [30,200]); title('Original');
subplot(1,3,3);
imshow(Img_filtered, [30,200]); title('Filtered');

%% Thresholding for segmentation
level = graythresh(Img_green);
int_threshold = level * threshold_factor;
Img_binary = imbinarize(Img_filtered, int_threshold);
Img_thres = Img_filtered .* double(Img_binary);
Img_bg = Img_filtered .* double(~Img_binary);

figure(2);
subplot(1,2,1);
imshow(Img_thres); title('Thresholded');

%% Active contour segmentation
mask_init = Img_thres > 0;
Img_segmented = activecontour(Img_filtered, mask_init, 15, 'Chan-Vese');

subplot(1,2,2);
imshow(Img_segmented); title('Active Contour');

%% Small object removal
Img_labeled = bwlabel(Img_segmented, 8);
Labels = unique(Img_labeled);

for i = 1:max(Labels)
    [xv, yv] = find(Img_labeled == i);
    selec_label = double(Img_labeled == i);
    props = regionprops(selec_label, 'Area');
    if props.Area <= 5
        for j = 1:length(xv)
            Img_segmented(xv(j), yv(j)) = 0;
        end
    end
end

%% Morphological smoothing
figure(3);
subplot(1,4,1);
imshow(Img_segmented); title('After Small Obj Remove');

se1 = strel('line', 5, 45);
Img_segmented = imdilate(Img_segmented, se1);

subplot(1,4,2);
imshow(Img_segmented); title('Dilated');

se2 = strel('line', 2, 45);
Img_segmented = imerode(Img_segmented, se2);

subplot(1,4,3);
imshow(Img_segmented); title('Eroded');

%% Final labeling and quantitative analysis
Img_labeled = bwlabel(Img_segmented, 8);
Img_masked = double(Img_green) .* Img_segmented;

Labels = unique(Img_labeled);
num_cond = max(Labels);
size_cond = [];
peri_cond = [];
circ_cond = [];
avg_int_cond = [];
sum_int_cond = [];

for i = 1:num_cond
    mask = Img_labeled == i;
    props = regionprops(mask, 'Area', 'Perimeter', 'Circularity');

    size_cond(i) = props.Area;
    peri_cond(i) = props.Perimeter;
    circ_cond(i) = props.Circularity;

    intensities = Img_masked(mask);
    avg_int_cond(i) = mean(intensities);
    sum_int_cond(i) = sum(intensities);
end
end