clc 
close all 
clear variables

%%%% Similarity Measure-Based Possibilistic FCM With Label Information %%%% 
%   define the dataset's path in "filepath" variable. If using IBSR dataset
%   just define the folder path (e.g. ".../20Normals_T1_brain/") otherwise
%   the path should be complete (folder + file) 
%   if you are going to segment an image set "is_image" to 1 otherwise set
%   it to 0
%   if you are using IBSR dataset, set "is_IBSR" to 1 otherwise set it to 0
%   if you are using IBSR dataset, set "seg_file" parameter to the
%   path of segmented images folder i.e. ".../20Normals_T1_seg/"
%   define the number of clusters in "cluster_n"
%   set the algorithm parameters in "options" defined whithin the code
%   to check the parameter each option shows refer to sim_pfcm_l.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = './20Normals_T1_brain/';
%filepath = 'ring.txt';
is_image = 1;
is_IBSR = 1;
cluster_n = 4;

if is_image == 1
    if is_IBSR == 1
        n = 17;  % image number from 1 to 20
        seg_file = './20Normals_T1_seg/';
        [seg_raw, img_raw] = IBSRread(seg_file, filepath, n);
        img = img_raw(51:end-50, 41:end-51, 20);    % eliminating image's boder
        seg = seg_raw(51:end-50, 41:end-51, 20);    % eliminating image's boder
        [m, n] = size(img);
        img_vec = reshape(im2double(img), [], 1);
        options = [NaN, NaN, 3, NaN, NaN, NaN, 20, NaN, m, n];  % algorithm parameters
        [U, T, obj_fcn] = sim_pfcm_l(img_vec, cluster_n, options);
        [~, i] = max(U, [], 1);
        figure
        subplot(1, 3, 1)
        imshow(img,[])
        title('Original Image')
        subplot(1, 3, 2)
        imshow(seg,[])
        title('Grand Truth Segmented Image')
        subplot(1, 3, 3)
        imshow(reshape(i, m, n), [])
        title('Proposed method output')
        fprintf('Dice Similarity Coeffcient For Each Cluster\n')
        dice_coeff(i, reshape(seg, [], 1)', cluster_n)
    else
        img = imread(filepath);
        img = imresize(img, 0.5);
        [m, n] = size(img);
        img_vec = reshape(im2double(img), [], 1);
        options = [NaN, NaN, NaN, NaN, NaN, NaN, 5, NaN, m, n]; % algorithm parameters
        [U, T, obj_fcn] = sim_pfcm_l(img_vec, cluster_n, options);
        [~, i] = max(U, [], 1);
        figure
        imshow(reshape(i, m, n), [])
        title('Proposed method output')
    end
else
    mat = txt2mat(filepath);
    options = [NaN, NaN, NaN, NaN, NaN, NaN, 100, NaN, 0, 0];   % algorithm parameters
    [U, T, obj_fcn] = sim_pfcm_l(mat(:, 1:2), cluster_n, options);
    [~, i] = max(U, [], 1);
    figure
    subplot(1, 2, 1)
    gscatter(mat(:,1), mat(:,2), mat(:,3))
    title('Ground Truth Segmented Points')
    subplot(1, 2, 2)
    gscatter(mat(:,1), mat(:,2), i')
    title('Proposed Method Output')
    fprintf('Dice Similarity Coeffcient For Each Cluster\n')
    dice_coeff(i, mat(:, 3)', cluster_n)
end

