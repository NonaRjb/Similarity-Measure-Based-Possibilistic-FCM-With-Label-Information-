clc 
close all 
clear variables

filepath = 'Mri2.bmp';
is_image = 1;
is_IBSR = 0;
cluster_n = 5;

if is_image == 1
    if is_IBSR == 1
        n = 10;
        seg_file = './20Normals_T1_seg/';
        [img_raw, seg_raw] = IBSRread(seg_file, filepath, n);
        img = img_raw(:, :, 20);
        seg = seg_raw(:, :, 20); 
        [m, n] = size(img);
        img_vec = reshape(im2double(img), [], 1);
        options = [NaN, NaN, NaN, NaN, NaN, NaN, 5, NaN, m, n];
        [U, T, obj_fcn] = sim_pfcm_l(img_vec, cluster_n, options);
        [~, i] = max(U, [], 1);
        figure
        imshow(reshape(i, m, n), [])
    else
        img = imread(filepath);
        img = imresize(img, 0.5);
        [m, n] = size(img);
        img_vec = reshape(im2double(img), [], 1);
        options = [NaN, NaN, NaN, NaN, NaN, NaN, 5, NaN, m, n];
        [U, T, obj_fcn] = sim_pfcm_l(img_vec, cluster_n, options);
        [~, i] = max(U, [], 1);
        figure
        imshow(reshape(i, m, n), [])
    end
else
    mat = txt2mat(filepath);
    options = [NaN, NaN, NaN, NaN, NaN, NaN, 5, NaN, 0, 0];
    [U, T, obj_fcn] = sim_pfcm_l(mat(:, 1:2), cluster_n, options);
    [~, i] = max(U, [], 1);
    figure
    gscatter(mat(:,1), mat(:,2), i')
end

