clc 
close all 
clear variables

img = imread('Mri2.bmp');
[m, n] = size(img);
img_vec = reshape(im2double(img), [], 1);

options = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, m, n]; 
cluster_n = 5;
[U, T, obj_fcn] = sim_pfcm_l(img_vec, cluster_n, options);
