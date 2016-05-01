% pca experiment
% clear all
clear all, close all
% load data
load ../data/MNIST.mat


% read training data and sample
tic
Train = Train_images';
Label = Train_labels';
toc

% read test data
tic
Test = Test_images';
GroundTruth = Test_labels';
toc

kernel_type = '+gauss';
kernel_width = 0.5;
PHI = sbl_kernelFunction(Train, Train, kernel_type, kernel_width);

p = length(unique(GroundTruth));
tdata = 
for i = 0:1:p-1
    