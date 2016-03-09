%% script_extract_images
% The source code are from 
% http://ufldl.stanford.edu/wiki/index.php/Using_the_MNIST_Dataset
clc; clear; close all;

Train_images = loadMNISTImages('./../../data/train-images-idx3-ubyte');
Train_labels = loadMNISTLabels('./../../data/train-labels-idx1-ubyte');

Test_images = loadMNISTImages('./../../data/t10k-images-idx3-ubyte');
Test_labels = loadMNISTLabels('./../../data/t10k-labels-idx1-ubyte');

save('./../../data/MNIST.mat');