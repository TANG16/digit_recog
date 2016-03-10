% script_LLE

load('./../../data/MNIST.mat');

k = 40;
d = 40;
svds_flag = true;
dist = 'euclidean';
[Y, Memb, ~, sigmaMemb, ~] = lle([Train_images, Test_images], 'k', k, ...
    'distance', dist, 'dim', d, 'svds', svds_flag);