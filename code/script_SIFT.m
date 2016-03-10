% Script_SIFT

%% load data
load('./../data/MNIST.mat');

numTrain = size(Train_images, 2);
numTest = size(Test_images, 2);
sizeImage = [28 28];

for i = 1: numTrain
    [f, d] = siftfeature(Train_images(:, i), sizeImage);
    feature{i} = f;
    descriptor{i} = d;
end


