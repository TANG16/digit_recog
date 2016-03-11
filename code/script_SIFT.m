% Script_SIFT

%% load data
load('./../data/MNIST.mat');

[~, I1] =sort(Train_labels);
[~, I2] =sort(Test_labels);
Train = Train_images(:,I1);
Test = Test_images(:,I2);
Label = Train_labels(I1);
GroundTruth = Test_labels(I2);

numTrain = size(Train, 2);
numTest = size(Test, 2);
sizeImage = [28 28];

for i = 1: numTrain
    [f, d] = siftfeature(Train(:, i), sizeImage);
    feature{i} = f;
    descriptor{i} = d;
end


