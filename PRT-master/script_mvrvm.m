% pca experiment
% clear all
clear all, close all
% load data
load ../data/MNIST.mat

sigma = 2;
Test_images_sample = Test_images(:,1:1000);
Train_images_sample = Train_images(:,1:6000);
Test_labels_sample = Test_labels(1:1000);
Train_labels_sample = Train_labels(1:6000);
tic
% Train = prtDataSetClass(Train_images', Train_labels);
% Test = prtDataSetClass(Test_images', Test_labels);
Train = prtDataSetClass(Train_images_sample', Train_labels_sample);
Test = prtDataSetClass(Test_images_sample', Test_labels_sample);
toc

classifier = prtClassBinaryToMaryOneVsAll;
classifier.baseClassifier = prtClassRvm;
classifier.internalDecider = prtDecisionMap;

kernel = prtKernelRbf('sigma', sigma);
classifier.baseClassifier.kernels = kernel;
classifier = classifier.train(Train);
classified = run(classifier, Test);