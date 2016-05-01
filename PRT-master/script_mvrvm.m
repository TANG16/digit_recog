% pca experiment
% clear all
clear all, close all
% load data
load ../data/MNIST.mat

sigma = 3;

tic
Train = prtDataGenUnimodal;
Test = prtDataGenUnimodal;
toc

classifier = prtClassBinaryToMaryOneVsAll;
classifier.baseClassifier = prtClassRvm;
classifier.internalDecider = prtDecisionMap;

kernel = prtKernelRbf('sigma', sigma);
classifier.baseClassifier.kernels = kernel;
classifier = classifier.train(Train);
classified = run(classifier, Test);