clear;clc;
addpath classifiers/liblinear-2.1/matlab;
addpath classifiers/libsvm-3.21/matlab;
addpath training/;
load ../data/tinyset_results.mat;

global config
config.C           = 10;               % C parameter for SVM  
config.LIBLINEAR_B = 10;               % LIBLINEAR (-B)
config.LIBSVM_d    = 7;                % LIBSVM poly degree
config.LIBSVM_r    = 1;                % LIBSVM r (coefficient)
config.LIBSVM_g    = (1/1)^2;                % LIBSVM gamma for rbf kernel  
config.BLOCKS      = [14 7 4; 14 7 4]; % block sizes for histogramming
config.DO_OVERLAP  = true;             % have overlapping blocks
config.NORI        = 12;               % number of orientations  
config.PATCH_W     = 28;               % patch width (do not change)
config.PATCH_H     = 28;               % patch height (do not change)
config.NORM_TYPE   = 'l2';             % or l1 (total pixels are normalized)
config.GRAD_TYPE   = 2;                % 0 - tap, 1-sobel, 2 - gaussian filters 
config.GRAD_SIGMA  = 2;                % sigma of the gaussian filter

tr_labels = Label';
tr_feats = (double(Train))';
te_labels = Truth';
te_feats = (double(Test))';

% fprintf('\n\n --- LINEAR SVM --- \n');
% config.KERNEL_TYPE = 0;    % linear SVM
% models = train_models(tr_labels,tr_feats);
% [acc_linear, pl_linear] = predict_labels(models,te_labels,te_feats);
% 
% fprintf('\n\n --- POLYNOMIAL KERNEL SVM --- \n');
% config.KERNEL_TYPE = 1;
% models = train_models(tr_labels,tr_feats);
% [acc_poly, pl_poly] = predict_labels(models,te_labels,te_feats);

fprintf('\n\n --- RBF KERNEL SVM --- \n');
config.KERNEL_TYPE = 2;
models = train_models(tr_labels,tr_feats);
[acc_rbf, pl_rbf] = predict_labels(models,te_labels,te_feats);

fprintf('\t-------------------------\n');
fprintf('\t Method\t Acc(%%)\t Err(%%)\n');
fprintf('\t-------------------------\n');
% fprintf('\t LINEAR\t%.2f%%\t %.2f%%\n',acc_linear,100-acc_linear);
% fprintf('\t POLY\t%.2f%%\t %.2f%%\n',acc_poly,100-acc_poly);
fprintf('\t RBF\t%.2f%%\t %.2f%%\n',acc_rbf,100-acc_rbf);

% figure;
% bar([acc_linear; acc_poly; acc_rbf]);
% set(gca,'XTickLabel',{'Linear SVM','POLY SVM','RBF SVM'});
% ylabel('Accuracy(%)'); colormap summer; grid on;
