% pca experiment
% clear all
clear all, close all
% load data
load ../../data/MNIST.mat

k = 40;% knn classifier

% read training data and sample
tic
Train = Train_images';
Label = Train_labels';
toc

% read test data
tic
Test = Test_images;
GroundTruth = Test_labels;
toc

%% Feature selection
%% HOG
tic
Train = Hog_digit(Train,28,28);
Test = Hog_digit(Test,28,28);
toc
%% SIFT




%% PCA
fprintf('PCA process');
tic
dim = 40;
[V, ~, ~] = pca(Train', 'Centered', true);
Train_low = compress(Train', V, mean(Train', 1), dim);
toc

tic
Test_low = compress(Test', V, mean(Test', 1), dim);
toc


%% LLE 
lle_k = 7;
lle_dim = 40;

fprintf('LLE process');
tic
[Y, Memb, ~, Sigma, ~] = lle([Train, Test], 'k', lle_k, 'dim', lle_dim, 'svds', true);
toc

Train_low_lle = Y(:, 1:6000);
Test_low_lle = Y(:, 6001:end);
% %% KPCA
% tic
% %simple kernel
% [Train_low_simple, V, ~] = kPCA(Train', k, 'simple', 0);
% Test_low_simple = kPCA_NewData(Test',Train', V,'simple', 0);
% toc
% 
% tic
% %poly kernel
% [Train_low_poly, V, ~] = kPCA(Train', k, 'poly', 3);
% Test_low_poly = kPCA_NewData(Test',Train', V,'poly', 3);
% toc
% 
% tic
% %gaussian kernel
% DIST=distanceMatrix(Train');
% DIST(DIST==0)=inf;
% DIST=min(DIST);
% para=5*mean(DIST);
% [Train_low_gaussian, V, ~] = kPCA(Train', k, 'gaussian',para); 
% Test_low_gaussian = kPCA_NewData(Test',Train', V,'gaussian', para);
% toc

%% ICA
% tic
% [Train_low_ica, A, T, mu] = myICA(Train, k, false);
% toc
% 
% tic
% Test_low_ica = A*T*(Test-repmat(mu,1,size(Test,2)));
% toc

%% kNN - NULL
fprintf('kNN in HOG feature space');
% training
tic
k = 5;
KNNMdl = fitcknn(Train', Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_origin, score_origin, cost] = predict(KNNMdl, Test');
toc


%% kNN - PCA
fprintf('kNN after PCA');
% training
tic
%k = 7;
KNNMdl = fitcknn(Train_low, Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_pca, score, cost] = predict(KNNMdl, Test_low);
toc

%% kNN - LLE
fprintf('kNN after LLE');
% training
tic
%k = 7;
KNNMdl_lle = fitcknn(Train_low_lle', Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_lle, score_lle, cost_lle] = predict(KNNMdl_lle, Test_low_lle');
toc

% %% KNN
% 
% % training
% tic
% k = 5;
% KNNMdl = fitcknn(Train_low_simple, Label', 'Numneighbors', k, 'Standardize', 1);
% toc
% 
% % testing
% tic
% [predictlabel_simple, score, cost] = predict(KNNMdl, Test_low_simple);
% toc
% 
% %% KNN
% 
% % training
% tic
% k = 5;
% KNNMdl = fitcknn(Train_low_poly, Label', 'Numneighbors', k, 'Standardize', 1);
% toc
% 
% % testing
% tic
% [predictlabel_poly, score, cost] = predict(KNNMdl, Test_low_poly);
% toc
% 
% %% KNN
% 
% % training
% tic
% k = 5;
% KNNMdl = fitcknn(Train_low_gaussian, Label', 'Numneighbors', k, 'Standardize', 1);
% toc
% 
% % testing
% tic
% [predictlabel_gaussian, score, cost] = predict(KNNMdl, Test_low_gaussian);
% toc


