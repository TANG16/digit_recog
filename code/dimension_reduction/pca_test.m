% pca experiment
% clear all
clear all, close all
% load data
load image10k.mat
load labels10k.mat
load trainimage.mat
load trainlabel.mat

k = 40;

% read training data and sample
tic
I1 = [];
I2 = [];
for i = 0:9
I1  = [I1; find(trainlabel == i, 600)];
I2  = [I2; find(labels10k == i, 100)];
end
Train = trainimage(:,I1);
Label = trainlabel(I1);
toc

% read test data
tic
Test = image10k(:,I2);
GroundTruth = labels10k(I2);
toc



%% PCA
tic
k = 40;
[V, ~, ~] = pca(Train', 'Centered', true);
Train_low = compress(Train', V, mean(Train', 1), k);
toc

tic
Test_low = compress(Test', V, mean(Test', 1), k);
toc


%% LLE 
tic
[Y, Memb, ~, Sigma, ~] = lle([Train, Test], 'k', 7, 'dim', k);
toc

Train_low_lle = Y(:, 1:6000);
Test_low_lle = Y(:, 6001:end);

%% KPCA
tic
%simple kernel
[Train_low_simple, V, ~] = kPCA(Train', k, 'simple', 0);
Test_low_simple = kPCA_NewData(Test',Train', V,'simple', 0);
toc

tic
%poly kernel
[Train_low_poly, V, ~] = kPCA(Train', k, 'poly', 3);
Test_low_poly = kPCA_NewData(Test',Train', V,'poly', 3);
toc

tic
%gaussian kernel
DIST=distanceMatrix(Train');
DIST(DIST==0)=inf;
DIST=min(DIST);
para=5*mean(DIST);
[Train_low_gaussian, V, ~] = kPCA(Train', k, 'gaussian',para); 
Test_low_gaussian = kPCA_NewData(Test',Train', V,'gaussian', para);
toc
%% ICA
% tic
% [Train_low_ica, A, T, mu] = myICA(Train, k, false);
% toc
% 
% tic
% Test_low_ica = A*T*(Test-repmat(mu,1,size(Test,2)));
% toc

%% kNN
% training
tic
k = 5;
KNNMdl = fitcknn(Train', Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_origin, score_origin, cost] = predict(KNNMdl, Test');
toc


%% kNN
% training
tic
%k = 7;
KNNMdl = fitcknn(Train_low, Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_pca, score, cost] = predict(KNNMdl, Test_low);
toc

%% kNN
% training
tic
%k = 7;
KNNMdl_lle = fitcknn(Train_low_lle', Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_lle, score_lle, cost_lle] = predict(KNNMdl_lle, Test_low_lle');
toc
%% KNN

% training
tic
k = 5;
KNNMdl = fitcknn(Train_low_simple, Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_simple, score, cost] = predict(KNNMdl, Test_low_simple);
toc

%% KNN

% training
tic
k = 5;
KNNMdl = fitcknn(Train_low_poly, Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_poly, score, cost] = predict(KNNMdl, Test_low_poly);
toc

%% KNN

% training
tic
k = 5;
KNNMdl = fitcknn(Train_low_gaussian, Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel_gaussian, score, cost] = predict(KNNMdl, Test_low_gaussian);
toc



% %% kNN
% % training
% tic
% k = 5;
% KNNMdl = fitcknn(Train_low, Label', 'Numneighbors', k, 'Standardize', 1);
% toc
% 
% % testing
% tic
% [predictlabel, score, cost] = predict(KNNMdl, Test_low);
% toc
% 
% %% kNN
% % training
% tic
% k = 5;
% KNNMdl_lle = fitcknn(Train_low_lle', Label', 'Numneighbors', k, 'Standardize', 1);
% toc
% 
% % testing
% tic
% [predictlabel_lle, score_lle, cost_lle] = predict(KNNMdl_lle, Test_low_lle');
% toc