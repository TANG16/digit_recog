function digit_LLE(outpath)
% load data
load ../../data/MNIST.mat

k = 40;% knn classifier

% read training data and sample
tic
Train = Train_images;
Label = Train_labels;
toc

% read test data
tic
Test = Test_images;
GroundTruth = Test_labels;
toc

k = 40;
d = 40;
svds_flag = true;
dist = 'euclidean';
fprintf('LLE Processing... \n')
tic
[Y, Memb, ~, sigmaMemb, ~] = lle([Train, Test], 'k', k, ...
    'distance', dist, 'dim', d, 'svds', svds_flag);
toc

Train_low_lle = Y(:, 1:60000);
Test_low_lle = Y(:, 60001:end);

%% kNN - LLE
fprintf('kNN after LLE... \n');
% training
fprintf('Training... \n');
tic
%k = 7;
KNNMdl_lle = fitcknn(Train_low_lle', Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
fprintf('Testing... \n');
tic
[predictlabel_lle, score_lle, cost_lle] = predict(KNNMdl_lle, Test_low_lle');
toc

fprintf( '\n***** WORKSPACE STORAGE *****\n\n' )

allvariables = who;
outfile = savesafe( outpath, allvariables{:} );

fprintf( '- Stored workspace in ''%s''\n\n', outfile );
end