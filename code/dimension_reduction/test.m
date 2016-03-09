% training
tic
k = 5;
KNNMdl = fitcknn(Train', Label', 'Numneighbors', k, 'Standardize', 1);
toc

% testing
tic
[predictlabel, score, cost] = predict(KNNMdl, Test');
toc

