function script_mvrvm(outpath)

% load data
load ../data/MNIST.mat

[~, I1] =sort(Train_labels);
[~, I2] =sort(Test_labels);
Train = Train_images(:,I1);
Test = Test_images(:,I2);
Label = Train_labels(I1);
GroundTruth = Test_labels(I2);

for i = 0:9
Traindata(:, (i*600+1):(i*600+600)) = Train(:, (i*6000+1):(i*6000+600));
Testdata(:,(i*100+1):(i*100+100)) = Test(:, (i*1000+1):(i*1000+100));
Labels((i*600+1):(i*600+600)) = Label((i*6000+1):(i*6000+600));
Truth((i*100+1):(i*100+100)) = GroundTruth((i*1000+1):(i*1000+100));
end

sigma = 3;

tic
Train = prtDataSetClass(Traindata', Labels');
Test = prtDataSetClass(Testdata', Truth');
toc

classifier = prtClassBinaryToMaryOneVsAll;
classifier.baseClassifier = prtClassRvm;
classifier.internalDecider = prtDecisionMap;

kernel = prtKernelRbf('sigma', sigma);
classifier.baseClassifier.kernels = kernel;

fprintf('***BEGIN TRAINING***\n\n');
classifier = classifier.train(Train);
fprintf('***BEGIN TESTING***\n\n');
classified = run(classifier, Test);

fprintf( '\n***** WORKSPACE STORAGE *****\n\n' )

allvariables = who;
outfile = savesafe( outpath, allvariables{:} );

fprintf( '- Stored workspace in ''%s''\n\n', outfile );
