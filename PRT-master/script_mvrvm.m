function script_mvrvm(outpath)

% load data
prtPath;
load ../data/MNIST.mat

[~, I1] =sort(Train_labels);
[~, I2] =sort(Test_labels);
Train = Train_images(:,I1);
Test = Test_images(:,I2);
Label = Train_labels(I1);
GroundTruth = Test_labels(I2);

for i = 0:9
Traindata(:, (i*60+1):(i*60+60)) = Train(:, (i*6000+1):(i*6000+60));
Testdata(:,(i*10+1):(i*10+10)) = Test(:, (i*1000+1):(i*1000+10));
Labels((i*60+1):(i*60+60)) = Label((i*6000+1):(i*6000+60));
Truth((i*10+1):(i*10+10)) = GroundTruth((i*1000+1):(i*1000+10));
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
