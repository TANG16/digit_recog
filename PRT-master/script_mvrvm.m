function script_mvrvm(outpath)

% load data
load ../data/MNIST.mat

sigma = 2;

tic
Train = prtDataSetClass(Train_images', Train_labels);
Test = prtDataSetClass(Test_images', Test_labels);
toc

classifier = prtClassBinaryToMaryOneVsAll;
classifier.baseClassifier = prtClassRvm;
classifier.internalDecider = prtDecisionMap;

kernel = prtKernelRbf('sigma', sigma);
classifier.baseClassifier.kernels = kernel;
classifier = classifier.train(Train);
classified = run(classifier, Test);

fprintf( '\n***** WORKSPACE STORAGE *****\n\n' )

allvariables = who;
outfile = savesafe( outpath, allvariables{:} );

fprintf( '- Stored workspace in ''%s''\n\n', outfile );
