function script_nn(outpath)


%% DIGIT-RECOGNITION : NEURAL NETWORK

filename = './MNIST/t10k-images-idx3-ubyte';


[I,labels,I_test,labels_test] = readMNIST(60000); 

%% Neural Network Structure

%Define the structure according to [2]
%Total number of layers
numLayers = 8; 
%Number of subsampling layers
numSLayers = 3; 
%Number of convolutional layers
numCLayers = 3; 
%Number of fully-connected layers
numFLayers = 2;
%Number of input images (simultaneously processed). Need for future
%releases, now only 1 is possible
numInputs = 1; 
%Image width
InputWidth = 32; 
%Image height
InputHeight = 32;
%Number of outputs
numOutputs = 10; 
%Create an empty convolutional neural network with deined structure
sinet = cnn(numLayers,numFLayers,numInputs,InputWidth,InputHeight,numOutputs);

%Now define the network parameters


%Due to implementation specifics layers are always in pairs. First must be
%subsampling and last (before fulli connected) is convolutional layer.
%That's why here first layer is dummy.
sinet.SLayer{1}.SRate = 1;
sinet.SLayer{1}.WS{1} = ones(size(sinet.SLayer{1}.WS{1}));
sinet.SLayer{1}.BS{1} = zeros(size(sinet.SLayer{1}.BS{1}));
sinet.SLayer{1}.TransfFunc = 'purelin';
%Weights 1
%Biases 1


%Second layer - 6 convolution kernels with 5x5 size 
sinet.CLayer{2}.numKernels = 6;
sinet.CLayer{2}.KernWidth = 5;
sinet.CLayer{2}.KernHeight = 5;
%Weights 150
%Biases 6

%Third layer
%Subsampling rate
sinet.SLayer{3}.SRate = 2;
%Weights 6
%Biases 6

%Forth layer - 16 kernels with 5x5 size 
sinet.CLayer{4}.numKernels = 16;
sinet.CLayer{4}.KernWidth = 5;
sinet.CLayer{4}.KernHeight = 5;
%Weights 150
%Biases 6

%Fifth layer
%Subsampling rate
sinet.SLayer{5}.SRate = 2;
%Weights 6
%Biases 6

%Sixth layer - outputs 120 feature maps 1x1 size
sinet.CLayer{6}.numKernels = 120;
sinet.CLayer{6}.KernWidth = 5;
sinet.CLayer{6}.KernHeight = 5;
%Weights 3000
%�������� 120

%Seventh layer - fully connected, 84 neurons
sinet.FLayer{7}.numNeurons = 84;
%Weights 10080
%Biases 84

%Eight layer - fully connected, 10 output neurons
sinet.FLayer{8}.numNeurons = 10;
%Weights 840
%Biases 10

%Initialize the network
sinet = init(sinet);

%According to [2] the generalisation is better if there's unsymmetry in
%layers connections. Yann LeCun uses this kind of connection map:
sinet.CLayer{4}.ConMap = ...
[1 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1;
 1 1 0 0 0 1 1 1 0 0 1 1 1 1 0 1;
 1 1 1 0 0 0 1 1 1 0 0 1 0 1 1 1;
 0 1 1 1 0 0 1 1 1 1 0 0 1 0 1 1;
 0 0 1 1 1 0 0 1 1 1 1 0 1 1 0 1; 
 0 0 0 1 1 1 0 0 1 1 1 1 0 1 1 1; 
]';
%but some papers proposes to randomly generate the connection map. So you
%can try it:
%sinet.CLayer{6}.ConMap = round(rand(size(sinet.CLayer{6}.ConMap))-0.1);
sinet.SLayer{1}.WS{1} = ones(size(sinet.SLayer{1}.WS{1}));
sinet.SLayer{1}.BS{1} = zeros(size(sinet.SLayer{1}.BS{1}));
%In my impementation output layer is ordinary tansig layer as opposed to
%[1,2], but I plan to implement the radial basis output layer

%sinet.FLayer{8}.TransfFunc = 'radbas';


%%
%Now the final preparations
%Number of epochs
sinet.epochs = 3;
%Mu coefficient for stochastic Levenberg-Markvardt
sinet.mu = 0.001;
%Training coefficient
%sinet.teta =  [50 50 20 20 20 10 10 10 5 5 5 5 1]/100000;
sinet.teta =  0.0005;
%0 - Hessian running estimate is calculated every iteration
%1 - Hessian approximation is recalculated every cnet.Hrecomp iterations
%2 - No Hessian calculations are made. Pure stochastic gradient
sinet.HcalcMode = 0;    
sinet.Hrecalc = 300; %Number of iterations to passs for Hessian recalculation
sinet.HrecalcSamplesNum = 50; %Number of samples for Hessian recalculation

%Teta decrease coefficient
sinet.teta_dec = 0.4;

%Images preprocessing. Resulting images have 0 mean and 1 standard
%deviation. Go inside the preproc_data for details
[Ip, labtrn] = preproc_data(I,60000,labels,0);
[I_testp, labtst] = preproc_data(I_test,10000,labels_test,0);
%Actualy training
sinet = train(sinet,Ip,labtrn,I_testp,labtst);

%% TEST
for i = 1: length(I_testp)
    [out, ~] = sim(sinet, I_testp{i});
    test_label(i) = find(out == max(out), 1) - 1;
end

accuracy = numel(find(test_label == labtst)) / numel(labtst);

fprintf( '\n***** WORKSPACE STORAGE *****\n\n' )

allvariables = who;
outfile = savesafe( outpath, allvariables{:} );

fprintf( '- Stored workspace in ''%s''\n\n', outfile );