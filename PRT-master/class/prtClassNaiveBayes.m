classdef prtClassNaiveBayes < prtClass
 %prtClassNaiveBayes Naive Bayes Classifier
 % 
 %    CLASSIFIER = prtClassNaiveBayes returns a Naive Bayes Classifier.
 %
 %    CLASSIFIER = prtClassFld(PROPERTY1, VALUE1, ...) constructs a
 %    prtClassNaiveBayes object CLASSIFIER with properties as specified by
 %    PROPERTY/VALUE pairs.
 %
 %    A prtClassNaiveBayes object inherits all properties from the abstract class
 %    prtClass. In addition is has the following properties:
 %
 %    baseRv             - The base type of random variable to be used in
 %                         training the model; baseRv is of type prtRv.
 %                         By default baseRv is a prtRvMvn.
 %
 %    A naive Bayes classification algorithm learns a distribution for the
 %    data under each hypothesis and assumes independence between the data
 %    features (columns) to simplify inference.  
 %
 %    A prtClassNaiveBayes object inherits the TRAIN, RUN, CROSSVALIDATE and
 %    KFOLDS methods from prtAction. It also inherits the PLOT method from
 %    prtClass.
 %
 %    Example:
 %
 %    TestDataSet = prtDataGenUniModal;       % Create some test and
 %    TrainingDataSet = prtDataGenUniModal;   % training data
 %    classifier = prtClassNaiveBayes;           % Create a classifier
 %    classifier = classifier.train(TrainingDataSet);    % Train
 %    classified = run(classifier, TestDataSet);         % Test
 %    subplot(2,1,1);
 %    classifier.plot;
 %    subplot(2,1,2);
 %    [pf,pd] = prtScoreRoc(classified,TestDataSet);
 %    h = plot(pf,pd,'linewidth',3);
 %    title('ROC'); xlabel('Pf'); ylabel('Pd');
 %  
 %   See also prtClass, prtClassLogisticDiscriminant, prtClassBagging,
 %   prtClassMap, prtClassCap, prtClassBinaryToMaryOneVsAll, prtClassDlrt,
 %   prtClassPlsda, prtClassKnn, prtClassRvm, prtClassGlrt,  prtClassSvm,
 %   prtClassTreeBaggingCap, prtClassKmsd, prtClassKnn







 
    properties
        baseRv = prtRvMvn;   % The base randon variable
    end
    properties (SetAccess=private)
        
        name = 'Naive Bayes Classifier' % Naive Bayes Classifier
        nameAbbreviation = 'NBC'        % NBC
        isNativeMary = true;            % true
        naiveRv = prtRvIndependent;     % The naive random variable
    end
    
    methods
        function Obj = prtClassNaiveBayes(varargin)
            Obj = prtUtilAssignStringValuePairs(Obj,varargin{:});
        end
    end
    methods %set methods
        function Obj = set.baseRv(Obj,val)
            if ~isa(val,'prtRv')
                error('prtClassNaiveBayes:baseRv','baseRv must be a prtRv; provided value was a %s',class(val));
            end
            Obj.baseRv = val;
        end
    end
    methods (Access=protected, Hidden = true)
        
        function Obj = trainAction(Obj,DataSet)
            
            %for each class in DataSet, builde an Independent RV; and store
            %class label
            tempIndepRv = repmat(prtRvIndependent('baseRv',Obj.baseRv),DataSet.nClasses,1);
            for classInd = 1:DataSet.nClasses
                x = DataSet.getObservationsByClassInd(classInd);
                tempIndepRv(classInd) = tempIndepRv(classInd).mle(x);
            end
            Obj.naiveRv = tempIndepRv;
        end
        
        function DataSet = runAction(Obj,DataSet)
            
            yOut = zeros(DataSet.nObservations,length(Obj.naiveRv));
            for i = 1:length(Obj.naiveRv)
                yOut(:,i) = Obj.naiveRv(i).logPdf(DataSet.getObservations);
            end
            %Replace poorly scaled samples
            logSum = prtUtilSumExp(yOut')';
            yOut = exp(bsxfun(@minus,yOut,logSum));
            DataSet = DataSet.setObservations(yOut);
        end
    end
end
