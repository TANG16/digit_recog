function [nidx, ndist] = neighborsSpectral (X, Y, varargin)
%
% NEIGHBORSSPECTRAL - Spectral neighbors of HSI pixels
%   
% SYNTAX
%
%   NIDX = NEIGHBORSSPECTRAL( X )
%   NIDX = NEIGHBORSSPECTRAL( X, Y )
%   NIDX = NEIGHBORSSPECTRAL( X, Y, 'k', K, ... )
%   NIDX = NEIGHBORSSPECTRAL( X, Y, 'distance', DIST, ... )
%   NIDX = NEIGHBORSSPECTRAL( X, Y, 'meanshift', MASKMS, ... )
%   NIDX = NEIGHBORSSPECTRAL( X, [], 'knnoffset', W, ... )
%   NIDX = NEIGHBORSSPECTRAL( X, [], 'knnoffset', W, 'dims', DIMX, .... )
%   NIDX = NEIGHBORSSPECTRAL( X, [], 'knnoffset', W, 'knnblock', BLK, ... )
%   NIDX = NEIGHBORSSPECTRAL( X, [], 'knnoffset', W, 'parallel', PARFLAG, ... )
%   [NIDX, NDIST] = NEIGHBORSSPECTRAL( X, ... )
%
% INPUT
%
%   X           Source feature vectors array            [D-by-NX]
%   Y           Target feature vectors array            [D-by-NY]
%               {default: X}
%   K           Number of spectral neighbors            [scalar]
%               {default: 10}
%   DIST        Distance metric                         [string]
%               (see KNNSEARCH)
%               {default: 'euclidean'}
%   MASKMS      Mask that indicates which feature       [D-by-1; logical]
%               dimensions to shift to zero mean
%               {default: true(D,1)}
%   W           Maximum spatial offset between returned [scalar]
%               neighbors (i.e. spectral neighbors are only
%               searched for in a (2W+1)-by-(2W+1) sliding window
%               {default: []}
%   DIMX        Dimensions of underlying image for X    [2-vector]
%               {default: []}
%   BLK         Block size (1D) for processing of       [scalar]
%               spatially bounded all-kNN search (refers to spatial
%               as well as spectral block-wise decomposition)
%               {default: []}
%   PARFLAG     Use the Parallel Computing Toolbox      [boolean]
%               within ALLKNNBOUNDED?
%               {default: []}
%
% OUTPUT
%
%   NIDX        Spectral neighbors' indices             [K-by-NX]
%   NDIST       Spectral neighbors' distances           [K-by-NX]
%
% DESCRIPTION
%
%   NIDX = NEIGHBORSSPECTRAL(X) returns the nearest spectral neighbors
%   of the feature vectors in X.
%
%   NDIX = NEIGHBORSSPECTRAL(X,Y) returns, for each feature vector
%   in X, its nearest neighbors in Y.
%
%   NIDX = NEIGHBORSSPECTRAL(X,Y,'k',K,...) specifies the desired number
%   of spectral neighbors.
%
%   NIDX = NEIGHBORSSPECTRAL(X,Y,'distance',DIST,...) also specifies
%   the distance metric (in a format recognized by KNNSEARCH).
%
%   NIDX = NEIGHBORSSPECTRAL(X,Y,'meanshift',MASKMS,...) also
%   specifies which feature dimensions must be shifted to zero mean
%   prior to distance metric calculations. Both X and Y are shifted
%   by mean(X(MASKMS,:),2).
%
%   NIDX = NEIGHBORSSPECTRAL(X,[],'knnoffset',W,...) only looks for
%   the spectral neighbors of each pixel within a (2W+1)-by-(2W+1)
%   window centered around it. This way, all-kNN searches may be
%   spatially bounded, thus allowing one to deal with larger HSIs. If
%   'knnoffset' is set and W is non-empty, then NEIGHBORSSPECTRAL
%   calls ALLKNNBOUNDED internally to compute the k-nearest spectral
%   neighbors.
%   
%   NIDX = NEIGHBORSSPECTRAL(X,[],'knnoffset',W,'dims',DIMX,...)
%   specifies the dimensions of the underlying HSI. This is necessary
%   in order for ALLKNNBOUNDED to assign spatial locations to the NX
%   input feature vectors. It is assumed that the feature vector
%   indices in X correspond to a column-major ordering of the HSI
%   pixel domain.
%
%   NIDX = NEIGHBORSSPECTRAL(X,[],'knnoffset',W,'knnblock',BLK,...)
%   specifies the blocking parameter for block-wise decomposition of
%   the ALLKNNBOUNDED processing. Blocking is carried out because
%   ALLKNNBOUNDED is (a) memory-intensive, and (b) parallelizable
%   across blocks.
%
%   [NIDX,NDIST] = NEIGHBORSSPECTRAL(X,...) also returns the distance
%   measurements between each point and its returned neighbors.
%
% NOTE
%
%    NEIGHBORSSPECTRAL currently supports shifting to zero mean
%    w.r.t. to the mean of X for both X and Y. This supports the case
%    where Y is either the same as X or a subset (e.g. reference
%    samples) of X. However, it might be more reasonable to shift both
%    with the mean of Y in the case where Y is an a priori available
%    dictionary of spectral features. This feature is not supported
%    yet.
%
%
% See also      knnsearch, hsiSegmentation, allknnBounded
%
    
    
    %% DEFAULTS
    
    % target feature vectors array
    if ~exist( 'Y', 'var' ) || isempty( Y )
        Y = X;
    end
    
    % spatial neighborhood size
    k = 10;
    
    % distance calculation method
    distMethod = 'euclidean';
    
    % mean-shift mask
    maskMS = true( size(X,1), 1 );
    
    % parameters for spatially bounded all-kNN search
    w       = [];
    dimImg  = sqrt( size( X, 2 ) );
    blk     = [];
    parflag = [];
    
    % parameter for kNN search in reference points
    LrefMask = [];

    
    %% NAME-VALUE INPUT PARSING
    
    [k, distMethod, LrefMask, maskMS, w, dimImg, blk, parflag] = ...
        parseOptArgs( k, distMethod, LrefMask, maskMS, w, dimImg, blk, ...
                      parflag, varargin{:} );
    
    if isscalar( maskMS )
        maskMS = repmat( maskMS, size(X,1), 1 );
    end
    
    
    %% INITIALIZATION
      
    % shift the mean of specified feature dimensions to zero
    Xmean = mean( X(maskMS,:), 2 );
    X(maskMS,:) = bsxfun( @minus, X(maskMS,:), Xmean );
    Y(maskMS,:) = bsxfun( @minus, Y(maskMS,:), Xmean );
    
    
    %% SPECTRAL NEIGHBORS
    
    % spatially-bounded all-kNN search
    if ~isempty( w ) && isequal( X, Y )
        
        % allknnBounded expects image-shaped input
        X = reshape( X.', [dimImg, size(X,1)] );
        
        % find k-nearest neighbors within bounded neighborhoods
        [nidx, ndist] = allknnBounded( X, k, w, distMethod, ...
                                       blk, parflag );
                                   
        % reshape allknnBounded output to (MN)-by-k
        nidx  = reshape( nidx , prod(dimImg), k );
        ndist = reshape( ndist, prod(dimImg), k );
        
        % kNN search in reference points
        if ~isempty( LrefMask )
            % find k-nearest neighbors within reference point set
            [nidx_ref, ndist_ref] = refknnSearch( X, LrefMask, k, distMethod);
            
            %permute the output to (MN)-by-k
            nidx_ref  = permute( nidx_ref, [2 1] );
            ndist_ref = permute( ndist_ref, [2 1] );
            
            % merge 2 kNN results
            [nidx_res, ndist_res] = mergekNN( k, nidx, ndist, nidx_ref, ndist_ref);
            
            %result
            nidx = nidx_res;
            ndist = ndist_res;
        end
        
        % permute output to k-by-(MN)
        nidx  = permute( nidx, [2 1] );
        ndist = permute( ndist, [2 1] );
        
        
    % generic kNN search
    else
        
        % if X=Y, nearest neighbor is trivially one's own self
        if isequal( X, Y )
            k     = k + 1;
            rmidx = 1;
        else
            rmidx = [];
        end
        
        % (KNNSEARCH works with N-by-D arrays)
        [nidx, ndist] = knnsearch( Y.', X.', ...
                                   'K', k, 'Distance', distMethod );
        
        % remove trivial neighbor and reshape results to K-by-N
        nidx( :,rmidx) = [];
        ndist(:,rmidx) = [];
        nidx  = nidx.';
        ndist = ndist.';
        
    end  % if (spatially bounded all-kNN vs. generic kNN selection)
    
    
end

%% LOCAL FUNCTION: OPTIONAL ARGUMENT PARSING

function [k, distMethod, LrefMask, maskMS, w, dimImg, blk, parflag] = ...
        parseOptArgs( k, distMethod, LrefMask, maskMS, w, dimImg, blk, ...
                      parflag, varargin )
    
    % inputParser allows for MATLAB-native parsing of name-value inputs
    ip = inputParser;
    
    % do not throw error for unrecognized name-value pairs
    ip.KeepUnmatched = true;
    
    % expand structs as name-value pairs
    ip.StructExpand = true;
    
    % string validation function
    validateString = @(x) assert( ischar(x), 'expected string'  );
    
    % specify recognized parameters and default values
    addParameter( ip, 'k'         , k );
    addParameter( ip, 'distance'  , distMethod, @(x) validateString(x) );
    addParameter( ip, 'LrefMask'  , LrefMask );
    addParameter( ip, 'meanshift' , maskMS );
    addParameter( ip, 'knnoffset' , w );
    addParameter( ip, 'dims'      , dimImg );
    addParameter( ip, 'knnblock'  , blk );
    addParameter( ip, 'parallel'  , parflag );
    
    % parse input arguments
    parse( ip, varargin{:} );
    
    % set output
    k          = ip.Results.k;
    distMethod = ip.Results.distance;
    LrefMask   = ip.Results.LrefMask;
    maskMS     = ip.Results.meanshift;
    w          = ip.Results.knnoffset;
    dimImg     = ip.Results.dims;
    blk        = ip.Results.knnblock;
    parflag    = ip.Results.parallel;
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Tiancheng Liu                       tiancheng.liu@duke.edu
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% VERSION
%
%   0.9 - April 16, 2015
%
% CHANGELOG
%
%   0.9 (Apr 16, 2015) - Alexandros
%       * shifting to zero mean is now done using the mean of X for
%         both X and Y (was with both); this is because in our HSI
%         experiments so far Y is either X or a subset of X
%
%   0.8 (Apr 10, 2015) - Alexandros
%       * added interface and optional arguments for spatially
%         bounded all-kNN search (integration with ALLKNNBOUNDED)
%
%   0.7 (Apr 02, 2015) - Alexandros
%       * replaced name-value arguments of version 0.6 with a mask
%         indicating the feature dimensions that require shifting
%         to zero mean; this allows flexible use of arbitrary features
%
%   0.6 (Mar 22, 2015) - Tiancheng
%       * add another 3 name-value arguments to determine whether we
%         do the zero-mean shifting on original part of the
%         features/derivative part of the features
%
%   0.5 (Feb 24, 2015) - Alexandros
%       * changed neighborhood size and distance method to optional
%         name-value arguments
%
%   0.4 (Feb 03, 2015) - Alexandros
%       * added neighbor distance scores as output array
%
%   0.3 (Feb 02, 2015) - Alexandros
%       * allowed different source and target feature-vector arrays
%       * changed expected input/output dimensions to (D/K)-by-N
%
%   0.2 (Jan 30, 2015) - Alexandros
%       * created NEIGHBORSSPECTRAAL as an independent module
%       * input|output arrays correspond to HSI dimensions
%         (M-by-N-by-(D|K))
%       * ensured zero-mean-shifted feature vectors for cosine
%         calculations 
%
%   0.1 (Jan 20, 2015) - Tiancheng
%       * initial implementation as part of LLE function
%
% ------------------------------------------------------------
