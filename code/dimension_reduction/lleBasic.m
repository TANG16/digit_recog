function [Y, Memb, idxNhood, sigmaMemb] = lleBasic (X, k, dimEnc)
%
% LLEBASIC - Low-dimensional manifold encoding via Locally Linear
%            Embeddings 
% 
% SYNTAX
%
%   [Y, MEMB, IDXNHOOD, SIGMAMEMB] = LLEBASIC( X, K, DIMENC )
%
% INPUT
%
%   X           Feature vectors array in ambient space  [D-by-N]
%   K           Size of local neighborhoods             [scalar]
%   DIMENC      LLE-based encoding dimensionality       [scalar]
%
% OUTPUT
%
%   Y           Low-dimensional encoding of X vectors   [DIMENC-by-N]
%   MEMB        LLE matrix (MEMB = I - W)               [N-by-N; sparse]
%   IDXNHOOD    Ambient-space neighborhood indices      [K-by-N]
%   SIGMAMEMB   Singular values of LLE matrix           [DIMENC-by-1]
%
% DESCRIPTION
%
%   [Y,MEMB,IDXNHOOD] = LLEBASIC(X,K,DENC) returns the LLE matrix and
%   a DENC-dimensional encoding of the feature vectors in X [1]. It
%   also returns the local neighborhood indices that were used for LLE.
%
% REFERENCES
%
%   [1] S. T. Roweis and L. K. Saul, "Nonlinear Dimensionality
%   Reduction by Locally Linear Embedding," Science, vol. 290,
%   no. 5500, pp. 2323-2326, Dec. 2000.
%
    
    
    %% LOCAL NEIGHBORHOODS
    
    idxNhood = localNeighborhoods( X, k );
    
    
    %% LLE MATRIX
    
    Memb = lleMatrixComputation( X, idxNhood );
    
    
    %% LOW-DIMENSIONAL ENCODING
    
    [Y, sigmaMemb] = lleEncoding( Memb, dimEnc );
    
    
end



%% FUNCTION: LOCAL NEIGHBORHOODS

function nidx = localNeighborhoods( X, k )
    
    % center data around origin
    X = bsxfun( @minus, X, mean( X, 2 ) );
    
    % get the k-nearest spectral neighbors for each point (removing each
    % point's own self)
    % [transpose input because KNNSEARCH works with N-by-D arrays]
    distMetric = 'cosine';
    nidx       = knnsearch( X.', X.', 'K', k+1, ...
                      'Distance', distMetric );
    nidx       = nidx.';
    nidx(1,:)  = [];     % remove self-neighbor
    
end



%% FUNCTION: LLE MATRIX COMPUTATION

function M = lleMatrixComputation( X, idxNhood )
    
    % regularization factor
    regfactor = 1e-3;
    
    % problem dimensions
    n = size( X, 2 );           % # of data vectors
    d = size( X, 1 );           % # of dimensions
    k = size( idxNhood, 1 );    % size of local neighborhoods
    
    % pre-allocate space for local solution vectors
    W = zeros( k, n );

    % iterate across all local systems (equivalently, data vectors)
    for i = 1 : n
        
        % create local neighborhood array and center neighbors
        % around reference (i-th) vector
        NXi = bsxfun( @minus, X(:, idxNhood(:,i)), X(:,i) );
        
        % create quadratic form of local neighborhood vectors array
        Ci = NXi' * NXi;
        
        % regularization (eigenspectrum lift)
        Ci = Ci + regfactor * trace(Ci) * eye(k);
        
        % solve local system and scale to unit algebraic sum
        W(:,i) = Ci \ ones(k,1);
        W(:,i) = W(:,i) ./ sum( W(:,i) );
        
    end  % for (i)
    
    % form global weights matrix in sparse format
    W = sparse( idxNhood(:), ...
                kron( (1 : n).', ones( k, 1 ) ), ...
                W(:), ...
                n, n );
    
    % form embedding matrix, M = I-W
    M = speye(n) - W;
    
end



%% FUNCTION: LOW-DIMENSIONAL ENCODING VIA LLE

function [Y, sigmaMemb] = lleEncoding( Memb, dimEnc )
    
    % number of embedding dimension (equal to the number of data points) 
    n = size( Memb, 1 );
    
    % compute the (d+1) smallest singular values and left singular vectors
    % of the embedding matrix 
    [Y, sigmaMemb, ~] = svds( Memb, dimEnc+1, 0 );
    
    % vectorize the singular values
    sigmaMemb = diag( sigmaMemb );
    
    % sort the singular values and vectors in asceding order
    [sigmaMemb, idx] = sort( sigmaMemb, 'ascend' );
    Y = Y(:,idx);
    
    % smallest singular value is 0; remove it and the corresponding vector
    Y(:,1)       = [];
    sigmaMemb(1) = [];
    
    % transpose to [D-by-N]
    Y = Y.';
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Tiancheng Liu                       tl137@duke.edu
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%   Xiaobai Sun                         xiaobai@cs.duke.edu
%
% VERSION
%
%   0.1 - February 23, 2015
%
% CHANGELOG
%
%   0.1 (Feb. 23, 2015) - Alexandros
%       * adapted from LLE function and its callees, by removing
%         any advanced/optional features
%
% ------------------------------------------------------------

