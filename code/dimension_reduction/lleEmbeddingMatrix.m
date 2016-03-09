function Memb = lleEmbeddingMatrix (X, idxNhood, varargin)
%
% LLEEMBEDDINGMATRIX - Locally Linear Embedding: Compute embedding matrix
%   
% SYNTAX
%
%   M = LLEEMBEDDINGMATRIX( X, NHOODS )
%   M = LLEEMBEDDINGMATRIX( X, NHOODS, 'constraint', CONSTRAINT, ... )
%   M = LLEEMBEDDINGMATRIX( X, NHOODS, 'options', OPTS, ... )
%
% INPUT
%
%   X           Feature vectors array                   [D-by-N]
%   NHOODS      Feature-space k-neighborhoods for each  [K-by-N]
%               point in X
%   CONSTRAINT  Local embedding vectors constraint for  [string]
%               unique solution
%               * 'algsum'      [ sum_{j} (w_{ij}) = 1 ]
%                 (see lclsAlgsum)
%               * 'norm2'       [ ||w_i||_2 = 1 ]
%                 (see lclsNorm2)
%               {default: 'algsum'}
%   OPTS        Additional optional arguments for the   [struct]
%               underlying solution function specified by 
%               CONSTRAINT)
%
% OUTPUT
%
%   M           Embedding matrix (I-W)                  [N-by-N; sparse]
%
% DESCRIPTION
%
%   M = LLEEMBEDDINGMATRIX(X,NHOODS) computes the embedding matrix
%   (I-W) for the LLE method [1].
%
%   M = LLEEMBEDDINGMATRIX(X,NHOODS,'constraint',CONSTRAINT,...) specifies
%   the constraint that is used for the solution of the local LS
%   minimization systems that yield the columns of W.
%
%   M = LLEEMBEDDINGMATRIX(X,NHOODS,'options',OPTS) specifies any
%   optional arguments that may be used by the function that solves
%   the locally constrained least squares (LCLS) problems. These
%   function are listed as dependencies below; refer to their
%   documentation for the corresponding lists of optional arguments.
%
% ALGORITHM
%
%   Each column (i) of the embedding matrix (I-W) is specified by
%   way of the local constrained LS minimization problem:
%
%     min_{w_i} = || X(:,i) - X(:,NHOODS(:,i)) * w_i ||_2
%
%     s.t. <constraint on w_i>
%
%   M = LLEEMBEDDINGMATRIX(...,'constraint','algsum') uses the
%   original LLE formulation [1], which sets the following constraint:
%
%                          sum(w_i) = 1
%
%   M = LLEEMBEDDINGMATRIX(...,'constraint','norm2') uses the
%   2-norm-based constraint [2]:
%
%                          ||w_i||_2 = 1
%
% REFERENCES
%
%   [1] S. T. Roweis and L. K. Saul, “Nonlinear Dimensionality
%   Reduction by Locally Linear Embedding,” Science, vol. 290,
%   no. 5500, pp. 2323–2326, Dec. 2000.
%
%   [2] T. Liu, A. S. Iliopoulos, and X. Sun, "Rank-adaptive,
%   Spatially Coherent Locally Linear Embeddings for Dimensionality
%   Reduction in Hyperspectral Image Analysis," Technical Report
%   Manuscript (in preparation).
%
% DEPENDENCIES
%
%   lclsAlgsum
%   lclsNorm2
%
%
% See also      neighborhoodsSpectral, lleLowdimEncoding, 
%               lclsAlgsum, lclsNorm2, lle
%
    
    
    %% DEFAULTS
    
    % costraint type
    constraint = 'algsum';
    
    % LCLS solution-method options
    options = struct;
    
    
    %% NAME-VALUE INPUT PARSING
    
    [constraint, options] = parseOptArgs( constraint, options, varargin{:} );
    
    
    %% INITIALIZATION
    
    % problem dimensions
    n = size( X       , 2 );    % # of data vectors
    k = size( idxNhood, 1 );    % size of local neighborhoods
    
    
    %% EMBEDDING MATRIX COMPUTATION
    
    % choose constraint type for local weights vector
    switch lower( constraint )
        
      % unit algebraic sum
      case {'algsum'}
        hsoln = @lclsAlgsum;
        
      % unit 2-norm
      case {'norm2'}
        error( 'lleEmbeddingMatrix:UnavailableNorm2', ...         % [TEMP]
               '''norm2''-constraint solver under development' ); % [TEMP]
        hsoln = @lclsNorm2;
        
      % error
      otherwise
        error( 'lleEmbeddingMatrix:InvalidConstraint', ...
               ['Unknown constraint type ''' constraint ''''] );
        
    end  % switch
        
    % compute the local weight vectors
    W = hsoln( X, idxNhood, options );
    
    % form global weights matrix in sparse format
    W = sparse( idxNhood(:), ...
                kron( (1 : n).', ones( k, 1 ) ), ...
                W(:), ...
                n, n );
    
    % form embedding matrix, M = I-W
    Memb = speye(n) - W;
    
    
end



%% LOCAL FUNCTION: OPTIONAL ARGUMENT PARSING

function [constraint, options] = parseOptArgs( constraint, options, varargin )
    
    % inputParser allows for MATLAB-native parsing of name-value inputs
    ip = inputParser;
    
    % do not throw error for unrecognized name-value pairs
    ip.KeepUnmatched = true;
    
    % do not expand structs
    ip.StructExpand = false;
    
    % string validation function
    validateString = @(x) assert( ischar(x), 'expected string'  );
    
    % specify recognized parameters and default values
    addParameter( ip, 'constraint', constraint, @(x) validateString(x) );
    addParameter( ip, 'options'   , options );
    
    % parse input arguments
    parse( ip, varargin{:} );
    
    % set output
    constraint = ip.Results.constraint;
    options    = ip.Results.options;
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% VERSION
%
%   0.2 - February 24, 2015
%
% CHANGELOG
%
%   0.2 (Feb. 24, 2015) - Alexandros
%       * optional arguments are now input as name-value pairs
%       * optional arguments to underlying solver must be packaged
%         in name-value-like struct
%
%   0.1 (Feb. 02, 2015) - Alexandros
%       * initial implementation
%       * global form of embedding matrix obtained directly from
%         sparse indexing of the local embedding vectors
%       * interface to local-weights-part of LLE by Tiancheng Liu
%       * two options for LCLS constraint (algebraic sum, 2-norm)
%
% ------------------------------------------------------------
