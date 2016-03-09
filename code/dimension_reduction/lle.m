function [Y, Memb, idxNhood, sigmaMemb, convFlag] = ...
        lle (X, varargin)
%
% LLE - Low-dimensional manifold encoding via Locally Linear Embeddings 
%   
% SYNTAX
%
%   [Y, MEMB] = LLE( X, K, DIMENC )
%   [Y, MEMB] = LLE( X, K, DIMENC, 'Name', Value, ... )
%   [Y, MEMB, NHOOD] = LLE( ... )
%   [Y, MEMB, NHOOD, SIGMAMEMB, CONVFLAG] = LLE( ... )
%
% INPUT
%
%   X           Feature vectors array in ambient space  [D-by-N]
%   K           Size of local neighborhoods             [scalar]
%   DIMENC      LLE-based encoding dimensionality       [scalar]
%   NHOODDIST   Distance metric for neighbor-finding    [string]
%               (see knnsearch)
%   NAME/VALUE  Name-value pairs for optional arguments of the
%               sub-functions of LLE:
%               * neighborsSpectral
%               * lleEmbeddingMatrix
%               * lleLowdimEncoding
%
% OUTPUT
%
%   Y           Low-dimensional encoding of X vectors   [DENC-by-N]
%   MEMB        LLE matrix (MEMB = I - W)               [N-by-N; sparse]
%   NHOOD       Ambient-space neighborhood indices      [K-by-N]
%   SIGMAMEMB   Singular values of MEMB                 [DENC-by-1]
%   CONVFLAG    Eigensolver convergence flag            [boolean]
%
% DESCRIPTION
%
%   [Y,MEMB,...] = LLE(X,K,DIMENC) returns as MEMB the K-neighbors LLE
%   matrix and as Y a DIMENC-dimensional encoding of the feature
%   vectors in X [1].
%
%   [Y,MEMB,...] = LLE(X,K,DIMENC,'Name',Value,...) specifies
%   name-value pairs for the optional arguments of the constituent
%   functions (lister as dependencies below) that LLE serves as a
%   wrapper for. The corresponding optional arguments must be
%   unambiguously named.
%
% REFERENCES
%
%   [1] S. T. Roweis and L. K. Saul, ???Nonlinear Dimensionality
%   Reduction by Locally Linear Embedding,??? Science, vol. 290,
%   no. 5500, pp. 2323???2326, Dec. 2000.
%
% DEPENDENCIES
%
%   neighborsSpectral
%   lleEmbeddingMatrix
%   lleLowdimEncoding
%
%
% See also      neighborsSpectral, lleEmbeddingMatrix,
%               lleLowdimEncoding, knnsearch
%
    
    
    %% ENSURE OPTIONAL ARGUMENTS ARE INPUT AS NAME-VALUE PAIRS
    
    if (nargin == 2)
        optargs = struct2namevalue( varargin{1} );
    else
        optargs = varargin;
    end
    
    
    %% LOCAL NEIGHBORHOODS
    
    idxNhood = neighborsSpectral( X, [], optargs{:} );
    
    
    %% LLE MATRIX
    
    Memb = lleEmbeddingMatrix( X, idxNhood, optargs{:} );
    
    
    %% LOW-DIMENSIONAL ENCODING
    
    [Y, sigmaMemb, convFlag] = lleLowdimEncoding( Memb, optargs{:} );
    
    
end
  


%% LOCAL FUNCTION: STRUCT --> NAME-VALUE PAIRS

function nvcell = struct2namevalue( strct )
    
    % struct field names
    fnames = fieldnames( strct );
    
    % struct field contents
    fvals  = struct2cell( strct );
    
    % interleave the field names and corresponding contents of the struct
    nvcell = cellfun( @(x,y) {x, y}, fnames, fvals, ...
                      'UniformOutput', false );
    
    % flatten the cell array into a cell-vector
    nvcell = cat( 2, nvcell{:} );
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Tiancheng Liu                       tl137@duke.edu
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% VERSION
%
%   0.5 - February 24, 2015
%
% CHANGELOG
%
%   0.5 (Feb. 24, 2015) - Alexandros
%       * renamed the function from LLE_NEW to LLE
%       * all optional arguments are expected as unambiguous
%         name-value pairs 
%
%   0.4 (Feb. 09, 2015) - Alexandros
%       * added optional input for selection of neighbor-finding
%         distance metric and for eigensolver selection (internal to
%         LLELOWDIMENCODING) 
%
%   0.3 (Feb. 06, 2015) - Alexandros
%       * LLE_NEW now also returns the neighborhood indices used in
%       the LLE computations
%
%   0.2 (Feb. 02, 2015) - Alexandros
%       * created LLE_NEW as a top-level, modular wrapper
%
%   0.1 (Jan. 20, 2015) - Tiancheng
%       * initial implementation
%
% ------------------------------------------------------------
