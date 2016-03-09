function [Y, sigmaMemb, convFlag] = lleLowdimEncoding (Memb, d, varargin)
%
% LLELOWDIMENCODING - Low-dimensional encoding of locally linear
%                     embedding (manifold) structure 
%   
% SYNTAX
%
%   Y = LLELOWDIMENCODING( M, 'dim', D )
%   Y = LLELOWDIMENCODING( M, 'dim', D, 'svds', SVDSFLAG, ... )
%   Y = LLELOWDIMENCODING( M, 'dim', D, 'zerothres', TAUZERO, ... )
%   [Y, SIGMAMEMB] = LLELOWDIMENCODING( M, 'dim', D, ... )
%   [Y, SIGMAMEMB, CONVFLAG] = LLELOWDIMENCODING( M, 'dim', D, ... )
%
% INPUT
%
%   M           LLE matrix (M = I-W)                    [N-by-N]
%   D           Manifold encoding dimensionality        [scalar]
%   SVDSFLAG    Use SVDS (otherwise EIGS) internally    [boolean]
%               {default: true}
%   TAUZERO     Threshold for removing numerical-zero   [scalar]
%               singular values and corresponding vectors
%               {default: 0 (no removal)}
%
%
% OUTPUT
%
%   Y           D-dimensional encoding of points on     [D-by-N]
%               the manifold described by M
%   SIGMAMEMB   Squared singular values of the          [1-by-N]
%               embedding matrix
%   CONVFLAG    SVDS/EIGS function convergence flag     [boolean]
%               (0: converged -- 1: did not converge)
%
% DESCRIPTION
%
%   Y = LLELOWDIMENCODING(M,D) returns a collection of N D-dimensional
%   vectors that preserve (in the 2-norm sense) the locally linear
%   embedding structure of the manifold described by M. The solution
%   is the D+1 smallest left singular vectors of M. (Singular values
%   and thei corresponding vectors are sorted in ascending order.)
%
%   The zero-eigenvalue vectors are not removed. if the
%   zero-eigenvalue is simple, then Y(1,:) is the constant vector.
%
%   Y = LLELOWDIMENCODING(M,D,'svds',false,...) explicitly forms the
%   Gramian of the embedding matrix (M*M') and computes the truncated
%   EVD of it, using the EIGS solver. This method may be slow for
%   large M matrices. Note, also, that the EIGS output may contain
%   non-zero imaginary parts.
%
%   Y = LLELOWDIMENCODING(M,D,'zerothres',TAUZERO) removes any
%   singular values that are below the specified threshold, as well as
%   the corresponding singular vectors. It is expected that there will
%   be a single zero singular value, hence the internal call to
%   SVDS/EIGS requests (D+1) modes.
%
%   [Y,SIGMAMEMB] = LLELOWDIMENCODING(M,D,...) also returns the
%   singular values of the embedding matrix.
%
%   Note that if EIGS is used instead of SVDS (default), then the
%   singular values are computed as the square roots of eigenvalues of
%   Gram(M); therefore, the magnitude of numerical error in the EIGS
%   computations will be halved.
%
%   [Y,SIGMAMEMB, CONVFLAG] = LLELOWDIMENCODING(M,D,...) also returns
%   a flag that specifies whether the SVDS/EIGS function (called
%   internally) converged to the returned eigenvalues.
%
% NOTES
%
%   The returned data (Y), which are in the low-dimensional (D) space,
%   are decorrelated; that is, the following should hold: Y * Y' = I.
%
%   It is assumed that the input LLE matrix is real.
%
%   There is an issue with the direction of the returned singular
%   vectors; that is, the signs of the SVDS/EIGS-returned vectors are
%   non-deterministic.
%
%   Also, LLELOWDIMENCODING currently suffers from an inconsistent
%   scaling in the encoding space.
%
% REFERENCES
%
%   [1] S. T. Roweis and L. K. Saul, ???Nonlinear Dimensionality
%   Reduction by Locally Linear Embedding,??? Science, vol. 290,
%   no. 5500, pp. 2323???2326, Dec. 2000.
%
%
% See also      lleEmbeddingMatrix, lle, eigs
%
    
    
    %% DEFAULTS
    
    % SVDS vs EIGS selection
    svdsFlag = true;
    
    % threshold for zero singular value removal
    tauZero = 0;
    
    
    %% NAME-VALUE INPUT PARSING
    
    [d, svdsFlag, tauZero, esOpts] = parseOptArgs( svdsFlag, tauZero, ...
                                                   varargin{:} );
    
    
    %% INITIALIZATION
    
    % number of data vectors
    n = size( Memb, 1 );
    
    
    %% LOW-DIMENSIONAL ENCODING
    
    % compute the d smallest singular values and left singular vectors of
    % the embedding matrix
    if svdsFlag
        % [Y, sigmaMemb, ~, convFlag] = svds( Memb, d+1, 0, esOpts );
        [Y, sigmaMemb, ~, convFlag] = svds_Laplace( Memb, d, false );
    else
        MembSq = Memb * Memb';
        [Y, lambdaMemb2, convFlag] = eigs( MembSq, d, 0, esOpts );
        sigmaMemb = sqrt( lambdaMemb2 );
    end
    
    % vectorize the singular values
    sigmaMemb = diag( sigmaMemb );
    
    % sort the singular values and vectors in asceding order
    [sigmaMemb, idx] = sort( sigmaMemb, 'ascend' );
    Y = Y(:,idx);
    
    % remove numerical-zero singular values
    maskZero = (sigmaMemb < tauZero);
    Y(:,maskZero)       = [];
    sigmaMemb(maskZero) = [];
    
    % transpose to [D-by-N]
    Y = Y.';
    
    
end



%% LOCAL FUNCTION: OPTIONAL ARGUMENT PARSING

function [d, svdsFlag, tauZero, esOpts] = ...
        parseOptArgs( svdsFlag, tauZero, varargin )
    
    % inputParser allows for MATLAB-native parsing of name-value inputs
    ip = inputParser;
    
    % do not throw error for unrecognized name-value pairs
    ip.KeepUnmatched = true;
    
    % specify mandatory arguments
    addRequired( ip, 'dim' );
    
    % specify recognized parameters and default values
    addParameter( ip, 'svds'     , svdsFlag );
    addParameter( ip, 'zerothres', tauZero  );
    
    % parse input arguments
    parse( ip, varargin{:} );
    
    % set output
    d        = ip.Results.dim;
    svdsFlag = ip.Results.svds;
    tauZero  = ip.Results.zerothres;
    if svdsFlag
        esOpts.isreal = 1;
    else
        esOpts.isreal = 1;
        esOpts.issym  = 1;
    end
    
    % if zero-threshold is non-zero, increment spectral dimension by 1
    if tauZero > eps
        d = d + 1;
    end
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Tiancheng Liu                       tl137@duke.edu
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%   Nikos Pitsianis                     nikos@cs.duke.edu
%
% VERSION
%
%   0.5 - April 16, 2015
%
% CHANGELOG
%
%   0.5 (Apr 17, 2015) - Alexandros
%       * integrated Xiaobai's SVDS_LLEW (stabilized SVDS for the LLE
%         embedding matrix) function into LLELOWDIMENCODING
%       * added optional argument for controlling removal of
%         numerical-zero singular values and corresponding vectors
%
%   0.4 (Feb 25, 2015) - Alexandros
%       * optional arguments are now input as name-value pairs
%       * eigensolver flag now is true for SVDS (used to be true
%         for EIGS)
%
%   0.3 (Feb 09, 2015) - Nikos & Alexandros
%       * replaced default internal call to EIGS (using the Gramiam
%         of the embedding matrix) to call to SVDS (using the
%         original embedding matrix)
%       * added option to choose between SVDS and EIGS
%       * explicitly sorted singular values, just in case
%       * removed singular vector scaling by sqrt(N)
%
%   0.2 (Feb 02, 2015) - Alexandros
%       * created LLELOWDIMENCODING as an independent module
%
%   0.1 (Jan 20, 2015) - Tiancheng
%       * initial implementation as part of the LLE function
%
% ------------------------------------------------------------
