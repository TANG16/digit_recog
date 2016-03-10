function X = lclsAlgsum (A, nidx, opts)
%
% LCLSALGSUM - Local-neighborhood, constrained least squares
%              fitting minimization subjecto to unit algebraic sum
%              constraint 
%   
% SYNTAX
%
%   X = LCLSALGSUM( A, NHOODS )
%   X = LCLSALGSUM( A, NHOODS, OPTS )
%
% INPUT
%
%   A           Data matrix                             [D-by-N]
%   NHOODS      Local K-neighborhood indices            [K-by-N]
%   OPTS        Optional arguments struct               [struct]
%   OPTS.REG    Regularization (eigenspectrum           [scalar]
%               lifting) factor
%               {default: 1e-3}
%   OPTS.TAU    Numerical threshold for singular value  [scalar]
%               magnitude
%               {default: sqrt(eps)}
%
% OUTPUT
%
%   X           Local solution vectors                  [K-by-N]
%
% DESCRIPTION
%
%   X = LCLSALGSUM(A,NHOODS,...) solves the following set of
%   local-neighborhood, constrained least squares fitting problems:
%
%     min_{x_i}  || A(:,i) - sum_j{ A(:,NHOODS(:,i)) * x_ij } ||_2
%     s.t.       sum_j{ x_ij } = 1
%
%   for all i (data vector index). 
%
% ALGORITHM
%
%   The X vectors are obtained by minimizing the corresponding
%   Lagrangian form; the latter easily leads to the following
%   (regularized) system:
%
%                   Ci + (r*tr(Ci)) * I = lambda_i * e
%
%   where C = Gram{A(:,i) - A(:,NHOODS(:,i))} is the "local,
%   Ai-centered correlation matrix", and r (OPTS.REG) is a
%   regularization factor, set to 0 if K <= D. More details may be
%   found in [1].
%
% REFERENCES
%
%   [1] S. T. Roweis and L. K. Saul, “Nonlinear Dimensionality
%   Reduction by Locally Linear Embedding,” Science, vol. 290,
%   no. 5500, pp. 2323–2326, Dec. 2000.
%
%
% See also      lclsNorm2, lleEmbeddingMatrix, lle
%
    
    
    %% PARAMETERS
    
    % eigenspectrum lift based solution for singular matrices?
    flagEigenlift = true;
    
    
    %% DEFAULTS
    
    % set emtpy optional-arguments struct if missing
    if ~exist( 'opts', 'var' ) || isempty( opts )
        opts = struct;
    end
    
    % regularization factor
    reg = 1e-3;
    
    % singular magnitude threshold
    tau = 1e-12; %sqrt( eps );
    
    
    %% OPTIONAL ARGUMENT PARSING
    
    [reg, tau] = parseOptArgs( reg, tau, opts );
    
    
    %% INITIALIZATION
    
    % problem dimensions
    n = size( A   , 2 );        % # of data vectors
    d = size( A   , 1 );        % # of dimensions
    k = size( nidx, 1 );        % size of local neighborhoods
    
    % pre-allocate space for local solution vectors
    X = zeros( k, n );
    
    
    %% LOCAL SYSTEM SOLUTIONS (SINGULAR OR NON-SINGULAR?)
    
    % iterate across all local systems (equivalently, data vectors)
    for i = 1 : n
        
        % create local neighborhood array and center neighbors
        % around reference (i-th) vector
        Ain = bsxfun( @minus, A(:,nidx(:,i)), A(:,i) );
        
        % check norm of local matrix row-sum
        normLocal = norm( sum( Ain, 2 ) );
        
        % if norm is too low, use solution for numerically singular system
        if (k > d) || (normLocal <= tau) || (rank( Ain, tau ) < k)
            if flagEigenlift
                hsoln = @(B) solveQPnonsingular( B, reg );
            else
                hsoln = @(B) solveQPsingular( B, sqrt(eps) );
            end
        else
            reg = 0;
            hsoln = @(B) solveQPnonsingular( B, reg );
        end
        
        % solve local system
        X(:,i) = hsoln( Ain );
        
    end  % for (i)
    
    
end



%% LOCAL FUNCTION: NON-SINGULAR LOCAL QP SOLUTION

function X = solveQPnonsingular( A, reg )
    
    % size of local neighborhood
    k = size( A, 2 );
    
    % create quadratic form of local matrix
    G = A' * A;
    
    % regularization (eigenspectrum lift)
    G = G + reg * trace(G) * eye(k);
    
    % solve local system and scale to unit algebraic sum
    X = G \ ones(k,1);
    X = X ./ sum( X );
    
end



%% LOCAL FUNCTION: SINGULAR LOCAL QP SOLUTION

function X = solveQPsingular( A, tau )
    
    % size of local neighborhood
    k = size( A, 2 );
    
    % set infinity norm of A to be equal to k (i.e. same as row-sum of
    % affine constraint vector)
    A = A / norm( A, Inf );
    
    % form expanded equality constraint matrix
    % (affine constraint & local homogeneous condition)
    C = vertcat( ones(1,k), A'*A );
    
    % solution vector is the first column of the pseudo-inverse of
    % the expanded equality constraint matrix
    [U,sigma,V] = svd( C );
    sigma( sigma <= tau ) = min( sigma( sigma > tau ) );
    X = V' * (1./sigma)' * U;
    
    X = X(:,1);
    
end



%% LOCAL FUNCTION: OPTIONAL ARGUMENT PARSING

function [reg, tau] = parseOptArgs( reg, tau, opts )
    
    % inputParser allows for MATLAB-native parsing of name-value inputs
    ip = inputParser;
    
    % expand struct to name-value pairs
    ip.StructExpand = true;
    
    % specify recognized parameters and default values
    addParameter( ip, 'reg', reg );
    addParameter( ip, 'tau', tau );
    
    % parse input arguments
    parse( ip, opts );
    
    % set output
    reg = ip.Results.reg;
    tau = ip.Results.tau;
    
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
%   0.4 - April 04, 2015
%
% CHANGELOG
%
%   0.4 (Apr 04, 2015) - Alexandros
%       * separate solution approaches for numerically singular and
%         non-singular local neighborhood matrices
%
%   0.3 (Feb 24, 2015) - Alexandros
%       * optional arguments are now input as struct fields
%
%   0.2 (Feb 02, 2015) - Alexandros
%       * created LCLSALGSUM as an independent module
%       * allowed for variable regularization factors
%
%   0.1 (Jan 20, 2015) - Tiancheng
%       * initial implementation as part of LLE function
%
% ------------------------------------------------------------

