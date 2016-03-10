function [ U, S, V , dflag] = svds_Laplace( B , dminSigma, flagDebug ) 
% 
%      
% 
%  [ U, S, V , dflag] = svds_Laplace( B , dminSigma, flagDebug ) ;
% 
%  
%  Compute DMINSIGMA SVD triplets, at the low end, of a large sparse 
%  matrix B in a more stable way than that via EIGS or SVDs 
%  
%  see test script : TEST_SVDS_LAPALACE  
% 
%  The left singular vector subspace U is given more attention and 
%  treatment in this function ; the right one V is used as a check 
%  the residual for the left subspace U' * B - S * V' 
%  
% Input 
% ------ 
%  B   : sparse and square matrix, real-valued, not necessarily symmetrc 
%        size( B,1) > 10 
%        sum( B, 1 ) = 0 
%        
%        When size(B) is large, use an utility routine to predict the 
%        sparsity of A = B * B' ; 
%        Alternatively, we may develop a version that takes 
%        customized routine for computing  v = ( BB')u 
% 
%  dminSigma: positive integer, to specify the number of the smallest 
%             singular values and their corresponding vectors to obtain 
% 
% Output 
% -------
% S :    the vector of the smallest singular vales as required, 
%        in non-decreasing order
% U, V : the left and right singular vectors corresponding to S 
% dflag : success-or-failure flag from EIGS 
% 
% flagDebug : for testing at every function revision 
% 
% dynamic rendering : 
% figures   
% messages  
% 
% Dependence 
% ------------------
% Callee functions 
% ----------------
% EIGS                                  % built-in sparse eigensolver 
% 
%                  The use of EIGS instead of SVDs is due to the 
%                  problem with the symmetric embedding with SVDs 
%                  which is good for the large singular values but 
%                  not for the small ones as required here 
%                  
%                  In particular, the symmetric embedding introduces 
%                  'shost' zero singular value triplets. 
% 
% QR               ise econmic version of it
% SVD              dense version for low-dim subpsace 
% 


%% Purpose of the Method  
% --------
%  To return the correct number of singular values at the low end and 
%    numerically correct values and vectors 
% 
%  To stabilize the singular vector subspaces 
% 
%  To prevent random flipping in the sign of singular vectors 
% 
%  To avoid the numerical problem with built-in SVDS for low-end 
%  singular values 
% 
%  While still utilize the built-in sparse eigensolvers 
% 

%% Two-steps method 
% --------
%  A. get the low-end of eigen-subspace of A = B * B'
%          numerically squared but not doublly mixed 
%          obtain the eigenspace A of higher dimension, 
%                2 * dminSigma 
% 
%  B. obtain the sparse SVD components of B via 
%          no longer suared, single sided 
%          gurantee the analytical zero singular value and the left vector 
%          fix the sign of vectors without random flipping 
% 
% *** No distribution without a permission from the programmer(s) 
%     listed at the end of the function 

%% ... Remaining numerical issues 
% 
%   numerical nodal domains of Fiedler function and higher harmonics 
%   in terms of numbers and shapes 
%   
%   the residual for the right singular space is relatively high, 
%   especially with V(:,1) ; need more detailed checking 
% 


%%  ===================== function body ================================
% 

%% ... get and set parameters 

fun_name = 'svds_Laplace' ;

fprintf( '\n   ... in %s \n\n', fun_name );

n  = size( B, 1 ); 

k  = 2 * dminSigma  ;         % double the subspace dimension
                              % for stabilization of the subspace 
                                      
kk = 5 ;                      % an internal parameter for determine 
                              % the sign of an eigenvector



%%  ... no change below unless revision is entailed ....  

if flagDebug 
    fprintf( '\n   in DEBUG mode ... ' );
end

%% ... Step A : ================================= 

fprintf( '\n ... %s in stage A \n', fun_name ) ;

%      get the least eigen-subspace of the sparse gram matrix 

A  = sparse( B * B' ) ;               % non-negative 
                           % sum(A,1) = 0, sum(A,2) = 0 
                           % no deflation because 
                           % the deflated matrix (A - 1/n) is full/dense  
                           
% A = A + A' ;               % to impose numerical symmetry as a record  
%                            % as implied in the 'sm' option 

esopts.issym  = 1;         % specify the symmetric and real-value properties 
esopts.isreal = 1;
esopts.maxit  = 60;        % reduce the default by 5X 

esopts.v0     = ones(n,1)/sqrt(n);
                           % override the default (random rhs)  
                           % use the theoretically constant vector as tne 
                           % first vector for the Lanczos process in EIGS 

                           
%% %  --- call EIGS solver 
%                          the choice of EIGS over SVDS is made 
%                          by the symmetric embedding used in SVDs 
%                          which may double the number of numerical 
%                          zero eigenvalues 

fprintf( '\n    %s calling built-in  EIGS \n', fun_name ) ;

[U, D2, dflag]  = eigs( A, k, 'sm' , esopts); 

D2      = real( diag(D2) ) ;
[D2, p] = sort( D2, 'ascend' ); 
U       = U(:,p); 

if dflag ~= 0 
    fprintf( [ fun_name, not all wanted eigenvalues converged in EIGS' ]); 
    disp(D2(1:kk));
end

fprintf( '\n    %s returning from EIGS \n', fun_name ) ;

if flagDebug                %% ... with display in figures  

 fprintf( '\n ... see pop-up figures \n' ) ; 
 
 figure 
 subplot(2,1,1) 
 plot( D2, 'o')             
                    % no log scale is used here because 
                    % numerically negative values may occur 
 xlabel( 'index to the squared singular values (in non-decreasing order) ' )

 u1   = U(:,1);
 err1 =  u1 - mean(u1) ; 

 subplot(2,1,2) 
 plot( log10 (abs( err1 ) + eps) ,'m.')
 xlabel( ' deviation in the first/leastt eigenvector from its mean ' ); 

 % ... subspace evaluation and refinement 
 %      introduce analytical eigenpair  
 %        (the zero-value-and-constant-vector pair) 
 %      re-orthogonalize the subspace 

 resEigen = norm( A * U - U*diag(D2), 'inf' ); 
 fprintf( '\n    the residual of AxU - UxD = %e \n ', resEigen ); 

end 

%% ... enforce the constant vector as the first one 

U  = [ -2*ones(n,1)/sqrt(n), U ] ;     
                      % insert and scale-by-2 to ensure the pivot position  
D3 = [0 ; D2 ];       % one precise zero, one numerial zero  
                      % the numerical zero to be pixoted out       


[ Q, R, p ] = qr( U, 0 ) ;             
                      % economic version with column pivoting 

if ( p(1) ~= 1 ) 
       error( [ fun_name, ': wrong padding' ] );
end 
     
U  = U(:,p);
D3 = D3(p) ;

if flagDebug 
    
  resEigen = norm( A * Q(:,1:k)*R(1:k,1:k) - Q(:,1:k)*R(1:k,1:k)*diag(D3(1:k)), 'inf' ); 
  fprintf( '\n    the residual of AxQxR - QxRxD = %e \n ', resEigen ); 

  figure
  imagesc( log10( abs(R) + eps )  ) 
  axis equal ; 
  colorbar
  title( 'Factor R in log scale' ); 
  ylabel( 'R(1,1) = 2, and R(:,end) = e1 in theory' ) ; 
  xlabel( 'row 1 and column 1 show numerical issues' ) ; 
end 

U = Q(:,1:k);                            % remove the redundant last column 
R = R(1:k,1:k); 
% R(1,1) = 1;                            % a check point, by theory 



%% ... Step B : ==================================
%               get the least singular value subspace of B : using B' 

fprintf( '\n ... %s in stage B \n', fun_name ) ;

B2 = B' * U;                              % B2(:,1) is numerically zero  

if flagDebug 
    fprintf( '\n   inf-norm of B2(:,1) = %e ', norm( B2(:, 1 ), 'inf') ) ; 
end

B2(:,1)  = zeros(n,1) ;                   % set B2(:,1) = e_1
B2(1,1)  = 1 + max( abs( B2(:) ) );       % in order to maintain domincant 
                                          % in pivoting 
[Q, R, p ]  = qr( B2, 0  );               % econmic version 

if p(1) ~= 1 | R(1,1) ~= B2(1,1) 
    error( [ fun_name, ' : wrong pivoting' ] ) ; 
end
B2(1,1) = 0 ; 
B2      = B2(:,p) ; 
U       = U(:,p) ; 
U(:,1) = ones(n,1)/sqrt(n);

% norm( Q(2:k,1), 'inf')   is numerical zero 

Q(:,1) = zeros(n,1);
Q(1,1) = 1; 

R = triu( R ); 
R(1,1) = 0; 

if flagDebug 
    
 fprintf( '\n   orhtogonality checking before the svd of R ' ) 
  [    norm( sum( U(:, 2:k), 1 ), 'inf' ) ... 
       norm( speye(k) - U'*U, 'inf' ),    ... 
       norm( speye(k) - Q'*Q, 'inf' ) ] 
   
 fprintf( '\n   Q^T B2 - R in inf norm = %e \n', norm( Q' * B2 - R, 'inf' ) ) ; 

end 



%% ... get the svd of B via the svd of R 

flagPad  = 0 ;

if flagPad 
    [ Vk, Sk, Uk ] = svd( R(:,2:k) ); 
else
    [ Vk, Sk, Uk ] = svd( R );      
 
end
% note the sides of Uk and Vk 
                                         
Sk        = diag( Sk ); 
[ Sk, p ] = sort( Sk, 'ascend' ) ;

Vk        = Vk(:,p);
Uk        = Uk(:,p); 

% ... complete Uk, Sk and Vk 

 
if flagPad    % R2(:,2:k) is used for the dense svd instead of R 
              % to complete Uk, Sk and Vk with one more dimension     
     
  Uk = [1, zeros(1,k-1) ; ... 
        zeros(k-1,1), Uk ] ;     
  Sk = [ 0 ; Sk ] ;

  P1 = eye(k) - Vk*Vk';                % projection to the complementary
  [ dmax, jmax] = max( diag(P1) ) ;    % CHOL fails to do the factorization
  v1 = P1(:, jmax );
  v1 = v1/ norm(v1) ; 

  Vk = [ v1, Vk ] ;
end 

if flagDebug 
   fprintf( '\n   the smallest singular value = %e \n',  Sk(1) ); 
   fprintf( '\n   Vk(:,1)^T R in inf norm = %e \n',  norm( Vk(:,1)'*R  , 'inf' ) ) ; 
   fprintf( '\n   Vk(:,1)^T x Vk(:,2:k) in inf norm = %e \n', norm( Vk(:,1)'*Vk(:,2:k) ) );
end

    

%% ... output rendering, 

S  =  Sk( 1:dminSigma ) ; 

U  = U * Uk ;  
U  = U(:, 1:dminSigma ); 

V  = Q * Vk ;                    % Q is the right-side space of B' 
V  = V(:, 1:dminSigma );


% ... error checking 

if flagDebug 
  fprintf( '\n\n   errors after subspace refinements : U-V orthogonalities \n\n' ); 

  [   norm( speye( dminSigma ) - U'*U, 'inf' ),    ... 
      norm( speye( dminSigma ) - V'*V, 'inf' ) ] 

  fprintf( '\n\n   The difference in the smalles  %d singular values', kk ) ; 
  format long 
  D1 = sqrt( D2 ); 
  [ S(1:kk) - D1(1:kk)  ] 

  svdRes = norm( U' * B * V - diag(S)  , 'inf' );
  fprintf( '\n\n   The residual : norm(  U^T B V - S , inf) = %e ', svdRes ); 
  
  svdRes = norm( U' * B - diag(S) * V' , 'inf' );
  fprintf( '\n\n   The residual : norm(  U^T B - S  V^T, inf) = %e ', svdRes ); 
  
  svdRes = norm(  B * V(:,2:end)- U(:,2:end)* diag(S(2:end)) , 'inf' );
  fprintf( '\n\n   The residual : norm(  B  V - U  S, inf) = %e ', svdRes ); 
  
  
  fprintf( '\n\n   ... see pop-up figures \n\n' );  
  
  
 figure

 subplot(2,1,1)
 semilogy( S+eps, 'o' ) ;  
 xlabel( 'the refined singular values' ) ; 

 subplot(2,1,2) 
 D1 = D1(1:dminSigma) ; 
 semilogy( abs( S - D1) , 'mo' ) ; 
 xlabel( ' the difference in singular values' ) ;

 figure 

 subplot(2,1,1) 
 plot( log10( abs( U(:,1) - mean(U(:,1)) )   ) , 'm.') ; 
 xlabel( 'Deviation of U(:,1), in magnitude at log10 scale, from mean ' ) 
 % there seems some randomness in the least bit 


 subplot(2,1,2) 
 plot( log10( abs( V(:,1) ) )  , '.' ) ; 
 xlabel( 'V(:,1) in magnitude at log10 scale' ) 

end 


S = diag( S );  % vector --> diagonal matrix (for seamless
                % replacement of SVDS)

%%  ... stablize the vector signs for visual rendering 
%      with vectors signs determined by the algebraic 
%      sum of the kk-largest elements in magnitude 

for jk = 1 : dminSigma  
    
   utemp      = U(:, jk ) ; 
   [umax, jmax] = max( abs( utemp ) ) ; 
   jknn = knnsearch ( utemp, utemp( jmax(1) ), 'k', kk ) ;
   
   jknnsum    = sum( utemp( jknn ) )  ;              % guaranteer nonzero 
   U(:, jk)   = utemp   * sign( jknnsum ) ;
   V(:, jk)   = V(:,jk) * sign( jknnsum ) ;
   
end


return 

%%  ======================================================================
%  Programmers 
% 
%   Xiaobai Sun ( initial draft )  
% 
%   Modification and feedback comments by 
%   Alexandros S. Iliopoulos 
%   Tiancheng Liu 
%   Nikos Pitsianis 
%   
%  Duke University 
%  Department of Computer Science 
% 
% 

% ======================================================================
%  Initial Draft April 10, 2014 
% 
%  ----------------------- revision history ----------------------------
% 
%  Revision on Aug. 23, 2015 
%       recheck the initial vector option, add additional test 
%       on singular invariant subspace, test up to 160,000 size, 
%       complete documentation 
% 
%  Revision on Aug. 20, 2015 
%       verified the correct file name suffix 
%       reordered a couple of lines related to zero value 
%                  and constant values 
%       to perform indepdent test (with a different interface) 
%
%  Revision on Aug. 14, 2015 
%       let the rhs generated randomly 
%       make the firt vector idealy constant 
%       remaining numerical issues with 2nd and so on 
% 
%  Revision April 12, 2015 
%    The remaining issue is resolved : the analytical zero 
%    and the corresponding constant left s-vector can now 
%    be used to descriminate from the other numerical zeros 
% 
%    Technically, the full svd does not return the orthogonal 
%    vector matrices as it claims; complementary completion 
%    is done .
% 
%    CHOL fails to factor a rank-1 orthogonal projector 
%    
%    A couple of other analytical-numerical issues are resolved 
%    as well, see the details in code 
% 
%  Revision on April 11, 2015 
%    -- bug in sparse casting mechanism 
%    -- add options to EIGS for real, symmetric matrix 
%    -- change 'SA' to 'SM' to EIGS call ( 'SA' is BAD ) 
%    -- modify the input from W to I-W 
%    -- stabilize the computatioon within the subspace 
%         
%  *** Revision in Dec.  2014 
%    -- as a converged version between sparse and non-sparse 
%        version, via kNN graph embedding 
% 
%    -- distinguished from the correlation matrix approach 
%        which is dense, linear, and lack of theoretical 
%        footing, also different from the so called sparsification 
% 
%  Remaining issue : 
%    The rank-1 deflation technique can not be used 
%    for the sparse case -- > the left singular 
%    vectors suffer from clusters at the zero as a result ; 
% 

% 


% 
