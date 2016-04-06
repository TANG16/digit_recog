function varargout = savesafe (filename, varargin)
%
% SAVESAFE - Save variables in fail-safe manner
%   
% SYNTAX
%
%   SAVESAFE( FILENAME, VAR1, VAR2, ... )
%   OUTFILE = SAVESAFE( FILENAME, VAR1, VAR2, ... )
%
% INPUT
%
%   FILENAME            Desired filename                [string]
%   VAR1, VAR2, ...     Variable names to be stored     [strings]
%
% OUTPUT
%
%   OUTFILE             Saved .mat filename             [string]
%
% DESCRIPTION
%
%   SAVESAFE(FILENAME,VAR1,VAR2,...) is essentially a wrapper for
%   MATLAB's SAVE function, with the following two features: If the
%   input FILENAME exists, then the suffix '_X' is appended to it,
%   where X is a unique serial number; and if the input variable set
%   is larger than 1GB, then that .mat file is stored using the
%   '-v7.3' option.
%   
%   OUTFILE = SAVESAFE(FILENAME,VAR1,VAR2,...) returns the filename of
%   the file where the variables where saved in OUTFILE.
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      save
%

% ------------------------------------------------------------
%
% AUTHOR
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% VERSION
%
%   0.2 - Novemember 19, 2013
%
% COPYRIGHT
%
%   Department of Computer Science, Duke University
%
% ------------------------------------------------------------
    
    
    %% PARAMETERS
    
    % large file option
    largenum = 2^30;    % 1GB
    largeopt = '-v7.3';
    
    
    %% INITIALIZATION 
    
    % convert input variable names from a cell array of strings into a
    % string of comma-delimited string literals
    strvars = cellfun( @(s) ['''' s ''', '], varargin, ...
                       'UniformOutput', false );
    strvars = [strvars{:}];
    strvars(end-1:end) = [];
    % strvars = ['''' strjoin(varargin, ''',''') ''''];
    
    
    %% SIZE OF VARIABLES
    
    % calculate the total number of bytes for the input variable names in
    % the caller function's workspace
    varstruct = evalin( 'caller', ['whos(' strvars ')'] );
    numbytes  = sum( cat( 1, varstruct.bytes ) );
    
    % use large .mat file compatibility mode?
    if numbytes >= largenum
        saveopt = largeopt;
    else
        saveopt = '';
    end
    
    
    %% SANITIZATION
    
    % remove extension from filename, if present
    [part_path, part_name, ~] = fileparts( filename );
    filename = [part_path filesep part_name];
    
    
    %% STORAGE
    
    % function to match filenames like 'filename.mat' or
    % 'filename_XXX.mat', where XXX is any number
    hMatchFiles = @(fn,fns) ...
        cellfun( @(s) regexp( s, ['^' fn '(|_[\d]*)[.]mat$'] ), ...
                 fns, 'UniformOutput', false );
    
    % function to count non-empty cell array elements (i.e. file matches)
    hCountFiles = @(c) nnz( ~cellfun( @isempty, c ) );
    
    % construct the output file name
    similarFiles = dir( [filename '*'] );
    similarFiles = { similarFiles.name };
    if isempty( similarFiles )
        n = 0;
    else
        n = hCountFiles( hMatchFiles( part_name, similarFiles ) );
    end
    outfile = [filename '_' num2str(n+1) '.mat'];
    
    % save the variables (which are to be found in the caller workspace)
    evalin( 'caller', ...
            ['save(''' outfile ''',' strvars ',''' saveopt ''')'] );
    
    % output out-filename?
    if nargout > 0
        varargout{1} = outfile;
    end
    
    
end
