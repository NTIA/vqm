function [result] = pars_find_par(par_struct, query, info, varargin)
% PARS_FIND_PAR
%  Search a parameter structure for a set of parameters, identified by
%  partial parameter names, and return the requested information (e.g., 
%  a parameter structure with fewer parameters, parameter index, etc.).
% SYNTAX
%  [result] = pars_find_par(par_struct, query, info);
%  [result] = pars_find_par(...,'PropertyName', ...);
% DESCRIPTION
%  "par_struct" is a parameter structure, of the same format as GPars.
%  "query" contains the partial parameter-name string.  query can be
%  either a string, or a cell array of multiple strings.
%  "info" is the type of information to extract from par_struct & return in result. 
%  info may be one of the following:
%       'par_name'
%       'index'
%       'par_struct'  
%
%  One of the following optional parameters may be specified:
%   'and'        Only find parameters matching all of the query strings.
%   'or'           Default.  Find parameters matching any of the query strings.
%   'nor'         Find parameters matching none of the query strings.


[num_pars, num_pts] = size(par_struct.data);
if ~iscell(query),
    query = {query};
end

if nargin > 4,
    error('pars_find_par:  Too many arguments specified\n');
end

% search for requested parameters.
% Find parameteres matching ANY of the strings.
want = zeros(num_pars,1);
if nargin <= 3 || strcmp(lower(varargin{1}),'or'),
    want = zeros(num_pars,1);
	for cnt=1:length(query),
        for i=1:num_pars,
            if size(strfind(par_struct.par_name{i},query{cnt}),2),
                want(i) = 1;
            end
        end
	end
% Find parameters matching ALL of the strings
elseif strcmp(lower(varargin{1}),'and'),
    want = ones(num_pars,1);
	for cnt=1:length(query),
        for i=1:num_pars,
            if ~size(strfind(par_struct.par_name{i},query{cnt}),2),
                want(i) = 0;
            end
        end
	end
% Find parameters matching NONE of the strings
elseif strcmp(lower(varargin{1}),'nor'),
    want = ones(num_pars,1);
	for cnt=1:length(query),
        for i=1:num_pars,
            if size(strfind(par_struct.par_name{i},query{cnt}),2),
                want(i) = 0;
            end
        end
	end
end

p_ind = [];
for i=1:num_pars,
    if want(i),
        p_ind = [p_ind i];
    end
end

if (size(p_ind,2) == 0)  % test for empty matrix
    disp('Requested data not found')
    result = 0;
    return
end

% pick out & return requested information.
switch lower(info)
    case {'par_name'},
        result = par_struct.par_name(p_ind);
    case {'index'}
        result = p_ind;
    case {'par_struct'}
        result.par_name = par_struct.par_name(p_ind);
        result.clip_name = par_struct.clip_name;
        result.inlsa_mos = par_struct.inlsa_mos;
        result.mos = par_struct.mos;
        result.data = par_struct.data(p_ind,:);
    otherwise  % user wants parameter data
        error('Request not recognized');
end



