function [pars] = join_parameters(varargin)
% JOIN_PARAMETERS
%  Take two or more parameter structures and combine them into a single
%  parameter structure. 
% SYNTAX
%  [pars] = join_parameters(par1, par2);
%  [pars] = join_parameters(par1, par2, ... parN);
% DESCRIPTION
%  This function takes two or more parameter structures as input, and 
%  combines them into a single parameter.
%  Input and returned arguments are all of the same format as
%  GPars, having elements 'clip_name', 'par_name', 'inlsa_mos', and 'data'.
%  
%  If input pars do not have the exact same list of clips, or the inlsa_mos
%  values do not match, then warnings will occur and this function will
%  fail.  See also join_parameters_robust.
%
% Note:  Clips in each parameter argument must be in the same order.  That
% clip ordering will be preserved by 'pars'.


% Match up list of clips.  Abort if they don't match.
for cnt = 2:nargin,
    if length(varargin{1}.clip_name) ~= length(varargin{cnt}.clip_name),
        error('ERROR: all parameters must have the same set of clips.');
    end
    for loop = 1:length(varargin{1}.clip_name),
        if varargin{1}.clip_name{loop} ~= varargin{cnt}.clip_name{loop},
            error('ERROR: all parameters must have the same set of clips.');
        end
        if varargin{1}.inlsa_mos(loop) ~= varargin{cnt}.inlsa_mos(loop),
            if ~isnan(varargin{1}.inlsa_mos(loop)) & ~isnan(varargin{cnt}.inlsa_mos(loop)),
                error('ERROR: all parameters must have the same inlsa_mos values.');
            end
        end
    end
end

% find lengths
num_pars = 0;
for cnt = 1:nargin,
    num_pars = length(varargin{cnt}.par_name);
end
num_clips = length(varargin{1}.clip_name);

% initialize space for return variable.
try
	pars.par_name = cell(1,num_pars);
	pars.clip_name = varargin{1}.clip_name;
	pars.inlsa_mos = varargin{1}.inlsa_mos;
	pars.mos = varargin{1}.mos;
	pars.data = zeros(num_pars,num_clips);
catch
    error('join_parameters is out of memory.  This many parameters simply cannot be joined on your computer.');
end

% combine the parameter data.
total = 1;
for cnt = 1:nargin,
    this_length = length(varargin{cnt}.par_name);
    
    pars.par_name(total:(total+this_length-1)) = varargin{cnt}.par_name;
    pars.data(total:(total+this_length-1), :) = varargin{cnt}.data;
    
    total = total + this_length;
end
