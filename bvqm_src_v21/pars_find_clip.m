function [result, status] = pars_find_clip(par_struct, query_name, query_value, info)
% PARS_FIND_CLIP
%  Search a parameter structure for a set of clips (identified by test,
%  scene or hrc) and return the requested information (e.g., a parameter
%  structure with fewer clips, array with MOS, clip index, etc.)
% SYNTAX
%  [result, status] = pars_find_clip(par_struct, query_name, query_value, info);
% DESCRIPTION
%  "par_struct" is a parameter structure, of the same format as GPars.
%  "query_name" is one of the following:
%       'test'
%       'scene'
%       'hrc'
%       'clip'
%  "query_value" is a cell array listing the specific tests, scenes, hrcs, or clips whose 
%  information is to be extracted.  If only one value is wanted, then this
%  can be a string, instead.  To list all clip names in GClips, use:
%               strcat([GClips.test],'_',[GClips.scene],'_',[GClips.hrc])
%  "info" is the type of information to extract from par_struct & return in result. 
%  info may be one of the following:
%       'clip_name'
%       'inlsa_mos'
%       'mos'
%       'index'
%       'par_struct'  
%       '<parameter_name>' The actual name of the parameter to be extracted. 
%  The order of clips will not be modified by this function (i.e., order of
%  clips in par_struct and info will be identical).
%
%  'status' returns 0 if operation succeeds; 1 if "query_name" is not
%  recognized; 2 if requested clips cannot be found.
%
% Example call for pars_find_clip:
%       [result, status] = pars_find_clip(GPars, 'test', 'nist1', 'Y_si13_8x8_0.2s_std_12_ratio_loss_below5%_10%');
%       extracts all of the given parameter's data for test nist1.
%
%       [result, status] = pars_find_clip(GPars, 'test', { 'nist1' 's3' t1a1'} , 'par_struct');
%       returns a new parameter structure containing only clips from
%       nist1, s3, and t1a1.

status = 0;

[num_pars, num_pts] = size(par_struct.data);

%  Parse the clip names into three cell arrays: test, scene, hrc.
%  Eliminate the underscores that separate these quantities.
tests = cell(1,num_pts);
scenes = cell(1,num_pts);
hrcs = cell(1,num_pts);

% findstr does not work for padded matrices, must use loop.
for i = 1:num_pts
    this_name = par_struct.clip_name{i};
    last_ind = size(this_name,2);  % the number of chars in this string
	under_ind = findstr('_',this_name);  
    tests{i} = this_name(1:under_ind(1)-1);
    scenes{i} = this_name(under_ind(1)+1:under_ind(2)-1);
    hrcs{i} = this_name(under_ind(2)+1:last_ind);
end

if ~iscell(query_value),
    query_value = {query_value};
end

% pick off the requested clips
switch lower(query_name)
    case {'test'}
        q_ind = [];
        for i=1:length(query_value),
            q_ind = [q_ind find(strcmp(query_value{i},tests))];
        end
    case {'scene'}
        q_ind = [];
        for i=1:length(query_value),
            q_ind = [q_ind find(strcmp(query_value{i},scenes))];
        end
    case {'hrc'}
        q_ind = [];
        for i=1:length(query_value),
            q_ind = [q_ind find(strcmp(query_value{i},hrcs))];
        end
    case {'clip'}
        q_ind = [];
        for i=1:length(query_value),
            q_ind = [q_ind find(strcmp(query_value{i},par_struct.clip_name))];
        end
    otherwise
        % disp('Unknown query_name')
        status = 1;
        result = 0;
        return;
end
if (size(q_ind,2) == 0)  % test for empty matrix
    % disp('Requested data not found')
    status = 2;
    result = 0;
    return
end

% sort clips, so that order will be retained.
q_ind = sort(unique(q_ind));

% return result requested
switch lower(info)
    case {'clip_name'}
        result = par_struct.clip_name(q_ind);
    case {'inlsa_mos'}
        result = par_struct.inlsa_mos(q_ind);
    case {'mos'}
        result = par_struct.mos(q_ind);
    case {'index'}
        result = q_ind;
    case {'par_struct'}
        result.par_name = par_struct.par_name;
        result.clip_name = par_struct.clip_name(q_ind);
        result.inlsa_mos = par_struct.inlsa_mos(q_ind);
        result.mos = par_struct.mos(q_ind);
        result.data = par_struct.data(:,q_ind);
    otherwise  % user wants parameter data
        p_ind = find(strcmp(info,par_struct.par_name));
        if (size(p_ind,2) == 0)  % test for empty matrix
            % disp('Info name not found')
            status = 2;
            result = 0;
        else
            result = par_struct.data(p_ind,q_ind);
        end
end
