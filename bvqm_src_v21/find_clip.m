function offset = find_clip(clip_structs, test, scene, hrc, varargin)
% FIND_CLIP
%  Find a clip's offset, given the test, scene, and hrc.
% SYNTAX
%  offsets = find_clip(clip_structs, test, scene, hrc);
%  offsets = find_clip(...,'PropertyName', ...);
% DESCRIPTION
%  'clip_structs' is an array of clips, all of the same format as GClips.
%  The return value, 'offsets', is an array of offsets of the requested clips.
%
%  'test', 'scene', and 'hrc' are the test, scene, and hrc names of the
%  clip you want to find.  Each of these can be one string, a cell array
%  listing multiple strings, or the string '*' for a wildcard matching of
%  all names.  
%
%  If a specified clip cannot be found, and error will occur.
%
% The following optional properties may be requested:
%   'not'   Return the complement of the clips meeting the above criteria,
%
%  Note:  'offsets' will be sorted according to the ordering of clip_structs.

% detect complement
do_complement = 0;
if nargin > 4,
    if nargin == 5 & strcmp(lower(varargin{1}),'not'),
        do_complement = 1;
    else
        error('Function argument problem, optional property not recognized');
    end
end


% split into individual requests, and concatonate return results.
if iscell(test),
    offset = [];
    for cnt = 1:length(test),
        offset = [offset find_clip(clip_structs, test{cnt}, scene, hrc)];
    end
    offset = sort(offset);
	% handle complement
	if do_complement,
            offset = setdiff(1:length(clip_structs), offset);
	end
    return;
elseif iscell(scene),
    offset = [];
    for cnt = 1:length(scene),
        offset = [offset find_clip(clip_structs, test, scene{cnt}, hrc)];
    end
    offset = sort(offset);
 	% handle complement
	if do_complement,
            offset = setdiff(1:length(clip_structs), offset);
	end
   return;
elseif iscell(hrc),
    offset = [];
    for cnt = 1:length(hrc),
        offset = [offset find_clip(clip_structs, test, scene, hrc{cnt})];
    end
    offset = sort(offset);
	% handle complement
	if do_complement,
            offset = setdiff(1:length(clip_structs), offset);
	end
    return;
end


% pick out test, scene & hrc names.
all_test = [clip_structs.test];
all_scene = [clip_structs.scene];
all_hrc = [clip_structs.hrc];

% find a specific clip's offset
if ~strcmp('*',test) & ~strcmp('*',scene) & ~strcmp('*',hrc),
	offset = find(strcmp(test,all_test) & strcmp(scene,all_scene)  & strcmp(hrc,all_hrc));
    
% one wild card
elseif strcmp('*',test) & ~strcmp('*',scene) & ~strcmp('*',hrc),
	offset = find(strcmp(scene,all_scene) & strcmp(hrc,all_hrc));
elseif ~strcmp('*',test) & strcmp('*',scene) & ~strcmp('*',hrc),
	offset = find(strcmp(test,all_test) & strcmp(hrc,all_hrc));
elseif ~strcmp('*',test) & ~strcmp('*',scene) & strcmp('*',hrc),
	offset = find(strcmp(test,all_test) & strcmp(scene,all_scene) );
    
% two wild cards
elseif ~strcmp('*',hrc),
	offset = find(strcmp(hrc,all_hrc));
elseif ~strcmp('*',scene),
	offset = find(strcmp(scene,all_scene));
elseif ~strcmp('*',test),
	offset = find(strcmp(test,all_test) );
    
% can't have all wild cards.
else
    error('Test, scene, and hrc cannot all be wildcards.  This would be meaningless\n');
end


% handle complement
if do_complement,
        offset = setdiff(1:length(clip_structs), offset);
end

% sort offsets
offset = sort(offset);
