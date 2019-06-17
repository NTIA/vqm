function [clip_test, clip_scene, clip_hrc] = pars_parse(par_struct);
% PARS_PARSE
%  Parse parameter's clip names (test_scene_hrc) to get the name of the
%  test, scene, and hrc for each clip.
% SYNTAX
%   [clip_test, clip_scene, clip_hrc] = pars_parse(par_struct);
% DEFINITION
%   Take a parameter structure, 'par_struct', of the same format as GPars,
%   and parse into test, scene, and hrc names for each clip.

num_clips = length(par_struct.clip_name);

%  Parse the clip names into three cell arrays: test, scene, hrc.
%  Eliminate the underscores that separate these quantities.
clip_test = cell(1,num_clips);
clip_scene = cell(1,num_clips);
clip_hrc = cell(1,num_clips);
% findstr does not work for padded matrices, must use loop.
for i = 1:num_clips
    this_name = par_struct.clip_name{i};
    last_ind = size(this_name,2);  % the number of chars in this string
	under_ind = findstr('_',this_name);  
    clip_test{i} = this_name(1:under_ind(1)-1);
    clip_scene{i} = this_name(under_ind(1)+1:under_ind(2)-1);
    clip_hrc{i} = this_name(under_ind(2)+1:last_ind);
end
