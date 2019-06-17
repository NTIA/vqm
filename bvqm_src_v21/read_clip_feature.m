function [data] = read_clip_feature(clip_struct, feature_base_dir,feature_name);
% READ_CLIP_FEATURE
%   Read one clip's feature stream.
% SYNTAX
%   [data] = read_clip_feature(clip_struct, feature_base_dir,feature_name);
% DESCRIPTION
%   Given one clip, 'clip_struct' (in the same format as GClips), a
%   directory path ('feature_base_dir') and the name of the feature 
%   (i.e., subdirectory) within that directory, read & return the specified
%   clip's features. 


% error checking.
if length(clip_struct) ~= 1,
    error('clip_struct must contain exactly one clip');
end

% load features. error if not computed.
load( sprintf('%s/%s/%s_%s_%s.mat',feature_base_dir,feature_name, ...
    clip_struct.test{1}, clip_struct.scene{1}, clip_struct.hrc{1}) );

% check whether feature is available.
[r,c,t] = size(data);
if r*c*t == 0;
    error('Data missing for %s:%s(%s).', ...
        clip_structs(cnt).test{1}, clip_structs(cnt).scene{1}, clip_struct.hrc{1});
end

