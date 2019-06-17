function orig_offset = find_original(clip_structs, proc_offset)
% FIND_ORIGINAL
%  Given a processed clip, find the original clip.
% SYNTAX
%  orig_offset = find_original(clip_structs, proc_offset);
% DESCRIPTION
%  'clip_structs' is an array of clips, all of the same format as GClips.
%  'proc_offset' is the address of a clip in 'clip_structs', a processed
%  clip.  Find and return the array address within 'clip_structs' of the
%  original clip associated with clip_structs(proc_offset).

if proc_offset > length(clip_structs),
    error('Requested clip number past end of array');
end
if strcmp(clip_structs(proc_offset).hrc{1},'original'),
    error('This clip is an original; functionality undefined');
end

% find the original clip's offset.
test = [clip_structs.test];
scene = [clip_structs.scene];
hrc = [clip_structs.hrc];
orig_offset = find(strcmp(test,test(proc_offset)) & strcmp(scene,scene(proc_offset)) ...
    & strcmp(hrc,'original'));