function [number] = search_test_list(test_structs, clip_struct)
% SEARCH_TEST_LIST
%  Find the test structure associated with one clip.
% SYNTAX
%  [number] = search_test_list(test_structs, clip_struct)
% DESCRIPTION
%  Given the array offset within 'test_structs' of the test to which 
%  clip 'clip_struct' belongs.  Variable 'clip_struct' is a structure with
%  the same format at SClips, and 'test_structs' is an array with the same
%  structure as GTests.

% loop through and find this clip.
[rows] = max(size(test_structs));
number = -1;
for cnt = 1:rows,
    if strcmp(test_structs(cnt).name{1},clip_struct.test),
        number = cnt;
        break;
    end
end

% check for missing entry.
if number < 0,
    error('Requested test, ''%s'', does not exist', clip_struct.test{1});
end

% check for multiple entries.
for cnt = number+1:rows,
    if strcmp(test_structs(cnt).name{1},clip_struct.test),
        error('Requested test, ''%s'', occurs twice:  #%d and #%d', clip_struct.test{1}, number, cnt);
    end
end
