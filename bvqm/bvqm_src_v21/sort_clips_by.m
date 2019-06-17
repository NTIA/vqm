function [offsets] = sort_clips_by (separate, clip_struct, test_struct)
% SORT_CLIPS_BY
%  Sort clips alphabetically & optionally separate into groups of
%  indexes with the same test, scene, or hrc.  
% SYNTAX
%  [offsets] = sort_clips_by (separate, clip_struct, test_struct)
% DESCRIPTION
%  This function sorts the list of clips in clip_struct (of the same format
%  as GClips) with the associated test information in test_struct (of the
%  same format as GTests).  Clips will first be sorted alphabetically by
%  tests, then scenes, then HRCs.  Next, clips will be separated by the 
%  criteria given in the variable 'separate'.
%  
%  Unless otherwise specified, the return value 'offsets' will be a cell
%  array, where each cell contains an array listing all clips with an
%  identical value for the criteria specified in 'separate' (e.g., all
%  clips from one test.)
%
%  Variable 'separate' can have the following values: 
%   'none'      Do not separate clips (i.e., just sort alphabetically).
%               In this case, the return value will be an array of offsets.
%   'test'      Separate clips by test. 
%   'scene'     Separate clips by scene, but treat scenes with the same
%               name as different scenes if they are from different tests.
%               The original clip (if it exists) will be first in the cell
%               array's list (note that this violates alphabetical ordering).
%   'hrc'       Separate clips by HRC, but treat HRCs with the same name as
%               different HRCs if they are from different tests.

% sort test_struct alphabetically.
if length(test_struct) > 1,
    [temp,order1] = sort([test_struct.name]);
else
    temp = test_struct(1).name;
    order1 = 1;
end
test_struct = test_struct(order1);
for i=1:length(test_struct),
	[temp,order1] = sort([test_struct(i).scenes]);
	test_struct(i).scenes = test_struct(i).scenes(order1);

    [temp,order1] = sort([test_struct(i).hrc]);
	test_struct(i).hrc = test_struct(i).hrc(order1);
end

% sort sort clipset alphabetically.
[temp, order1] = sort([clip_struct.hrc]);
new_clip_struct = clip_struct(order1);

[temp, order2] = sort([new_clip_struct.scene]);
new_clip_struct = new_clip_struct(order2);
order = order1(order2);

[temp, order3] = sort([new_clip_struct.test]);
new_clip_struct = new_clip_struct(order3);
order = order(order3);

clear order1 order2 order3 temp;

% 'order' tells how to sort the original clip_struct into alphabetcial order.
clip_struct = new_clip_struct;

% Just sort alphabetically.  Don't do anything else.
if strcmp(lower(separate),'none'),
    offsets = order;
    return;
elseif strcmp(lower(separate),'test'),
    hold_test = [clip_struct.test];
    curr = 1;
    for loop = 1:length(test_struct),
        test_name = test_struct(loop).name{1};
        temp = find(strcmp(hold_test,test_name));
        if length(temp),
            offsets{curr} = temp;
            curr = curr + 1;
        end
    end
elseif strcmp(lower(separate),'scene'),
    hold_test = [clip_struct.test];
    hold_scene = [clip_struct.scene];
    curr = 1;
    % loop through all possible tests
    for loop = 1:length(test_struct),
        test_name = test_struct(loop).name{1};
        test_scenes = test_struct(loop).scenes;
        % loop through all scenes for this test
        for cnt = 1:length(test_scenes);
            temp = find(strcmp(hold_test,test_name) & ...
                strcmp(hold_scene,test_scenes{cnt}));
            if length(temp),
                % if there is an original, put it first.
                for index=1:length(temp),
                    if strcmp(clip_struct(temp(index)).hrc{1},'original'),
                        temp = [ temp(index) temp(1:index-1) temp(index+1:length(temp))];
                    end
                end
                % record
                offsets{curr} = temp;
                curr = curr + 1;
            end
        end
    end
elseif strcmp(lower(separate),'hrc'),
    hold_test = [clip_struct.test];
    hold_hrc = [clip_struct.hrc];
    curr = 1;
    % loop through each test
    for loop = 1:length(test_struct),
        test_name = test_struct(loop).name{1};
        test_hrcs = test_struct(loop).hrc;
        % loop through all HRCs within that test
        for cnt = 1:length(test_hrcs);
            temp = find(strcmp(hold_test,test_name) & ...
                strcmp(hold_hrc,test_hrcs{cnt}));
            % if at least 1 clip corresponds to that HRC, record it.
            if length(temp),
                offsets{curr} = temp;
                curr = curr + 1;
            end
        end
    end
else
    error('Separation type not recognized\n');
end

% now, roll back the offsets to the numbering prior to the alphabetical
% sorting.
for loop = 1:length(offsets),
    hold = offsets{loop};
    for cnt = 1:length(hold);
        hold(cnt) = order(hold(cnt));
    end
    offsets{loop} = hold;
end

