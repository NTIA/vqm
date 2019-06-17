function [jclips_mc] = bvqm_pc_calimport(jclips,jtests,c_rpt, need_temporal)
     
jclips = jclips( sort_clips_by('none',jclips,jtests) );
jclips_mc = jclips;
status_report = 0;

% load HRC specific data everywhere
[hrc_offsets] = sort_clips_by('hrc',jclips,jtests);

c_rpt_sheet3 = sprintf('%s_sheet3.csv', c_rpt);
c_rpt_sheet2 = sprintf('%s_sheet2.csv', c_rpt);
c_rpt_sheet1 = sprintf('%s_sheet1.csv', c_rpt);

fid = fopen(c_rpt_sheet1);
tmp = textscan(fid,'%s %n %n %n %n %n %n','HeaderLines',2, 'delimiter',',');
fclose(fid);
hrc_data = [tmp{2}, tmp{3}, tmp{4}, tmp{5}, tmp{6}, tmp{7}];

[row,col]=size(hrc_data);
if row ~= length(hrc_offsets) || col ~= 6,
    temp = nan(length(hrc_offsets),6);
    temp(1:row,1:col) = hrc_data;
    hrc_data = temp;
end

loading_cnt = 1;
for cnt = 1:size(hrc_offsets,2),
    if strcmpi('original',jclips(hrc_offsets{cnt}(1)).hrc{1}),
        continue;
    end
    for jc = hrc_offsets{cnt},
        jclips_mc(jc).scale.horizontal = hrc_data(loading_cnt,1);
        jclips_mc(jc).scale.vertical = hrc_data(loading_cnt,2);
        jclips_mc(jc).spatial.horizontal = hrc_data(loading_cnt,3);
        jclips_mc(jc).spatial.vertical = hrc_data(loading_cnt,4);
        jclips_mc(jc).luminance_gain = hrc_data(loading_cnt,5);
        jclips_mc(jc).luminance_offset = hrc_data(loading_cnt,6);
    end
    loading_cnt = loading_cnt + 1;
end

        
% load Scene specific data.
[scene_offsets] = sort_clips_by('scene',jclips,jtests);

fid = fopen(c_rpt_sheet2);
tmp = textscan(fid,'%s %n %n %n %n','HeaderLines',2, 'delimiter',',');
fclose(fid);
scene_data = [tmp{2}, tmp{3}, tmp{4}, tmp{5}];

[row,col]=size(scene_data);
if row ~= length(scene_offsets) || col ~= 4,
    temp = nan(length(scene_offsets),4);
    temp(1:row,1:col) = scene_data;
    scene_data = temp;
end

loading_cnt = 1;
for cnt = 1:size(scene_offsets,2),
    for jc = scene_offsets{cnt},
        jclips_mc(jc).cvr.top = scene_data(cnt,1);
        jclips_mc(jc).cvr.left = scene_data(cnt,2);
        jclips_mc(jc).cvr.bottom = scene_data(cnt,3);
        jclips_mc(jc).cvr.right = scene_data(cnt,4);
    end
end



% load Clip specific data.
[clip_offsets] = sort_clips_by('none',jclips,jtests);

fid = fopen(c_rpt_sheet3);
tmp = textscan(fid,'%s %s %n %n %n %n %n','HeaderLines',2, 'delimiter',',');
fclose(fid);
clip_data = [ tmp{3}, tmp{4}, tmp{5}, tmp{6}, tmp{7}];

[row,col]=size(clip_data);
if row ~= length(clip_offsets) || col ~= 5,
    temp = nan(length(clip_offsets),5);
    temp(1:row,1:col) = clip_data;
    clip_data = temp;
end

% check temporal registration against 'nan'.  Want 'nan' if don't need temporal
% registration; can't have 'nan' otherwise.
temporal_status_report = 0;
for cnt = 1:length(clip_offsets),
    jclips_mc(cnt).loc_start = clip_data(cnt,1);
    jclips_mc(cnt).loc_stop = clip_data(cnt,2);
    jclips_mc(cnt).mos = clip_data(cnt,5);
    
    if need_temporal,
        jclips_mc(cnt).align_start = clip_data(cnt,3);
        jclips_mc(cnt).align_stop = clip_data(cnt,4);
        if isnan(jclips_mc(cnt).align_start) || isnan(jclips_mc(cnt).align_stop),
            temporal_status_report = 1;
        end
    else
        jclips_mc(cnt).align_start = nan;
        jclips_mc(cnt).align_stop = nan;
    end
end

% check whether the user's data is valid . . . or not!
[status_report, message] = check_clips(jclips_mc, jtests, 'quiet');

if status_report == 1,
    
    button = questdlg( ...
        sprintf('%s\nIf you want to try again, edit the manual calibration file, save it, and press "Try Again".  Press "No Calibration" to continue without manual calibration.', message), ...
        'Fatal Error in File', 'Try Again', 'No Calibration', 'Try Again');
    if strcmpi(button,'No Calibration'),
        jclips_mc = jclips;
        return;
    elseif strcmpi(button,'try again'),
        [jclips_mc] = bvqm_pc_calimport(jclips,jtests,c_rpt,need_temporal);
        return;
    end

elseif status_report == 2,
    button = questdlg( ...
        sprintf('%s\nThe above errors cannot be ignored.  If you want to try again, edit the manual calibration file, save it, and press "Try Again".  Press "No Calibration" to continue without manual calibration.', message), ...
        'Invalid Data', 'Try Again', 'No Calibration', 'Try Again');
    if strcmpi(button,'No Calibration'),
        jclips_mc = jclips;
        return;
    elseif strcmpi(button,'try again'),
        [jclips_mc] = bvqm_pc_calimport(jclips,jtests,c_rpt,need_temporal);
        return;
    end

elseif status_report == 3,
    button = questdlg( ...
        sprintf('%s\nThe above warnings should be addressed.  Ignoring these warnings might cause VQM inaccuracies or errors.  If you want to try again, edit the manual calibration file, save it, and press "Try Again".  Press "continue" to use current values.', message), ...
        'Suspicous Data', 'Try Again', 'Continue', 'Try Again');
    if strcmpi(button,'continue'),
        ;
    elseif strcmpi(button,'try again'),
        [jclips_mc] = bvqm_pc_calimport(jclips,jtests,c_rpt,need_temporal);
        return;
    end
end

% now, deal with missing temporal registration, if appropriate.  Only get
% here with status_report==3 & decide to continue; or status_report==0.
if temporal_status_report == 4,
    button = questdlg( ...
        'Temporal alignment missing for some clips.  This problem cannot be ignored.  If you want to try again, edit the manual calibration file, save it, and press "Try Again".  Press "No Calibration" to continue without manual calibration.', ...
        'Missing Data', 'Try Again', 'No Calibration', 'Try Again');
    if strcmpi(button,'No Calibration'),
        jclips_mc = jclips;
        return;
    elseif strcmpi(button,'try again'),
        [jclips_mc] = bvqm_pc_calimport(jclips,jtests,c_rpt,need_temporal);
        return;
    end
end

function data = lame_comma_read(filename, top, left, bottom, right)

% Open the file for writing
fid = fopen(filename, 'w');

% write out the data
for (i = 1:size(data,1))
    for (j = 1:size(data,2))

        % write out one data element and place a comma
        a=data{i,j};
        if ischar(a)
            fprintf(fid, [a, ',']);
        else
            fprintf(fid, [mat2str(a), ',']);
        end
    end

    % erase the last comma and go to the next line
    fseek(fid, -1, 'eof');
    fprintf(fid, '\n');
end;

% erase the last comma and go to the next line
fseek(fid, -1, 'eof');
fprintf(fid, '\n');

fclose(fid);


