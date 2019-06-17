function error_status = bvqm_pc_calexport(jclips,jtests,c_rpt)
        
jclips = jclips( sort_clips_by('none',jclips,jtests) );

% print out HRC specific data
hrc_data{1,1} = ' ';
hrc_data{1,2} = 'Scale';
hrc_data{1,3} = 'Scale';
hrc_data{1,4} = 'Shift';
hrc_data{1,5} = 'Shift';
hrc_data{1,6} = 'Luminance';
hrc_data{1,7} = 'Luminance';

hrc_data{2,1} = 'HRC';
hrc_data{2,2} = 'Horizontal';
hrc_data{2,3} = 'Vertical';
hrc_data{2,4} = 'Horizontal';
hrc_data{2,5} = 'Vertical';
hrc_data{2,6} = 'Gain';
hrc_data{2,7} = 'Offset';

write_cnt = 3;

[hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
for cnt = 1:length(hrc_offsets),
    jc = (hrc_offsets{cnt}(1));
    if strcmpi('original',jclips(jc).hrc{1}),
        continue;
    end
    hrc_data{write_cnt,1} = jclips(jc).hrc{1};
    hrc_data{write_cnt,2} = jclips(jc).scale.horizontal;
    hrc_data{write_cnt,3} = jclips(jc).scale.vertical;
    hrc_data{write_cnt,4} = jclips(jc).spatial.horizontal;
    hrc_data{write_cnt,5} = jclips(jc).spatial.vertical;
    hrc_data{write_cnt,6} = jclips(jc).luminance_gain;
    hrc_data{write_cnt,7} = jclips(jc).luminance_offset;
    write_cnt = write_cnt + 1;
end

        
% print out Scene specific data.
scene_data{1,1} = 'Scene';
scene_data{1,2} = 'Valid';
scene_data{1,3} = 'Valid';
scene_data{1,4} = 'Valid';
scene_data{1,5} = 'Valid';

scene_data{2,1} = 'Scene';
scene_data{2,2} = 'Top';
scene_data{2,3} = 'Left';
scene_data{2,4} = 'Bottom';
scene_data{2,5} = 'Right';

write_cnt = 3;

[scene_offsets] = sort_clips_by('scene',jclips,jtests);
for cnt = 1:length(scene_offsets),
    jc = (scene_offsets{cnt}(1));
    scene_data{write_cnt,1} = jclips(jc).scene{1};
    scene_data{write_cnt,2} = jclips(jc).cvr.top;
    scene_data{write_cnt,3} = jclips(jc).cvr.left;
    scene_data{write_cnt,4} = jclips(jc).cvr.bottom;
    scene_data{write_cnt,5} = jclips(jc).cvr.right;
    write_cnt = write_cnt + 1;
end


% print out clip specific data 
clip_data{1,1} = '';
clip_data{1,2} = '';
clip_data{1,3} = 'File';
clip_data{1,4} = 'File';
clip_data{1,5} = 'Temporal';
clip_data{1,6} = 'Temporal';
clip_data{1,7} = '';


clip_data{2,1} = 'Scene';
clip_data{2,2} = 'HRC';
clip_data{2,3} = 'First Valid Frame';
clip_data{2,4} = 'Last Valid Frame';
clip_data{2,5} = 'Start Frame';
clip_data{2,6} = 'Stop Frame';
clip_data{2,7} = 'MOS';

write_cnt = 3;

[clip_offsets] = sort_clips_by('none',jclips,jtests);
for cnt = 1:length(clip_offsets),
    jc = clip_offsets(cnt);
    clip_data{write_cnt,1} = jclips(jc).scene{1};
    clip_data{write_cnt,2} = jclips(jc).hrc{1};
    clip_data{write_cnt,3} = jclips(jc).loc_start;
    clip_data{write_cnt,4} = jclips(jc).loc_stop;
    clip_data{write_cnt,5} = jclips(jc).align_start;
    clip_data{write_cnt,6} = jclips(jc).align_stop;
    clip_data{write_cnt,7} = jclips(jc).mos;
    write_cnt = write_cnt + 1;
end

% export to CSV files
c_rpt_sheet3 = sprintf('%s_sheet3.csv', c_rpt);
c_rpt_sheet2 = sprintf('%s_sheet2.csv', c_rpt);
c_rpt_sheet1 = sprintf('%s_sheet1.csv', c_rpt);
if exist(c_rpt_sheet3,'file') || exist(c_rpt_sheet2,'file') || exist(c_rpt_sheet1,'file'),
    uiwait(msgbox(sprintf('Files "%s", "%s", and "%s" will be overwritten.  Press "OK" to continue.',...
        c_rpt_sheet3, c_rpt_sheet2, c_rpt_sheet1), ...
        'Warning! Deleting File'));  
    pause(0.2);
end
try
    warning off;
    delete(c_rpt_sheet3);
    lame_comma_write(c_rpt_sheet3, clip_data);
    delete(c_rpt_sheet2);
    lame_comma_write(c_rpt_sheet2, scene_data);
    delete(c_rpt_sheet1);
    lame_comma_write(c_rpt_sheet1, hrc_data);
    error_status = 1;
catch
    error_status = 0;
end
warning on;


function lame_comma_write(filename, data)

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

%disp('Done writing CSV');

