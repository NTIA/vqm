function bvqm_pc_shortrpt(model_selected,jclips,jtests,model_op,cal_type,working_dir,model_to_run,s_rpt,ssf, uncert, status, lgo_status, jclips_vr_unfilt, jclips_sr_unfilt, jclips_lgo_unfilt)

% Generate the short report, to be displayed on the screen and written to a
% file.

if ispc ; path_sep = '\'; else; path_sep='/'; end

% clk=fix(clock);
% yr=sprintf('%4.0f',clk(1)) ; mo=sprintf('%02.0f',clk(2));
% da=sprintf('%02.0f',clk(3));
% hr=sprintf('%02.0f',clk(4)); mn=sprintf('%02.0f',clk(5)); sc=sprintf('%02.0f',clk(6));
%
% date=strcat(yr,'-',mo,'-',da);
% time=strcat(hr,':',mn);

%report_file=strcat(working_dir, path_sep, 'bvqm_summary_report-', model_to_run, '.txt');
%report_file=strcat (working_dir,path_sep,'bvqm_summary_report-',model_to_run,'-',date,'@',time, '.txt');
report_file = s_rpt;
fid = fopen(report_file, 'wt');
num_clips = size(model_op.clip_name, 2);
sz_jclips = size(jclips,2); % size of jclips is different that the size of model_op.clip_name!


fprintf(fid, 'BVQM_v2.0    S U M M A R Y   R E P O R T   P A G E    %s\n\n\n', date);
if isfield(status,'algorithm'),
    fprintf(fid, 'Calibration: %s -- %s\n', cal_type, status.algorithm);
else
    fprintf(fid, 'Calibration: %s \n', cal_type);
end
fprintf(fid, 'Model: %s \n', cell2mat(model_selected));
fprintf(fid, 'Video Standard: %s \n', jclips(1).video_standard);
fprintf(fid, 'Clip Path: %s \n', cell2mat(jtests.path));
fprintf(fid, 'Results Path: %s \n', strcat(working_dir,path_sep));  %cell2mat(jtests.path));
fprintf(fid, '\n\n\n');

fprintf(fid, '--- VQM ---\n\n');

if exist('status') && isfield(status,'model') && status.model == 1,
    fprintf(fid, 'Model calculation failed\n');
    fprintf(fid, 'Error:  %s\n', lasterr);
    
else

    line1 = ['Clip Name                       ' ' VQM' '\n'];
    line2 = ['--------------                  ' '------' '\n'];
    fprintf(fid, line1);
    fprintf(fid, line2);


    for count = 1:num_clips
        clip_name = cell2mat(model_op.clip_name(count));
        clip_data = model_op.data(1,count);
        %     fprintf(fid, '%s%1.4f\n', clip_name, clip_data);
        fprintf(fid, '%-30s  %1.4f\n', clip_name, clip_data);
    end


    line2 = ['--------------------------' '\n'];
    fprintf(fid, '\n\n');
    fprintf(fid, 'HRC Name                      Average VQM \n');
    fprintf(fid, '--------                      ----------- \n');

    hrc_stats=ave_par_values(model_op,'scene');
    num_hrcs=length(hrc_stats.data(1,:)); % These are avg vqm values
    hrc_stats.clip_name; % These are the hrc names
    for hrc_index=1:num_hrcs
        prt_hrc=char(hrc_stats.clip_name(hrc_index));
        [pa, pb]=strtok(prt_hrc,'_');
        [pc, pd]=strtok(pb, '_');
        [pe, pf] = strtok(pd, '_');
    %     fprintf(fid, '%-30s    %7.4f\n', char(hrc_stats.clip_name(hrc_index)), hrc_stats.data(1,hrc_index));
        fprintf(fid, '%-30s    %7.4f\n', pe, hrc_stats.data(1,hrc_index));

    %     fprintf(fid, '%-30s  %1.4f\n', char(hrc_stats.clip_name(hrc_index)), hrc_stats.data(1,hrc_index));
    end
end
fprintf(fid, '\n\n\n--- OPTIONS AND CALIBRATION ---\n\n');

fprintf(fid, '.......... Parsing ..........\n');

if status.is_parsed == 1,
    fprintf(fid,'Files longer than 15 seconds truncated to 15 seconds\n\n');
    
elseif status.is_parsed == 2,
    fprintf(fid,'All files parsed into %d second segments, shifting by %d seconds\n', ...
        status.parse_length, status.parse_overlap);
    
else
    fprintf(fid, 'No Parsing Performed.\n');
end


% BVQM_PC_CAL SUMMARY:-
ct_T='Temporal Registration and Valid Region Only';
ct_TM='Temporal Registration after Manual Calibration';
ct_NC='No Calibration';
ct_MC='Manual Calibration';
ct_RR2='Reduced Reference Calibration Version 2';

if isequal(cal_type, ct_NC),
    fprintf(fid,'\n.......... Calibration ..........\n\n');

    fprintf(fid, 'Calibration metrics neither estimated nor removed.\n');
    fprintf(fid, 'Assume processed video matchs the original video perfectly with respect to:\n');
    fprintf(fid, '- spatial scaling\n');
    fprintf(fid, '- spatial shift\n');
    fprintf(fid, '- temporal shift\n');
    fprintf(fid, '- luminance gain offset\n');
    fprintf(fid, 'Discard overscan (if any), and assume picture inside the standard overscan \n');
    fprintf(fid, 'contains only valid video.\n');
    
    if strcmpi(cell2mat(model_selected),'PSNR with Variable Frame Delay'),
        fprintf(fid,'\n.......... Calibration Overriden by Model ..........\n\n');

        line2 = ['Clip Name                          Gain    Offset\n'];
        line3 = ['--------------                    ------   ------\n'];

        fprintf(fid, line2); fprintf(fid, line3);
        for count = 1:num_clips
            clip_name = cell2mat(model_op.clip_name(count));
            for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
                jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
                if strcmp(jclips_name, clip_name)
                    y_gain = jclips(jc).luminance_gain;
                    y_offset = jclips(jc).luminance_offset;
                    fprintf (fid, '%-30s    %7.4f  %7.3f\n' ,clip_name, y_gain, y_offset);
                end
            end
        end
    end

    fclose(fid);
    return;
end

if isequal(cal_type, ct_MC),
    fprintf(fid,'\n.......... Calibration ..........\n\n');

    fprintf(fid, 'Manually entered calibration values used.  See detailed report.\n');

    fclose(fid);
    return;
end


% Spatial scaling option:
ssf_txt='\n.......... Spatial Scaling ..........\n';
fprintf(fid,ssf_txt);
% fprintf(fid, 'Spatial Scaling Option:\n');

if isequal(cal_type, ct_TM),
    fprintf(fid, 'Values manually entered for spatial scaling.\n');
    fprintf(fid, 'See detailed report for values\n');
elseif ssf==0
    fprintf(fid, 'Assume processed video has not been spatially scaled (e.g., image not stretched).\n\n');

else  % create table:

    if exist('jclips_sr_unfilt') ,
        
        if isfield(status,'spatial') && isfield(status.spatial,'error') && status.spatial.error > 1,
            fprintf(fid, 'ERROR\nFatal error encountered during Spatial Scaling.  \n');
            fprintf(fid, 'Assume the process video sequences do not contain spatial scaling.\n');
        else
            if isfield(status.spatial,'scale') && status.spatial.scale > 0,
                fprintf(fid, '%d clips appear to have scaling in excess of what this algorithm can predict\n\n', status.spatial.scale);
            end

            line1 = ['HRC Name                        Horizontal Vertical\n'];
            line2 = ['--------------                  ---------- --------\n'];
            fprintf(fid, line1);fprintf(fid, line2);

            % now, assign values to table
            [hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
            for count = 1:length(hrc_offsets),
                jc = (hrc_offsets{count}(1));
                jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
                [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
                if search_status == 0,
                    hscale = jclips(jc).scale.horizontal;
                    vscale = jclips(jc).scale.vertical;
                    fprintf (fid, '%-30s     %4d      %4d\n' ,jclips(jc).hrc{1}, hscale, vscale);
                end

            end
        end
    else
        line1 = ['Clip Name                       Horizontal Vertical\n'];
        line2 = ['--------------                  ---------- --------\n'];
        fprintf(fid, line1);fprintf(fid, line2);
        for count = 1:num_clips
            clip_name = cell2mat(model_op.clip_name(count));
            for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
                jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
                if strcmp(jclips_name, clip_name)
                    hscale = jclips(jc).scale.horizontal;
                    vscale = jclips(jc).scale.vertical;

                    fprintf (fid, '%-30s     %4d      %4d\n' ,clip_name, hscale, vscale);
                end
            end
        end
    end


end



fprintf (fid, '\n\n.......... Spatial Shift and Luminance Gain & Offset ..........\n');

if isfield(status, 'gain_offset') && status.gain_offset,
    % error, 
    fprintf(fid, '\nERROR\nFatal error computing Gain & Offet.\n');
    fprintf(fid, 'Assume no luminance gain & offset unchanged (i.e., gain = 1.0 & offset = 0.0)\n\n');
    
end
if isfield(status,'spatial'),
    if isfield(status.spatial,'error') && status.spatial.error >= 1,
        fprintf(fid, '\nERROR\nFatal error encountered during Spatial Shift. \n');
        fprintf(fid, 'Assume processed video sequence is not spatially shifted  with respect to original.\n\n');
    else
        if isfield(status.spatial,'shift') && status.spatial.shift > 0,
            fprintf(fid, '%d clips appear to have spatial shift in excess of what this algorithm can predict\n\n', status.spatial.scale);
        end
    end
    if isfield(status.spatial,'failed') && status.spatial.failed >= 1,
        fprintf(fid, '\nWARNING\nSpatial Shift failed for %d clips.\n', status.spatial.failed);
        fprintf(fid, 'Assume processed video sequence is not spatially shifted  with respect to original.\n\n');
    end
end
if strcmpi(cell2mat(model_selected),'PSNR with Variable Frame Delay'),
    fprintf(fid, '\nWARNING\nLuminance gain & offset overriden by model. Model values reported.\n\n');
end

if isequal(cal_type, ct_T),
    fprintf(fid, 'Assume the processed video has not been shifted spatially.\n');
    fprintf(fid, 'Assume luminance gain & offset of processed video matches original.\n\n');
elseif isequal(cal_type, ct_TM),
    fprintf(fid, 'Values manually entered for spatial shift, luminance gain and luminance offset.\n');
    fprintf(fid, 'See detailed report for values\n');
else
    if exist('jclips_sr_unfilt') && exist('jclips_lgo_unfilt'),
        line2 = ['HRC Name                           H   V       Gain    Offset\n'];
        line3 = ['--------------                    --- ---     ------   ------\n'];

        fprintf(fid, line2); fprintf(fid, line3);

        % now, assign values to table
        [hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
        for count = 1:length(hrc_offsets),
            jc = (hrc_offsets{count}(1));
            jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
            [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
            if search_status == 0,
                spatial_h = jclips(jc).spatial.horizontal;
                spatial_v = jclips(jc).spatial.vertical;
                y_gain = jclips(jc).luminance_gain;
                y_offset = jclips(jc).luminance_offset;
                fprintf (fid, '%-30s    %3d %3d    %7.4f  %7.3f\n', jclips(jc).hrc{1}, spatial_h, spatial_v, y_gain, y_offset);
            end

        end
    else
    % no filtering
        for count = 1:num_clips
            clip_name = cell2mat(model_op.clip_name(count));
            for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
                jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
                if strcmp(jclips_name, clip_name)
                    spatial_h = jclips(jc).spatial.horizontal;
                    spatial_v = jclips(jc).spatial.vertical;
                    y_gain = jclips(jc).luminance_gain;
                    y_offset = jclips(jc).luminance_offset;
                    fprintf (fid, '%-30s    %3d %3d    %7.4f  %7.3f\n' ,clip_name, spatial_h, spatial_v, y_gain, y_offset);
                end
            end
        end
    end
end


fprintf (fid, '\n\n.......... Color Gain & Offset ..........\n');

if isequal(cal_type, ct_RR2 ),
    if isfield(status, 'gain_offset') && status.gain_offset,
        % error, 
        fprintf(fid, '\nERROR\nFatal error computing Color Gain & Offet.\n\n');
    else
    fprintf(fid, 'Color Gain & Offset were estimated but not removed.\n\n');

        line2 = ['HRC Name                          Cb Gain  Cb Offset  Cr Gain  Cr Offset\n'];
        line3 = ['--------------                    -------  ---------  -------  ---------\n'];

        fprintf(fid, line2); fprintf(fid, line3);

        % now, assign values to table
        [hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
        for count = 1:length(hrc_offsets),
            jc = (hrc_offsets{count}(1));
            jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
            [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
            if search_status == 0,
                cb_gain = jclips(jc).cb_gain;
                cb_offset = jclips(jc).cb_offset;
                cr_gain = jclips(jc).cr_gain;
                cr_offset = jclips(jc).cr_offset;
                fprintf (fid, '%-30s    %7.4f %7.3f    %7.4f  %7.3f\n', jclips(jc).hrc{1}, cb_gain, cb_offset, cr_gain, cr_offset);
            end

        end
    end
else
    fprintf(fid, 'Color Gain & Offset were not estimated.');
end


% COMMON VALID REGION SUMMARY:
fprintf (fid, '\n\n.......... Valid Region ..........\n');
if isfield(status, 'valid') && status.valid,
    % error, 
    fprintf(fid,'\nERROR\nFatal error computing valid region.  Discard overscan region, if any.\n\n');
elseif isequal(cal_type, ct_TM),
    fprintf(fid, 'Values manually entered for valid region.  See detailed report for values\n\n');
elseif exist('jclips_vr_unfilt') ,

    line2 = ['Scene Name                      top  bottom    left right\n'];
    line3 = ['--------------                  ---- ------    ---- -----\n'];
    fprintf(fid, line2);fprintf(fid, line3);

    % now, assign values to table
    [scene_offsets] = sort_clips_by('scene',jclips,jtests);
    for count = 1:length(scene_offsets),
        jc = (scene_offsets{count}(2));
        jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
        [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
        if search_status == 0,
            top  = jclips(jc).cvr.top;
            bot  = jclips(jc).cvr.bottom;
            left  = jclips(jc).cvr.left;
            right = jclips(jc).cvr.right;

            fprintf (fid, '%-30s  %4d  %4d     %4d  %4d\n',jclips(jc).scene{1}, top, bot, left, right);
        end

    end
end


% TEMPORAL REG SUMMARY:
fprintf (fid, '\n\n.......... Temporal Registration ..........\n');

if isfield(status,'temporal') && isfield(status.temporal, 'error') && status.temporal.error == 1,
    fprintf(fid, '\nERROR\nFatal error encountered during temporal registration.\n');
    fprintf(fid, 'Assume that the first valid frame of the processed video file aligns in time with\n');
    fprintf(fid, 'the first valid frame of the matching original video file.\n');
elseif isfield(status,'temporal'),
    fprintf(fid,'\n');

    if isfield(status.temporal,'uncertainty') && status.temporal.uncertainty > 0,
        if isfield(status, 'temporal2'),
            fprintf(fid,'\nFirst run of temporal registration indicates:\n');
        end
        fprintf(fid, '%d clips should be re-run with a larger temporal uncertainty\n', status.temporal.uncertainty);
        fprintf(fid, 'Temporal registration uncertainty was +/- %f seconds\n', uncert);
    end
    fprintf(fid,'\n');
    if isfield(status, 'temporal2'),
        if isfield(status.temporal2, 'error') && status.temporal2.error == 1,
            fprintf(fid, '\nOn second run of temporal registration,\nFatal error encountered during temporal registration.\n');
            fprintf(fid, 'Assume that the first valid frame of the processed video file aligns in time with\n');
            fprintf(fid, 'the first valid frame of the matching original video file.\n');
        else
        
            if isfield(status.temporal2,'uncertainty') && status.temporal2.uncertainty > 0,
                fprintf(fid,'\nSecond run of temporal registration indicates:\n');
                fprintf(fid, '%d clips should be re-run with a larger temporal uncertainty\n', status.temporal2.uncertainty);
                fprintf(fid, 'Temporal registration uncertainty was +/- %f seconds\n', uncert);
            end
        end
    
    end
end

line1 = ['.                                      Original       Processed\n'];
line2 = ['Clip Name                       delay  start  stop    start stop\n'];
line3 = ['--------------                  -----  ----- -----    ----- -----\n'];

fprintf(fid, line1);fprintf(fid, line2); fprintf(fid, line3);

% now, assign values to table
for count = 1:num_clips
    clip_name = cell2mat(model_op.clip_name(count));
    for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
        jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
        orig_name=strcat(jclips(jc).test, '_', jclips(jc).scene, '_', 'original');
        % find original clip for current hrc, and extract the start/stop
        % frame values from the original clip:
        for count2=1:sz_jclips
            jclips_name2=strcat(jclips(count2).test, '_', jclips(count2).scene, '_', jclips(count2).hrc);
        end

        % get align start/stop values for processed clips:
        if strcmp(jclips_name, clip_name)
            align_start = jclips(jc).align_start;
            align_stop = jclips(jc).align_stop;
            orig_start = jclips( find_original(jclips, jc) ).align_start;
            orig_stop = jclips( find_original(jclips, jc) ).align_stop;
            delay = align_start - orig_start; 
            if ~strcmp(jclips(jc).video_standard,'progressive') & mod(jclips(jc).spatial.vertical, 2) == 1,
                delay = delay + 0.5;
            end
            fprintf (fid, '%-30s  %5.1f  %3d    %3d     %3d    %3d\n', clip_name, delay, orig_start, orig_stop, align_start, align_stop);
        end
    end
end



fclose(fid);
