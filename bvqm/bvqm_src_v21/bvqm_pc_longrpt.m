function detailed_report(model_selected,jclips,jtests,model_op,cal_type,working_dir,model_to_run,d_rpt,ssf,uncert, status, lgo_status, jclips_vr_unfilt, jclips_sr_unfilt, jclips_lgo_unfilt)

%
% WRITTEN BY : Mark A. McFarland, P.E.    303/497-4132
%

if ispc ; path_sep = '\'; else; path_sep='/'; end

clk=fix(clock);
yr=sprintf('%4.0f',clk(1)) ; mo=sprintf('%02.0f',clk(2)); da=sprintf('%02.0f',clk(3));
hr=sprintf('%02.0f',clk(4)); mn=sprintf('%02.0f',clk(5)); sc=sprintf('%02.0f',clk(6));

date=strcat(yr,'-',mo,'-',da);
time=strcat(hr,':',mn);

%report_file=strcat(working_dir, path_sep, 'bvqm_summary_report-', model_to_run, '.txt');
%report_file=strcat
%(working_dir,path_sep,'bvqm_detail_report-',model_to_run,'-',date,'@',time, '.txt');
report_file=d_rpt;
fid = fopen(report_file, 'wt');
num_clips = size(model_op.clip_name, 2);
sz_jclips = size(jclips,2); % size of jclips is different that the size of model_op.clip_name!


fid = fopen(report_file, 'wt');
fprintf(fid, 'BVQM_v2.0    D E T A I L E D   R E P O R T   P A G E    %s \n\n\n', date);
if isfield(status,'algorithm'),
    fprintf(fid, 'Calibration: %s -- %s\n', cal_type, status.algorithm);
else
    fprintf(fid, 'Calibration: %s \n', cal_type);
end
fprintf(fid, 'Model: %s \n', cell2mat(model_selected));
fprintf(fid, 'Video Standard: %s \n', jclips(1).video_standard);
fprintf(fid, 'Frames Per Second: %1.2f \n', jclips(1).fps);
fprintf(fid, 'Image Size, Columns: %1.0f  Rows: %1.0f  \n', jclips(1).image_size.cols, jclips(1).image_size.rows);
fprintf(fid, 'Clip Path: %s \n', cell2mat(jtests.path));
fprintf(fid, 'Results Path: %s \n', strcat(working_dir, path_sep));
fprintf(fid, '\n\n\n');


fprintf(fid, '--- Clip File Mappings ---\n');
fprintf(fid, '%-25s %-10s %-10s\n' ,'Clip Name', 'Scene Name', 'HRC Name');
fprintf(fid, '%-25s %-10s %-10s\n' ,'---------', '----------', '--------');
for jc = 1:sz_jclips
     jc_n = cell2mat(strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc));
     jc_s=cell2mat(jclips(jc).scene) ; jc_h=cell2mat(jclips(jc).hrc);
     jc_fn=cell2mat(jclips(jc).file_name);

     fprintf (fid, '%-25s %-10s %-10s\n', jc_n , jc_s, jc_h);
end 
fprintf(fid,'\n\n');


fprintf(fid, '%-25s %-35s \n' ,'Clip Name', 'File Name');
fprintf(fid, '%-25s %-35s \n' ,'---------', '---------');
for jc = 1:sz_jclips
     jc_n = cell2mat(strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc));
     jc_s=cell2mat(jclips(jc).scene) ; jc_h=cell2mat(jclips(jc).hrc);
     jc_fn=cell2mat(jclips(jc).file_name);

     fprintf (fid, '%-25s %-34s \n', jc_n , jc_fn);
end

%fprintf(fid, '\n\n\n--- Clip File Associations ---\n\n');
% Print out scene table
fprintf(fid, '\n%-15s %-10s \n','Scene Name', '# Clips');
fprintf(fid, '%-15s %-10s \n','----------', '-------');

uqs=unique([jclips.scene]);
nmc=0; % number of matching clips
for jc=1:size(unique([jclips.scene]),2)
    unique_name=cell2mat(uqs(jc));
    nmc=0;
    for ccount=1:sz_jclips
        jc_s=cell2mat(jclips(ccount).scene);
        if isequal(unique_name, jc_s)
            nmc=nmc+1;
        end
    end
     fprintf(fid, '%-15s %d \n', unique_name, nmc);   
end

% Print out HRC table
fprintf(fid, '\n%-15s %-10s \n','HRC   Name', '# Clips');
fprintf(fid, '%-15s %-10s \n','----------', '-------');

uqs=unique([jclips.hrc]);
nmc=0; % number of matching clips
for jc=1:size(unique([jclips.hrc]),2)
    unique_name=cell2mat(uqs(jc));
    nmc=0;
    for ccount=1:sz_jclips
        jc_h=cell2mat(jclips(ccount).hrc);
        if isequal(unique_name, jc_h)
            nmc=nmc+1;
        end
    end
     fprintf(fid, '%-15s %d \n', unique_name, nmc);   
end


fprintf(fid, '\n\n\n--- VQM ---\n\n');

if exist('status') && isfield(status,'model') && status.model == 1,
    fprintf(fid, 'Model calculation failed\n');
    fprintf(fid, 'Error:  %s\n', lasterr);
    
else

    total_clips=length(jclips); %total number of clips including originals
    cliplist=strcat([jclips.test],'_',[jclips.scene],'_',[jclips.hrc]); %make a list of clips
    for clip_index=1:total_clips

        if clip_index ~= 1,
            fprintf(fid,'\n..........\n\n');
        end
        fprintf(fid, 'File: %s\n\n', char( jclips(clip_index).file_name));


        %let's print out the clip-specific informataion from jclips
        % do we lack the mean opinion scores or standard deviations?
        lack_mos = isnan( jclips(clip_index).mos );
        lack_stdev =isnan( jclips(clip_index).stdev );
        fprintf(fid, 'Scene: %s\n',char(jclips(clip_index).scene));
        fprintf(fid, 'HRC: %s\n',char(jclips(clip_index).hrc));

        if (lack_mos == 0)
            fprintf(fid, 'Mean Opinion Score: %1.4f\n',jclips(clip_index).mos);
        end
        if (lack_stdev == 0 )
            fprintf(fid, 'Standard Deviation: %1.4f\n',jclips(clip_index).stdev);
        end


        %Extract the stuff we need from model_op that is specific to this clip
        %  NOTE: error: pars_find_clip outputs to screen!!! o/p = "Requested data not found"
        clip_specific_struct=pars_find_clip(model_op,'clip',cliplist(clip_index),'par_struct');

        %now, if we are looking at an 'original' clip than we won't have any
        %vqm data and thus clip_specific_struct won't actually be a struct
        %if this is the case than we need not output any parameters
        if isstruct(clip_specific_struct) == 1
            fprintf(fid, '\nParameter                    Value\n');
            line2 =       ['--------------              -------\n'];
            fprintf(fid,line2);
            num_param=length(model_op.par_name);
            %pad all of the parameter names to be the same length by adding spaces
            parameter_names=char([model_op.par_name, '                            ']');
            for par_index=1:num_param

                fprintf(fid,'%s%7.4f\n',char( parameter_names(par_index,:) ),clip_specific_struct.data(par_index));
            end
        end

    end

    fprintf(fid, '\n\n.......... Average Results by HRC ..........\n\n');
    fprintf(fid, 'HRC Name                        Average VQM \n');
    fprintf(fid, '--------                        ----------- \n');

    hrc_stats=ave_par_values(model_op,'scene');
    num_hrcs=length(hrc_stats.data(1,:)); % avg vqm values
    hrc_stats.clip_name; % hrc names
    for hrc_index=1:num_hrcs
        prt_hrc=char(hrc_stats.clip_name(hrc_index));
        [pa, pb]=strtok(prt_hrc,'_');
        [pc, pd]=strtok(pb, '_');
        [pe, pf] = strtok(pd, '_');
    %     fprintf(fid, '%-30s    %7.4f\n', char(hrc_stats.clip_name(hrc_index)), hrc_stats.data(1,hrc_index));
        fprintf(fid, '%-30s    %7.4f\n', pe, hrc_stats.data(1,hrc_index));
    end
end

fprintf(fid, '\n\n\n--- OPTIONS AND CALIBRATION ---\n\n');

fprintf(fid, '.......... Parsing ..........\n');

num_clips = size(model_op.clip_name, 2);

if status.is_parsed == 1,
    fprintf(fid,'Files longer than 15 seconds truncated to 15 seconds\n\n');
    
elseif status.is_parsed == 2,
    fprintf(fid, 'All files parsed into %d second segments, shifting by %d seconds\n', ...
        status.parse_length, status.parse_overlap);
    
else
    fprintf(fid, 'No Parsing Performed.\n');
end


ct_T='Temporal Registration and Valid Region Only';
ct_TM='Temporal Registration after Manual Calibration';
ct_NC='No Calibration';
ct_MC='Manual Calibration';
ct_RR1='Reduced Reference Calibration Version 1';
ct_RR2='Reduced Reference Calibration Version 2';
ct_FR='Full Reference Calibration';
ct_PSNR='Peak Signal to Noise Ratio';
ct_RR2_PSNR='ITU-T J.244 followed by ITU-T J.340';
ct_FR_PSNR='ITU-T J.144 followed by ITU-T J.340';

if(isequal(cal_type,ct_PSNR))
    %Clear variables that are not used in order to prevent crashes below!
    clear jclips_sr_unfilt;
    clear jclips_lgo_unfilt;
end
if(isequal(cal_type,ct_RR2_PSNR) || isequal(cal_type,ct_FR_PSNR))
    %Clear variables that are not used in order to prevent crashes below!
    clear jclips_lgo_unfilt;
end

%%%%%% No Calibration comments %%%%%%%%%

% Note frames discarded from file
fprintf(fid,'\n.......... First & Last Valid Frame ..........\n');

if isequal(cal_type,ct_T) || isequal(cal_type,ct_TM) || isequal(cal_type,ct_RR1) || ...
        isequal(cal_type,ct_RR2) || isequal(cal_type,ct_FR),
    fprintf(fid, 'Algorithms begin with the assumtion that the first valid frame of the processed\n');
    fprintf(fid, 'video file aligns in time with the first valid frame of the matching original\n');
    fprintf(fid, 'video file.\n');
    fprintf(fid,'\nFrames before first valid frame, and after last valid frame,\n');
    fprintf(fid,'are ignored by all algorithms.\n');
end

fprintf(fid,'\nClip Name                       first  last\n');
fprintf(fid,'--------------                  ----- -----\n');
for jc = 1:length(jclips),
    jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
    fprintf (fid, '%-30s  %3d    %3d\n', jclips_name{1}, jclips(jc).loc_start, jclips(jc).loc_stop);
end
fprintf(fid,'\n\n');


% Spatial scaling option:
% fprintf(fid, 'Spatial Scaling Option:\n');

fprintf(fid,'\n.......... Spatial Scaling ..........\n');

if isfield(status,'spatial')&& isfield(status.spatial,'error') && status.spatial.error > 1,
    fprintf(fid, '\nERROR\nFatal error encountered during Spatial Scaling.\n');
    fprintf(fid, 'Assume processed video has not been spatially scaled (e.g., image not stretched).\n\n');

elseif (isequal(cal_type,ct_RR1) || isequal(cal_type,ct_RR2)) && exist('jclips_sr_unfilt') && ssf ~= 0,
    if isfield(status,'spatial') && isfield(status.spatial,'scale') && status.spatial.scale > 0,
        fprintf(fid, '%d clips appear to have scaling in excess of what this algorithm can predict.\n\n', status.spatial.scale);
    end

    fprintf(fid,'\n   - - -   Un-Filtered    - - - \n');
    line1 = ['Clip Name                       Horizontal Vertical\n'];
    line2 = ['--------------                  ---------- --------\n'];
    fprintf(fid, line1);fprintf(fid, line2);
    for count = 1:num_clips
        clip_name = cell2mat(model_op.clip_name(count));
        for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
            jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
            if strcmp(jclips_name, clip_name)
                hscale = jclips_sr_unfilt(jc).scale.horizontal;
                vscale = jclips_sr_unfilt(jc).scale.vertical;

                fprintf (fid, '%-30s     %4d      %4d\n' ,clip_name, hscale, vscale);
            end
        end
    end

    ssf_txt='\n   - - -   Filtered    - - - \n';
    fprintf(fid,ssf_txt);
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
    
elseif isequal(cal_type, ct_MC) || isequal(cal_type,ct_TM),

    % manual calibration

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
    
else
    fprintf(fid, 'Assume processed video has not been spatially scaled (e.g., image not stretched).\n\n');

end

% check for warnings.
[hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
num_warn = 0;
for count = 1:length(hrc_offsets),
    jc = (hrc_offsets{count}(1));
    jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
    [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
    if search_status == 0,
        if jclips(jc).scale.horizontal ~= 1000 || jclips(jc).scale.vertical ~= 1000,
            if num_warn == 0;
                num_warn = 1;
                fprintf(fid,'\n');
            end
            fprintf (fid, 'Warning: Video scaling detected in %s; visual confirmation recommended.\n', jclips(jc).hrc{1});
        end
    end
end
if num_warn == 1,
    fprintf(fid,'\n');
end


fprintf (fid, '\n\n.......... Spatial Shift ..........\n');



if isfield(status,'spatial') && isfield(status.spatial,'error') && status.spatial.error >= 1,
    fprintf(fid, '\nERROR\nFatal error encountered during Spatial Shift.\n');
    fprintf(fid, 'Assume the processed video has not been shifted spatially.\n\n');

elseif isequal(cal_type, ct_NC) || isequal(cal_type, ct_T),
    fprintf(fid, 'Assume the processed video has not been shifted spatially.\n');
    
elseif exist('jclips_sr_unfilt') ,

    if isfield(status,'spatial') && isfield(status.spatial,'shift') && status.spatial.shift > 0,
        fprintf(fid, '%d clips appear to have spatial shift in excess of what this algorithm can predict\n\n', status.spatial.scale);
    end
    if isfield(status.spatial,'failed') && status.spatial.failed >= 1,
        fprintf(fid, '\nWARNING\nFinal (filtered) spatial shift failed for %d clips.\n', status.spatial.failed);
        fprintf(fid, 'Un-filtered spatial shift will indicate "nan" for clips where algorithm failed.\n');
        fprintf(fid, 'Assume no shift for these clips -- which is likely incorrect.\n\n');
    end

    ssf_txt='\n   - - -   Un-Filtered    - - - \n';
    fprintf(fid,ssf_txt);
    line2 = ['Clip Name                          H   V \n'];
    line3 = ['--------------                    --- --- \n'];

    fprintf(fid, line2); fprintf(fid, line3);

    % now, assign values to table
    for count = 1:num_clips
        clip_name = cell2mat(model_op.clip_name(count));
        for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
            jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
            if strcmp(jclips_name, clip_name),
                spatial_h = jclips_sr_unfilt(jc).spatial.horizontal;
                spatial_v = jclips_sr_unfilt(jc).spatial.vertical;
                fprintf (fid, '%-30s    %3d %3d \n' ,clip_name, spatial_h, spatial_v);
            end
        end
    end

    ssf_txt='\n   - - -   Filtered    - - - \n';
    fprintf(fid,ssf_txt);

    line2 = ['HRC Name                           H   V  \n'];
    line3 = ['--------------                    --- --- \n'];

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
            fprintf (fid, '%-30s    %3d %3d\n', jclips(jc).hrc{1}, spatial_h, spatial_v);
        end

    end

else
    % manual spatial shift

    line2 = ['HRC Name                           H   V  \n'];
    line3 = ['--------------                    --- --- \n'];

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
            fprintf (fid, '%-30s    %3d %3d\n', jclips(jc).hrc{1}, spatial_h, spatial_v);
        end

    end
    
end

% now, check for errors
[hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
num_warn = 0;
for count = 1:length(hrc_offsets),
    jc = (hrc_offsets{count}(1));
    jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
    [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
    if search_status == 0,
        if abs(jclips(jc).spatial.horizontal) >  8 || abs(jclips(jc).spatial.vertical) > 5,
            if num_warn == 0;
                num_warn = 1;
                fprintf(fid,'\n');
            end
            fprintf (fid, 'Warning: Extreme spatial shift detected in %s\n', jclips(jc).hrc{1});
        end
    end

end
if num_warn == 1,
    fprintf(fid,'\n');
end

fprintf (fid, '\n\n.......... Luminance Gain & Offset ..........\n');
if strcmpi(cell2mat(model_selected),'PSNR with Variable Frame Delay'),
    fprintf(fid, '\nWARNING\nLuminance gain & offset overriden by model. Model values reported.\n\n');
end
if isfield(status, 'gain_offset') && status.gain_offset ...
        && ~strcmpi(cell2mat(model_selected),'PSNR with Variable Frame Delay'),
    % error, 
    fprintf(fid, '\nERROR\nFatal error computing Gain & Offet.\n');
    fprintf(fid, 'Assume luminance gain & offset of processed video matches original.\n\n');

elseif (isequal(cal_type, ct_NC) || isequal(cal_type, ct_T)) ...
        && ~strcmpi(cell2mat(model_selected),'PSNR with Variable Frame Delay'),
    fprintf(fid, 'Assume luminance gain & offset of processed video matches original.\n\n');
    
elseif exist('jclips_lgo_unfilt') ...
        && ~strcmpi(cell2mat(model_selected),'PSNR with Variable Frame Delay'),


    ssf_txt='\n   - - -   Un-Filtered    - - - \n';
    fprintf(fid,ssf_txt);
    line2 = ['Clip Name                          Gain    Offset\n'];
    line3 = ['--------------                    ------   ------\n'];

    fprintf(fid, line2); fprintf(fid, line3);

    % now, assign values to table
    for count = 1:num_clips
        clip_name = cell2mat(model_op.clip_name(count));
        for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
            jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
            if strcmp(jclips_name, clip_name)
                y_gain = jclips_lgo_unfilt(jc).luminance_gain;
                y_offset = jclips_lgo_unfilt(jc).luminance_offset;
                fprintf (fid, '%-30s    %7.4f  %7.3f\n' ,clip_name, y_gain, y_offset);
            end
        end
    end

    ssf_txt='\n   - - -   Filtered    - - - \n';
    fprintf(fid,ssf_txt);

    line2 = ['HRC Name                           Gain    Offset\n'];
    line3 = ['--------------                    ------   ------\n'];

    fprintf(fid, line2); fprintf(fid, line3);

    % now, assign values to table
    [hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
    for count = 1:length(hrc_offsets),
        jc = (hrc_offsets{count}(1));
        jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
        [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
        if search_status == 0,
            y_gain = jclips(jc).luminance_gain;
            y_offset = jclips(jc).luminance_offset;
            fprintf (fid, '%-30s    %7.4f  %7.3f\n', jclips(jc).hrc{1}, y_gain, y_offset);
        end

    end

    % Report if output from luminance_gain_offset has extreme values:
    if ~isempty(lgo_status)
        fprintf(fid, '   Extreme values for luminance gain and/or offset detected.\n');
        fprintf(fid, '   Algorithm failure assumed -- results for impacted HRCs discarded.\n');
        fprintf(fid, '   Assume luminance gain & offset of processed video matches original.\n    ');
        for jc = 1:length(lgo_status),
            fprintf(fid, '        %s\n', lgo_status{jc});
        end
        fprintf('\n');
    end
else
    % manual 

    line2 = ['HRC Name                           Gain    Offset\n'];
    line3 = ['--------------                    ------   ------\n'];

    fprintf(fid, line2); fprintf(fid, line3);

    % now, assign values to table
    [hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
    for count = 1:length(hrc_offsets),
        jc = (hrc_offsets{count}(1));
        jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
        [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
        if search_status == 0,
            y_gain = jclips(jc).luminance_gain;
            y_offset = jclips(jc).luminance_offset;
            fprintf (fid, '%-30s    %7.4f  %7.3f\n', jclips(jc).hrc{1}, y_gain, y_offset);
        end

    end
end

% check for warnings
[hrc_offsets] = sort_clips_by('hrc',jclips,jtests);
num_warn = 0;
for count = 1:length(hrc_offsets),
    jc = (hrc_offsets{count}(1));
    jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
    [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
    if search_status == 0,
        if jclips(jc).luminance_gain < 0.9 || jclips(jc).luminance_gain > 1.1,
            if num_warn == 0;
                num_warn = 1;
                fprintf(fid,'\n');
            end
            fprintf(fid, 'Warning: extreme luminance gain detected in %s\n', jclips(jc).hrc{1});
        end
        if jclips(jc).luminance_offset < -20 || jclips(jc).luminance_offset > 20,
            if num_warn == 0;
                num_warn = 1;
                fprintf(fid,'\n');
            end
            fprintf(fid, 'Warning: extreme luminance offset detected in %s\n', jclips(jc).hrc{1});
        end
    end

end
if num_warn == 1,
    fprintf(fid,'\n');
end

    



fprintf (fid, '\n\n.......... Color Gain & Offset ..........\n');

if ~isequal(cal_type, ct_RR2),
    fprintf(fid, 'Color Gain & Offset were not estimated.\n');
    
elseif isfield(status, 'gain_offset') && status.gain_offset,
    % error, 
    fprintf(fid,'\nERROR\nFatal error computing Color Gain & Offet.\n');
    fprintf(fid, 'Color Gain & Offset were not estimated.\n\n');
    
else
    fprintf(fid, 'Color Gain & Offset were estimated but not removed.\n\n');

    ssf_txt='\n   - - -   Un-Filtered    - - - \n';
    fprintf(fid,ssf_txt);
    
    line2 = ['Clip Name                         Cb Gain  Cb Offset  Cr Gain  Cr Offset\n'];
    line3 = ['--------------                    -------  ---------  -------  ---------\n'];

    fprintf(fid, line2); fprintf(fid, line3);

    % now, assign values to table
    for count = 1:num_clips
        clip_name = cell2mat(model_op.clip_name(count));
        for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
            jclips_name = strcat(jclips_lgo_unfilt(jc).test, '_', jclips_lgo_unfilt(jc).scene, '_', jclips_lgo_unfilt(jc).hrc);
            if strcmp(jclips_name, clip_name)
                cb_gain = jclips_lgo_unfilt(jc).cb_gain;
                cb_offset = jclips_lgo_unfilt(jc).cb_offset;
                cr_gain = jclips_lgo_unfilt(jc).cr_gain;
                cr_offset = jclips_lgo_unfilt(jc).cr_offset;
                fprintf (fid, '%-30s    %7.4f %7.3f    %7.4f  %7.3f\n', jclips_name{1}, cb_gain, cb_offset, cr_gain, cr_offset);
            end
        end
    end

    ssf_txt='\n   - - -   Filtered    - - - \n';
    fprintf(fid,ssf_txt);

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


% COMMON VALID REGION SUMMARY:
fprintf (fid, '\n\n.......... Valid Region ..........\n');

if isfield(status, 'valid') && status.valid,
    % error, 
    fprintf(fid,'\nERROR\nFatal error computing valid region.\n');
    fprintf(fid,'Assume standard valid regions (e.g., discard overscan).\n\n');

elseif isequal(cal_type, ct_NC),
    fprintf(fid,'Assume standard valid regions (e.g., discard overscan).\n\n');

elseif exist('jclips_vr_unfilt') ,

    ssf_txt='\n   - - -   Un-Filtered    - - - \n';
    fprintf(fid,ssf_txt);

    line2 = ['Clip Name                       top  bottom    left right\n'];
    line3 = ['--------------                  ---- ------    ---- -----\n'];

    fprintf(fid, line2); fprintf(fid, line3);

    % now, assign values to table
    for count = 1:num_clips
        clip_name = cell2mat(model_op.clip_name(count));
        for jc = 1:sz_jclips % loop through all jclips items, and assign spatial values
            jclips_name = strcat(jclips_vr_unfilt(jc).test, '_', jclips_vr_unfilt(jc).scene, '_', jclips_vr_unfilt(jc).hrc);
            if strcmp(jclips_name, clip_name)
                top  = jclips_vr_unfilt(jc).cvr.top;
                bot  = jclips_vr_unfilt(jc).cvr.bottom;
                left  = jclips_vr_unfilt(jc).cvr.left;
                right = jclips_vr_unfilt(jc).cvr.right;

                fprintf (fid, '%-30s  %4d %4d      %4d %4d\n',clip_name, top, bot, left, right);        end
        end
    end
    
    ssf_txt='\n   - - -   Filtered    - - - \n';
    fprintf(fid,ssf_txt);

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
else

    line2 = ['Scene Name                      top  bottom    left right\n'];
    line3 = ['--------------                  ---- ------    ---- -----\n'];
    fprintf(fid, line2);fprintf(fid, line3);

    % now, assign values to table
    [scene_offsets] = sort_clips_by('scene',jclips,jtests);
    for count = 1:length(scene_offsets),
        if length(scene_offsets{count}) < 2,
            continue;
        end
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

% check for problems
[scene_offsets] = sort_clips_by('scene',jclips,jtests);
num_warn = 0;
for count = 1:length(scene_offsets),
    if length(scene_offsets{count}) < 2,
        continue;
    end
    jc = (scene_offsets{count}(2));
    jclips_name = strcat(jclips(jc).test, '_', jclips(jc).scene, '_', jclips(jc).hrc);
    [answer, search_status] = pars_find_clip(model_op, 'clip', jclips_name, 'index');
    if search_status == 0,
        if ( jclips(jc).cvr.bottom - jclips(jc).cvr.top + 1) / jclips(jc).image_size.rows < 0.55 || ...
                ( jclips(jc).cvr.right - jclips(jc).cvr.left + 1) / jclips(jc).image_size.cols < 0.80,
            if num_warn == 0;
                num_warn = 1;
                fprintf(fid,'\n');
            end
            fprintf(fid, 'Warning: Greatly reduced valid video region detected in %s.\n', jclips(jc).scene{1});
        end
        
    end
end
if num_warn == 1,
    fprintf(fid,'\n');
end

% TEMPORAL REG SUMMARY:
fprintf (fid, '\n\n.......... Temporal Registration ..........\n');
if isfield(status,'temporal') && isfield(status.temporal, 'error') && status.temporal.error == 1,
    fprintf(fid, '\nERROR\nFatal error encountered during temporal registration.\n');
    fprintf(fid, 'Assume that the first valid frame of the processed video file aligns in time with\n');
    fprintf(fid, 'the first valid frame of the matching original video file.\n');

elseif isequal(cal_type, ct_NC)
    fprintf(fid, 'Assume that the first valid frame of the processed video file aligns in time with\n');
    fprintf(fid, 'the first valid frame of the matching original video file.\n');

else

    if uncert >= 0,
        fprintf(fid, '\nTemporal registration uncertainty was +/- %f seconds\n\n', uncert);
    end
    
    if isfield(status,'temporal'),
        
        if isfield(status, 'temporal2'),
            fprintf(fid, 'First run of temporal registration indicates:\n');
        end
        if isfield(status.temporal,'uncertainty'),
            fprintf(fid, '%d clips should be re-run with a larger temporal uncertainty\n', status.temporal.uncertainty);
        end
        if isfield(status.temporal, 'still'),
            fprintf(fid, '%d clips came from still scenes (or extremely impaired scenes)\n', status.temporal.still);
        end
        if isfield(status.temporal, 'ambiguous'),
            fprintf(fid, '%d clips have an ambiguous temporal registration\n', status.temporal.ambiguous);
        end
        fprintf(fid,'\n');
    end
    
    if isfield(status, 'temporal2'),
        if isfield(status.temporal2, 'error') && status.temporal2.error == 1,
            fprintf(fid, '\nSecond run of temporal registration,\nFatal error encountered during temporal registration.\n');
            fprintf(fid, 'Default temporal registration used (see below).\n');
        else
        
            fprintf(fid,'\nSecond run of temporal registration indicates:\n');
            if isfield(status.temporal2,'uncertainty'),
                fprintf(fid, '%d clips should be re-run with a larger temporal uncertainty\n', status.temporal2.uncertainty);
            end
            if isfield(status.temporal2, 'still'),
                fprintf(fid, '%d clips came from still scenes (or extremely impaired scenes)\n', status.temporal2.still);
            end
            if isfield(status.temporal2, 'ambiguous'),
                fprintf(fid, '%d clips have an ambiguous temporal registration\n', status.temporal2.ambiguous);
            end
            fprintf(fid,'\n');
        end
    
    end

    line1 = ['.                                      Original       Processed\n'];
    line2 = ['Clip Name                       delay  start  stop    start  stop\n'];
    line3 = ['--------------                  -----  ----- -----    -----  -----\n'];

    fprintf(fid, line1); fprintf(fid, line2); fprintf(fid, line3);

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
                fprintf (fid, '%-30s  %5.1f  %3d    %3d     %3d     %3d\n', clip_name, delay, orig_start, orig_stop, align_start, align_stop);
            end
        end
    end
end


fclose(fid);
