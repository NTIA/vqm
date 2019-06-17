function status = feature_frame_color_ati_sec(test_structs, clip_structs,...
    feature_base_dir, tis_sec, varargin)
% FEATURE_FRAME_COLOR_ATI_SEC
%  Loop through a list of clips, compute color ATI features with more than
%  1F differences, and save each clip's features in a separate file.
%  Do frame-based / whole frame features.
% SYNTAX
%  feature_frame_color_ati_sec(test_structs, clip_structs,feature_base_dir, tis_sec)
%  feature_frame_color_ati_sec(...,'PropertyName',...);
% DESCRIPTION
%  Compute ATI feature
% 
%  This function takes variable 'test_structs' (of the same format as
%  GTests), which describes each video test and the location of the
%  associated video; and variable 'clip_structs' (of the same format as
%  GClips), which lists the clips to be processed. 'feature_base_dir' is
%  the path to the directory where the feature's sub-directory should be
%  created.  That sub-directory's name will be this feature's name.
%  'tis_sec' is the fraction of a second spanning images to be used for the TI
%  calculation (e.g., 2/30 for 2F, 0.2 for 6F NTSC / 5F Pal, etc).  Will be
%  rounded to the nearest integer.  Note that a 2F span will difference
%  image 1F appart (i.e., range must include end points).
%
%  The following optional parameters are also available.  
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).
%
%  Return 0 if operated correctly; 1 if an error was encountered.
% NOTES
%  This feature loop function holds one time-slice of images in memory all at
%  once; and also holds in memory the 3D feature matrix associated with the 
%  current video clip.
status = 0;

verbose = 0;
for cnt = 1:length(varargin),
    if strcmp(lower(varargin{cnt}),'verbose'),
        verbose = 1;
    elseif strcmp(lower(varargin{cnt}),'quiet'),
        verbose = 0;
    else
        error('optional property not recognized');
    end
end

tslice_length_sec = 1/30; % one image at a time, exactly

specific = sprintf('ati%.4gs', tis_sec);
feature_name_y_mean = standard_feature_names('Y',...
        'fullimage',...
        'block_stat','mean', 'specific',specific);
feature_name_y_std = standard_feature_names('Y',...
        'fullimage',...
        'block_stat','std', 'specific',specific);
feature_name_y_rms = standard_feature_names('Y',...
        'fullimage',...
        'block_stat','rms', 'specific',specific);

feature_name_cb_mean = standard_feature_names('Cb',...
        'fullimage',...
        'block_stat','mean', 'specific',specific);
feature_name_cb_std = standard_feature_names('Cb',...
        'fullimage',...
        'block_stat','std', 'specific',specific);
feature_name_cb_rms = standard_feature_names('Cb',...
        'fullimage',...
        'block_stat','rms', 'specific',specific);

feature_name_cr_mean = standard_feature_names('Cr',...
        'fullimage',...
        'block_stat','mean', 'specific',specific);
feature_name_cr_std = standard_feature_names('Cr',...
        'fullimage',...
        'block_stat','std', 'specific',specific);
feature_name_cr_rms = standard_feature_names('Cr',...
        'fullimage',...
        'block_stat','rms', 'specific',specific);
    
feature_name_ycbcr_mean = standard_feature_names('YCbCr',...
        'fullimage',...
        'block_stat','mean', 'specific',specific);
feature_name_ycbcr_std = standard_feature_names('YCbCr',...
        'fullimage',...
        'block_stat','std', 'specific',specific);
feature_name_ycbcr_rms = standard_feature_names('YCbCr',...
        'fullimage',...
        'block_stat','rms', 'specific',specific);

        
for loop = 1:max(size(clip_structs)),
    try,
        t = clock;
        [name_exists_mean] = write_feature_standard(0,feature_base_dir,feature_name_y_mean,clip_structs(loop),'exists');
        [name_exists_std] = write_feature_standard(0,feature_base_dir,feature_name_y_std,clip_structs(loop),'exists');
        [name_exists_rms] = write_feature_standard(0,feature_base_dir,feature_name_y_rms,clip_structs(loop),'exists');
        if name_exists_mean && name_exists_std && name_exists_rms,
            if verbose,
                fprintf('\tfeature already computed for clip %s:%s(%s)\n',...
                    clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
                    clip_structs(loop).hrc{1});
            end
            continue;
        end
        if verbose,
            fprintf('Clip %d of %d ==> %s:%s(%s) at %d:%d\n', loop, max(size(clip_structs)), ...
                clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
                clip_structs(loop).hrc{1}, t(4), t(5) );
        end
	
        % Find the offset of the test structure for this clip.
        tnum = search_test_list(test_structs, clip_structs(loop));
        
        % Compute the default adjusted SROI and number of blocks available, so
        % can allocate memory to hold all features for this clip.
        [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(loop));
        number_tslices = total_tslices(clip_structs(loop),tslice_length_sec, 'TIS', tis_sec);
        data_y_mean = zeros(1, 1, number_tslices);
        data_y_std = zeros(1, 1, number_tslices);
        data_y_rms = zeros(1, 1, number_tslices);
        data_cb_mean = zeros(1, 1, number_tslices);
        data_cb_std = zeros(1, 1, number_tslices);
        data_cb_rms = zeros(1, 1, number_tslices);
        data_cr_mean = zeros(1, 1, number_tslices);
        data_cr_std = zeros(1, 1, number_tslices);
        data_cr_rms = zeros(1, 1, number_tslices);
        data_ycbcr_mean = zeros(1, 1, number_tslices);
        data_ycbcr_std = zeros(1, 1, number_tslices);
        data_ycbcr_rms = zeros(1, 1, number_tslices);
        
        % For each time-slice in this clip, read the time-slice; perceptually filter
        % (none); and compute block statistic.
        [tis_frames] = tslice_conversion(tis_sec, clip_structs(loop).fps);
        for cnt = 1:number_tslices,
            [y,cb,cr] = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec,cnt, ...
                'TIS', tis_sec, ...
                'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
            [rows,cols,time] = size(y);

            ati_y = filter_ati(y, tis_frames);
            clear y;
            [data_y_mean(:,:,cnt), data_y_std(:,:,cnt), data_y_rms(:,:,cnt)]  = ...
                block_statistic (ati_y, rows, cols, 'mean','std','rms');

            ati_cb = filter_ati(cb, tis_frames);
            clear cb;
            [data_cb_mean(:,:,cnt), data_cb_std(:,:,cnt), data_cb_rms(:,:,cnt)]  = ...
                block_statistic (ati_cb, rows, cols, 'mean','std','rms');

            ati_cr = filter_ati(cr, tis_frames);
            clear cr;
            [data_cr_mean(:,:,cnt), data_cr_std(:,:,cnt), data_cr_rms(:,:,cnt)]  = ...
                block_statistic (ati_cr, rows, cols, 'mean','std','rms');
            
            ati_ycbcr = sqrt(ati_y.^2 + ati_cb.^2 + ati_cr.^2);
            [data_ycbcr_mean(:,:,cnt), data_ycbcr_std(:,:,cnt), data_ycbcr_rms(:,:,cnt)]  = ...
                block_statistic (ati_ycbcr, rows, cols, 'mean','std','rms');

            clear ati_y ati_cb ati_cr ati_ycbcr;
        end
        
        % Write features all at once to a single file.
        write_feature_standard(data_y_mean, feature_base_dir,feature_name_y_mean, clip_structs(loop));
        write_feature_standard(data_y_std, feature_base_dir,feature_name_y_std, clip_structs(loop));
        write_feature_standard(data_y_rms, feature_base_dir,feature_name_y_rms, clip_structs(loop));
        write_feature_standard(data_cb_mean, feature_base_dir,feature_name_cb_mean, clip_structs(loop));
        write_feature_standard(data_cb_std, feature_base_dir,feature_name_cb_std, clip_structs(loop));
        write_feature_standard(data_cb_rms, feature_base_dir,feature_name_cb_rms, clip_structs(loop));
        write_feature_standard(data_cr_mean, feature_base_dir,feature_name_cr_mean, clip_structs(loop));
        write_feature_standard(data_cr_std, feature_base_dir,feature_name_cr_std, clip_structs(loop));
        write_feature_standard(data_cr_rms, feature_base_dir,feature_name_cr_rms, clip_structs(loop));
        write_feature_standard(data_ycbcr_mean, feature_base_dir,feature_name_ycbcr_mean, clip_structs(loop));
        write_feature_standard(data_ycbcr_std, feature_base_dir,feature_name_ycbcr_std, clip_structs(loop));
        write_feature_standard(data_ycbcr_rms, feature_base_dir,feature_name_ycbcr_rms, clip_structs(loop));
        
        % Erase the large variable 'data' from memory to prevent problems.
        clear data_y_mean data_y_std data_y_rms;
        clear data_cb_mean data_cb_std data_cb_rms;
        clear data_cr_mean data_cr_std data_cr_rms;
        clear data_ycbcr_mean data_ycbcr_std data_ycbcr_rms;
    catch
        status = 1;
        if verbose,
            fprintf('\n%s\n', lasterr);
            fprintf('\tSkipping clip.\n'); 
        end
    end
end
