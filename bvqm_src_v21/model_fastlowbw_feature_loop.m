function status = model_fastlowbw_feature_loop(test_structs, clip_structs, feature_base_dir, varargin)
% model_fastlowbw_feature_loop
%  Loop through a list of clips; compute features for fast lowbw model.
% SYNTAX
%  [status] = model_fastlowbw_feature_loop(test_structs, clip_structs,feature_base_dir)
%  model_fastlowbw_feature_loop(...,'PropertyName',...);
% DESCRIPTION
%  This function computes lowbandwidth model features.  Store in standard
%  named directories, as if used feature_loop_si_hv_adapt, etc.
%  This function takes variable 'test_structs' (of the same format as
%  GTests), which describes each video test and the location of the
%  associated video; and variable 'clip_structs' (of the same format as
%  GClips), which lists the clips to be processed. 
%   
%  Filter size adjusts automatically for the image size, according to
%  function 'adaptive_filter'.
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
%
% These features could have been computed with the following functions:
%     feature_frame_color_ati_sec(test_structs,clip_structs_sroi,feature_base_dir,0.2);
%     feature_loop_si_hv_adapt(test_structs,clip_structs_pvr, 30,30,1,feature_base_dir);
%     feature_loop_coherent_color(test_structs,clip_structs_sroi, 30,30,1,feature_base_dir);
%     feature_loop_cont(test_structs,clip_structs_sroi, 30,30,1,feature_base_dir,'BlockMean');
% However, the above would not have the +/-1 pixel shifts.

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

try
    hsize = 30;
    vsize = 30;
    
    feature_fastlowbw_name = 'feature_fastlowbw_model/';
    
    feature_si_name = 'feature_avg1s_Y_siA_30x30_std/';
    feature_hv_name = 'feature_avg1s_Y_hvA_30x30_mean_ratio/';
    feature_y_name = 'feature_Y_cont_30x30_1s_mean/';
    feature_cb_name = 'feature_color_coher_cb_30x30_1s_mean/';
    feature_cr_name = 'feature_color_coher_cr_30x30_1s_mean/';
    feature_ati_name = 'feature_Y_randomati0.2s_image_rms/';
    
    for loop = 1:max(size(clip_structs)),
        t = clock;
        [exists, file_name] = write_feature_standard(0,feature_base_dir,...
            feature_fastlowbw_name,clip_structs(loop),'exists');
        
        if min(exists) > 0,
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
        
        sroi = clip_structs(loop).cvr;
        [filter_size, extra] = adaptive_filter(clip_structs(loop).image_size);
        [valid, cvr, sroi] = model_lowbw_sroi(extra, sroi.top, sroi.left, sroi.bottom, sroi.right);

        if strcmp('original', clip_structs(loop).hrc{1}),
            [si_orig hv_feat_orig luma_orig cb_orig cr_orig ati_orig] = ...
                model_fastlowbw_features ('clip', test_structs(tnum), clip_structs(loop), sroi);
            
            data.si_std = si_orig;
            data.hv_ratio = hv_feat_orig;
            data.y_mean = luma_orig; 
            data.cb_mean = cb_orig;
            data.cr_mean = cr_orig;
            data.ati_rms = ati_orig;
            center = 1;
        else
            [data] = model_fastlowbw_features_shift ('clip', test_structs(tnum), clip_structs(loop), sroi);
            center = 5;
        end
        
        % Write each features all at once to a single file.
        write_feature_standard(data, feature_base_dir,feature_fastlowbw_name, clip_structs(loop));

        % write center shift features to individual files, for separate
        % searches & examination if desired.
        write_feature_standard(data(center).si_std, feature_base_dir,feature_si_name, clip_structs(loop));
        write_feature_standard(data(center).hv_ratio, feature_base_dir,feature_hv_name, clip_structs(loop));
        write_feature_standard(data(center).y_mean, feature_base_dir,feature_y_name, clip_structs(loop));
        write_feature_standard(data(center).cb_mean, feature_base_dir,feature_cb_name, clip_structs(loop));
        write_feature_standard(data(center).cr_mean, feature_base_dir,feature_cr_name, clip_structs(loop));
        write_feature_standard(data(center).ati_rms, feature_base_dir,feature_ati_name, clip_structs(loop));

        % Erase the large variable 'data' from memory to prevent problems.
        clear data;
    end
catch
    status = 1;
    if verbose,
        fprintf('\n%s\n', lasterr);
        fprintf('\tSkipping clip.\n');
    end
end