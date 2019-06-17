function output = bvqm_model_vqm_vfd(testset, clipset, feature_base_dir, varargin)

%#function nn_8par
%#function network

feature_base_dir = [feature_base_dir '\'];

verbose = 0;
t_uncert = 30;
i = 1;
while i <= length(varargin)
    if strcmp(varargin{i}, 'verbose')
        verbose = 1;
        i = i + 1;
    elseif strcmp(varargin{i}, 't_uncert')
        t_uncert = ceil(varargin{i+1}*clipset(1).fps);
        i = i + 2;
    else
        error('Unknown optional input for bvqm_model_vqm_vfd');
    end
end

% Sort the clipset
clipset = clipset(sort_clips_by('none', clipset, testset));

% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This cell estimates and stores the VFD information for the clipset. This
% VFD information is referenced to the alignment information stored in 
% clipset. The VFD alignment information is stored in mat files called
% "test_vfd.mat" (in the current directory), where test is the video test.
% Thus, VFD results from all clips within the same test are stored
% together.

vfd_file = [feature_base_dir 'vfd_results.csv'];
fid_vfd = fopen(vfd_file,'w');
fprintf(fid_vfd,'Test, Scene, HRC, (Proc Orig)\n');
fclose(fid_vfd);

offsets = find_clip(clipset,'*','*','original','not'); % offsets of the processed clips
num_clips = length(offsets);
clear resultst results_rmset results_fuzzyt results_fuzzy_mset;
% Step through each processed clip and generate its VFD alignment results
curr = 1;
output.par_name{1} = 'vqm_vfd';
for i = 1:num_clips
    output.clip_name{curr} = [clipset(offsets(i)).test{:} '_' clipset(offsets(i)).scene{:} '_' clipset(offsets(i)).hrc{:}];
    output.inlsa_mos(curr) = clipset(offsets(i)).inlsa_mos;
    output.mos(curr) = clipset(offsets(i)).mos;
    curr = curr + 1;
    % Compute the default adjusted SROI and number of blocks available, so
    % can allocate memory to hold all features for this clip.
    [sroi] = adjust_requested_sroi (clipset(offsets(i)));
    
    % Modify SROI, because we want an even number of lines and top line to 
    % always be odd.
    if mod(sroi.top,2) == 0,
        sroi.top = sroi.top + 1;
    end
    if mod(sroi.bottom,2) == 1,
        sroi.bottom = sroi.bottom -1;
    end
    
    if ~exist(sprintf('vfd_%s',char(clipset(offsets(1)).test)),'var')
        try
            % Figure out number of processed frames available, aligned.
            tslice_length_sec = (clipset(offsets(i)).align_stop - ...
                clipset(offsets(i)).align_start + 1) / clipset(offsets(i)).fps;
            if is_reframing_indicated(clipset(offsets(i))),
                tslice_length_sec = tslice_length_sec - 1/clipset(offsets(i)).fps;
            end
            number_tslices = total_tslices(clipset(offsets(i)),tslice_length_sec);
            if number_tslices ~= 1,
                error('Computation of processed number tslices wrong');
            end

            % Read in processed clip:
            % if you want to read in the full image (with black border where
            % invalid), add 'full' to the list of arguments.
            yp = read_tslice(testset, clipset(offsets(i)), tslice_length_sec, 1, ...
                'sroi', sroi.top, sroi.left, sroi.bottom, sroi.right);

%             % Figure out number of original frames available, unaligned.
%             orig_offset = find_original(clipset, offsets(i));
%             tslice_length_sec = (clipset(orig_offset).loc_stop - ...
%                 clipset(orig_offset).loc_start + 1) / clipset(orig_offset).fps;
%             number_tslices = total_tslices(clipset(orig_offset),tslice_length_sec, 'unaligned');
%             if number_tslices ~= 1,
%                 error('Computation of original number tslices wrong');
%             end

            % compute the range of original video frames to be read.
            % Set the initial alignment point (in frames) of the original
            orig_offset = find_original(clipset, offsets(i));

            orig_start = max(clipset(orig_offset).align_start - t_uncert, clipset(orig_offset).loc_start);
            orig_stop = min(clipset(orig_offset).align_stop + t_uncert, clipset(orig_offset).loc_stop);
            first_align = clipset(orig_offset).align_start - orig_start + 1;

            % Read in original clip:  Need a different starting point, and a different
            % read length!
            yo = read_tslice(testset,clipset(orig_offset), 1, 1, ...
                'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right, 'all_frames', ...
                'align_start', orig_start, 'align_stop', orig_stop);
        catch
            fprintf('Error.  Video read failure!\n');
            return;
        end

        % Set the interlaced flag
        if (strcmpi(clipset(orig_offset).video_standard, 'interlace_lower_field_first') == 1)
            is_interlaced = 1;
            field_first = 1;
            first_align = 2*first_align-1;  % convert to fields
        elseif (strcmpi(clipset(orig_offset).video_standard, 'interlace_upper_field_first') == 1)
            is_interlaced = 1;
            field_first = 2;
            first_align = 2*first_align-1;  % convert to fields
        else
            is_interlaced = 0;  % progressive
        end

        % Estimate the variable frame delays of the processed clip
        if (is_interlaced)
            if verbose
                [resultst{i} results_rmset{i} results_fuzzyt{i} results_fuzzy_mset{i}] = ...
                    est_var_frame_delays(yp, yo, 'interlaced', field_first, 'reframe', ...
                    'first_align', first_align, 'causal', 'normalize', 't_uncert', t_uncert, 'verbose');
            else
                [resultst{i} results_rmset{i} results_fuzzyt{i} results_fuzzy_mset{i}] = ...
                    est_var_frame_delays(yp, yo, 'interlaced', field_first, 'reframe', ...
                    'first_align', first_align, 'causal', 'normalize', 't_uncert', t_uncert);
            end
        else
            if verbose
                [resultst{i} results_rmset{i} results_fuzzyt{i} results_fuzzy_mset{i}] = ...
                    est_var_frame_delays(yp, yo, 'first_align', first_align, 'causal', ...
                    'normalize', 'verbose', 't_uncert', t_uncert);
            else
                [resultst{i} results_rmset{i} results_fuzzyt{i} results_fuzzy_mset{i}] = ...
                    est_var_frame_delays(yp, yo, 'first_align', first_align, 'causal', ...
                    'normalize', 't_uncert', t_uncert);
            end
        end
        
        [~,~,nframes] = size(resultst{i});
        if (resultst{i} == 0)
            vfd_failed = 1;  % Set a logical variable to record that the VFD algorithm failed.
            if is_interlaced  % interlaced
                npts = nframes*2;
            else  % progressive
                npts = nframes;
            end
            resultst{i} = first_align:first_align+npts-1;
        else
            vfd_failed = 0;
            npts = length(resultst{i});
        end
        % note field or frame numbers for processed images examined.
        % also, compute matching original field or frame numbers. Note that
        % "resultst" starts at 1 where "orig_start" is located, thus the
        % subtraction of 1 from orig_start.
        if isequal(is_interlaced,0)
            proc_indices = (clipset(offsets(i)).align_start-1) + (1:length(resultst{i}));
            orig_indices = (orig_start - 1) + resultst{i}; 
        else % interlaced
            proc_indices = 2*(clipset(offsets(i)).align_start-1) + (1:length(resultst{i}));
            orig_indices = 2*(orig_start-1) + resultst{i};
        end
        
        fid_vfd = fopen(vfd_file,'a');
        do_reframing = is_reframing_indicated(clipset(offsets(i)));
        fprintf(fid_vfd,'%s, %s, %s,', clipset(offsets(i)).test{:}, clipset(offsets(i)).scene{:}, ...
            clipset(offsets(i)).hrc{:});
        for k = 1:npts-1
            fprintf(fid_vfd,'%d %d,',proc_indices(k)+do_reframing, orig_indices(k));
        end
        fprintf(fid_vfd,'%d %d\n',proc_indices(k+1)+do_reframing, orig_indices(k+1));
        fclose(fid_vfd);
    end
    
end
clear yp yo;
pause(0.1);

% % % modify alignment and VFD information. Discard frames off the beginning
% % % and end of every clip, so that the constant alignments are kept intact
% % % and no processed clip uses the first or last original frame.
% 
% % % variables to change: jclips, results, results_rmse, results_fuzzy,
% % % results_fuzzy_mse
% 
% okay this is rediculously complicated. I'm going to try to simplify.

clipset_proc = clipset(offsets);
all_offsets = sort_clips_by('scene', clipset_proc, testset);
for src=1:length(all_offsets),
    % curr_offsets is all PVS matching one SRC, with offsets that match
    % varaiable 'resultst'
    curr_offsets = all_offsets{src};
    
    % find the matching SRC in 'clipset'
    curr = find_clip(clipset, clipset_proc(curr_offsets(1)).test{1}, ...
        clipset_proc(curr_offsets(1)).scene{1}, 'original');
    o_clip = clipset(curr);

    % assume all of it is valid
    first_valid = 1;
    last_valid = length(resultst{curr_offsets(1)});

    % loop through all PVS.
    if is_interlaced,
        for pvs = 1:length(curr_offsets),
            % loop through resultst for this PVS frame by frame
            % find first frame that doesn't align to original frame #1
            tempR = resultst{curr_offsets(pvs)};
            for cnt = 1:length(tempR),
                if tempR(cnt) > 2,
                    break;
                end
            end
            first_valid = max(first_valid,cnt);

            % loop through resultst for this PVS frame by frame
            % find last frame that doesn't align to original frame #1
            for cnt = length(resultst{curr_offsets(pvs)}):-1:1,
                if tempR(cnt) < 2*o_clip.loc_stop-1,
                    break;
                end
            end
            last_valid = min(last_valid,cnt);
        end
        % convert first_valid and last_valid from fields to frames
        first_valid = ceil(first_valid/2);
        last_valid = floor(last_valid/2);
    else
        for pvs = 1:length(curr_offsets),
            % loop through resultst for this PVS frame by frame
            % find first frame that doesn't align to original frame #1
            tempR = resultst{curr_offsets(pvs)};
            for cnt = 1:length(tempR),
                if tempR(cnt) > 1,
                    break;
                end
            end
            first_valid = max(first_valid,cnt);

            % loop through resultst for this PVS frame by frame
            % find last frame that doesn't align to original frame #1
            for cnt = length(resultst{curr_offsets(pvs)}):-1:1,
                if tempR(cnt) < o_clip.loc_stop,
                    break;
                end
            end
            last_valid = min(last_valid,cnt);
        end
    end
    
    % loop through all clips. Move all of them inward.
    % need to change numbering from clipset with only PVSs (clipset_proc)
    % to the actual clipset used everywhere else (clipset).
    
    % also, print results to log file
    fid_vfd = fopen(vfd_file,'a');
    fprintf(fid_vfd,'\n');
    for pvs = 1:length(curr_offsets),
        curr = find_clip(clipset, clipset_proc(curr_offsets(pvs)).test{1}, ...
            clipset_proc(curr_offsets(pvs)).scene{1}, clipset_proc(curr_offsets(pvs)).hrc{1});
        
        fprintf(fid_vfd,'%s, %s, %s,', clipset(curr).test{1}, ...
            clipset(curr).scene{1}, clipset(curr).hrc{1});
        if is_interlaced,
            % figure out whether or not this clip was reframed.  If so, the
            % first and last field of the range won't be used.  Thus, the
            % effective range grows by 2!
            if is_reframing_indicated(clipset(curr)),
                clipset(curr).align_stop = clipset(curr).align_start + last_valid; 
                clipset(curr).align_start = clipset(curr).align_start + first_valid - 1; 

                fprintf(fid_vfd,' VQM-VFD model will use processed fields [%d..%d]\n',...
                    clipset(curr).align_start*2, ...
                    clipset(curr).align_stop*2 -1);
            else
                clipset(curr).align_stop = clipset(curr).align_start + last_valid - 1; 
                clipset(curr).align_start = clipset(curr).align_start + first_valid - 1; 

                fprintf(fid_vfd,' VQM-VFD model will use processed fields [%d..%d]\n',...
                    clipset(curr).align_start*2 - 1, ...
                    clipset(curr).align_stop*2);
            end

        else
            clipset(curr).align_stop = clipset(curr).align_start + last_valid - 1; 
            clipset(curr).align_start = clipset(curr).align_start + first_valid - 1; 

            fprintf(fid_vfd,' VQM-VFD model will use processed frames [%d..%d]\n',...
                clipset(curr).align_start, clipset(curr).align_stop);
        end
    end
    fclose(fid_vfd);

    % if interlaced, change range from frames to fields
    if is_interlaced,
        first_valid = first_valid * 2 - 1;
        last_valid = last_valid * 2;
    end
    % loop through all results. Move all of them inward.
    for cnt = 1:length(curr_offsets),
        resultst{curr_offsets(cnt)} = resultst{curr_offsets(cnt)}(first_valid:last_valid);
        results_fuzzyt{curr_offsets(cnt)} = results_fuzzyt{curr_offsets(cnt)}(:,first_valid:last_valid);
        results_fuzzy_mset{curr_offsets(cnt)} = results_fuzzy_mset{curr_offsets(cnt)}(:,first_valid:last_valid);
    end

 
end




% Save intermediate VFD results to mat files called vfd_test.mat. This is
% to maintain compatibility with existing code, where the VFD information
% is stored by test name.
clipset_proc = clipset(offsets); % clipset with only the proc clips, in the order they were computed
offsets_by_test = sort_clips_by('test', clipset_proc, testset); % Holds offsets of each test, proc clips only
offsets_by_test_all = sort_clips_by('test', clipset, testset); % Holds offsets of each test, proc + orig clips only
num_tests = length(offsets_by_test);
for i = 1:num_tests
    these_offsets = offsets_by_test{i};
    these_offsets_all = offsets_by_test_all{i};
    test_name = char(clipset_proc(these_offsets(1)).test);
    test_offset = find_test(testset, test_name);
    this_test = testset(test_offset);
    this_clipset = clipset(these_offsets_all); % Need to include proc and orig clips in this_clipset
    results = resultst(these_offsets);
    results_rmse = results_rmset(these_offsets);
    results_fuzzy = results_fuzzyt(these_offsets);
    results_fuzzy_mse = results_fuzzy_mset(these_offsets);
    save([feature_base_dir 'vfd_' test_name '.mat'], 'results','results_rmse','results_fuzzy','results_fuzzy_mse','this_test','this_clipset');
end


% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This cell performs the VFD feature extraction. It's inefficient in that
% the videos are re-read for every feature that is extracted. But the
% advantage is that existing code can be used "as is".


%  Define the spatial degrees and time extent of the ST blocks
deg_size = 0.4;  % in angular degrees, a square block spatially
time_size = 0.2;  % in seconds

% Extract SI and HV features
if verbose
    bvqm_vfd_feature_loop_si_hv_adapt(testset, clipset, deg_size, time_size, feature_base_dir, 'verbose');
else
    bvqm_vfd_feature_loop_si_hv_adapt(testset, clipset, deg_size, time_size, feature_base_dir);
end

% Extract Y (CONT) features
if verbose
    bvqm_vfd_feature_loop_cont(testset, clipset, deg_size, time_size, feature_base_dir, 'blockmean', 'verbose');
else
    bvqm_vfd_feature_loop_cont(testset, clipset, deg_size, time_size, feature_base_dir, 'blockmean');
end

% Extract TI features
if verbose
    bvqm_vfd_feature_loop_ti(testset, clipset, deg_size, time_size, feature_base_dir, 'blockmean', 'verbose');
else
    bvqm_vfd_feature_loop_ti(testset, clipset, deg_size, time_size, feature_base_dir, 'blockmean');
end

% Extract Mean Square Error (MSE) features
if verbose
    bvqm_vfd_feature_loop_mse(testset, clipset, deg_size, time_size, feature_base_dir, 'blockmean', 'verbose');
else
    bvqm_vfd_feature_loop_mse(testset, clipset, deg_size, time_size, feature_base_dir, 'blockmean');
end

% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This cell performs the parameter extraction. There are 8 parameters in
% the new VFD model. They will be computed in the following order:
% 
% 1. vfd_exp_hv_loss_below5%_mean_square_clip_0.06
% 2. Y_vfd_hvA_0.4deg_0.2s_mean_3_divide_log_gain_rms_rms
% 3. Y_vfd_siA_0.4deg_0.2s_std_12_ratio_loss_mean_above90%
% 4. Y_vfd_siA_0.4deg_0.2s_std_8_log_gain_above98%tail_rms
% 5. Y_vfd_ti_0.4deg_0.2s_rms_3_log_gain_STabove95%tail
% 6. Y_FR_vfd_0.4deg_0.2s_rms_minus_gain_STmean
% 7. vfd_par1, see NTIA TM-11-475
% 8. vfd_par1*psnr_vfd, see NTIA TM-11-475 for psnr_vfd description
%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hv_loss calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfd_hv_loss = vfd_clippar_loop_exp_hv_loss(testset, clipset, feature_base_dir);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hv_gain calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfd_hv_gain = vfd_parameter_dual_loop(testset, clipset, feature_base_dir, ...
         'feature_Y_vfd_hvA_0.4deg_0.2s_mean', 'feature_Y_vfd_hvbarA_0.4deg_0.2s_mean', ...
         'feature_Y_vfd_hvA_0.4deg_0.2s_mean', {'rms'}, {'rms'}, ...
         'MinThreshold', 3, 'divide', 'compare', 'log_gain');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% si_loss calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfd_si_loss = vfd_parameter_loop(testset, clipset, feature_base_dir, ...
        'feature_Y_vfd_siA_0.4deg_0.2s_std', 'ratio_loss', {'mean'}, {'above90%'}, 'MinThreshold', 12);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% si_gain calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfd_si_gain = vfd_parameter_loop(testset, clipset, feature_base_dir, ...
        'feature_Y_vfd_siA_0.4deg_0.2s_std', 'log_gain', {'above98%tail'}, {'rms'}, 'MinThreshold', 8);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ti_gain calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfd_ti_gain = vfd_parameter_loop(testset, clipset, feature_base_dir, ...
        'feature_Y_vfd_ti_0.4deg_0.2s_rms', 'log_gain', {'above95%tail'}, {'above95%tail'}, 'MinThreshold', 3, '3D');  % Use 3D collapsing

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rmse (FR root mean squared error).
% Like PSNR but just the average rmse
% (over ST blocks), but not log based.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfd_rmse = vfd_parameter_loop(testset, clipset, feature_base_dir, ...
        'feature_Y_FR_vfd_0.4deg_0.2s_rms', 'minus_gain', {'mean'}, {'mean'}, '3D');  % Use 3D collapsing

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vfd_par1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfd_par1 = vfd_clippar_loop_par1(testset, clipset, feature_base_dir);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vfd_par1*psnr_vfd (vfd_xpar)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psnr_vfd = vfd_clippar_loop_psnr(testset, clipset, feature_base_dir);
vfd_xpar = psnr_vfd;
vfd_xpar.par_name = {'vfd_par1*psnr_vfd'};
vfd_xpar.data = vfd_par1.data.*psnr_vfd.data;

% Combine all the parameters into one structure
vfd_8par = join_parameters(vfd_hv_loss, vfd_hv_gain, vfd_si_loss, vfd_si_gain, vfd_ti_gain, vfd_rmse, vfd_par1, vfd_xpar);

% Pull off the parameter data and format for Neural Network input
nn_obj = vfd_8par.data;


% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This cell computes the model by combining all the parameters using a
% previously trained Neural Network (NN). The trained NN is stored in
% nn_8par.mat, variable nn_8par. It requires an input matrix that is
% num_pars x num_clips, where the pars must be ordered as given in the
% above cell.


try
    load nn_8par;
catch err,
    h=msgbox('File nn_8par.mat must be in the same directory as BVQM.exe. Otherwise the VQM_VFD model cannot run.', 'Cannot load file nn_8par.mat', 'error');
    waitfor(h);
    error('Cannot load file nn_8par');
end
%#function nn_8par
%#function network
nn_model = nn_8par(nn_obj);
output.data = nn_model;

