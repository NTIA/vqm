function [output clipset] = bvqm_model_psnr_vfd(testset, clipset, feature_base_dir, varargin)
% returned clipset has the luminance gain & offset as calculated by the
% PSNR_VFD model.


% reorder clipset by SRC then HRC.
order = sort_clips_by('scene', clipset, testset);
all_order = [];
for i=1:length(order),
    all_order = [all_order, order{i}];
end
clipset = clipset(all_order);
clear all_order;
clear order;


% make sure path ends with separator.
feature_base_dir = [feature_base_dir '\'];

% Define the peak signal level
peak = 255.0;

verbose = 0;
quiet = 0;
t_uncert = 30;
i = 1;
while i <= length(varargin)
    if strcmp(varargin{i}, 't_uncert')
        t_uncert = ceil(varargin{i+1}*clipset(1).fps);
        i = i + 2;
    elseif strcmp(varargin{i}, 'verbose')
        verbose = 1;
        i = i + 1;
    elseif strcmp(varargin{i}, 'quiet')
        quiet = 1;
        i = i + 1;
    elseif strcmp(varargin{i}, 'peak')
        if strcmp(varargin{i+1}, '255')
            peak = 255.0;
        else
            peak = 235.0;
        end
        i = i + 2;
    else
        error('Invalid optional argument for bvqm_model_psnr_vfd');
    end
end

if quiet && verbose
    error('bvqm_model_psnr_vfd: Cannot choose both ''quiet'' and ''verbose'' as options');
end

% Define the sub-sampling factor on the pixels for the final gain and
% offset adjusting fit, which is performed right before calculation of
% PSNR_VFD.
fraction_sampled = 1;

dataset = clipset(1).test;
for i = 2:length(clipset)
    if ~strcmp(clipset(i).test, dataset)
        error('Only include clips from one data set');
    end
end

% sort clips by HRCs
offsets = sort_clips_by('scene', clipset, testset);

output.clip_name = [];
output.insla_mos = [];
output.mos = [];
output.data = [];
output.par_name{1} = sprintf('psnr_vfd_%d', peak);

vfd_file = [feature_base_dir 'vfd_results.csv'];
fid_vfd = fopen(vfd_file,'w');
fprintf(fid_vfd,'Test, Scene, HRC, (Proc Orig)\n');
fclose(fid_vfd);

vfd_text = sprintf('\n');
vfd_model_text = sprintf('\nTest, Scene, HRC, PSNR_VFD_%d, Par1, Par2\n', peak);


        
index = 1;
curr = 1;
num_scenes = length(offsets);
for i = 1:num_scenes
    curr_offsets = offsets{i};
    
    num_clips = length(curr_offsets) - 1;
    
    psnr_ave = 0;  % psnr_vfd average summer for this HRC
    par1_ave = 0;  % par1_vfd average summer for this HRC
    par2_ave = 0;  % par2_vfd average summer for this HRC
    
    % skip #1 -- the original
    for j = 2:(num_clips + 1)
        scan_type = clipset(curr_offsets(j)).video_standard;
   
        % Assign the scene information to the results_vfd structure
        results_vfd(index).test = clipset(curr_offsets(j)).test{:};
        results_vfd(index).scene = clipset(curr_offsets(j)).scene{:};
        results_vfd(index).hrc = clipset(curr_offsets(j)).hrc{:};
        
        clip_name = sprintf('%s_%s_%s', results_vfd(index).test, results_vfd(index).scene, results_vfd(index).hrc);
        
        % Compute the default adjusted SROI
        sroi = adjust_requested_sroi(clipset(curr_offsets(j)));

        
        tslice_length_sec = (clipset(curr_offsets(j)).align_stop - ...
            clipset(curr_offsets(j)).align_start + 1) / clipset(curr_offsets(j)).fps;
        if is_reframing_indicated(clipset(curr_offsets(j)))
            tslice_length_sec = tslice_length_sec - 1/clipset(curr_offsets(j)).fps;
        end
        number_tslices = total_tslices(clipset(curr_offsets(j)),tslice_length_sec);        
        if number_tslices ~= 1
            error('Computation of number tslices wrong');
        end        
        

        vfd_results_exists = write_feature_standard(0, feature_base_dir, 'vfd_results_psnr', clipset(curr_offsets(j)), 'exists');
        
        if (vfd_results_exists)
            if (verbose)
                fprintf('\tvfd info already computed for clip %s:%s(%s)\n', clipset(curr_offsets(j)).test{:}, clipset(curr_offsets(j)).scene{:}, clipset(curr_offsets(j)).hrc{:});
            end
            load([feature_base_dir '\vfd_results_psnr\' clip_name '.mat']);
            results = data;
            results_fuzzy = datao;
            clear data datao;
        else
            % Read in yp
            y_proc = read_tslice(testset, clipset(curr_offsets(j)), tslice_length_sec, 1, ...
                'sroi', sroi.top, sroi.left, sroi.bottom, sroi.right);

            % compute the range of original video frames to be read.
            % Set the initial alignment point (in frames) of the original
            orig_offset = find_original(clipset, curr_offsets(j));
            
            orig_start = max(clipset(orig_offset).align_start - t_uncert, clipset(orig_offset).loc_start);
            orig_stop = min(clipset(orig_offset).align_stop + t_uncert, clipset(orig_offset).loc_stop);
            first_align = clipset(orig_offset).align_start - orig_start + 1;

            % Read in yo.  Need a different starting point, and a different
            % read length!
            y_orig = read_tslice(testset, clipset(orig_offset), 1, 1, ...
                'sroi', sroi.top, sroi.left, sroi.bottom, sroi.right, ...
                'all_frames', 'align_start', orig_start, 'align_stop', orig_stop);

            % Convert everythign to double precision before any calculations
            % are performed.
            y_orig = double(y_orig);
            y_proc = double(y_proc);

            % Set the interlaced flag
            if (strcmpi(clipset(orig_offset).video_standard, 'interlace_lower_field_first') == 1)
                is_interlaced = 1;
                field_first = 1;
                first_align = 2*first_align - 1; % convert to fields
            elseif (strcmpi(clipset(orig_offset).video_standard, 'interlace_upper_field_first') == 1)
                is_interlaced = 1;
                field_first = 2;
                first_align = 2*first_align - 1; % convert to fields
            else 
                is_interlaced = 0; % progressive
            end

            if is_interlaced
                if ~verbose
                    [results results_rmse results_fuzzy results_fuzzy_mse] = est_var_frame_delays( ...
                        y_proc, y_orig, 'normalize', 'causal', 'first_align', first_align, ...
                        'interlaced', field_first, 'reframe', 't_uncert', t_uncert);
                else
                    [results results_rmse results_fuzzy results_fuzzy_mse] = est_var_frame_delays(...
                        y_proc, y_orig, 'normalize', 'causal', 'first_align', first_align, ...
                        'interlaced', field_first, 'reframe', 'verbose', 't_uncert', t_uncert);
                end
            else
                if ~verbose
                    [results results_rmse results_fuzzy results_fuzzy_mse] = est_var_frame_delays(...
                        y_proc, y_orig, 'normalize', 'causal', 'first_align', first_align, ...
                        't_uncert', t_uncert);
                else
                    [results results_rmse results_fuzzy results_fuzzy_mse] = est_var_frame_delays(...
                        y_proc, y_orig, 'normalize', 'causal', 'first_align', first_align, 'verbose', ...
                        't_uncert', t_uncert);
                end
            end

            write_feature_standard(results, feature_base_dir, 'vfd_results_psnr', clipset(curr_offsets(j)), 'vfd', results_fuzzy);
        end
           

        %  If the VFD alignment algorithm failed, then results==0.
        %  In that case, use the time alignment given by the psnr_file and
        %  generate psuedo results where the alignment just increases by
        %  one frame (or field) at a time.
        [nrows, ncols, nframes] = size(y_proc);
        if (results == 0)
            vfd_failed = 1;  % Set a logical variable to record that the VFD algorithm failed.
            if (~strcmpi(scan_type,'progressive'))  % interlaced
                npts = nframes*2;
            else  % progressive
                npts = nframes;
            end
            results = first_align:first_align+npts-1;
            if quiet == 0
                fprintf('WARNING: VFD algorithm failed for clip %s_%s_%s, using psnr_file time alignment.\n', ...
                    test, this_scene, this_hrc);
            end
        else
            vfd_failed = 0;
            npts = length(results);
        end

        % This code translates the VFD results to use the orig and proc FILE indexing
        if (strcmpi(scan_type,'progressive'))
            proc_indices = (clipset(curr_offsets(j)).align_start-1) + (1:length(results));
            orig_indices = (orig_start-1) + results;
        else % interlaced
            proc_indices = 2*(clipset(curr_offsets(j)).align_start-1) + (1:length(results));
            orig_indices = 2*(orig_start-1) + results;
        end
        results_vfd(index).orig_indices = orig_indices;
        results_vfd(index).proc_indices = proc_indices;
        
        % Write out the psnr_vfd results for this clip into the vfd_file
        do_reframing = is_reframing_indicated(clipset(curr_offsets(j)));
        fid_vfd = fopen(vfd_file,'a');
        fprintf(fid_vfd,'%s, %s, %s,', results_vfd(index).test, results_vfd(index).scene, ...
            results_vfd(index).hrc );
        for k = 1:npts
            fprintf(fid_vfd,'%d %d,',proc_indices(k)+do_reframing, orig_indices(k));
        end
        fprintf(fid_vfd,'\n');
        fclose(fid_vfd);

        % Apply the VFD correction to the original clip to make it look
        % like the processed clip.
        [y_orig, first_valid, last_valid] = vfd_match(y_orig, scan_type, results);
        
        % discard invalid frames from beginning and end.  Update variables.
        y_orig = y_orig(:,:,first_valid:last_valid);
        y_proc = y_proc(:,:,first_valid:last_valid);
        nframes = last_valid - first_valid + 1;
        if (~strcmpi(scan_type,'progressive'))  % interlaced
            npts = nframes*2;
            results = results((2*first_valid-1):(2*last_valid)); 
            results_fuzzy = results_fuzzy(:,(2*first_valid-1):(2*last_valid));
        else  % progressive
            npts = nframes;
            results = results(first_valid:last_valid); 
            results_fuzzy = results_fuzzy(:,first_valid:last_valid);
        end
        
        % print valid range to file
        vfd_text = sprintf('%s%s, %s, %s, ', vfd_text, clipset(curr_offsets(j)).test{1}, ...
            clipset(curr_offsets(j)).scene{1}, clipset(curr_offsets(j)).hrc{1});

        align_start = clipset(curr_offsets(j)).align_start;
        if ~strcmpi(scan_type,'progressive'),
            vfd_text = sprintf('%s PSNR-VFD model will use processed fields [%d..%d]\n',...
                vfd_text, ...
                (first_valid*2-1) + (align_start-1)*2 + do_reframing, ...
                last_valid*2 + (align_start-1)*2 + do_reframing);
        else
            vfd_text = sprintf('%s PSNR-VFD model will use processed frames [%d..%d]\n',...
                vfd_text, ...
                first_valid + (align_start-1), ...
                last_valid + (align_start-1));
        end


        % Reshape for gain/offset fit and PSNR calculation
        y_proc = reshape(y_proc,nrows*ncols*nframes,1);
        y_orig = reshape(y_orig,nrows*ncols*nframes,1);

        % Perform the final gain and offset fit using randomly sub-sampled pixels
        rand_nums = round(randperm((nrows*ncols*nframes))); %Randomizes numbers from 1 to nrows*ncols*nframes
        this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nframes*fraction_sampled))),...
                                y_orig(rand_nums(1:round(nrows*ncols*nframes*fraction_sampled))),1);

        clear rand_nums;
        results_vfd(index).gain_adjust = this_fit(1);
        results_vfd(index).offset_adjust = this_fit(2);
        
        clipset(curr_offsets(j)).luminance_gain = clipset(curr_offsets(j)).luminance_gain/this_fit(1);
        clipset(curr_offsets(j)).luminance_offset = clipset(curr_offsets(j)).luminance_offset - this_fit(2);

        %  Calculate the final PSNR_VFD
        this_psnr_vfd = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc+this_fit(2))-y_orig).^2)/(nrows*ncols*nframes)));
        results_vfd(index).psnr_vfd = this_psnr_vfd;
        clear y_orig;  % Done with y_orig        
        
        %  Reshape y_proc to 3D and calculate TI_RMS
        y_proc = reshape(y_proc,nrows,ncols,nframes);
        y_proc = cat(3,y_proc(:,:,1),y_proc);  % Add extra frame at the beginning for TI calculation
        y_proc = diff(y_proc,1,3);  % first order difference along 3rd dimension
        y_proc = reshape(y_proc,nrows*ncols,nframes);
        y_proc = y_proc.^2;
        y_proc = sqrt(sum(y_proc)./(nrows*ncols));
        if (~strcmpi(scan_type,'progressive'))  % Replicate every other sample for interlaced video
            y_proc = reshape(repmat(y_proc,2,1), 1, 2*nframes);
        end        
        
        %  Calculate the diff of the VFD information, which forms the basis
        %  for both par1_vfd and par2_vfd.  The VFD information is
        %  converted to a vector that gives Abnormal Frame Jumps (AFJs).
        % Subtracting 1 and maxing with zero produces a parameter that (1)
        % does not penalize for normal field/frame delivery (where the VFD
        % field/frame indices increase by one from one field/frame to the
        % next), (2) does not penalize for frame/field repeats (e.g., where
        % the VFD frame indices stay fixed from one frame to the next), and
        % (3) does not penalize for interlaced frame repeats (where the VFD
        % field indices jump back one in time from one field to the next).
        % A non-impairment value of 0 is used for the first field/frame
        % (which must be padded since it's diff is not available).  
        if (vfd_failed)
            
            % This equation assumes that both the early and late time sides
            % used for the frame jump estimates are absolutely correct
            % (i.e., no fuzzy alignments).  The fuzzy alignment information
            % is not available here as the VFD algorithm failed.  Abnormal
            % Frame Jumps (AFJ) is then calculated as:
            afj = max([0 abs(diff(results))-1], 0);
            
        else
            
            % This code assumes fuzzy uncertainty on both the early and
            % late time sides when calculating frame jumps (the uncertainty
            % is given by results_fuzzy array).  Here, frame jumps are only
            % penalized when they are absolutely certain to be correct
            % (i.e., no fuzzy overlapping). 
            fuzzy_max_early = min(max(results_fuzzy(:,1:npts-1)),results(2:npts));
            fuzzy_min_late = max(min(results_fuzzy(:,2:npts)),fuzzy_max_early);
            afj = max([0 fuzzy_min_late-fuzzy_max_early-1], 0);
            
        end
        
        %  Calculate the pure VFD parameter (par1_vfd) and the TI weighted
        %  variant VFD parameter (par2_vfd).
        par1_vfd = log10(sqrt(mean(afj.^2))+1);
        results_vfd(index).par1_vfd = par1_vfd;
        par2_vfd = log10(sqrt(mean((afj.*log10(1+y_proc)).^2))+1);
        results_vfd(index).par2_vfd = par2_vfd;
        
        %  Output the clip information, current time, and the psnr_vfd
        t = clock;
        if (verbose)
            fprintf('Clip %s_%s_%s at %d:%d, psnr_vfd = %5.4f, par1_vfd = %5.4f, par2_vfd = %5.4f\n', ...
                clipset(curr_offsets(j)).test{:}, clipset(curr_offsets(j)).scene{:}, clipset(curr_offsets(j)).hrc{:}, t(4), t(5), this_psnr_vfd, par1_vfd, par2_vfd);
        end
        
        % Write out the psnr_vfd results for this clip into the vfd_file
        do_reframing = is_reframing_indicated(clipset(curr_offsets(j)));
        vfd_model_text = sprintf('%s%s, %s, %s, %5.4f, %5.4f, %5.4f\n', vfd_model_text,...
            results_vfd(index).test, results_vfd(index).scene, ...
            results_vfd(index).hrc, results_vfd(index).psnr_vfd, results_vfd(index).par1_vfd, results_vfd(index).par2_vfd);
        
        output.clip_name{curr} = clip_name;
        output.data(curr) = results_vfd(index).psnr_vfd;
        output.inlsa_mos(curr) = clipset(curr_offsets(j)).inlsa_mos;
        output.mos(curr) = clipset(curr_offsets(j)).mos;
        curr = curr + 1;
        
        %  Add to the HRC summer
        psnr_ave = psnr_ave + this_psnr_vfd;
        par1_ave = par1_ave + par1_vfd;
        par2_ave = par2_ave + par2_vfd;
        
        %  Increment the results_vfd counter
        index = index+1;        
    end
    
    % Compute average psnr_vfd for this HRC
    psnr_ave = psnr_ave/(num_clips);
    par1_ave = par1_ave/(num_clips);
    par2_ave = par2_ave/(num_clips);
    
    if(verbose)
        fprintf('\nHRC = %s, psnr_vfd_ave = %5.4f, par1_vfd_ave = %5.4f, par2_vfd_ave = %5.4f\n\n', ...
            clipset(curr_offsets(j)).hrc{:}, psnr_ave, par1_ave, par2_ave);
    end        
    
end
fid_vfd = fopen(vfd_file,'a');
fprintf(fid_vfd, '%s', vfd_text);
fprintf(fid_vfd, '%s', vfd_model_text);
fclose(fid_vfd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yout, first_valid, last_valid] = vfd_match(yin, scan_type, results)
%  function [yout] = vfd_match(yin, scan_type, results)
%  This function converts 3D input video array 'yin' into another 3D video
%  array 'yout' with frames and/or fields ordered according to the Variable
%  Frame Delay (VFD) 'results'.  The scan_type must be either
%  'progressive', 'interlaced_lff', or 'interlaced_uff'.
%

[nrows, ncols, nframes_orig] = size(yin);

if (strcmpi(scan_type,'progressive'))
    nframes = length(results);
    is_interlaced = 0;
elseif (strcmpi(scan_type, 'interlace_lower_field_first'))  
    nframes = length(results)/2;  % These are field results, so must half to get nframes
    is_interlaced = 1;
    field_first = 1;  % Use same field_first definition as est_var_frame_delays
elseif (strcmpi(scan_type, 'interlace_upper_field_first'))
    nframes = length(results)/2;
    is_interlaced = 1;
    field_first = 2;
else
    error('Unsupported scan_type in function vfd_match.');
end

% yout will be the VFD-corrected original and will have the same number of
% frames as the processed clip.  Use the same 'single' or 'double' rule as
% read_tslice for the image precision of yout.
if (nrows > 650)
    yout = zeros(nrows,ncols,nframes,'single');
else
    yout = zeros(nrows,ncols,nframes,'double');
end

if(is_interlaced)
    
    for j = 1:nframes
        % Get matching original field for the early processed field
        orig_frame_num = ceil(results(2*j-1)/2);  % The frame number that contains the original field
        early_field = mod(results(2*j-1),2);  % =1 if early field, =0 if late field
        [yo1 yo2] = split_into_fields(squeeze(yin(:,:,orig_frame_num)));
        if (early_field)
            switch field_first
                case(1)
                    this_orig1 = yo1;
                case(2)
                    this_orig1 = yo2;
            end
        else % late field
            switch field_first
                case(1)
                    this_orig1 = yo2;
                case(2)
                    this_orig1 = yo1;
            end
        end
        % Get matching original field for the late processed field
        orig_frame_num = ceil(results(2*j)/2);  % The frame number that contains the original field
        early_field = mod(results(2*j),2);  % =1 if early field, =0 if late field
        [yo1 yo2] = split_into_fields(squeeze(yin(:,:,orig_frame_num)));
        if (early_field)
            switch field_first
                case(1)
                    this_orig2 = yo1;
                case(2)
                    this_orig2 = yo2;
            end
        else % late field
            switch field_first
                case(1)
                    this_orig2 = yo2;
                case(2)
                    this_orig2 = yo1;
            end
        end
        % Joint the two original fields into a frame
        switch field_first
            case(1)
                this_orig = join_into_frames(this_orig1,this_orig2);
            case(2)
                this_orig = join_into_frames(this_orig2,this_orig1);
        end
        yout(:,:,j) = this_orig;
    end
    clear this_orig this_orig1 this_orig2 yo1 yo2;

    % find the first valid frame (i.e., the processed frame does not align
    % to the first original frame, and so we can be confident that the
    % alignment isn't just smack into the search boundary).
    [~, ~, len_of_yin] = size(yin);
    for first_valid=1:nframes*2,
        if results(first_valid) > 2,
            break;
        end
    end
    first_valid = ceil(first_valid/2);
    % find the last valid frame (i.e., the processed frame does not align
    % to the last original frame, and so we can be confident that the
    % alignment isn't just smack into the search boundary).
    for last_valid=nframes*2:-1:1,
        if results(last_valid) < (len_of_yin*2-1),
            break;
        end
    end
    last_valid = floor(last_valid/2);
    
else  % progressive
    
    for j = 1:nframes
        yout(:,:,j) = yin(:,:,results(j));
    end
    
    % find the first valid frame (i.e., the processed frame does not align
    % to the first original frame, and so we can be confident that the
    % alignment isn't just smack into the search boundary).
    [~, ~, len_of_yin] = size(yin);
    for first_valid=1:nframes,
        if results(first_valid) > 1,
            break;
        end
    end
    % find the last valid frame (i.e., the processed frame does not align
    % to the last original frame, and so we can be confident that the
    % alignment isn't just smack into the search boundary).
    for last_valid=nframes:-1:1,
        if results(last_valid) < len_of_yin,
            break;
        end
    end
    
end

return

        
        


