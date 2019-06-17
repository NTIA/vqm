function results = dll_psnr_search(varargin)
% DLL_PSNR_SEARCH
%   Estimate the Y-channel PSNR (PSNR) of a designated processed clip
%   compared to an original reference clip. This program is intended to be
%   used by cvqm.m and will use dll_video and dll_calib_video calls to
%   obtain the video clip information.  These clips must be in the AVI
%   format and already stored as global variables.
%
%   A peak signal of 255 is used for calculation of PSNR.  Double precision
%   calculations are used everywhere.  A 64-bit operating system with at
%   least 4 GB of free memory is recommended since the entire double
%   precision versions of the original and processed sequences must be held
%   in memory.
%
% SYNTAX
%   results = dll_psnr_search(options)
%
%   The options syntax is based around cell arrays.  Each option is
%   represented by a cell array, each with the name of the option in the
%   first position and the specifications for that option following in the
%   next cells.  For example, calling the 'spacial_uncertainty' option
%   would look like
%                {'spatial_uncertainty',5, 5}
%
%   A full function call may look like
%           output = dll_psnr_search({'spatial_uncertainty', 5, 5}, ...
%               {'temporal_uncertainty', 2});
%
%   The results are output through a structure.  The structure has field
%   names yshift, xshift, tshift, gain, offset, and pnsr.
%
% DESCRIPTION
%   This function will process a video clip stored as a global variable by
%   dll_video and dll_calib_video, estimate the Y-channel PSNR of the clip,
%   and output the PSNR results into a structure.
%
%   The algorithm performs an exhaustive search for max PSNR over plus or
%   minus the spatial_uncertainty (in pixels) and plus or minus the
%   temporal_uncertainty (in frames).  The processed video segment is fixed
%   and the original video segment is shifted over the search range.  For
%   each spatial-temporal shift, a linear fit between the processed pixels
%   and the original pixels is performed such that the mean square error of
%   [original-gain*processed+offset] is minimized (hence maximizing PSNR).
%
%   Any or all of the following optional properties may be requested.
%
%
%   'sroi' top left bottom right    Only use the specified spatial region
%                                   of interest (sroi) for the PSNR
%                                   calculation.  This is the sroi of the
%                                   processed sequence, which remains fixed
%                                   over all spatial shifts.  By default,
%                                   sroi is the entire image reduced by the
%                                   spatial uncertainty.  If the user
%                                   inputs a sroi, allowance must be made
%                                   for the spatial search specified by
%                                   'spatial_uncertainty'.
%
%   'fraction_sampled' p   Specifies the fraction of pixels to be sub
%                        sampled.  This affects only the gain and offset
%                        calculation.  This is a deviation from ITU-T
%                        Recommendation J.340 because it uses a randomly
%                        sub-sampled version of the video sequence.  There
%                        is a performance gain in speed and considerable
%                        less memory is used and the results seem to agree
%                        closely with using the full video sequence.  The
%                        p value ranges from (0:1].  However the lowest
%                        recommended value is .001.  By default, this is
%                        set to 1 for 100 percent sampling.
%
%   'spatial_uncertainty' x y   Specifies the spatial uncertainty (plus
%                               or minus, in pixels) over which to
%                               search.  The processed remains fixed and
%                               the original is shifted.  By default,
%                               this is set to zero.  These x and y shifts
%                               are always with respect to the video FRAME.
%
%   'temporal_uncertainty' t  Specifies the temporal uncertainty
%                             (plus or minus, in frames) over which
%                             to search.  The processed remains fixed
%                             and the original is shifted.  By
%                             default, this is set to zero.
%
%
% RESULTS
%   Each field in the results structure contains the following information for
%   the processed clip:
%
%   yshift  The y shift for maximum PSNR.  This is in frames.
%   xshift  The x shift for maximum PSNR.  This is in frames.
%   tshift  The time shift for maximum PSNR.  This is in frames, however
%   when a .5 is present, for interlace, this means that there is a field
%   shift and not a whole frame shift.
%   gain    The gain*processed+offset for maximum PSNR.
%   offset
%   psnr    The maximum PSNR observed over the search range.
%
% EXAMPLE
%   This example illustrates the format for the function call.
%
%   dll_psnr_search({'sroi', 5, 5, 140, 172}, {'spatial_uncertainty', 1, 1}, {'temporal_uncertainty', 8});


% Define the peak signal level
peak = 255.0;

% Validate input arguments and set their defaults
x_uncert = 0;
y_uncert = 0;
t_uncert = 0;
is_whole_image = 1;
fraction_sampled = 1;
video_standard = dll_video('get_video_standard', 1);
cnt=1;

% Set the rewind point for dll_video to go back to at the end of the
% function.
dll_video('set_rewind', 1);
dll_video('set_rewind', 2);

% Change values based on options selected by user
while cnt <= length(varargin)
    if strcmpi(varargin{cnt}{1}, 'fraction_sampled') == 1
        if length(varargin{cnt}) ~= 2
            error('Too many or too few arguments for dll_psnr_search option ''fraction_sampled''');
        end
        fraction_sampled = varargin{cnt}{2};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt}{1}, 'spatial_uncertainty') == 1
        if length(varargin{cnt}) ~= 3
            error('Too many or too few arguments for dll_psnr_search option ''spatial_uncertainty''');
        end
        x_uncert = varargin{cnt}{2};
        y_uncert = varargin{cnt}{3};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt}{1}, 'temporal_uncertainty') == 1
        if length(varargin{cnt}) ~= 2
            error('Too many or too few arguments for dll_psnr_search option ''temporal_uncertainty''');
        end
        t_uncert = varargin{cnt}{2};
        cnt = cnt + 1;
    else
        error('Option name passed into dll_psnr_search not recognized.');
    end
end

% Validate 'fraction_sampled' value
if (fraction_sampled <= 0 || fraction_sampled > 1)
    error('The value of fraction_sampled is invalid.');
end
if (fraction_sampled < .001)
    warning('Fraction_sampled is below the recommended limit.');
end

% If not progressive and user inputs an SROI, they must have an odd top and
% an even bottom.  Otherwise the field ordering will reverse.
if (~strcmpi(video_standard,'progressive') && ~is_whole_image && (~mod(top,2) || mod(bottom,2)))
    error('SROI top must be odd and bottom must be even for interlaced video.');
end

% For interlaced video, multiple the y_uncert by two since it will later be
% divided by two (i.e., field line shifts are searched rather than frame
% lines).
if(~strcmpi(video_standard,'progressive'))
    y_uncert = 2*y_uncert;
end


%Results struct to store results, shifts are how much the original must be
%shifted with respect to the processed
results = struct('test', [], 'scene', [], 'hrc', [], 'yshift', [], ...
    'xshift', [], 'tshift', [], 'gain', [], 'offset', [], 'psnr', []);


% Read original and processed video files
[fps] = dll_video('fps', 1);
[fps_proc] = dll_video('fps', 2);

tframes = dll_video('total_frames', 1);
tframes_proc = dll_video('total_frames', 2);

% tframes = dur*fps;  % total frames in orig file 
% tframes_proc = dur_proc*fps_proc;
% Validate that orig and proc have the same number of frames
if (tframes ~= tframes_proc)
    tframes = min(tframes,tframes_proc);
end
% Set/Validate the time segment to use
% use whole time segment less uncertainty
fstart= 1+t_uncert;
fstop = tframes-t_uncert;

if(strcmpi(video_standard,'progressive')) %Reads in normal amount of frames
    % Read in video and clear color planes to free up memory
    dll_video('set_tslice', 1, ((fstop+t_uncert)-(fstart-t_uncert)+1)/fps); 
    y_orig = dll_calib_video('tslice', 1, 0);
    dll_video('discard', 2, (fstart-1)/fps_proc);
    dll_video('set_tslice', 2, (fstop - fstart + 1)/fps_proc);
    y_proc = dll_calib_video('tslice', 2);
else %Reads in one less processed frame for interlace
    % Read in video and clear color planes to free up memory
    dll_video('set_tslice', 1, ((fstop+t_uncert)-(fstart-t_uncert)+1)/fps); 
    y_orig = dll_calib_video('tslice', 1, 0);
    dll_video('discard', 2, (fstart-1)/fps_proc);
    dll_video('set_tslice', 2, (fstop - fstart)/fps_proc);
    y_proc = dll_calib_video('tslice', 2);
end

y_orig = double(y_orig);
y_proc = double(y_proc);

% Crop the processed video based on uncertainty
y_proc = y_proc((1+y_uncert):(end-y_uncert), (1+x_uncert):(end-x_uncert), :);
% If the video standard is interlaced, you must crop off an extra row off
% the top and bottom of each array.
if(~strcmpi(video_standard, 'progressive'))
    y_orig = y_orig(:, (1+x_uncert):(end-x_uncert),:);
    y_proc = y_proc(:, (1+x_uncert):(end-x_uncert),:);
end
[nrows, ncols, nsamps] = size(y_proc);

% Compute PSNR for each spatial-temporal shift
best_psnr = -inf;
best_xshift = 0;
best_yshift = 0;
best_tshift = 0;
best_gain = 1;
best_offset = 0;


if(fraction_sampled ~= 1)  %Random Sampling
    rand_nums = round(randperm((nrows*ncols*nsamps))); %Randomizes numbers from 1 to nrows*ncols*nsamps
end

% Interlace lower field first search
if (strcmpi(video_standard, 'interlace_lower_field_first'))
    
    %split_into_fields -- y_orig
    [row, col, time_orig] = size(y_orig);
    y_orig = reshape(y_orig,2,row/2,col,time_orig);
    late_field_orig_master = squeeze(y_orig(1,:,:,:));
    early_field_orig_master = squeeze(y_orig(2,:,:,:));
    clear y_orig;
    
    % Reshape y_proc for the PSNR calculation
    y_proc = reshape(y_proc,nrows*ncols*nsamps,1);
    
    y_uncert_field = (y_uncert/2);
    for curr_time = -t_uncert:t_uncert
        for curr_h = -x_uncert:x_uncert
            for curr_v = -y_uncert_field:y_uncert_field
                
                if(mod(curr_v,2)) % Odd - reframing
                    late_field_v = ceil(curr_v/2);
                    early_field_v = floor(curr_v/2);
                    
                    % This code assumes lower field first, so first
                    % lower field is discarded & need one extra
                    % field to make up for this at the end.
                    % no change to length in time.
                    
                    lower_field_orig = late_field_orig_master(...
                        (y_uncert_field + 1 + late_field_v):(y_uncert_field + (nrows/2) + late_field_v),...
                        (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                        (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                    
                    upper_field_orig = early_field_orig_master(...
                        (y_uncert_field + 1 + early_field_v):(y_uncert_field + (nrows/2) + early_field_v),...
                        (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                        (t_uncert + 2 + curr_time):(t_uncert + nsamps + curr_time + 1));
                    
                else % Even - no reframing
                    lower_field_orig = early_field_orig_master((y_uncert_field + 1 + curr_v/2):(y_uncert_field + (nrows/2) + curr_v/2),...
                        (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                        (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                    
                    upper_field_orig = late_field_orig_master(...
                        (y_uncert_field + 1 + curr_v/2):(y_uncert_field + (nrows/2) + curr_v/2),...
                        (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                        (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                    
                end
                
                % Join into frames and reshape into column vector
                [row1, col1, time1] = size(lower_field_orig);
                y_orig(1,:,:,:) = upper_field_orig;
                y_orig(2,:,:,:) = lower_field_orig;
                clear lower_field_orig upper_field_orig;
                y_orig = reshape(y_orig, row1*2, col1, time1);  % 3D
                y_orig = reshape(y_orig, nrows*ncols*nsamps, 1);  % Column Vector
                
                % Perform gain and level offset calculation
                if(fraction_sampled == 1)
                    this_fit = polyfit(y_proc,y_orig,1);
                else
                    this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),...
                        y_orig(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),1);
                end
                
                % Calculate the PSNR
                this_psnr = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc + this_fit(2))-...
                    y_orig).^2)/(nrows*ncols*nsamps)));
                clear y_orig;
                
                if(this_psnr > best_psnr)
                    best_psnr = this_psnr;
                    best_yshift = curr_v;
                    best_xshift = curr_h;
                    %  Correct curr_time for reframed video
                    if(mod(best_yshift,2))
                        best_tshift = curr_time + 0.5;
                    else
                        best_tshift = curr_time;
                    end
                    best_gain = this_fit(1);
                    best_offset = this_fit(2);
                end
                
            end
        end
    end
    
elseif (strcmpi(video_standard, 'interlace_upper_field_first'))  % Interlace upper field first search
    
    %split_into_fields -- y_orig
    [row, col, time] = size(y_orig);
    y_orig = reshape(y_orig,2,row/2,col,time);
    early_field_orig_master = squeeze(y_orig(1,:,:,:));
    late_field_orig_master = squeeze(y_orig(2,:,:,:));
    clear y_orig;
    
    % Reshape y_proc for the PSNR calculation
    y_proc = reshape(y_proc,nrows*ncols*nsamps,1);
    
    y_uncert_field = (y_uncert/2);
    for curr_time = -t_uncert:t_uncert
        for curr_h = -x_uncert:x_uncert
            for curr_v = -y_uncert_field:y_uncert_field
                
                if(mod(curr_v,2)) % Odd - reframing
                    early_field_v = ceil(curr_v/2);
                    late_field_v = floor(curr_v/2);
                    
                    % This code assumes upper field first, so first
                    % upper field is discarded & need one extra
                    % field to make up for this at the end.
                    % No change to length in time.
                    
                    lower_field_orig = early_field_orig_master(...
                        (y_uncert_field + 1 + early_field_v):(y_uncert_field + (nrows/2) + early_field_v),...
                        (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                        (t_uncert + 2 + curr_time):(t_uncert + nsamps + curr_time + 1));
                    
                    upper_field_orig = late_field_orig_master(...
                        (y_uncert_field + 1 + late_field_v):(y_uncert_field + (nrows/2) + late_field_v),...
                        (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                        (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                    
                else % Even - no reframing
                    lower_field_orig = late_field_orig_master(...
                        (y_uncert_field + 1 + curr_v/2):(y_uncert_field + (nrows/2) + curr_v/2),...
                        (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                        (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                    
                    upper_field_orig = early_field_orig_master((y_uncert_field + 1 + curr_v/2):(y_uncert_field + (nrows/2) + curr_v/2),...
                        (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                        (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                    
                end
                
                % Join into frames and reshape into column vector
                [row1, col1, time1] = size(upper_field_orig);
                y_orig(1,:,:,:) = upper_field_orig;
                y_orig(2,:,:,:) = lower_field_orig;
                clear lower_field_orig upper_field_orig;
                y_orig = reshape(y_orig, row1*2, col1, time1);  % 3D
                y_orig = reshape(y_orig, nrows*ncols*nsamps, 1);  % Column Vector
                
                % Perform gain and level offset calculation
                if(fraction_sampled == 1)
                    this_fit = polyfit(y_proc,y_orig,1);
                else
                    this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),...
                        y_orig(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),1);
                end
                
                % Calculate the PSNR
                this_psnr = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc + this_fit(2))-...
                    y_orig).^2)/(nrows*ncols*nsamps)));
                clear y_orig;
                
                if(this_psnr > best_psnr)
                    best_psnr = this_psnr;
                    best_yshift = curr_v;
                    best_xshift = curr_h;
                    %  Correct curr_time for reframed video
                    if(mod(best_yshift,2))
                        best_tshift = curr_time + 0.5;
                    else
                        best_tshift = curr_time;
                    end
                    best_gain = this_fit(1);
                    best_offset = this_fit(2);
                end
                
            end
        end
    end
    
else % progressive or interlace with no reframing
    y_proc = reshape(y_proc,nrows*ncols*nsamps,1);
    
    for curr_time = -t_uncert:t_uncert
        for curr_h = -x_uncert:x_uncert
            for curr_v = -y_uncert:y_uncert
                
                if(fraction_sampled == 1)
                    this_fit = polyfit(y_proc,reshape(y_orig(1+curr_v+y_uncert:curr_v+y_uncert+nrows,...
                        1+curr_h+x_uncert:curr_h+x_uncert+ncols,1+curr_time+t_uncert:curr_time+t_uncert+nsamps),...
                        nrows*ncols*nsamps,1),1);
                else
                    %Reshapes y_orig
                    reshape_y_orig = reshape(y_orig(1+curr_v+y_uncert:curr_v+y_uncert+nrows, 1+curr_h+x_uncert:curr_h+x_uncert+ncols,...
                        1+curr_time+t_uncert:curr_time+t_uncert+nsamps), nrows*ncols*nsamps, 1);
                    
                    % Perform gain and level offset calculation
                    this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),...
                        reshape_y_orig(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),1);
                end
                
                % Calculate the PSNR
                if(fraction_sampled == 1)
                    this_psnr = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc + this_fit(2))-...
                        reshape(y_orig(1+curr_v+y_uncert:curr_v+y_uncert+nrows,1+curr_h+x_uncert:curr_h+x_uncert+ncols,...
                        1+curr_time+t_uncert:curr_time+t_uncert+nsamps),nrows*ncols*nsamps,1)).^2)/(nrows*ncols*nsamps)));
                else
                    this_psnr = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc + this_fit(2))-...
                        reshape_y_orig).^2)/(nrows*ncols*nsamps)));
                    clear reshape_y_orig;
                end
                
                if(this_psnr > best_psnr)
                    best_psnr = this_psnr;
                    best_yshift = curr_v;
                    best_xshift = curr_h;
                    best_tshift = curr_time;
                    best_gain = this_fit(1);
                    best_offset = this_fit(2);
                end
            end
        end
    end
end

% Insert results into a structure to be returned.
results.yshift = best_yshift;
results.xshift = best_xshift;
results.tshift = best_tshift;
results.gain = best_gain;
results.offset = best_offset;
results.psnr = best_psnr;

dll_video('rewind', 1);
dll_video('rewind', 2);


    
