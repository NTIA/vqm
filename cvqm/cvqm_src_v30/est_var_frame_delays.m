function [results, results_rmse, results_fuzzy, results_fuzzy_mse] = est_var_frame_delays(proc, orig, varargin)
% EST_VAR_FRAME_DELAYS
%   Estimate the variable delays of each frame in a processed video clip
%   (i.e., proc) given an original video clip (i.e., orig).  The processed
%   and original clips are three dimensional (rows x cols x frames) Y (luma) 
%   matrices that may have a different number of frames.  The results row 
%   vector gives the best matching frame number in the original clip for each
%   frame in the processed clip, where original frame indices start from 1.
%   The user may optionally specify interlaced format, and then the results
%   vector gives the best matching field number (i.e., the results vector
%   will be twice the frame-length of the processed clip).  Note that this
%   routine does not perform spatial registration and/or gain and level
%   offset so the orig and proc clips are assumed to have been fully
%   calibrated.
%
% SYNTAX
%   [results, results_rmse, results_fuzzy, results_fuzzy_mse] = est_var_frame_delays (proc, orig, options)
%
% DESCRIPTION
%   For each frame in the processed clip, this algorithm finds the frame in
%   the original clip that minimizes the mean squared error, subject to the
%   constraints imposed by the optional inputs.  The following output
%   arguments can be requested:
%
%   results       A row vector of the same length as the processed clip in
%                 frames (or fields) that gives the best matching original
%                 frame (or field) for each processed frame (or field).
%                 The original frame (or field) indices are assumed to
%                 start at one.  If the 'causal' option is NOT specified, this
%                 is merely the original frame (or field) with the smallest
%                 Mean Squared Error (MSE) when compared to the processed
%                 frame.  If the 'causal' option is specified, then the
%                 results row vector contains the results from the causal
%                 filtering algorithm (see below).
%
%   results_rmse  If two output arguments are requested, the second one
%                 will contain the Root Mean Squared Error (RMSE) between
%                 the causal alignment and the unfiltered and possibly
%                 non-causal alignment (in frames or fields).  The higher
%                 the value, the more difference there is between the
%                 the two algorithms, which probably indicates that the
%                 scene is difficult to align (e.g., a scene with a small
%                 amount of motion or repetitive motion), or that not
%                 enough temporal uncertainty (t_uncert) was used.
%
%   results_fuzzy If three output arguments are requested, the third one
%                 will contain the fuzzy alignment results.  This is a
%                 matrix where each column gives a set of rank sorted fuzzy
%                 original frame (or field) alignments, sorted from most
%                 likely to least likely.  The first element of each column
%                 vector gives the most likely alignment for that processed
%                 frame (field) and is the same as the first output
%                 argument (results).  The number of rows in each column is
%                 equal to how many frames (or fields) were searched, which
%                 depends upon t_uncert.  But only likely alignments are
%                 included and the remainder of the rows are filled in by
%                 'NaN' (Not-a-Number).  The fuzzy matrix is influenced by
%                 the 'causal' option since the range of possible
%                 alignments may be expanded to include the causal
%                 alignments.
%
%   results_fuzzy_mse If four output arguments are requested, the fourth one
%                     will contain the Mean Squared Error of each fuzzy
%                     frame (or field) alignment, so this matrix is the
%                     same size as results_fuzzy.
%
% OPTIONS
%   Any or all of the following optional inputs may be requested to control
%   the behavior of the algorithm.
%
%   'interlaced',first_field  Specifies interlaced format, where first_field
%                             specifies which field is first in time, 
%                             first_field = 1 if the lower field is first,
%                             first_field = 2 if the upper field is first.
%                             The orig and proc clips must have an even
%                             number of rows for this option.
%
%   'sroi',top,left,bottom,right    Only use the specified spatial region 
%                                   of interest (sroi) for the frame delay
%                                   calculation.  sroi is referenced to the
%                                   original frame, where the (top, left)
%                                   corner of the image is (1, 1).  For
%                                   interlaced video, top must be odd and
%                                   bottom must be even. By default, sroi
%                                   is the entire image.
%
%   'first_align',a           Specifies the best guess for the original
%                             frame (or field) number that corresponds to
%                             the first frame (or field) in the processed
%                             clip.  Set to 1 by default.  Note that for
%                             interlaced video, this value must be in
%                             fields, not frames!
%
%   't_uncert',t              Specifies the temporal uncertainty (always
%                             plus or minus t frames, NOT fields) over which 
%                             to search.  The processed remains fixed and
%                             the original is shifted.  The center (zero 
%                             shift) point for the first frame (or field)
%                             is given by first_align. By default,
%                             temporal uncertainty is set to 30.  When the
%                             original cannot be shifted by the
%                             temporal uncertainty (e.g., near the ends of
%                             the sequence), the original will be shifted
%                             by the maximum possible.
%
%   'normalize'     Perform a normalization on the original and processed
%                   clips such that the processed and original clips have
%                   zero mean, unit variance.  All frames in the processed
%                   clip are used to estimate its mean and variance.  For
%                   the original clip, the first_align point is used to
%                   select an equal number of frames from the original clip
%                   upon which to base its mean and variance (provided
%                   these corresponding frames are available, if not, the
%                   number of frames for the original's mean and variance
%                   estimate will be reduced).  By default, no
%                   normalization is performed.  However, normalization is
%                   a recommended option unless the processed video is
%                   perfectly calibrated wrt gain and level offset.
%
%   'causal'   Impose causality constraint so that later frames (fields) in
%              the processed clip cannot align to original frames (fields)
%              that are earlier in time than found for the proceeding
%              processed frames (fields).  For interlaced video, a
%              one-field jump back in time is allowed since this is
%              indicative of a frozen frame.  By default, causality is
%              turned off (yes, codecs can output non-causal sequences).
%              But specifying the causal option is usually recommended. 
%
%   'reframe'  Allow for the possibility that the processed video clip has
%              been reframing.  Must also specify the 'interlaced' option.
%              Reframing can vary throughout the processed clip, although
%              this should be rare.  This option will increase the
%              runtime substantially since extra spatial shifts must be
%              examined.
%
%   'verbose'   Display output during processing.  verbose mode is turned
%               off by default.
%
% EXAMPLES
%   

% Return values of function if failure
results = 0;  
results_rmse = 0;
results_fuzzy = 0;
results_fuzzy_mse = 0;

% Set input arguments to their defaults and assign optional inputs
interlaced = 0;
is_whole_image = 1;
first_align = 1;
t_uncert = 30;
normalize = 0;
is_causal = 0;
reframe = 0;
verbose = 0;
% Assign optional inputs
cnt=1;
while cnt <= length(varargin),
    if strcmpi(varargin(cnt),'interlaced') == 1
        interlaced = 1;
        first_field = varargin{cnt+1};
        if (first_field~=1 && first_field~=2)
            error('first_field must be 1 or 2.');
        end
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'sroi') == 1
        is_whole_image = 0;
        top = varargin{cnt+1};
        left = varargin{cnt+2};
        bottom = varargin{cnt+3};
        right = varargin{cnt+4};
        cnt = cnt + 5;
    elseif strcmpi(varargin(cnt),'first_align') == 1
        first_align = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt), 't_uncert') == 1
        t_uncert = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt), 'normalize') == 1
        normalize = 1;
        cnt = cnt + 1;
    elseif strcmpi(varargin(cnt), 'causal') == 1
        is_causal = 1;
        cnt = cnt + 1;
    elseif strcmpi(varargin(cnt), 'reframe') == 1
        reframe = 1;
        cnt = cnt +1;
    elseif strcmpi(varargin(cnt),'verbose') == 1
        verbose = 1;
        cnt = cnt +1;
    else
        error('Property value passed into est_frame_delays not recognized');
    end
end

% Temporal uncertainty search bound check
if(t_uncert < 1)
    error('t_uncert must be at least 1 frame.');
end

%  Find image resolution and number of frames in orig and proc clips
[nr, nc, nf] = size(proc);
[nr_o, nc_o, nf_o] = size(orig);
if (nr ~= nr_o || nc ~= nc_o)
    error('Orig and proc clips must have the same number of rows and cols.');
end

% Validate the SROI
if (is_whole_image) % make ROI whole image
    top = 1;
    left = 1;
    bottom = nr;
    right = nc;
elseif (top<1 || left<1 || bottom>nr || right>nc || top>bottom || left>right)
    error('Requested SROI incompatible with image size.');
end

% Additional checks on SROI for interlaced video
if (interlaced && (mod(nr,2)~=0))
    error('Number of rows must be even.')
end
if (interlaced && (mod(top+1,2)~=0 || (mod(bottom,2)~=0)))
    error('Requested SROI invalid for interlaced video.');
end
if (reframe && ~interlaced)
    error('For the reframe option, you must also specify the interlaced option.')
end

% Split into fields from the get-go if this is interlaced video
if (interlaced)
    
    % This is the normal split (not reframing)
    if (first_field == 2)  % upper field is first
        origf = reshape(orig(top:bottom,left:right,:), 2, (bottom-top+1)/2, right-left+1, nf_o);
        orig = reshape(permute(origf, [2 3 1 4]), (bottom-top+1)/2, right-left+1, 2*nf_o);
        clear origf;
        procf = reshape(proc(top:bottom,left:right,:), 2, (bottom-top+1)/2, right-left+1, nf);
        proc = reshape(permute(procf, [2 3 1 4]), (bottom-top+1)/2, right-left+1, 2*nf);
        clear procf;
    else  % lower field is first
        origf = flipdim(reshape(orig(top:bottom,left:right,:), 2, (bottom-top+1)/2, right-left+1, nf_o),1);
        orig = reshape(permute(origf, [2 3 1 4]), (bottom-top+1)/2, right-left+1, 2*nf_o);
        clear origf;
        procf = flipdim(reshape(proc(top:bottom,left:right,:), 2, (bottom-top+1)/2, right-left+1, nf),1);
        proc = reshape(permute(procf, [2 3 1 4]), (bottom-top+1)/2, right-left+1, 2*nf);
        clear procf;
    end
    
    % Adjust t_uncert and number of fields in orig and proc
    t_uncert = 2*t_uncert;
    nf = 2*nf;
    nf_o = 2*nf_o;
    
    % Number of samples used for mse calculation, field based
    nrows = (bottom-top+1)/2;
    ncols = right-left+1;
    nsamps = nrows*ncols;
    
else  % Progressive sequence
    
    % Window out desired SROI
    orig = orig(top:bottom, left:right, :);
    proc = proc(top:bottom, left:right, :);
    
    % Number of samples used for mse calculation
    nrows = bottom-top+1;
    ncols = right-left+1;
    nsamps = nrows*ncols;
    
end

% Check on alignment of first frame (or field) before we start searching
if (first_align<1 || first_align>nf_o)
    error('Requested first_align is not valid for orig clip.');
end

%%%%%%%%%%
% Stage 1
% Normalize video sequences and compute the Mean Squared Error (MSE)
% between each processed frame (or field) and the set of original frames
% (or fields) within the temporal search window.
%%%%%%%%%%

% Perform the normalization on the orig and proc clips (if requested)
if (normalize)
    
    % Compute mean and stdev of the proc clip and normalize
    proc_mean = sum(reshape(proc, nsamps*nf, 1))/(nsamps*nf);
    proc_std = sqrt(sum(reshape(proc, nsamps*nf, 1).^2)/(nsamps*nf) - proc_mean^2);
    proc = (proc-proc_mean)/proc_std;
    
    % Compute the mean and stdev of the orig clip and normalize
    nf_o1 = first_align;  % The first frame (or field) in the orig to use
    nf_o2 = min(first_align + nf -1, nf_o);  % The final frame (or field) in the orig to use
    orig_mean = sum(reshape(orig(:,:,nf_o1:nf_o2), nsamps*(nf_o2-nf_o1+1), 1))/(nsamps*(nf_o2-nf_o1+1));
    orig_std = sqrt(sum(reshape(orig(:,:,nf_o1:nf_o2), nsamps*(nf_o2-nf_o1+1), 1).^2)/(nsamps*(nf_o2-nf_o1+1)) - orig_mean^2);
    orig = (orig-orig_mean)/orig_std;
    
end

%  Find the best matching orig frame (field) for every proc frame (Field)
if (verbose)
    fprintf('proc frame/field   orig frame/field   mse\n');
    fprintf('----------------   ----------------   ---\n');
end
max_size = 2*t_uncert+1;  % the max size of the alignment results
mse = NaN(max_size,nf);  % holds the mse alignment results
offsets = []; % holds the orig index offsets
for t = 1:nf  % Loop over all processed frames
    
    this_align = first_align+t-1;  % when t = 1, this is the reference alignment
    neg = max(this_align-t_uncert,1);  % the negative most orig index to search
    pos = min(this_align+t_uncert,nf_o);  % the positive most orig index to search
    this_proc = repmat(squeeze(proc(:,:,t)),[1 1 pos-neg+1]);  % Create replicas of this processed frame

    if(~reframe)  % Progressive or interlaced with no reframing
        
        this_mse = sum(reshape((orig(:,:,neg:pos)-this_proc).^2, nsamps, pos-neg+1), 1) / nsamps;
        [best_mse best_ind] = min(this_mse);
        best_ind = best_ind + neg - 1;

    else  % Must do extra comparisons shifted by 1 line, but all comparisons use the same lines from orig
        
        origt = orig(2:nrows-1,:,neg:pos);
        
        %  First comparison, proc is not shifted
        this_mse1 = sum(reshape((origt-this_proc(2:nrows-1,:,:)).^2, nsamps-2*ncols, pos-neg+1), 1) / (nsamps-2*ncols);
        
        %  Second comparison, proc shifted down by one line wrt orig
        this_mse2 = sum(reshape((origt-this_proc(1:nrows-2,:,:)).^2, nsamps-2*ncols, pos-neg+1), 1) / (nsamps-2*ncols);
        
        %  Third comparison, proc shifted up by one line wrt orig
        this_mse3 = sum(reshape((origt-this_proc(3:nrows,:,:)).^2, nsamps-2*ncols, pos-neg+1), 1) / (nsamps-2*ncols);
        
        %  Combine comparisons and take min overall
        this_mse = min([this_mse1;this_mse2;this_mse3]);
        [best_mse best_ind] = min(this_mse);
        best_ind = best_ind + neg - 1;

    end
  
    % Save the offset indices; the mse row indices plus this offset
    % is the original frame (field) indices
    offsets(t) = neg-1;

    % save mse vector time history for later fuzzy processing 
    mse(1:pos-neg+1,t) = this_mse';
    
    if (verbose)
        
        % Plot the correlation function for this processed frame
        figure(1)
        plot(offsets(t)+(1:pos-neg+1),this_mse,'LineWidth',2);
        hold on
        set(gca,'LineWidth',2)
        set(gca,'FontName','Ariel')
        set(gca,'fontsize',12)
        xlabel('Original Frame or Field');
        ylabel('Mean Squared Error (MSE)');
        this_title = ['Processed Frame or Field ' int2str(t)];
        title(this_title);
        grid on
        hold off
        pause(0.01);
        
        % Print out the best matching point
        fprintf('     %4i               %4i          %5.4e\n', t, best_ind, best_mse);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Process the correlation results to compute the four output arguments.
%  WARNING:
%  This algorithm is a highly complicated heuristic multi-stage algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[range num_frames] = size(mse);

% Generate matrix that provides the matching original frame/field for
% every MSE value in the matrix
orig_index = repmat((1:range)',1,num_frames) + repmat(offsets,range,1);

% Sort the MSE values and their orig indices for each proc frame/field
[mse_sort index_sort] = sort(mse);
orig_index_sort = [];
for j = 1:num_frames
    orig_index_sort(:,j) = orig_index(index_sort(:,j),j);
end

% Define a Boolean array that is 1 when the best alignment has hit the
% edge of the search range and zero otherwise (i.e., not enough
% temporal uncertainty or perhaps the correlation function has gone
% berserk).  This array will be used as a filter later to de-weight
% these points in the causal alignment estimation.  When the
% best aligned frame (first element of index_sort) is equal to 1 or
% range, you have hit the edge.
edge = zeros(1, num_frames);
edge(index_sort(1,:) == 1) = 1;  % left search edge
edge(index_sort(1,:) == range) = 1;  % right search edge

%%%%%%%%%%
% Stage 2
% Step thru the frames for this clip and find an initial set of fuzzy alignments
%%%%%%%%%%
final_fuzzy_index = NaN(range, num_frames); % The final fuzzy alignment from most to least likely at end of stage 1
final_fuzzy_mse = NaN(range, num_frames); % Their corresponding MSEs
for j = 1:num_frames
    
    % Since subtractive correlation is being used, the correlation function
    % will be a minimum at the best match and increase from there.  This
    % code uses a thresholding scheme that sets a threshold at the minimum
    % correlation value + thres * maximum correlation value, but the
    % maximum correlation value used is the top (95%) rank sorted value
    % for robustness.  All alignments within this threshold from the best
    % correlation value are considered possible alignments for later
    % processing by the algorithm.
    thres = 0.005;  % Fraction of the maximum correlation level above the minimum for valid fuzzy
    top = 0.95;  % Rank sorted fraction used to determine the max correlation level, for robustness
    this_mse = mse(:,j);
    this_mse_valid = find(~isnan(this_mse));  % array can contain NaN if orig frames did not exist
    range_valid = length(this_mse_valid);
    this_mse_sort = sort(this_mse);
    this_index = orig_index(:,j);
    this_index_sort = orig_index_sort(:,j);
    
    max_corr = this_mse_sort(floor(top*range_valid));
    min_corr = this_mse_sort(1);
    fuzzy = find(this_mse >= min_corr & this_mse <= min_corr+thres*max_corr);
    num_fuzzy = length(fuzzy);  % will always be at least one
    fuzzy_sort = [];
    for k = 1:num_fuzzy
        fuzzy_sort(k) = find(this_index_sort == this_index(fuzzy(k)));
    end
    
    % Expand the fuzzy alignments to include the range of original indices
    % that are covered by the minimum MSE points.  This will include extra
    % points that were not in the minimum MSE fuzzy set.
    bigger_fuzzy = (fuzzy(1):fuzzy(length(fuzzy)))';
    bigger_fuzzy_mse = this_mse(bigger_fuzzy);
    bigger_fuzzy_index = this_index(bigger_fuzzy);
    
    % Added fuzzy points that were not included in the original set
    added_fuzzy = setxor(fuzzy,bigger_fuzzy);
    added_fuzzy_mse = this_mse(added_fuzzy);
    added_fuzzy_index = this_index(added_fuzzy);
    num_added_fuzzy = length(added_fuzzy);
    added_fuzzy_sort = [];
    for k = 1:num_added_fuzzy
        added_fuzzy_sort(k) = find(this_index_sort == this_index(added_fuzzy(k)));
    end
    
    % Sort the bigger fuzzy information from the most likely alignment
    % to the least likely alignment, based on MSE.  The array
    % final_fuzzy_index gives the best fuzzy alignments at this stage of
    % the algorithm without causal processing.  So final_fuzzy_index and
    % final_fuzzy_mse will be assigned to results_fuzzy and
    % results_fuzzy_mse if causal output was not requested by the user.
    [bigger_fuzzy_mse_sort mse_order] = sort(bigger_fuzzy_mse);
    bigger_fuzzy_index_sort = bigger_fuzzy_index(mse_order);
    final_fuzzy_index(1:length(bigger_fuzzy_index_sort),j) = bigger_fuzzy_index_sort;
    final_fuzzy_mse(1:length(bigger_fuzzy_index_sort),j) = bigger_fuzzy_mse_sort;
    
end

%  The array best holds the best alignments based on min MSE
best = final_fuzzy_index(1,:);

%%%%%%%%%%
% Stage 3
% Find normal causal segments
%%%%%%%%%%
% Find normal causal segments using only the best alignment point.  A
% normal causal segment is defined as follows:
%
% For interlaced video: A segment where field n jumps forward in
% alignment by 0 to 2*cjump fields (a field repeating system will cause a
% 2 field jump) with respect to field n-1.  A jump back of 1 field during
% the segment is allowed for the frame repeating case.
%
% For progressive video: A segment where frame n jumps forward in
% alignment by 0 to cjump frames with respect to frame n-1.
%
% A frame alignment that is on the edge will not be included in any
% normal causal segment.
run_beg = [];  % Holds the beginning index of a run
run_length = [];  % Holds the length of the run
ri = 1;  % Index counter for the causal runs
cjump = 2;  % Specifies the jump forward (in frames) that is allowed for 'normal' causal

k = 1;  % begin search at the first frame
while(edge(k))
    k = k+1;
    if (k > num_frames)
        break;
    end
end

if(k <= num_frames)  % Found at least one causal alignment
    run_beg(ri) = k;  % Beginning of first causal run
    run_length(ri) = 1;
    k = k+1;
    
    while(k <= num_frames)  % Complete search for causal runs
        
        if(~edge(k))
            if(~interlaced)  % Progressive algorithm
                if( (best(k)-best(k-1)>=0) && (best(k)-best(k-1)<=cjump) )
                    run_length(ri) = run_length(ri)+1;
                    k = k +1;
                else
                    % Skip ahead to the next valid causal alignment
                    while(edge(k))
                        k = k+1;
                        if (k > num_frames)
                            break;
                        end
                    end
                    if (k <= num_frames)
                        ri = ri + 1;
                        run_beg(ri) = k;
                        run_length(ri) = 1;
                        k = k + 1;
                    end
                end
            else  % Interlaced algorithm
                if( (min(repmat(best(k),1,run_length(ri))-best(k-run_length(ri):k-1))>=-1) && (best(k)-best(k-1)<=2*cjump) )
                    run_length(ri) = run_length(ri)+1;
                    k = k +1;
                else
                    % Skip ahead to the next valid causal alignment
                    while(edge(k))
                        k = k+1;
                        if (k > num_frames)
                            break;
                        end
                    end
                    if (k <= num_frames)
                        ri = ri + 1;
                        run_beg(ri) = k;
                        run_length(ri) = 1;
                        k = k +1;
                    end
                end
            end
        else % hit an edge alignment, close off this causal segment
            % Skip ahead to the next valid causal alignment
            while(edge(k))
                k = k+1;
                if (k > num_frames)
                    break;
                end
            end
            if (k <= num_frames)
                ri = ri + 1;
                run_beg(ri) = k;
                run_length(ri) = 1;
                k = k +1;
            end
        end
        
    end
    
    % Sort the causal runs according to their length
    [run_length_sort length_order] = sort(run_length,'descend');
    run_beg_sort = run_beg(length_order);
    
end

%%%%%%%%%%
% Stage 4
% Fill normal causal segments from the longest to the shortest
%%%%%%%%%%
causal = zeros(1,num_frames);  % Holds the non-fuzzy causal alignment

%  This code will fill in the normal causal segments from the longest
%  to the shortest.  To qualify for a fill, a normal causal segment
%  must have at least min_length points.
if(~isempty(run_length))
    
    min_length = 2;
    for k = 1:length(run_length_sort)
        
        if (run_length_sort(k) >= min_length)
            
            % See if there are any runs before or after the current run
            if (run_beg_sort(k)+run_length_sort(k) <= num_frames)
                after_index = find(causal(run_beg_sort(k)+run_length_sort(k):num_frames) ~= 0);
            else
                after_index = [];
            end
            if (run_beg_sort(k)-1 >= 1)
                before_index = find(causal(1:run_beg_sort(k)-1) ~= 0);
            else
                before_index = [];
            end
            
            % Determine the index of the nearest after and before for
            % progressive algorithm
            if(~isempty(before_index))
                before_orig = before_index(length(before_index));
            end
            if (~isempty(after_index))
                after_orig = run_beg_sort(k)+run_length_sort(k)+after_index(1)-1;
            end
            
            % If causal, fit this run's alignment into the existing
            if (~interlaced) % Progressive algorithm
                if ( isempty(before_index) && isempty(after_index) )
                    causal(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1) = ...
                        best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1);
                elseif ( isempty(before_index) )
                    if (best(after_orig)>=best(run_beg_sort(k)+run_length_sort(k)-1))
                        causal(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1) = ...
                            best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1);
                    end
                elseif ( isempty(after_index) )
                    if (best(before_orig)<=best(run_beg_sort(k)))
                        causal(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1) = ...
                            best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1);
                    end
                elseif ((best(before_orig)<=best(run_beg_sort(k))) && (best(after_orig)>=best(run_beg_sort(k)+run_length_sort(k)-1)))
                    causal(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1) = ...
                        best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1);
                end
            else  % Interlaced algorithm allows for up to 1 field jump back in time for the segment you are adding into timeline
                if ( isempty(before_index) && isempty(after_index) )
                    causal(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1) = ...
                        best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1);
                elseif ( isempty(before_index) )
                    if (min(best(after_index+run_beg_sort(k)+run_length_sort(k)-1)) >= ...
                            max(best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1))-1)
                        causal(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1) = ...
                            best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1);
                    end
                elseif ( isempty(after_index) )
                    if (max(best(before_index))<=min(best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1))+1)
                        causal(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1) = ...
                            best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1);
                    end
                elseif ( (min(best(after_index+run_beg_sort(k)+run_length_sort(k)-1)) >= ...
                        max(best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1))-1) && ...
                        (max(best(before_index))<=min(best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1))+1) )
                    causal(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1) = ...
                        best(run_beg_sort(k):run_beg_sort(k)+run_length_sort(k)-1);
                end
            end
            
        end
        
    end
    
end

%%%%%%%%%%
% Stage 5
% Fill in the missing holes in the causal array - see algorithm description 
%%%%%%%%%%

%  At this point the causal array may still contain holes (zeros).  The 
%  algorithm for filling in these segments is as follows: The earlier
%  time point is extended later in time by the value in the best array
%  if this is allowed by the causality rules, then the later time point
%  is extended earlier in time by the value in the best array if this
%  is allowed by the causality rules, alternating back and forth until
%  no other extensions are possible.  If causality cannot be achieved,
%  and the final_fuzzy_index has valid next best alignments (i.e.,
%  something other than 'NaN'), these are examined to see if causality
%  can be achieved using them rather than the best array.  If
%  insertions are made from the final_fuzzy_index array, these are
%  compared to what interpolation would have yielded to fill in the
%  segment hole - and the choice that produces the minimum RMSE with
%  respect to the best array is chosen.
%
%  In this manner, missing segments are filled in as they are
%  encountered from early to late time. When the missing segment occurs
%  at the beginning or end of the time sequence, then the holes are
%  filled by the best array alignments (if they are causal), or the
%  optimum choice between final_fuzzy_index alignments (if available
%  and if causal) and the last good causal alignment which is
%  replicated/extended. The last good causal alignment is
%  replicated/extended if no other options are available.

holes = find(causal==0);
while(~isempty(holes))  % keep filling hole segments from early to late time
    num_holes = length(holes);
    
    % Find the length of the first hole
    hole_start = holes(1);
    hole_stop = hole_start;
    if (num_holes > 1)
        for k = 2:num_holes
            if (holes(k)-holes(k-1) == 1)
                hole_stop = hole_stop+1;
            else
                break;
            end
        end
    end
    
    % Fill the hole segment - there are four cases
    % 1. Segment is at beginning of clip but is followed by valid causal
    % 2. Segment is at end of clip but is preceded by valid causal
    % 3. Segment is preceded and followed by valid causal (most likely)
    % 4. Segment is not proceeded or followed by valid causal (function aborts) 
    
    % Do progressive filling algorithm first
    if(~interlaced)
        
        % Case 3, probably the most likely case
        if (hole_start > 1 && hole_stop < num_frames)
            
            % This code toggles from beg to end, starting at beg
            hole_pts = hole_stop-hole_start+1;  % number of points in this hole segment
            valid_early = hole_start-1;  % last valid causal point before time segment
            valid_late = hole_stop+1; % last valid causal point after time segment
            for k = 1:hole_pts
                if (mod(k,2)==1) % try to fill from the beginning first and then from the end if that fails
                    if (best(valid_early+1)>=causal(valid_early) && best(valid_early+1)<=causal(valid_late))
                        causal(valid_early+1) = best(valid_early+1);
                        valid_early = valid_early + 1;
                    elseif (best(valid_late-1)<=causal(valid_late) && best(valid_late-1)>=causal(valid_early))
                        causal(valid_late-1) = best(valid_late-1);
                        valid_late = valid_late - 1;
                    else % see if other possible alignments exist besides the best array
                        early_insert = 0;  % tells if this is an early or late insert
                        succeed = 0;  % set to 1 when a causal substitution is made
                        depth = 2;  % the final_fuzzy_index level, depth = 1 is the best array
                        while (~succeed && (~isnan(final_fuzzy_index(depth,valid_early+1)) || ...
                                ~isnan(final_fuzzy_index(depth,valid_late-1))) )
                            if (~isnan(final_fuzzy_index(depth,valid_early+1)))
                                if (final_fuzzy_index(depth,valid_early+1) >= causal(valid_early) && ...
                                        final_fuzzy_index(depth,valid_early+1)<=causal(valid_late))
                                    causal(valid_early+1) = final_fuzzy_index(depth,valid_early+1);
                                    valid_early = valid_early + 1;
                                    succeed = 1;
                                    early_insert = 1;
                                end
                            else
                                if (final_fuzzy_index(depth,valid_late-1)<=causal(valid_late) && ...
                                        final_fuzzy_index(depth,valid_late-1)>=causal(valid_early))
                                    causal(valid_late-1) = final_fuzzy_index(depth,valid_late-1);
                                    valid_late = valid_late - 1;
                                    succeed = 1;
                                end
                            end
                            depth = depth +1;
                            if (depth > range)
                                break;
                            end
                        end
                        if (~succeed) % linear interpolation at beginning as last resort after all possible alignments exhausted
                            causal(valid_early+1) = causal(valid_early) + ...
                                round((causal(valid_late)-causal(valid_early))/(valid_late-valid_early));
                            valid_early = valid_early + 1;
                        else % see if interpolation yields better RMSE than the alignment from the final_fuzzy_index
                            if (early_insert)
                                interp_value = causal(valid_early-1) + ...
                                    round((causal(valid_late)-causal(valid_early-1))/(valid_late-(valid_early-1)));
                                if (abs(causal(valid_early)-best(valid_early)) > abs(interp_value-best(valid_early)))
                                    causal(valid_early) = interp_value;
                                end
                            else % late insert
                                interp_value = causal(valid_late+1) - ...
                                    round((causal(valid_late+1)-causal(valid_early))/(valid_late+1-valid_early));
                                if (abs(causal(valid_late)-best(valid_late)) > abs(interp_value-best(valid_late)))
                                    causal(valid_late) = interp_value;
                                end
                            end
                        end
                    end
                else % try to fill from the end first and then from the beginning if that fails
                    if (best(valid_late-1)<=causal(valid_late) && best(valid_late-1)>=causal(valid_early))
                        causal(valid_late-1) = best(valid_late-1);
                        valid_late = valid_late - 1;
                    elseif (best(valid_early+1)>=causal(valid_early) && best(valid_early+1)<=causal(valid_late))
                        causal(valid_early+1) = best(valid_early+1);
                        valid_early = valid_early + 1;
                    else
                        early_insert = 0; % tells if this is an early or late insert
                        succeed = 0;  % set to 1 when a causal substitution is made
                        depth = 2;  % the final_fuzzy_index level, depth = 1 is the best array
                        while (~succeed && (~isnan(final_fuzzy_index(depth,valid_early+1)) || ...
                                ~isnan(final_fuzzy_index(depth,valid_late-1))) )
                            if (~isnan(final_fuzzy_index(depth,valid_late-1)))
                                if (final_fuzzy_index(depth,valid_late-1)<=causal(valid_late) && ...
                                        final_fuzzy_index(depth,valid_late-1)>=causal(valid_early))
                                    causal(valid_late-1) = final_fuzzy_index(depth,valid_late-1);
                                    valid_late = valid_late - 1;
                                    succeed = 1;
                                end
                            else
                                if (final_fuzzy_index(depth,valid_early+1)>=causal(valid_early) && ...
                                        final_fuzzy_index(depth,valid_early+1)<=causal(valid_late))
                                    causal(valid_early+1) = final_fuzzy_index(depth,valid_early+1);
                                    valid_early = valid_early + 1;
                                    succeed = 1;
                                    early_insert = 1;
                                end
                            end
                            depth = depth +1;
                            if (depth > range)
                                break;
                            end
                        end
                        if (~succeed) % perform linear interpolation at end as a last resort after all possible alignments exhausted
                            causal(valid_late-1) = causal(valid_late) - ...
                                round((causal(valid_late)-causal(valid_early))/(valid_late-valid_early));
                            valid_late = valid_late - 1;
                        else % see if interpolation yields better RMSE than the alignment from the final_fuzzy_index
                            if (early_insert)
                                interp_value = causal(valid_early-1) + ...
                                    round((causal(valid_late)-causal(valid_early-1))/(valid_late-(valid_early-1)));
                                if (abs(causal(valid_early)-best(valid_early)) > abs(interp_value-best(valid_early)))
                                    causal(valid_early) = interp_value;
                                end
                            else % late insert
                                interp_value = causal(valid_late+1) - ...
                                    round((causal(valid_late+1)-causal(valid_early))/(valid_late+1-valid_early));
                                if (abs(causal(valid_late)-best(valid_late)) > abs(interp_value-best(valid_late)))
                                    causal(valid_late) = interp_value;
                                end
                            end
                        end
                    end
                end
            end
            
            % Case 1, always fill from the end
        elseif (hole_start == 1 && hole_stop < num_frames)
            hole_pts = hole_stop-hole_start+1;  % number of points in this hole segment
            valid_late = hole_stop+1; % last valid causal point after time segment
            for k = 1:hole_pts
                % try to fill from the end
                if (best(valid_late-1)<=causal(valid_late))
                    causal(valid_late-1) = best(valid_late-1);
                    valid_late = valid_late - 1;
                else
                    succeed = 0;  % set to 1 when a causal substitution is made
                    depth = 2;  % the final_fuzzy_index level, depth = 1 is the best array
                    while (~succeed && ~isnan(final_fuzzy_index(depth,valid_late-1)))
                        if (final_fuzzy_index(depth,valid_late-1)<=causal(valid_late))
                            causal(valid_late-1) = final_fuzzy_index(depth,valid_late-1);
                            valid_late = valid_late - 1;
                            succeed = 1;
                        end
                        depth = depth +1;
                        if (depth > range)
                            break;
                        end
                    end
                    if(~succeed) % extend last good alignment
                        causal(valid_late-1) = causal(valid_late);
                        valid_late = valid_late - 1;
                    else % see if extension yields better RMSE than the alignment from the final_fuzzy_index
                        interp_value = causal(valid_late+1);
                        if (abs(causal(valid_late)-best(valid_late)) > abs(interp_value-best(valid_late)))
                            causal(valid_late) = interp_value;
                        end
                    end
                end
            end
            
            % Case 2, always fill from the beginning
        elseif (hole_start > 1 && hole_stop == num_frames)
            hole_pts = hole_stop-hole_start+1;  % number of points in this hole segment
            valid_early = hole_start-1;  % last valid causal point before time segment
            for k = 1:hole_pts
                % try to fill from the beginning
                if (best(valid_early+1)>=causal(valid_early))
                    causal(valid_early+1) = best(valid_early+1);
                    valid_early = valid_early + 1;
                else
                    succeed = 0;
                    depth = 2;
                    while (~succeed && ~isnan(final_fuzzy_index(depth,valid_early+1)))
                        if (final_fuzzy_index(depth,valid_early+1)>=causal(valid_early))
                            causal(valid_early+1) = final_fuzzy_index(depth,valid_early+1);
                            valid_early = valid_early + 1;
                            succeed = 1;
                        end
                        depth = depth +1;
                        if (depth > range)
                            break;
                        end
                    end
                    if (~succeed) % extend last good alignment
                        causal(valid_early+1) = causal(valid_early);
                        valid_early = valid_early + 1;
                    else % see if extension yields better RMSE than the alignment from the final_fuzzy_index
                        interp_value = causal(valid_early-1);
                        if (abs(causal(valid_early)-best(valid_early)) > abs(interp_value-best(valid_early)))
                            causal(valid_early) = interp_value;
                        end
                    end
                end
            end
            
            % Case 4, super rare case might be possible, no valid causal segments anywhere in the entire clip
        else
            warning('Warning: No valid 2-pt causal alignment segments found in clip - aborting est_var_frame_delays.');
            return;
            
        end
        
    else  % Interlaced algorithm allows for 1 field jump back in time for segment you are adding
        
        % Case 3, probably the most likely case
        if (hole_start > 1 && hole_stop < num_frames)
            
            % This code toggles from beg to end, starting at beg
            hole_pts = hole_stop-hole_start+1;  % number of points in this hole segment
            valid_early = hole_start-1;  % last valid causal point before time segment
            valid_late = hole_stop+1; % last valid causal point after time segment
            for k = 1:hole_pts
                causal_temp = causal(valid_late:num_frames);
                min_valid_late = min(causal_temp(causal_temp~=0));  % don't include causal=0 points
                max_valid_early = max(causal(1:valid_early)); % causal=0 points don't affect this calculation
                if (mod(k,2)==1) % try to fill from the beginning first and then from the end if that fails
                    if (best(valid_early+1)>=max_valid_early-1 && best(valid_early+1)<=min_valid_late+1)
                        causal(valid_early+1) = best(valid_early+1);
                        valid_early = valid_early + 1;
                    elseif (best(valid_late-1)<=min_valid_late+1 && best(valid_late-1)>=max_valid_early-1)
                        causal(valid_late-1) = best(valid_late-1);
                        valid_late = valid_late - 1;
                    else % see if other possible alignments exist besides the best array
                        early_insert = 0; % tells if this is an early or late insert
                        succeed = 0;
                        depth = 2;
                        while (~succeed && (~isnan(final_fuzzy_index(depth,valid_early+1)) || ...
                                ~isnan(final_fuzzy_index(depth,valid_late-1))) )
                            if (~isnan(final_fuzzy_index(depth,valid_early+1)))
                                if (final_fuzzy_index(depth,valid_early+1) >= ...
                                        max_valid_early-1 && final_fuzzy_index(depth,valid_early+1)<=min_valid_late+1)
                                    causal(valid_early+1) = final_fuzzy_index(depth,valid_early+1);
                                    valid_early = valid_early + 1;
                                    succeed = 1;
                                    early_insert = 1;
                                end
                            else
                                if (final_fuzzy_index(depth,valid_late-1)<=min_valid_late+1 && ...
                                        final_fuzzy_index(depth,valid_late-1)>=max_valid_early-1)
                                    causal(valid_late-1) = final_fuzzy_index(depth,valid_late-1);
                                    valid_late = valid_late - 1;
                                    succeed = 1;
                                end
                            end
                            depth = depth +1;
                            if (depth > range)
                                break;
                            end
                        end
                        if (~succeed) % linear interpolation at beginning as last resort after all possible alignments exhausted
                            causal(valid_early+1) = causal(valid_early) + ...
                                round((causal(valid_late)-causal(valid_early))/(valid_late-valid_early));
                            valid_early = valid_early + 1;
                        else % see if interpolation yields better RMSE than the alignment from the final_fuzzy_index
                            if (early_insert)
                                interp_value = causal(valid_early-1) + ...
                                    round((causal(valid_late)-causal(valid_early-1))/(valid_late-(valid_early-1)));
                                if (abs(causal(valid_early)-best(valid_early)) > abs(interp_value-best(valid_early)))
                                    causal(valid_early) = interp_value;
                                end
                            else % late insert
                                interp_value = causal(valid_late+1) - ...
                                    round((causal(valid_late+1)-causal(valid_early))/(valid_late+1-valid_early));
                                if (abs(causal(valid_late)-best(valid_late)) > abs(interp_value-best(valid_late)))
                                    causal(valid_late) = interp_value;
                                end
                            end
                        end
                    end
                else % try to fill from the end first and then from the beginning if that fails
                    if (best(valid_late-1)<=min_valid_late+1 && best(valid_late-1)>=max_valid_early-1)
                        causal(valid_late-1) = best(valid_late-1);
                        valid_late = valid_late - 1;
                    elseif (best(valid_early+1)>=max_valid_early-1 && best(valid_early+1)<=min_valid_late+1)
                        causal(valid_early+1) = best(valid_early+1);
                        valid_early = valid_early + 1;
                    else
                        early_insert = 0; % tells if this is an early or late insert
                        succeed = 0;  % set to 1 when a causal substitution is made
                        depth = 2;  % the final_fuzzy_index level, depth = 1 is the best array
                        while (~succeed && (~isnan(final_fuzzy_index(depth,valid_early+1)) || ...
                                ~isnan(final_fuzzy_index(depth,valid_late-1))) )
                            if (~isnan(final_fuzzy_index(depth,valid_late-1)))
                                if (final_fuzzy_index(depth,valid_late-1)<=min_valid_late+1 && ...
                                        final_fuzzy_index(depth,valid_late-1)>=max_valid_early-1)
                                    causal(valid_late-1) = final_fuzzy_index(depth,valid_late-1);
                                    valid_late = valid_late - 1;
                                    succeed = 1;
                                end
                            else
                                if (final_fuzzy_index(depth,valid_early+1)>=max_valid_early-1 && ...
                                        final_fuzzy_index(depth,valid_early+1)<=min_valid_late+1)
                                    causal(valid_early+1) = final_fuzzy_index(depth,valid_early+1);
                                    valid_early = valid_early + 1;
                                    succeed = 1;
                                    early_insert = 1;
                                end
                            end
                            depth = depth +1;
                            if (depth > range)
                                break;
                            end
                        end
                        if (~succeed) % perform linear interpolation at end
                            causal(valid_late-1) = causal(valid_late) - ...
                                round((causal(valid_late)-causal(valid_early))/(valid_late-valid_early));
                            valid_late = valid_late - 1;
                        else % see if interpolation yields better RMSE than the alignment from the final_fuzzy_index
                            if (early_insert)
                                interp_value = causal(valid_early-1) + ...
                                    round((causal(valid_late)-causal(valid_early-1))/(valid_late-(valid_early-1)));
                                if (abs(causal(valid_early)-best(valid_early)) > abs(interp_value-best(valid_early)))
                                    causal(valid_early) = interp_value;
                                end
                            else % late insert
                                interp_value = causal(valid_late+1) - ...
                                    round((causal(valid_late+1)-causal(valid_early))/(valid_late+1-valid_early));
                                if (abs(causal(valid_late)-best(valid_late)) > abs(interp_value-best(valid_late)))
                                    causal(valid_late) = interp_value;
                                end
                            end
                        end
                    end
                end
            end
            
            % Case 1, always fill from the end
        elseif (hole_start == 1 && hole_stop < num_frames)
            hole_pts = hole_stop-hole_start+1;  % number of points in this hole segment
            valid_late = hole_stop+1; % last valid causal point after time segment
            for k = 1:hole_pts
                causal_temp = causal(valid_late:num_frames);
                min_valid_late = min(causal_temp(causal_temp~=0));  % don't include causal=0 points
                % try to fill from the end
                if (best(valid_late-1) <= min_valid_late+1)
                    causal(valid_late-1) = best(valid_late-1);
                    valid_late = valid_late - 1;
                else
                    succeed = 0;  % set to 1 when a causal substitution is made
                    depth = 2;  % the final_fuzzy_index level, depth = 1 is the best array
                    while (~succeed && ~isnan(final_fuzzy_index(depth,valid_late-1)))
                        if (final_fuzzy_index(depth,valid_late-1)<=min_valid_late+1)
                            causal(valid_late-1) = final_fuzzy_index(depth,valid_late-1);
                            valid_late = valid_late - 1;
                            succeed = 1;
                        end
                        depth = depth +1;
                        if (depth > range)
                            break;
                        end
                    end
                    if (~succeed) % extend last good alignment
                        causal(valid_late-1) = causal(valid_late);
                        valid_late = valid_late - 1;
                    else % see if extension yields better RMSE than the alignment from the final_fuzzy_index
                        interp_value = causal(valid_late+1);
                        if (abs(causal(valid_late)-best(valid_late)) > abs(interp_value-best(valid_late)))
                            causal(valid_late) = interp_value;
                        end
                    end
                end
            end
            
            % Case 2, always fill from the beginning
        elseif (hole_start > 1 && hole_stop == num_frames)
            hole_pts = hole_stop-hole_start+1;  % number of points in this hole segment
            valid_early = hole_start-1;  % last valid causal point before time segment
            for k = 1:hole_pts
                max_valid_early = max(causal(1:valid_early));  % causal=0 points don't affect this calculation
                % try to fill from the beginning
                if (best(valid_early+1) >= max_valid_early-1)
                    causal(valid_early+1) = best(valid_early+1);
                    valid_early = valid_early + 1;
                else
                    succeed = 0;
                    depth = 2;
                    while (~succeed && ~isnan(final_fuzzy_index(depth,valid_early+1)))
                        if (final_fuzzy_index(depth,valid_early+1)>=max_valid_early-1)
                            causal(valid_early+1) = final_fuzzy_index(depth,valid_early+1);
                            valid_early = valid_early + 1;
                            succeed = 1;
                        end
                        depth = depth +1;
                        if (depth > range)
                            break;
                        end
                    end
                    if (~succeed) % extend last good alignment
                        causal(valid_early+1) = causal(valid_early);
                        valid_early = valid_early + 1;
                    else % see if extension yields better RMSE than the alignment from the final_fuzzy_index
                        interp_value = causal(valid_early-1);
                        if (abs(causal(valid_early)-best(valid_early)) > abs(interp_value-best(valid_early)))
                            causal(valid_early) = interp_value;
                        end
                    end
                end
            end
            
            % Case 4, very rare case might be possible, no valid 2-pt causal segments anywhere in the entire clip
        else
            warning('Warning: No valid 2-pt causal alignment segments found in clip - aborting est_var_frame_delays.');
            return;
            
        end
        
    end
    
    % Update the holes that are remaining
    holes = find(causal==0);
    
end

% Compute the average alignment RMSE between the causal alignment and the
% minimum MSE alignment, this is return argument number 2.
results_rmse = sqrt(mean((causal-best).^2));

if (verbose)
    
    % Plot final causal alignment
    figure(2)
    plot(best,'LineWidth',1.25)
    hold on
    set(gca,'LineWidth',2)
    set(gca,'FontName','Ariel')
    set(gca,'fontsize',12)
    plot(causal,'r','LineWidth',3)
    xlabel('Processed Frame or Field');
    ylabel('Original Frame or Field');
    title('Processed Clip Alignment Results');
    grid on
    hold off;
    pause(0.01);
    
    %  Print out results for verbose mode.
    fprintf('results_rmse %f\n',results_rmse);
    
    pause(5);

end

% Compute the causal index and its MSE and put into an array for later use
causal_index = zeros(1,num_frames);
causal_mse = zeros(1,num_frames);
for j = 1:num_frames
    this_mse = mse(:,j);
    this_index = orig_index(:,j);
    causal_index(j) = find(this_index == causal(j));
    causal_mse(j) = this_mse(causal_index(j));
end

if (verbose)
    
    % Plot the 3d correlation function, x = proc frame, y = orig frame, z =
    % mse, where x is the second index and the y is the first index
    orig_last = max(max(orig_index));  % last orig frame index that has data
    [x, y] = meshgrid(1:num_frames, 1:orig_last);
    [ysize, xsize] = size(x);
    z = NaN(ysize, xsize);  % missing values will be NaN
    for j = 1:num_frames  % fill column by column
        for i = 1:range  % row element
            z(orig_index(i,j),j) = mse(i,j);
        end
    end
    
    figure(3)
    mesh(x,y,z, 'EdgeColor', 'black');
    hold on
    hidden on
    set(gca,'LineWidth',2)
    set(gca,'FontName','Ariel')
    set(gca,'fontsize',12)
    xlabel('Processed Frame or Field');
    ylabel('Original Frame or Field');
    zlabel('MSE');
    title('3D MSE Function');
    shading interp

    % Now overlay the causal function MSE using a red line
    x_causal = 1:num_frames;
    y_causal = causal;
    z_causal = causal_mse;
    warning off
    alpha(0.9);  % Sets some face transparency to see points
    scatter3(x_causal,y_causal,z_causal,'r','filled');
    hold off
    pause(5);  % pause 5 sec to display plots properly

end

% Assign results array as requested by the user
if (~is_causal)
    
    results = best;
    results_fuzzy = final_fuzzy_index;
    results_fuzzy_mse = final_fuzzy_mse;

else % further processing is required for fuzzy alignments
    
    results = causal;
    
    %%%%%%%%%%
    % Stage 6
    % Expand the fuzzy alignments to include the causal alignments.
    % Some of this code is replicated from stage 1.
    %%%%%%%%%%
    
    % Step thru the frames for this clip and find the causal fuzzy
    % alignments.  The causal fuzzy alignments expands the set of final
    % fuzzy alignments found previously to include the causal alignment
    % point and all points in-between this point and the prior set of final
    % fuzzy alignments.
    final_fuzzy_causal_index = NaN(range, num_frames); % The final fuzzy causal alignment from most to least likely
    final_fuzzy_causal_mse = NaN(range, num_frames); % The corresponding MSEs
    for j = 1:num_frames
        
%         % Plot the correlation function
%         figure(1)
%         plot(orig_index(:,j),mse(:,j),'LineWidth',2);
%         hold on
%         grid on
%         set(gca,'LineWidth',2)
%         set(gca,'FontName','Ariel')
%         set(gca,'fontsize',12)
%         xlabel('Original Frame or Field');
%         ylabel('Mean Squared Error (MSE)');
%         this_title = ['Processed Frame or Field ' int2str(j)];
%         title(this_title);
%         hold off
        
        % Since subtractive correlation is being used, the correlation function
        % will be a minimum at the best match and increase from there.  This
        % code uses a thresholding scheme that sets a threshold at the minimum
        % correlation value + thres * maximum correlation value, but the
        % maximum correlation value used is the top (95%) rank sorted value
        % for robustness.  All alignments within this threshold from the best
        % correlation value are considered possible alignments for later
        % processing by the algorithm.
        this_mse = mse(:,j);
        this_mse_valid = find(~isnan(this_mse));  % array can contain NaN if orig frames did not exist
        range_valid = length(this_mse_valid);
        this_mse_sort = sort(this_mse);
        this_index = orig_index(:,j);
        this_index_sort = orig_index_sort(:,j);
        
        max_corr = this_mse_sort(floor(top*range_valid));
        min_corr = this_mse_sort(1);
        fuzzy = find(this_mse >= min_corr & this_mse <= min_corr+thres*max_corr);
        num_fuzzy = length(fuzzy);  % will always be at least one
        fuzzy_sort = [];
        for k = 1:num_fuzzy
            fuzzy_sort(k) = find(this_index_sort == this_index(fuzzy(k)));
        end
        
%         % Add the fuzzy alignments to the plot, red points meet threshold
%         figure(1)
%         hold on
%         plot(this_index(fuzzy),this_mse(fuzzy),'r.','MarkerSize',15);
%         hold off
        
        % Expand the fuzzy alignments to include the range of original indices
        % that are covered by the minimum MSE points.  This will include
        % extra points that were not in the minimum MSE fuzzy set.
        bigger_fuzzy = (fuzzy(1):fuzzy(length(fuzzy)))';
        
        % Added fuzzy points that were not included in the original set
        added_fuzzy = setxor(fuzzy,bigger_fuzzy); % just the added points, for plotting
        added_fuzzy_mse = this_mse(added_fuzzy);
        added_fuzzy_index = this_index(added_fuzzy);
        num_added_fuzzy = length(added_fuzzy);
        added_fuzzy_sort = [];
        for k = 1:num_added_fuzzy
            added_fuzzy_sort(k) = find(this_index_sort == this_index(added_fuzzy(k)));
        end
        
%         % Add the expanded fuzzy points to the plots in green squares.
%         figure(1)
%         hold on
%         plot(this_index(added_fuzzy),this_mse(added_fuzzy),'gs','MarkerSize',5,'MarkerFaceColor','g','MarkerEdgeColor','g');
%         hold off
        
%         % Overlay the causal point with a large black diamond
%         figure(1)
%         hold on
%         plot(this_index(causal_index(j)),causal_mse(j),'kd','MarkerSize',10,'LineWidth',1.25);
%         hold off
        
        % Now expand the fuzzy set to include the causal point
        added_causal = [];
        added_causal_mse = [];
        cplot = 0;
        if (this_index(causal_index(j)) < this_index(fuzzy(1))) % Must expand on left
            added_causal = (causal_index(j):(fuzzy(1)-1))';
            added_causal_mse = this_mse(added_causal);
            bigger_fuzzy = (causal_index(j):fuzzy(length(fuzzy)))';
            cplot = 1;
        elseif (this_index(causal_index(j)) > this_index(fuzzy(length(fuzzy)))) % Must expand on right
            added_causal = ((fuzzy(length(fuzzy))+1):causal_index(j))';
            added_causal_mse = this_mse(added_causal);
            bigger_fuzzy = (fuzzy(1):causal_index(j))';
            cplot = 1;
        end
        
%         % Overlay the added causal points in solid black diamonds
%         if(cplot)
%             figure(1)
%             hold on
%             plot(this_index(added_causal),added_causal_mse,'kd','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
%             hold off
%         end
        
        % Sort the bigger fuzzy information from the most likely alignment
        % to the least likely alignment, based on MSE.  Then move the causal
        % alignment point to the first element if it is not there already.
        % The array final_fuzzy_causal_index gives the fuzzy alignments
        % with causal processing.
        bigger_fuzzy_mse = this_mse(bigger_fuzzy);
        bigger_fuzzy_index = this_index(bigger_fuzzy);
        [bigger_fuzzy_mse_sort mse_order] = sort(bigger_fuzzy_mse);
        bigger_fuzzy_index_sort = bigger_fuzzy_index(mse_order);
        tlen = length(bigger_fuzzy_index_sort);
        tind = find(bigger_fuzzy_index_sort == this_index(causal_index(j)));
        if (tind ~= 1)  % reorder so causal alignment is first
            if (tind ~= tlen)
                bigger_fuzzy_index_sort = cat(1, bigger_fuzzy_index_sort(tind), bigger_fuzzy_index_sort(1:tind-1), ...
                    bigger_fuzzy_index_sort(tind+1:tlen));
                bigger_fuzzy_mse_sort = cat(1, bigger_fuzzy_mse_sort(tind), bigger_fuzzy_mse_sort(1:tind-1), ...
                    bigger_fuzzy_mse_sort(tind+1:tlen));
            else
                bigger_fuzzy_index_sort = cat(1, bigger_fuzzy_index_sort(tind), bigger_fuzzy_index_sort(1:tind-1));
                bigger_fuzzy_mse_sort = cat(1, bigger_fuzzy_mse_sort(tind), bigger_fuzzy_mse_sort(1:tind-1));
            end
        end
        final_fuzzy_causal_index(1:tlen,j) = bigger_fuzzy_index_sort;
        final_fuzzy_causal_mse(1:tlen,j) = bigger_fuzzy_mse_sort;
        
    end
    
    % Assign the fuzzy results for a causal user request
    results_fuzzy = final_fuzzy_causal_index;
    results_fuzzy_mse = final_fuzzy_causal_mse;
    
end

return;

end

