function this_par = dll_vfd_clippar_loop_par1(clip_structs)
% VFD_CLIPPAR_LOOP_PAR1
%  Computes the vfd_par1 parameter of the VQM_VFD model. The vfd_par1
%  parameter quantifies the perceptual motion distortion of a video stream
%  with variable frame delay (VFD).  This parameter is necessary to capture
%  information that is lost by the VFD feature-based clip parameters, which
%  use an original sequence that is changed to match the processed
%  sequence. The vfd_par1 parameter is described in detail in NTIA
%  TM-11-475. The function takes the variable 'clip_structs' (of thes ame
%  format as GClips), which specifies the set of video clips and their
%  properties.  The returned par will be the value of par1. This function 
%  works in conjunction with other dll functions for aditional information.
% SYNTAX
%  this_par = dll_vfd_clippar_loop_par1(clip_structs)
% DESCRIPTION
%  See NTIA TM-11-475 for a technical description of the vfd_par1
%  parameter. This function loads the VFD information from dll_calib_video
%  and uses the information to exctract the vfd_par1 parameter for the pair
%  of video clips.
%

% Set the initial alignment point (in frames) for the VFD
% correction.  first_align is considered frame or field
% number 1 in the VFD results file.
first_align = clip_structs(2).align_start - clip_structs(2).loc_start + 1;

% Set the interlaced flag for the VFD correction.
if (strcmpi(clip_structs(1).video_standard, 'interlace_lower_field_first') == 1)
    is_interlaced = 1;
    field_first = 1;
    first_align = 2*first_align-1;  % convert to fields
elseif (strcmpi(clip_structs(1).video_standard, 'interlace_upper_field_first') == 1)
    is_interlaced = 1;
    field_first = 2;
    first_align = 2*first_align-1;  % convert to fields
else
    is_interlaced = 0;  % progressive
end

% obtain the vfd information
[this_result this_result_fuzzy] = dll_calib_video('get_vfd');

% Determine the total number of aligned frames that are available
% (in frames) for the processed clip.  This will be used to
% generate a gclips alignment result if the VFD algorithm failed.
nframes = clip_structs(2).align_stop-clip_structs(2).align_start+1;
% Reduce by one if the processed clip is reframed
if (is_interlaced && mod(clip_structs(2).spatial.vertical,2))
    nframes = nframes-1;
end

npts = length(this_result);

% Convert the VFD information to a vector that gives abnormal
% forward field/frame jumps.  Subtracting 1 from the abs(diff)
% and maxing with zero produces a parameter that (1) does not
% penalize for normal field/frame delivery (where the VFD
% field/frame indices increase by one from one field/frame to
% the next), (2) does not penalize for progressive frame
% repeats (where the VFD frame indices stay fixed from one
% frame to the next), and (3) does not penalize for interlaced
% frame repeats (where the VFD field indices jump back in time
% from one field to the next). A non-impairment value is used
% for the first field/frame (which must be padded since it's
% diff is not available).  If one assumes that both the
% early and late time sides used for the frame jump estimates
% are absolutely correct (i.e., no fuzzy alignments), then this
% equation would result:
%
% this_result = max([0 abs(diff(this_result))-1], 0);
%
% One could also use a modification of the above diff that
% subtracts the current frame position from the max fuzzy
% causal alignment of the prior frame.  This is thus a more
% conservative estimate of the frame jumps (i.e., smaller frame
% jumps) at any given position since it allows uncertainty on
% the early time side. However, this allows no uncertainty on
% the late time side, since that alignment is considered
% absolutely correct. This equation would result:
%
% fuzzy_max = min(max(this_result_fuzzy(:,1:npts-1)),this_result(2:npts));
% this_result = max([0 abs(this_result(2:npts)-fuzzy_max)-1], 0);
%
% A final modification assumes uncertainty on both the early and
% late time sides so it is the most forgiving of all three.
% Thus, frame jumps are only penalized when they are absolutely
% certain to be correct (i.e., no fuzzy overlapping).  Even
% then, when analyzing this algorithm for the test scene
% t1y1_kiel12_A1pass, it is clear that the VFD estimation
% algorithm is overconfident for some alignments (e.g., not
% enough fuzzy alignments are included for field 511). The
% final equation is:
fuzzy_max_early = min(max(this_result_fuzzy(:,1:npts-1)),this_result(2:npts));
fuzzy_min_late = max(min(this_result_fuzzy(:,2:npts)),fuzzy_max_early);
this_result = max([0 abs(fuzzy_min_late-fuzzy_max_early)-1], 0);

% Calculate vfd_par1 for this clip using log10 of the RMS error.
this_par = log10(sqrt(mean(this_result.^2))+1);






