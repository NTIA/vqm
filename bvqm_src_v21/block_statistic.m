function [s1,s2,s3] = block_statistic(y,vsize,hsize,varargin)
% BLOCK_STATISTIC
%  Extract feature from each spatial-temporal (S-T) region, producing one number 
%  for each S-T block.  Takes a block of perceptually filtered images 
%  and produces features.  
% SYNTAX
%  [s1] = block_statistic(y,vsize,hsize,'Stat1');
%  [s1,s2] = block_statistic(y,vsize,hsize,'Stat1','Stat2');
%  [s1,s2,s3] = block_statistic(y,vsize,hsize,'Stat1','Stat2','Stat3');
% DEFINITION
%  [s1] = block_statistic(y,'Stat1'); divides time-slice of images 'y' into
%  abutting Spatial-Temporal (S-T) regions that contain 'vsize' pixels
%  vertically and 'hsize' pixels horizontally.  Then, statstic 'Stat1' is
%  computed over each S-T region, and the results are returned in 's1'.
%  When called with the names of two statistics, two statistics are 
%  computed and returned, and so forth.  Available statistics are:
%  'mean' ==> compute the mean over each S-T region.
%  'std' ==> compute the standard deviation over each S-T region.
%  'rms' ==> compute the RMS over each S-T region.
%  'fraction' ==> compute fraction of pixels that are greater than or
%                 equal to 1.0
% REMARKS
%  Functionality tested pretty well.

want_mean = 0;
want_std = 0;
want_rms = 0;
want_fraction = 0;

for cnt = 1:nargin-3,
    if strcmp(lower(varargin(cnt)),'mean'),
        want_mean = 1;
        want_mean_at = cnt;
    elseif strcmp(lower(varargin(cnt)),'std'),
        want_std = 1;
        want_std_at = cnt;
    elseif strcmp(lower(varargin(cnt)),'rms'),
        want_rms = 1;
        want_rms_at = cnt;
    elseif strcmp(lower(varargin(cnt)),'fraction'),
        want_fraction = 1;
        want_fraction_at = cnt;
    end
end

if want_mean + want_std + want_rms + want_fraction ~= nargin - 3,
    error('block_statistic ''Stat'' not recognized or repeated.');
end

% check block size request.
[row,col,time] = size(y);
if mod(row,vsize) ~= 0,
    error('vertical size of block must evenly divide the SROI.');
end

if mod(col,hsize) ~= 0,
    error('horizontal size of block must evenly divide the SROI.');
end

if want_mean | want_std,
    temp = sum(y,3);  % sum over time
    temp = sum( reshape(temp,vsize,row*col/vsize)); % sum block vertically
    temp = reshape(temp,row/vsize,col)';
    temp = sum(reshape(temp,hsize,row*col/(vsize*hsize))); % sum block horizontally
    y_sum = reshape(temp,col/hsize,row/vsize)' ./ (hsize*vsize*time); % reshape
end

if want_std | want_rms,
    temp = sum(y.^2,3);
    temp = sum(reshape(temp,vsize,row*col/vsize));
    temp = reshape(temp,row/vsize,col)';
    temp = sum(reshape(temp,hsize,row*col/(vsize*hsize)));
    y_square = reshape(temp,col/hsize,row/vsize)' ./ (hsize*vsize*time);
end

if want_fraction,
    y(find(y >= 1.0)) = 1.0;
    y(find(y < 1.0)) = 0.0;
    
    temp = sum(y,3);  % sum over time
    temp = sum( reshape(temp,vsize,row*col/vsize)); % sum block vertically
    temp = reshape(temp,row/vsize,col)';
    temp = sum(reshape(temp,hsize,row*col/(vsize*hsize))); % sum block horizontally
    y_fraction = reshape(temp,col/hsize,row/vsize)' ./ (hsize*vsize*time); % reshape
end


if want_mean,
    switch want_mean_at,
        case 1,
            s1 = y_sum;
        case 2,
            s2 = y_sum;
        case 3,
            s3 = y_sum;
        otherwise
            error('Code defect')
    end
end

if want_std,
    y_std = sqrt( max(y_square - y_sum .^ 2,0));
    switch want_std_at,
        case 1,
            s1 = y_std;
        case 2,
            s2 = y_std;
        case 3,
            s3 = y_std;
        otherwise
            error('Code defect')
    end
end

if want_rms,
    y_rms = sqrt( y_square );
    switch want_rms_at,
        case 1,
            s1 = y_rms;
        case 2,
            s2 = y_rms;
        case 3,
            s3 = y_rms;
        otherwise
            error('Code defect')
    end
end


if want_fraction,
    switch want_fraction_at,
        case 1,
            s1 = y_fraction;
        case 2,
            s2 = y_fraction;
        case 3,
            s3 = y_fraction;
        otherwise
            error('Code defect')
    end
end

