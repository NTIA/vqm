function [data] = block_statistic_shift(y,vsize,hsize,varargin)
% BLOCK_STATISTIC_SHIFT
%  Extract feature from each spatial-temporal (S-T) region.  Produce one set 
%  of features for each of 9 shifts:  +/- 1 horizontally & vertically.
%  See also function block_statistic.
% SYNTAX
%  [data] = block_statistic_shift(y,vsize,hsize,'Stat1');
%  [data] = block_statistic_shift(y,vsize,hsize,'Stat1','Stat2');
%  [data] = block_statistic_shift(y,vsize,hsize,'Stat1','Stat2','Stat3');
% DEFINITION
%  [s1] = block_statistic_shift(y,'Stat1'); divides time-slice of images 'y' into
%  abutting Spatial-Temporal (S-T) regions that contain 'vsize' pixels
%  vertically and 'hsize' pixels horizontally.  'y' must contain an extra 1
%  pixels on all sides, to be used for shifts.  Statstic 'Stat1' is
%  computed over each S-T region, and the results are returned in 's1'.
%  When called with the names of two statistics, two statistics are 
%  computed and returned, and so forth.  Available statistics are:
%  'mean' ==> compute the mean over each S-T region, in data(:).mean
%  'std' ==> compute the standard deviation over each S-T region, in data(:).std.
%  'rms' ==> compute the RMS over each S-T region, in data(:).rms
%  'fraction' ==> compute fraction of pixels that are greater than or
%                 equal to 1.0, , in data(:).fraction
%  Return value 'data' is an array length 9, with one element for each
%  shift.

want_mean = 0;
want_std = 0;
want_rms = 0;
want_fraction = 0;

for cnt = 1:nargin-3,
    if strcmp(lower(varargin(cnt)),'mean'),
        want_mean = 1;
    elseif strcmp(lower(varargin(cnt)),'std'),
        want_std = 1;
    elseif strcmp(lower(varargin(cnt)),'rms'),
        want_rms = 1;
    elseif strcmp(lower(varargin(cnt)),'fraction'),
        want_fraction = 1;
    else
        error('block_statistic_shift ''Stat'' not recognized.');
    end
end


% check block size request.
[row,col,time] = size(y);
row = row - 2;
col = col - 2;
if mod(row,vsize) ~= 0,
    error('vertical size of block must evenly divide the SROI.');
end

if mod(col,hsize) ~= 0,
    error('horizontal size of block must evenly divide the SROI.');
end

if want_mean | want_std,
    tempT = sum(y,3);  % sum over time
    loop = 1;
    for cnt1=1:3,
        rng1=(cnt1+row-1);
        for cnt2=1:3,
            rng2=(cnt2+col-1);
            temp = sum( reshape(tempT(cnt1:rng1,cnt2:rng2),vsize,row*col/vsize)); % sum block vertically
            temp = reshape(temp,row/vsize,col)';
            temp = sum(reshape(temp,hsize,row*col/(vsize*hsize))); % sum block horizontally
            data(loop).mean = reshape(temp,col/hsize,row/vsize)' ./ (hsize*vsize*time); % reshape
            loop = loop + 1;
        end
    end
end

if want_std | want_rms,
    tempT = sum(y.^2,3);
    loop = 1;
    for cnt1=1:3,
        for cnt2=1:3,
            rng1=(cnt1+row-1);
            rng2=(cnt2+col-1);
            temp = sum(reshape(tempT(cnt1:rng1,cnt2:rng2),vsize,row*col/vsize));
            temp = reshape(temp,row/vsize,col)';
            temp = sum(reshape(temp,hsize,row*col/(vsize*hsize)));
            data(loop).rms = reshape(temp,col/hsize,row/vsize)' ./ (hsize*vsize*time);
            loop = loop + 1;
        end
    end
end

if want_fraction,
    y(find(y >= 1.0)) = 1.0;
    y(find(y < 1.0)) = 0.0;
    
    tempT = sum(y,3);  % sum over time
    loop = 1;
    for cnt1=1:3,
        for cnt2=1:3,
            rng1=(cnt1+row-1);
            rng2=(cnt2+col-1);
            temp = sum( reshape(tempT(cnt1:rng1,cnt2:rng2),vsize,row*col/vsize)); % sum block vertically
            temp = reshape(temp,row/vsize,col)';
            temp = sum(reshape(temp,hsize,row*col/(vsize*hsize))); % sum block horizontally
            data(loop).fraction = reshape(temp,col/hsize,row/vsize)' ./ (hsize*vsize*time); % reshape
            loop = loop + 1;
        end
    end
end

if want_std,
    for loop=1:9,
        data(loop).std = sqrt( max(data(loop).rms - data(loop).mean .^ 2,0));
    end
end

if want_rms | want_std,
    for loop=1:9,
        data(loop).rms = sqrt( data(loop).rms );
    end
end
