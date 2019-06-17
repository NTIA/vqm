function [si, hv, hvb] = filter_si_hv(y, varargin)
% FILTER_SI_HV
%
%  Filters Y with the 13x13 gradient (c=2) filters described in SPIE 1999 paper
%
% SYNTAX
%
%  [SI] = filter_si_hv(Y)
%  [SI] = filter_si_hv(Y, rmin, theta)
%  [SI, HV, HVB] = filter_si_hv(...)
%
% DESCRIPTION
%
%  [SI] = filter_si_hv(Y)  Perceptually fiters liminence image Y using 
%  the 13x13 Horizontal and Vertical gradient filters in a fashion similar 
%  to the sobel filter.
%
%  If Y is a 3 dimensional matrix, Y will be presumed to contain multiple
%  images as follows:  (row, col, time).  No execution time penalties occur.
%
%  [SI] = filter_si_hv(Y, rmin, theta)  allows the user to over ride
%  the default values for rmin and theta. 
%
%  [SI, HV, HVB] = filter_si_hv(...)  returns three perceptually fitered
%  versions of image Y:  the SI filtered image, the HV filtered image
%  (containing horiziontal & vertical edges) and the HVB image (containing
%  diagonal edges.)
%
% REMARKS
%
%  rmin defaults to 20, where pixels with a radius (i.e., SI value) less 
%  than rmin are set to zero in HV and HVB images.  
%
%  Theta defaults to 0.225 radians.  Theta is the maximum angle deviation 
%  from the H and V axis for pixels to be considered HV pixels.  
%
%  Returned images (SI, HV, and HVB) are the same size as Y; except that a
%  border of 6 pixels around the edge of each image is invalid.
%

% if pass in a time-slice of 2+ images, reshape into 2-D.
if ndims(y) == 3, 
    must_reshape = 1;
    [row_size, col_size, time_size] = size(y);
    y = reshape(y, row_size, col_size * time_size);   
elseif ndims(y) == 2,
    must_reshape = 0;
    [row_size, col_size, time_size] = size(y);
else
    error('Function ''filter_si_hv'' requires Y to be a 2-D or 3-D image');
end

%  Assign defaults
[row_size, col_size] = size(y);
rmin = 20;
theta = .225;

if (length(varargin) == 2);
    rmin = varargin{1};
    theta = varargin{2};
end

%
if row_size < 13 | col_size < 13,
    error('Function ''filter_si_hv'' requires images to be at least 13x13');
end

% compute angle as a ratio of HV and HVbar.
ratio_threshold = tan(theta);

%  The weights for a single row of the H filter 
%  is given by: w(x) = k*(x/c)*exp{-(1/2)*(x/c)^2}, where x = {-6, -5, ..., 5, 6}, 
%  and k is a normalization constant chosen such that this filter produces the same 
%  amplitude response on an H V edge as the Sobel filter.

%  Generate the 13 long filter mask, in one dimension.
c=2;
filter_mask = zeros(1,13);
for x = -6:1:6
   filter_mask(x+7) = (x/c)*exp(-(1/2)*(x/c)^2);
end
filter_mask = (filter_mask./(13 * sum(filter_mask(1:6)))) * 4;

%  Convolve 13x13 mask with y in horizontal & vertical direction.
%  do two 1x13 convolutions instead of one 13x13, for speed.
horiz = conv2(y,filter_mask,'same');
horiz = conv2(horiz,ones(13,1),'same');

vert = conv2(y,filter_mask','same');
vert = conv2(vert,ones(1,13),'same');

%% for debugging, comment in the below lines.
% h_in = horiz;
% v_in = vert;

% Construct SI image
si = sqrt(horiz.^2 + vert.^2);

% If use only wants to compute SI, skip HV & HVB.  If need be, reshape back
% into 3-D
if nargout == 1
    if must_reshape == 1,
        si = reshape(si,row_size,col_size/time_size,time_size);
    end
    si = si(7:row_size-6,7:col_size/time_size-6,:);
    return;
end

% Start calculation of HV.
% We don't want to use atan2 (because that is slow) so we are going to
% compute the ratio between h & v, putting the smaller value on top and the
% larger value on the bottom.  Ignore divide by zero, because later code
% checking against rmin will catch that.  Essentially, fold angle into pi/4.
horiz = abs(horiz);
vert = abs(vert);
warning off MATLAB:divideByZero;
ratio = min(horiz,vert) ./ max(horiz,vert);
warning on MATLAB:divideByZero;

clear horiz vert;

% Split image into small values (set to 0) and HV versus HVbar areas.
find_below = find(ratio < ratio_threshold);
find_zeros = find(si <= rmin);

% Start generating HVbar image.  Zero out areas where SI is too small.
hvb = si;
hvb(find_zeros) = 0;

% Generate HV image.  Use HVbar image, so don't have to repeat the zeroing
% out of small SI values.  Then, zero out HV area.
hv = zeros(row_size,col_size);
hv(find_below) = hvb(find_below);

% Finnish generating HVbar image.  Zero out HVbar area.
hvb(find_below) = 0;

% if needed, reshape back into 3-D
if must_reshape == 1,
    si = reshape(si,row_size,col_size/time_size,time_size);
    hv = reshape(hv,row_size,col_size/time_size,time_size);
    hvb = reshape(hvb,row_size,col_size/time_size,time_size);
end

% take off invalid border around the edge.
    si = si(7:row_size-6,7:col_size/time_size-6,:);
    hv = hv(7:row_size-6,7:col_size/time_size-6,:);
    hvb = hvb(7:row_size-6,7:col_size/time_size-6,:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % alternate implementation, for debugging purposes.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Construct angle image
% a = abs(atan2(v_in, h_in));  % values from -pi to pi, folded 0 to pi
% temp = find(a > pi/2);
% a(temp) = pi - a(temp);  % folded 0 to pi/2
% 
% % Find locations of all pixels to zero in the hv image
% zeros_hv = find((si <= rmin) | ((a > theta) & (a < pi/2-theta)));
% 
% % Find locations of all pixels to zero in the hvb image
% zeros_hvb = find((si <= rmin) | (a <= theta) | (a >= pi/2-theta));
% 
% % Generate hv and hvb images
% hv_alt = si;
% hv_alt(zeros_hv) = 0;
% 
% hvb_alt = si;
% hvb_alt(zeros_hvb) = 0;
% 
% %%% check
% [i,j] = find(hv ~= hv_alt | hvb ~= hvb_alt);
% fprintf('%d points do not match:\n', size(i)); 
% % will say '0 do not, 1 do not' if 0 points mis-match.
% for cnt=1:size(i),
%     fprintf('(%d,%d) H=%f V=%f HV=%f aHV=%f HVb=%f aHVb=%f\n', ...
%         i(cnt),j(cnt),h_in(i(cnt),j(cnt)),v_in(i(cnt),j(cnt)),...
%         hv(i(cnt),j(cnt)),hv_alt(i(cnt),j(cnt)), ...
%         hvb(i(cnt),j(cnt)),hvb_alt(i(cnt),j(cnt)));
% end
% 
