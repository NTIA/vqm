function [si, hv, hvb] = filter_si_hv_adapt(y, filter_size, extra, varargin)
% FILTER_SI_HV_ADAPT
%
%  Filters Y with the NxN gradient (c=2) filters described in SPIE 1999
%  paper.  Like funciton 'filter_si_hv', but can choose other than 13x13.
%
% SYNTAX
%
%  [SI] = filter_si_hv_adapt(Y, N, EXTRA)
%  [SI] = filter_si_hv_adapt(Y, N, EXTRA, rmin, theta)
%  [SI, HV, HVB] = filter_si_hv_adapt(...)
%
% DESCRIPTION
%
%  [SI] = filter_si_hv_adapt(Y, N, EXTRA)  Perceptually fiters liminence image Y using 
%  the NxN Horizontal and Vertical gradient filters in a fashion similar 
%  to the sobel filter.  Discard the filtered border, of width EXTRA, where
%  ( EXTRA >= floor(N/2) ).  Filter size, N, must be an odd number.
%
%  If Y is a 3 dimensional matrix, Y will be presumed to contain multiple
%  images as follows:  (row, col, time).  No execution time penalties occur.
%
%  [SI] = filter_si_hv(Y, N, EXTRA, rmin, theta)  allows the user to over ride
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
%  border of EXTRA pixels around the edge of each image is invalid.
%

if mod(filter_size,2) == 0,
    error('SI filter size must be an odd number');
end

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
if row_size < filter_size | col_size < filter_size,
    error(sprintf('Function ''filter_si_hv'' requires images to be at least %dx%d', filter_size, filter_size));
end

% compute angle as a ratio of HV and HVbar.
ratio_threshold = tan(theta);

%  The weights for a single row of the H filter 
%  is given by: w(x) = k*(x/c)*exp{-(1/2)*(x/c)^2}, where x = {-6, -5, ..., 5, 6}, 
%  and k is a normalization constant chosen such that this filter produces the same 
%  amplitude response on an H V edge as the Sobel filter.

%  Generate the filter_size long filter mask, in one dimension.
filter_mask = filter_h(filter_size);
filter_mask = filter_mask(1,:);

%  Convolve mask with y in horizontal & vertical direction.
%  do two 1xfilter_size convolutions instead of one filter_sizexfilter_size, for speed.
horiz = conv2(y,filter_mask,'same');
horiz = conv2(horiz,ones(filter_size,1),'same');

vert = conv2(y,filter_mask','same');
vert = conv2(vert,ones(1,filter_size),'same');

% for debugging, comment in the below lines.
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
    si = si((extra+1):row_size-extra,(extra+1):col_size/time_size-extra,:);
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
    si = si((extra+1):row_size-extra,(extra+1):col_size/time_size-extra,:);
    hv = hv((extra+1):row_size-extra,(extra+1):col_size/time_size-extra,:);
    hvb = hvb((extra+1):row_size-extra,(extra+1):col_size/time_size-extra,:);

    
function h = filter_h(l)
%  H = FILTER_H(L)
%  Returns a horizontal bandpass filter H (like Figure 27 left, in NTIA
%  Report 02-392) given the filter width L as input.  L must be an odd
%  positive integer greater than or equal to 3.  The vertical bandpass 
%  filter is the transpose of this H filter.
%
%  The weights for a single row of the H filter are given by: 
%  w(x) = k*(x/c)*exp{-(1/2)*(x/c)^2}, where x = {-m, -(m-1), ..., 0, ..., m-1, m}, 
%  m = (L-1)/2, and k is a normalization constant chosen such that the filter produces 
%  the same amplitude response on a vertical edge as the Sobel filter.
%
%  For the reference 13 x 13 filter, L = 13 and C = 2.  For other filter sizes, C
%  is scaled appropriately, e.g.,
%  for the scaled 9 x 9 filter, L = 9, C = 4/3.
%  for the corresponding 5 x 5 filter, L = 5, C = 2/3.
%

%  Generate the H filter masks
m = (l-1)/2;  % half filter width
c = 2.0*m/6;  % c = 2 for m = 6, the reference filter
x = -m:m;
h1 = (x/c).*exp(-(1/2)*(x/c).^2);
h = repmat(h1,l,1);

%  Normalize for Sobel filter energy
h = h./(sum(sum(abs(h)))/8);

    
    