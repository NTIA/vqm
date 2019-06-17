function [rgb1, rgb2, rgb3] = ycbcr2rgb_double(one, two, three, four)
% YCBCR2RGB_DOUBLE
%   Convert image from YCbCr space into RGB space
% SYNTAX
%   [rgb] = ycbcr2rgb_double(ycbcr);
%   [r,g,b] = ycbcr2rgb_double(y, cb, cr);
%   [...] = ycbcr2rgb_double(...,'128');
% DESCRIPTION
%  Takes 'ycbcr'  -- an nr x nc x 3 YCbCr double precision image,
%  Converts 'ycbcr' into an RGB image, 'rgb'
%
%  Alternately, each image plane may be passed separately, in 'y', 'cb' and
%  'cr' input arguments.  In this case, the RGB image will be returned in
%  separate image planes, 'r', 'g', and 'b'.
%
%  Nominal input Y values are on [16,235] and nominal input Cb and Cr values 
%  are on [16,240].  RGB values are on [0,255].  This routine does not
%  round final RGB values to the nearest integer.
%  
%  When optional argument '128' is present, Cb and Cr values will be
%  presumed to be on a range from -128 to 127.
%
%   Reference: 
%     Charles Poynton ColorFAQ.pdf (page 15), available from www.poynton.com.
%

if nargin == 1 | nargin == 2,
    
    if nargin == 2 & strcmp(two,'128'),
        one(:,:,2:3) = one(:,:,2:3) + 128;
    end
    
    [nr,nc,np] = size(one);
    if (np ~= 3)
        disp('Must have three image planes (Y, Cb, Cr) for third dimension');
        return
    end
    
    % Transformation for each pixel is given by:
    % [Y Cb Cr]' = a0 + a1 * [R G B]'; RGB on [0, 255]
    a0 = [16;128;128];  % Offset
    a1 = [65.481 128.553 24.966; -37.797 -74.203 112; 112 -93.786 -18.214]/255;  % Matrix

    rgb = inv(a1) * (reshape(one,nr*nc,np)' - repmat(a0,1,nr*nc));
    rgb = reshape(rgb',nr,nc,np);

    % Clip at 0 and 255
    rgb1 = max(0, min(rgb, 255));
    
    
elseif nargin == 3 | nargin == 4,
    
    if nargin == 4 & strcmp(four,'128'),
        two = two + 128;
        three = three + 128;
    end
    
    % error checking
    
    [nr,nc,np] = size(one);
    if (np ~= 1)
        error('arguments must be single image planes');
    end
    if size(one) ~= size(two) | size(two) ~= size(three),
        error('images must be the same size');
    end
    
    % Transformation for each pixel is given by:
    % [Y Cb Cr]' = a0 + a1 * [R G B]'; RGB on [0, 255]
%     a0 = [16;128;128];  % Offset
    a1 = [65.481 128.553 24.966; -37.797 -74.203 112; 112 -93.786 -18.214]/255;  % Matrix

    temp = [ reshape(one - 16,1,nr*nc) ; reshape(two - 128,1,nr*nc); reshape(three - 128,1,nr*nc)];
    rgb = inv(a1) * (temp);
    
    % Clip at 0 and 255
    rgb = max(0, min(rgb, 255));

    % reshape & return
    rgb1 = reshape(rgb(1,:),nr,nc);
    rgb2 = reshape(rgb(2,:),nr,nc);
    rgb3 = reshape(rgb(3,:),nr,nc);
end

    





