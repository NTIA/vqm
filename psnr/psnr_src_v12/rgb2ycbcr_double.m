function [ycbcr1,ycbcr2,ycbcr3] = rgb2ycbcr_double(one, two, three, four)
% RGB2YCBCR_DOUBLE
%   Convert image from RGB space into YCbCr space
% SYNTAX
%   [ycbcr] = rgb2ycbcr_double(rgb);
%   [y, cb, cr] = rgb2ycbcr_double(r,g,b);
%   [...] = rgb2ycbcr_double(...,'128');
% DESCRIPTION
%  Takes 'rgb'  -- an nr x nc x 3 RGB double precision image,
%  Converts 'rgb' into an YCbCr image, 'ycbcr'
%
%  Alternately, each image plane may be passed separately, in 'r', 'g' and
%  'b' input arguments.  In this case, the YCbCr image will be returned in
%  separate image planes, 'y', 'cb', and 'cr'.
%
%  Nominal input Y values are on [16,235] and nominal input Cb and Cr values 
%  are on [16,240].  RGB values are on [0,255].  This routine does not
%  round final RGB values to the nearest integer.
%  
%  When optional argument '128' is present, Cb and Cr values will be
%  returned on a range from -128 to 127.
%
%   Reference: 
%     Charles Poynton ColorFAQ.pdf (page 15), available from www.poynton.com.
%


if nargin == 1 || nargin == 2,
    [nr,nc,np] = size(one);
    if (np ~= 3)
        disp('Must have three image planes (R, G, B) for third dimension');
        return
    end

    % Transformation for each pixel is given by:
    % [Y Cb Cr]' = a0 + a1 * [R G B]'; RGB on [0, 255]
    a0 = [16;128;128];  % Offset
    a1 = [65.481 128.553 24.966; -37.797 -74.203 112; 112 -93.786 -18.214]/255;  % Matrix

    ycbcr = repmat(a0,1,nr*nc) + a1*reshape(one,nr*nc,np)';
    ycbcr = reshape(ycbcr',nr,nc,np);

    % Clip at 0 and 255
    ycbcr1 = max(0, min(ycbcr, 255));
    
    if nargin == 2 && strcmp(two,'128'),
        ycbcr1(:,:,2:3) = ycbcr1(:,:,2:3) - 128;
    end

elseif nargin == 3 || nargin == 4
    
    [nr,nc,np] = size(one);
    if (np ~= 1)
        error('arguments must be single image planes');
    end
    if size(one) ~= size(two) | size(two) ~= size(three),
        error('images must be the same size');
    end

    % Transformation for each pixel is given by:
    % [Y Cb Cr]' = a0 + a1 * [R G B]'; RGB on [0, 255]
    

    ycbcr1 = 16.0 + ((65.481/255.0)*one) + ((128.553/255.0)*two) + ((24.966/255.0)*three);
    ycbcr2 = 128.0 + ((-37.797/255.0)*one) + ((-74.203/255.0)*two) + ((112.0/255.0)*three);
    ycbcr3 = 128.0 + ((112.0/255.0)*one) + ((-93.786/255.0)*two) + ((-18.214/255.0)*three);
    
    % An alternate implementation of the above three lines is:
%     a0 = [16;128;128];  % Offset
%     a1 = [65.481 128.553 24.966; -37.797 -74.203 112; 112 -93.786 -18.214]/255;  % Matrix
%     temp = [ reshape(one,1,nr*nc) ; reshape(two,1,nr*nc); reshape(three,1,nr*nc)];
%     ycbcr = a1 * (temp);
%     % reshape & return
%     ycbcr1 = reshape(ycbcr(1,:),nr,nc) + 16;
%     ycbcr2 = reshape(ycbcr(2,:),nr,nc) + 128;
%     ycbcr3 = reshape(ycbcr(3,:),nr,nc) + 128;
    
    
%     % Clip at 0 and 255 -- Not needed, this cannot happen when RGB is
%     % within legal 0 to 255 range.
%     ycbcr1 = max(0, min(ycbcr1, 255));
%     ycbcr2 = max(0, min(ycbcr2, 255));
%     ycbcr3 = max(0, min(ycbcr3, 255));
    
    if nargin == 4 && strcmp(four,'128'),
        ycbcr2 = ycbcr2 - 128;
        ycbcr3 = ycbcr3 - 128;
    end
end

