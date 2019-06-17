function [image] = resample_image(image, v, h, varargin);
% RESAMPLE_IMAGE
%  stretch or shrink an image
% SYNTAX
%  [scaled_image] = resample_image(image, v, h);
%  [...] = resample_image(...'PropertyName',...);
% DESCRIPTION
%  This function applies horizontal scaling factor (h / 1000) 
%  and vertical scaling factor (v / 1000).  The returned image,
%  'scaled_image', will be of the same size as the input image.
%
%  The following optional properties may be requested.  Fast, Linear, and
%  Quadratic are mutually exclusive (i.e., only one of these may
%  be selected):
%
%   'Fast'  When this option is selected, the function will use a very fast
%           but significantly less accurate resampling algorithm.  The
%           nearest neighbor pixel value will be used (1-point).  This 
%           approach appears to be sufficient for color planes (Cb and Cr)
%           but not luminance (Y). 
%
%   'Linear' When this option is selected, the function will use linear
%           interpolation (2-point).
%
%   'Quadratic' When this option is selected, the function will use
%           quadratic interpolation (3-point).  This is the default.
%
%   'Interlace' The image is interlaced, so vertical scaling needs to be
%           performed on each field separately. 


do_fast = 0;
do_linear = 0;
do_quadratic = 1;
do_interlace = 0;
is_type = 'quadratic';

cnt = 1;
while cnt <= nargin - 3,
    if strcmpi(varargin(cnt),'fast') == 1,
        do_fast = 1;
        do_linear = 0;
        do_quadratic = 0;
        is_type = 'fast';
        cnt = cnt + 1;
    elseif strcmpi(varargin(cnt),'linear') == 1,
        do_fast = 0;
        do_linear = 1;
        do_quadratic = 0;
        is_type = 'linear';
        cnt = cnt + 1;
    elseif strcmpi(varargin(cnt),'quadratic') == 1,
        do_fast = 0;
        do_linear = 0;
        do_quadratic = 1;
        is_type = 'quadratic';
        cnt = cnt + 1;
    elseif strcmpi(varargin(cnt),'interlace') == 1,
        do_interlace = 1;
        cnt = cnt + 1;
    else
        error('optional argument not recognized.');

    end
end

% vertical scaling on interlaced images must be done on fields.
% handle this here, by recursing (calling this function on each field).
if do_interlace && v ~= 1000,
    % split into fields
    [image1, image2] = split_into_fields(image);
    
    % call this routine for each field
    [image1] = resample_image(image1, v, h, is_type);
    [image2] = resample_image(image2, v, h, is_type);
    
    % join into frames
    image = join_into_frames(image1, image2);
    return;
end

% split into frames, if needed, and run each separately
if ndims(image) == 3,
    [rows,cols,time] = size(image);
    if time > 1,
        for i=1:time,
            image(:,:,i) = resample_image(image(:,:,i), v, h, is_type);
        end
    end
    return;
elseif ndims(image) > 3,
    error('Function resample_image cannot work on arrays of images with 4 or more dimensions');
end


if do_linear,
    
    [rows,cols] = size(image);

    if v ~= 1000,
        % scale coordinates by vertical factor
        offset = (1:rows)' * 1000 / v + (rows/2 - (1000/v) * (rows/2));
        
        % find pixels invalidated by this operation
        invalid = find(offset < 1 | offset > rows);

        % if exceed boundary of image, clip at edge.
        offset = max(min(offset, rows), 1);

        % find closest pixels and distance (alpha)
        offset_before = floor(offset);
        offset_after = ceil(offset);
        alpha = offset - offset_before;

        % change from vectors to matrixes, to apply this to the entire image
        alpha = repmat(alpha, 1, cols);

        % apply linear scaling
        image = (1 - alpha) .* image(offset_before, :) + alpha .* image(offset_after, :);
        
        % fill invalid portion with zeros
        image(invalid, :) = 0;
    end
    
    if h ~= 1000,
        % scale coordinates by vertical factor
        offset = (1:cols) * 1000 / h + (cols/2 - (1000/h) * (cols/2));

        % find pixels invalidated by this operation
        invalid = find(offset < 1 | offset > cols);

        % if exceed boundary of image, clip at edge.
        offset = max(min(offset, cols), 1);

        % find closest pixels and distance (alpha)
        offset_before = floor(offset);
        offset_after = ceil(offset);
        alpha = offset - offset_before;

        % change from vectors to matrixes, to apply this to the entire image
        alpha = repmat(alpha, rows, 1);

        % apply linear scaling
        image = (1 - alpha) .* image(:, offset_before) + alpha .* image(:, offset_after);

        % fill invalid portion with zeros
        image(:, invalid) = 0;
    end
    

elseif do_quadratic,
    
    [rows,cols] = size(image);

    if v ~= 1000,
        % scale coordinates by vertical factor
        offset = (1:rows)' * 1000 / v + (rows/2 - (1000/v) * (rows/2));
        
        % find pixels invalidated by this operation
        invalid = find(offset < 1 | offset > rows);

        % find closest pixels and distance (alpha)
        y1 = round(offset);
        y0 = y1 - 1;
        y2 = y1 + 1;
        
        % find distance
        d = offset - y0;

        % change d from vector to matrix, to apply this to the entire image
        d = repmat(d, 1, cols);
        d2 = d.^2;
        
        % if exceed boundary of image, clip at edge (i.e., take last pixel value).
        y0 = max(min(y0, rows), 1);
        y1 = max(min(y1, rows), 1);
        y2 = max(min(y2, rows), 1);
        
        % apply to image
        image_y0 = image(y0,:);
        image_y1 = image(y1,:);
        image_y2 = image(y2,:) ./ 2;

        % calculate weights a, b & c
        % where c = image_y0;
        % where a = y(2) / 2 - y(1) + y(0) / 2
        % where b = -y(2) / 2 + 2 * y(1) - 3/2 * y(0)
        a = image_y2 - image_y1 + image_y0 ./ 2;
        b = -image_y2 + 2 .* image_y1 - (3/2) .* image_y0;
        
        % apply quadratic scaling: a .* d2 + b .* d + c
        image = a .*d2 + b .* d + image_y0;
        
        % fill invalid portion with zeros
        image(invalid, :) = 0;
    end
    
    if h ~= 1000,
        % scale coordinates by vertical factor
        offset = (1:cols)' * 1000 / h + (cols/2 - (1000/h) * (cols/2));
        
        % find pixels invalidated by this operation
        invalid = find(offset < 1 | offset > cols);

        % find closest pixels and distance (alpha)
        y1 = round(offset);
        y0 = y1 - 1;
        y2 = y1 + 1;
        
        % find distance
        d = offset - y0;

        % change d from vector to matrix, to apply this to the entire image
        d = repmat(d', rows,1);
        d2 = d.^2;
        
        % if exceed boundary of image, clip at edge (i.e., take last pixel value).
        y0 = max(min(y0, cols), 1);
        y1 = max(min(y1, cols), 1);
        y2 = max(min(y2, cols), 1);
        
        % apply to image
        image_y0 = image(:,y0);
        image_y1 = image(:,y1);
        image_y2 = image(:,y2) ./ 2;

        % calculate weights a, b & c
        % where c = image_y0;
        % where a = y(2) / 2 - y(1) + y(0) / 2
        % where b = -y(2) / 2 + 2 * y(1) - 3/2 * y(0)
        a = image_y2 - image_y1 + image_y0 ./ 2;
        b = -image_y2 + 2 .* image_y1 - (3/2) .* image_y0;
        
        % apply quadratic scaling: a .* d2 + b .* d + c
        image = a .*d2 + b .* d + image_y0;
        
        % fill invalid portion with zeros
        image(:,invalid) = 0;
    end


elseif do_fast,
    
    [rows,cols] = size(image);

    if v ~= 1000,
        % scale coordinates by vertical factor
        offset = (1:rows)' * 1000 / v + (rows/2 - (1000/v) * (rows/2));
        
        % find pixels invalidated by this operation
        invalid = find(offset < 1 | offset > rows);

        % if exceed boundary of image, clip at edge.
        offset = max(min(offset, rows), 1);

        % find closest pixels
        offset = round(offset);

        % apply linear scaling
        image = image(offset, :);
        
        % fill invalid portion with zeros
        image(invalid, :) = 0;
    end
    
    if h ~= 1000,
        % scale coordinates by vertical factor
        offset = (1:cols) * 1000 / h + (cols/2 - (1000/h) * (cols/2));

        % find pixels invalidated by this operation
        invalid = find(offset < 1 | offset > cols);

        % if exceed boundary of image, clip at edge.
        offset = max(min(offset, cols), 1);

        % find closest pixels
        offset = round(offset);

        % apply linear scaling
        image = image(:, offset);

        % fill invalid portion with zeros
        image(:, invalid) = 0;
    end
    

else
    error('argument not recognized. "resample" no longer available');
end
