function this_par = dll_vfd_clippar_loop_psnr(this_par)
% VFD_CLIPPAR_LOOP_PSNR
%  Computes the psnr_vfd parameter of the VQM_VFD model. The only input is
%  a feature generated by dll_ feature_loop_mse.  The returned par will be
%  the PSNR value calculated by the function.
% SYNTAX
%  par = dll_vfd_clippar_loop_psnr(this_par)
% DESCRIPTION
%  This function uses the feature generated by dll_feature_loop_mse to
%  compute the psnr_vfd described in NTIA TM-11-475.  This feature is
%  converted to PSNR using a peak signal level of 255.  PSNR is clipped at
%  an upper level of 48 dB.
%

peak_signal = 255;  % Defines the peak signal level for PSNR
psnr_limit = 48;  % limit on PSNR max value, otherwise PSNR could correlate poorly

% Convert RMSE to MSE of each block by squaring.  Then compute mean
% over all space and time (i.e., the MSE of the video sequence).  Then
% convert to PSNR.
this_par = this_par.^2;  % convert rmse to mse

% Calculate mean over ST
this_par = st_collapse('mean', this_par, '3D');

% Convert to PSNR
if (this_par ~= 0)  % Might be zero MSE in rare cases
    this_par = 10*log10(peak_signal*peak_signal/this_par);
    this_par = min(this_par, psnr_limit);
else
    this_par = psnr_limit;
end

