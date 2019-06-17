function [aba] = ave_best_aligned(results_fuzzy_mse)
% AVE_BEST_ALIGNED (results_fuzzy_mse)
%  This function takes a results_fuzzy_mse matrix for one video clip,
%  which is obtained from results_fuzzy_mse{index}, and calculates the
%  average number of best aligned fields/frames with a Mean Squared Error
%  (MSE) that is less than or equal to the MSE of the best aligned frame.
%  This is obtained by comparing the MSE of the first row of the  
%  matrix with the MSE of the subsequent rows of the matrix (e.g., the
%  causal alignment option may pick an alignment that is not the minimum
%  MSE) and counting the number of fields/frames that have MSEs less than
%  or equal to the chosen alignment.  The results_fuzzy_mse cell array is
%  obtained from loading the VFD results file for the given test (e.g.,
%  vfd_test.mat).  The average number of best aligned fields/frames is
%  computed over all columns (e.g., output fields/frames).
%
%  For a perfectly aligned sequence where the chosen alignment has the
%  minimum MSE for all frames, this routine should return a 1.
%

[nrows,ncols] = size(results_fuzzy_mse);
chosen_mse = repmat(results_fuzzy_mse(1,:),nrows,1);  % The MSEs of the chosen alignments
error_mat = results_fuzzy_mse-chosen_mse;  % This will be less than or equal to zero for the best aligned frames
aba = length(find(error_mat<=0))/ncols;

end