function [  ] = Second(subject,score,MCT )
% Method for post experiment sceening of subjects from ITU-R BT.1788
% Pearson and Spearman correlation
% This code is intended to be called from function 'screening.m'
%
% Margaret H. Pinson
% mpinson@its.bldrdoc.gov
% U.S. Department of Commerce, NTIA/ITS
  
     numberOfViewers = max(subject);
     numberOfPVS = length(subject)/numberOfViewers;
     total = 0;
      
     %pre allocate size for speed
     x = zeros(1,numberOfPVS);
     r1 = zeros(numberOfViewers,1);
     r2 = zeros(numberOfViewers,1);
     r = zeros(numberOfViewers,1);
    
     %calculate x_i = means score of all observers per PVS
     for i = 1:numberOfPVS,
        for j =0:numberOfViewers-1,
            total = total + score(i + j*numberOfPVS);   
        end
        x(i) = total/numberOfViewers;
        total = 0;
     end
     
     %calculate y = individual score of one observer per PVS
     %r1 = Pearson correlation
     %r2 = Spearman correlation
     %r = min of Pearson and Spearman correlation
     for i = 0:numberOfViewers-1,
        y = score(i*numberOfPVS +1:(i+1)*numberOfPVS);
        r1(i+1) = corr(x',y');
        r2(i+1) = corr(x',y', 'type','spearman');
        r(i+1) = min(r1(i+1),r2(i+1));
     end

     %average of the correlations of all observers
     mean_r = mean(r);
     %standard deviation of all observers
     std_r = std(r);

     %set reject threshold
     if(mean_r - std_r) > MCT
         reject_threshold = MCT;
     else
         reject_threshold = (mean_r-std_r);
     end 
     
     %decide which observers are discarded and print results
     for i = 1:numberOfViewers,
         if r(i) < reject_threshold
            disp(['Exclude viewers: ' num2str(i)]);
         end
     end

end

