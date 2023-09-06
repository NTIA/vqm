function [ ] = Third( subject,score,hrc )
% Method for post experiment sceening of subjects from VQEG MM Testplan Annex VI
% Pearson and Pearson per HRC
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
     
     
     %calculate x_i = means score of all observers per PVS
     for i = 1:numberOfPVS,
        for j =0:numberOfViewers-1,
            total = total + score(i + j*numberOfPVS);   
        end
        x(i) = total/numberOfViewers;
        total = 0;
     end
     
     %calculate y = individual score of one viewer for corresponding PVSs
     %r_1 = Pearson correlation per observer
     for i = 0:numberOfViewers-1,
        y = score(i*numberOfPVS +1:(i+1)*numberOfPVS);
        r1(i+1) = corr(x',y'); 
     end
    
     %remove \r from end of hrc's
     for i=1:length(hrc),
        tmp = hrc{i};
        want = strfind(tmp,'\r');
        if want,
            tmp = tmp(1:length(tmp)-2);
            hrc{i} = tmp;
        end
     end
     HRCS = unique(hrc);
     
     %pre allocate for speed
     hrc_values=zeros(length(HRCS),1);
     hrc_values2=zeros(length(HRCS),1);
     x2=zeros(1,length(HRCS));
     y2=zeros(1,length(HRCS));

     %sort scores of all viewers by HRC 
     for i=1:length(score)
         for j = 1:length(HRCS)
             if strcmp(hrc(i), HRCS(j))
                 hrc_values(j,i+1) = score(i);
                 hrc_values(j,1) = hrc_values(j,1) + 1;
             end
         end
     end
     
     %find average score of each HRC across all PVSs
     for j = 1:length(HRCS)
         x2(j)=sum(hrc_values(j,2:length(score)))/hrc_values(j,1);
     end
     

     for i = 0:numberOfViewers-1,
         l = 2; %variable to start in second column since first is used as counter
         
         %sort scores of one viewer by HRC 
         for k = i*numberOfPVS+1:(i+1)*numberOfPVS
              for j = 1:length(HRCS)
                 if strcmp(hrc(k), HRCS(j))
                     hrc_values2(j,l) = score(k);
                     %keep track of number of PVS in certain HRC, used in
                     %average
                     hrc_values2(j,1) = hrc_values2(j,1) + 1;
                 end
              end
              l = l+1;
         end
         
         %y2 = individual MOS of one viewer for corresponding HRC
         %r2 = Pearson correlation per HRC
         for j = 1:length(HRCS)
             y2(j)=sum(hrc_values2(j,2:length(hrc_values2)))/hrc_values2(j,1);
         end
         r2(i+1) = corr(x2',y2'); 
         hrc_values2=0;
       
     end
     
     %decide which observers discarded and print results
     for i = 1:numberOfViewers,
         if ((r1(i) < 0.75) && (r2(i) < 0.8))
            disp(['Exclude viewers: ' num2str(i)]);
         end
     end

end

