function [  ] = First( subject,score )
% Method for post experiment sceening of subjects from VQEG HDTV Annex I 
% Pearson correlation applied per-clip
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
    r = zeros(numberOfViewers,1);
    
    %calculate x_i = MOS of all viewers per PVS
    for i = 1:numberOfPVS,
        for j =0:numberOfViewers-1,
            total = total + score(i + j*numberOfPVS);   
        end
        x(i) = total/numberOfViewers;
        total = 0;
    end
    
    %calculate y = individual score of one viewer for corresponding PVSs
    %r = Pearson correlation per observer
    for i = 0:numberOfViewers-1,
        y = score(i*numberOfPVS +1:(i+1)*numberOfPVS);
        r(i+1) = corr(x',y');
    end
    
    %decide which observers discarded and print results
    for i = 1:numberOfViewers,
        if r(i) < 0.75
            disp(['Exclude viewers: ' num2str(i)]);
        end
    end
    
end

