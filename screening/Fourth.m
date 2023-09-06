function [ ] = Fourth( subject,score )
% Method for post experiment sceening of subjects from ITU-R BT.500-12
% This code is intended to be called from function 'screening.m'
%
% Margaret H. Pinson
% mpinson@its.bldrdoc.gov
% U.S. Department of Commerce, NTIA/ITS 

    numberOfViewers = max(subject);
    numberOfPVS = length(subject)/numberOfViewers;
    
    m_4sum = 0;
    m_2sum = 0;
    
    %pre allocate size for speed
    s = zeros(1,numberOfViewers);
    beta_2 = zeros(1,numberOfPVS);
    mean_1 = zeros(1,numberOfPVS);
    std_1 = zeros(1,numberOfPVS);
    P = zeros(1,numberOfViewers);
    Q = zeros(1,numberOfViewers);

    %calculate mean and standard deviation of each PVS
    for i = 1:numberOfPVS,
        total = 0;
        for j =0:numberOfViewers-1,
            total = total + score(i + j*numberOfPVS); 
            s(j+1) = score(i+j*numberOfPVS);
        end
        mean_1(i) = total/numberOfViewers;
        std_1(i) = std(s);
    end
    
    %m_4 = fourth order moment
    %m_2 = second order moment
    %beta_2 = kurtosis coefficient
    for i = 1:numberOfPVS,
        for j =0:numberOfViewers-1, 
            m_4sum = m_4sum + power((score(i+j*numberOfPVS)- mean_1(i)),4);
            m_2sum = m_2sum + power((score(i+j*numberOfPVS)- mean_1(i)),4);
        end
        m_4 = m_4sum/numberOfViewers;
        m_2 = m_2sum/numberOfViewers;
        beta_2(i) = m_4/(power(m_2,2));
    end
    
    %determine if distribution is normal or not
    %Increment P and Q counters
    for i = 1:numberOfPVS,
        if beta_2(i)<=4 || beta_2(i) >=2
            for j = 0:numberOfViewers-1
                if score(i+j*numberOfPVS) >= (mean_1(i) + 2*std_1(i))
                    P(j+1) = P(j+1) +1;
                end
                if score(i+j*numberOfPVS) <= (mean_1(i) - 2*std_1(i))
                    Q(j+1) = Q(j+1) +1;
                end
            end         
        else
            for j = 0:numberOfViewers-1
                if score(i+j*numberOfPVS) >= (mean_1(i) + sqrt(20)*std_1(i))
                    P(j+1) = P(j+1) +1;
                end
                if score(i+j*numberOfPVS) <= (mean_1(i) - sqrt(20)*std_1(i))
                    Q(j+1) = Q(j+1) +1;
                end
            end  
        end
    end
    
    %decide which observers discarded and print results
    for i = 1:numberOfViewers
        if ((P(i)+Q(i))/numberOfPVS) > 0.05 && abs((P(i)-Q(i))/(P(i)+Q(i)))
            disp(['Exclude viewers: ' num2str(i)]);
        end
    end
    

end

