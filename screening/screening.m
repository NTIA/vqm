function [] = screening(subject, score, hrc) 
% SYNTAX 
%     screening(subject, score)
%     screening(subject, score, hrc)
% SEMANTICS
% This function takes three input vectors, each with the same length. Each
% element of the vectors contains one subjective rating, as follows: 
% 'subject' contains the ID number associated with each score
% 'score' contains one raw subjective score associated 
% 'hrc' (optional) contains the ID number associated with each score's
% hypothetical reference circuit (HRC). This variable may be left out if
% the algorithm chosen does not use it.
%
% User will be prompted for the techique to be used and the rating scale.
% Subjects recommended for elimination will be printed.
% if no subjects should be discarded, then nothing will print.
%
% Assumptions:
% 1) The clips are listed in the same order for each subject.
% 2) The subjects are numbered [1..N]
%
% Margaret H. Pinson
% mpinson@its.bldrdoc.gov
% U.S. Department of Commerce, NTIA/ITS

prompt = 'Pick the screening method and enter the corresponding number [1..5]:\n(1) VQEG HDTV Annex I\n(2) ITU-R BT.1788\n(3) VQEG MM Annex VI\n(4) ITU-R Rec. BT.500 Annex 2\n(5) compare all methods\n: ';
response = input(prompt); 

%Method 1 Pearson correlation
if response == 1
    disp('VQEG HDTV Annex I screening: ');
    First(subject,score);
    
%Method 2 Pearson and Spearman correlation
elseif response == 2 
    %Get MCT value and pass to function
    while(1)
         prompt1 = 'Pick which method for second function is being used and enter corresponding number: SAMVIQ(1), DSCQS(1), SS(2), DSIS(2),or none of the above(3): ';
         method = input(prompt1);
         if method == 1
             MCT = 0.85;
             break
         elseif method == 2
             MCT = 0.7;
             break
         elseif method == 3
             prompt2 = 'Enter threshold to be used: ';
             MCT = input(prompt2);
             break
         else
             disp('Incorrect input')
         end
    end
     disp('ITU-R Rec. BT.1788 screening: ');
    Second(subject,score,MCT);
     
%Method 3 Pearson and Pearson per HRC
elseif response == 3
    disp('VQEG MM Annex VI screening: ');
    Third( subject,score,hrc );
    
%Method 4 Kurtosis Coefficient    
elseif response == 4
    disp('ITU-R Rec. BT.500 screening: ');
    Fourth( subject,score );
   
%comparison of all the screening methods    
elseif response == 5
    
    %Get MCT value for second method
     while(1)
         prompt1 = 'Pick which method for second function is being used and enter corresponding number: SAMVIQ(1), DSCQS(1), SS(2), DSIS(2),or none of the above(3): ';
         method = input(prompt1);
         if method == 1
             MCT = 0.85;
             break
         elseif method == 2
             MCT = 0.7;
             break
         elseif method == 3
             prompt2 = 'Enter threshold to be used: ';
             MCT = input(prompt2);
             break
         else
             disp('Incorrect input')
         end
     end
     
    %output results from all four methods
    disp('Results from all four methods');
    
    fprintf('\n');
    disp('VQEG HDTV Annex I screening: ');
    First(subject,score);
    fprintf('\n');
    disp('ITU-R Rec. BT.1788 screening: ');
    Second(subject,score,MCT);
    fprintf('\n');
    disp('VQEG MM Annex VI screening: ');
    Third(subject,score,hrc);
    fprintf('\n');
    disp('ITU-R Rec. BT.500 screening: ');
    Fourth(subject,score);
     
else
    %have user input a valid option if input is invalid
    disp('Pick a valid option');
    screening(subject,score,hrc);
end

   