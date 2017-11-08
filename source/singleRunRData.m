function [X,A] = singleRunRData(Subj)
%% Perform single run of go_psiicos for subject Subj
%
% FORMAT:
%    function [X,A] = singleRunRData(Subj)
% INPUTS: 
%   Subj            - subject name
% OUTPUTS:
%   X               - {Nactive_pairs x Time} recovered timecourses
%   A               - recovered indices of connected pairs
% __________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    path(path,'../input');
    [ConData, G] = PrepRealData(Subj);
    CT = CrossSpectralTimeseries(ConData{3}.Trials );
    [X,A] = go_psiicos(0.06, 500, G, CT, ConData{1}.CrossSpecTime, Subj);
end
