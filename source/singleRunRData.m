% and perform single run of go_psiicos
% function [X,A] = singleRunRData(Subj)
% X - recovered timecourses
% A - recovered interaction pairs

function [X,A] = singleRunRData(Subj)
	path(path,'../input');
	[ConData, G] = PrepRealData(Subj);
	CT = CrossSpectralTimeseries(ConData{2}.Trials );
	[X,A] = go_psiicos(0.08, 25, G, CT, ConData{1}.CrossSpecTime);
