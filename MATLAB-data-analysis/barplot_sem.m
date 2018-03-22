function [barplot] = barplot_sem(DATA)
%BARPLOT_SEM Bar plot with mean and sem of input data
%
%   giobarbera@neuralmapper.com
%
% DATA(cell(M)): input data (M variables)

M = size(DATA,1);



DATA_MEAN = zeros(M,1);
DATA_SEM = zeros(M,1);


for k = 1:M
    curr_vect = DATA{k};
    curr_vect(isnan(curr_vect)) = [];
    DATA_MEAN(k) = mean(curr_vect);
    DATA_SEM(k) = std(curr_vect)/sqrt(length(curr_vect));
end

figure; barplot = bar([1:M],DATA_MEAN); hold on;
errorbar([1:M],DATA_MEAN,DATA_SEM,'Linestyle','none');
xlim([0 M+1]);