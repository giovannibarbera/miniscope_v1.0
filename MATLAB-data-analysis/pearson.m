function [r] = pearson(X, Y)
%PEARSON calculate Pearson correlation between vector X and Y
%
%
%   giobarbera@neuralmapper.com

m_x = mean(X);
m_y = mean(Y);

num = sum((X-m_x).*(Y-m_y));
den = (sqrt(sum((X-m_x).^2)))*(sqrt(sum((Y-m_y).^2)))+0.00000000001;
r = num/den;
