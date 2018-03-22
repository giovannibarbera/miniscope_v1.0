function [FILTER] = gauss_filt(filt_size, sigma)
%GAUSS_FILT Create two dimensional Gauss filter of arbitrary size and
%standard deviation
%   giobarbera@neuralmapper.com


k = 1/(2*pi*sigma^2);
FILTER = zeros(filt_size, filt_size);

for m = 1:filt_size
    for n = 1:filt_size
       % FILTER(m,n) = (k*exp(-((m-1-floor(filt_size/2))^2+((n-1-floor(filt_size/2)))^2)/(2*sigma^2)));
       % FILTER(m,n) = (k*exp(-(((m-ceil(filt_size/2))/ceil(filt_size/(2*3)))^2+(((n-ceil(filt_size/2)))/ceil(filt_size/(2*3)))^2)/(2*sigma^2))); 
       FILTER(m,n) = (k*exp(-((m-1-floor(filt_size/2))^2+((n-1-floor(filt_size/2)))^2)/(2*sigma^2)));
    end
end

