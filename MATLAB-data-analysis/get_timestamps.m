function [t time_str] = get_timestamps(filePath)
%GET_TIMESTAMPS Get time stamps and number of experiments in the specified
%folder
%
%   giobarbera@neuralmapper.com



curr_timestamp = '';
k = 0;

dir_content = dir(filePath);
prefix = '';
for kk = 1:length(dir_content)
    if strcmp(dir_content(kk).name(max(1,length(dir_content(kk).name)-3):end),'.txt') && ~isempty(strfind(dir_content(kk).name,'timestamp'))
        prefix = dir_content(kk).name;
        tmp = strfind(prefix,'timestamp');
        prefix = prefix(1:tmp(1)-2);
    end
end

if isempty(prefix)
    error('No NeuralMapper txt files found!');
end

curr_name = strcat(filePath,prefix,'_timestamp_exp_',num2str(k),'.txt');%'img_data_timestamp_exp_'

while (exist(curr_name, 'file') ~= 2)  && k < 10000
    % for k = 1:10
    k = k+1;
    % curr_name = strcat(filePath,prefix,'_exp_',num2str(k),'.txt');%'img_data_timestamp_exp_'
    curr_name = strcat(filePath,prefix,'_timestamp_exp_',num2str(k),'.txt');%'img_data_timestamp_exp_'
    if exist(curr_name, 'file') == 2
        curr_timestamp = curr_name;
        curr_img_data = strcat(filePath,prefix,'image_exp_',num2str(k),'.txt');%'img_data_exp_'
    end
end

delimiterIn = '\t';
A = importdata(curr_timestamp,delimiterIn);


% t = (A.data(:,2)-A.data(1,2))./1000;
t = A.data(:,2)./1000;


time_str = A.textdata{1}(end-24:end-1);




