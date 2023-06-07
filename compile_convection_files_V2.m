%% Compile convection files
clear all
close all
%% Find compiled file

fullfilename=uigetfile('*');

k = strfind(fullfilename,'_');

filename=fullfilename(1:k(end-2)-1);

%% Read in data from compiled file
readfile=lvm_import(fullfilename,0);
grab_data=readfile.Segment1.data(:,:);

prompt1={'Enter number of files'};
name1 = 'File number';
defaultans1 = {'2'};
options.Interpreter = 'tex';
answer1 = inputdlg(prompt1,name1,[1 40],defaultans1,options);

if str2double(answer1{1})~=1
for i=2:str2double(answer1{1})
    DialogTitle=strcat('Select the file:  ',filename,'...',num2str(i));
    nextfilename=uigetfile('*',DialogTitle);
    %datatmp=dlmread(nextfilename,'\t',23,0);
    readfile_tmp=lvm_import(nextfilename,0);
    grab_data_tmp=readfile_tmp.Segment1.data(:,:);
    grab_data=[grab_data;grab_data_tmp];
    indie=i; 
    
    clear nextfilename
    clear grab_data_tmp
end
end

figure(1)
plot(grab_data(:,1),grab_data(:,2:13))

time=grab_data(:,1);
time_start=grab_data(1,1);
time_end=grab_data(end,1)
delta_t=grab_data(3,1)-grab_data(2,1);
samp_freq=1/delta_t;

prompt2={'Enter start time'};
name2 = 'Start Time';
defaultans2 = {'1000'};
start_time = inputdlg(prompt2,name2,[1 40],defaultans2,options);
start_time_index = floor(str2double(start_time{1})*samp_freq + 1);

prompt3={'Enter end time'};
name3 = 'End Time';
defaultans3 = {num2str(grab_data(end:1))};
end_time = inputdlg(prompt3,name3,[1 40],defaultans3,options);
end_time_index = floor(str2double(end_time{1})*samp_freq + 1);


% if str2num(end_time{1})==1e9
%     end_time_index=char('end');
% else
%     end_time_index = str2num(end_time{1})*10 + 1;
% end


% figure(2)
plot(grab_data(start_time_index:end_time_index),grab_data(start_time_index:end_time_index,2:7))

datasave=grab_data(start_time_index:end_time_index,:);


writefile=strcat(filename,'.txt');
% cd ..

dlmwrite(writefile,datasave,'Delimiter','\t','precision','%.6f');

% pause(1);
% 
% datacheck=dlmread(writefile,'\t',0,0);

% figure(3)
% plot(grab_data(:,1),grab_data(:,2:7),'-')
% hold on
% hold all
% plot(datacheck(:,1),datacheck(:,2:7),'o')

load train
sound(y,Fs)

    