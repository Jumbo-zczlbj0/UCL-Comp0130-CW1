clc
clear
%DR_GNSS_Integration = csvread('zzz_DR_GNSS_output.csv');
aaa = readmatrix('zzz_Sensor_output.csv');
epoch_output = csvread('zzz_single_GNSS_output.csv');




% plot(DR_GNSS_Integration(:,2),DR_GNSS_Integration(:,3))

% % sensor
plot(aaa(:,2), aaa(:,3))

% % single
% plot(epoch_output(:,2),epoch_output(:,3))

%plot(DR_GNSS_Integration(:,2),DR_GNSS_Integration(:,3),'r',L_lambda_v(:,2), L_lambda_v(:,3),'b','LineWidth', 3)

%plot(DR_GNSS_Integration(:,2),DR_GNSS_Integration(:,3),L_lambda_v(:,2), L_lambda_v(:,3),epoch_output(:,2),epoch_output(:,3))