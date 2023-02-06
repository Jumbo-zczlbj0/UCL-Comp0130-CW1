clc
clear
Define_Constants

% time position velocity (NED)(851,7)
GNSS_Pos_Vel = csvread('GNSS_Pos_Vel_NED.csv');
% time/L/lambda/v_N/v_E (851,5)
DR_Pos_Vel = csvread('single_DR.csv');

f = matfile("GNSS_P_matrix.mat");
GNSS_P_matrix = f.P_matrix; %(8,8,851)

time = DR_Pos_Vel(:,1);
num_time = length(time);
% x initilisation
x_0_pos = zeros(4,1);

% P initilisation (sigma_v/sigma_r)
sigma_v = 0.03+0.02+0.05;
sigma_r = 0; % TODO
L_0 = deg2rad(GNSS_Pos_Vel(1,2));
h_0 = GNSS_Pos_Vel(1,4);
[R_N,R_E] = Radii_of_curvature(L_0);
P_0_pos = [sigma_v^2, 0, 0, 0;
           0, sigma_v^2, 0, 0;
           0, 0, (sigma_r^2)/((R_N+h_0)^2), 0;
           0, 0, 0, (sigma_r^2)/((R_E+h_0)^2*(cos(L_0))^2)];

tau_s = 0.5;
% initilisation
% S_DR: coursework indicated (0.01m^2/s^3)
S_DR = 0.2;

% Corrected DR Solution (Degree)
DR_GNSS_Integration = zeros(num_time,5);
x_k_1_pos = x_0_pos; % rad
P_k_1_pos = P_0_pos;
% loop
for i=1:num_time
    if i == 1
        L_k_1 = L_0;
        h_k_1 = h_0;
        [R_N,R_E]= Radii_of_curvature(L_k_1);
    else
        L_k_1 = deg2rad(GNSS_Pos_Vel(i-1,2));
        h_k_1 = GNSS_Pos_Vel(i-1,4);
        [R_N,R_E]= Radii_of_curvature(L_k_1);
    end
    
    phi_k_1 = [1, 0, 0, 0;
               0, 1, 0, 0;
               tau_s/(R_N+h_k_1), 0, 1, 0;
               0, tau_s/((R_E+h_k_1)*cos(L_k_1)), 0, 1];
    
    Q_k_1 = [S_DR*tau_s, 0, (1/2)*((S_DR*tau_s^2)/(R_N+h_k_1)), 0;
             0, S_DR*tau_s, 0, (1/2)*((S_DR*tau_s^2)/((R_E+h_k_1)*cos(L_k_1)));
             (1/2)*((S_DR*tau_s^2)/(R_N+h_k_1)), 0, (1/3)*((S_DR*tau_s^3)/((R_N+h_k_1)^2)), 0;
             0, (1/2)*((S_DR*tau_s^2)/((R_E+h_k_1)*cos(L_k_1))), 0, (1/3)*((S_DR*tau_s^3)/((R_E+h_k_1)^2*(cos(L_k_1)^2)))];
    
    x_k_neg = phi_k_1*x_k_1_pos;
    
    P_k_neg = phi_k_1*P_k_1_pos*phi_k_1.'+Q_k_1;
    
    H_k = [ 0,  0, -1, 0;
            0,  0, 0, -1;
           -1,  0, 0,  0;
            0, -1, 0,  0;];
    
    L_k = deg2rad(GNSS_Pos_Vel(i,2));
    h_k = GNSS_Pos_Vel(i,4);
    [R_N,R_E]= Radii_of_curvature(L_k);
    
    % measurement error std from GNSS P
    GNSS_P = GNSS_P_matrix(:,:,i);
    sigma_G_L = GNSS_P(1,1);
    sigma_G_lambda = GNSS_P(2,2);
    sigma_G_vn = GNSS_P(4,4);
    sigma_G_ve = GNSS_P(5,5);
    R_k = [(sigma_G_L^2)/(R_N+h_k)^2, 0, 0, 0;
           0, (sigma_G_lambda^2)/((R_E+h_k)^2*(cos(L_k))^2), 0, 0;
           0, 0, sigma_G_vn^2, 0;
           0, 0, 0, sigma_G_ve^2];
    
    
    K_k = P_k_neg*H_k.'*inv(H_k*P_k_neg*H_k.'+R_k);
    
    temp = [deg2rad(GNSS_Pos_Vel(i,2)) - deg2rad(DR_Pos_Vel(i,2));
            deg2rad(GNSS_Pos_Vel(i,3)) - deg2rad(DR_Pos_Vel(i,3));
            GNSS_Pos_Vel(i,5) - DR_Pos_Vel(i,4);
            GNSS_Pos_Vel(i,6) - DR_Pos_Vel(i,5);];
    
    delta_z = temp - H_k*x_k_neg;
    
    x_k_pos = x_k_neg + K_k*delta_z;
    P_k_pos = (eye(4)-K_k*H_k)*P_k_neg;
    
    x_k_1_pos = x_k_pos;
    P_k_1_pos = P_k_pos;
    
    % DR solution minus difference
    x_k_C = [DR_Pos_Vel(i,2:3)*deg_to_rad, DR_Pos_Vel(i,4:5)]-[x_k_pos(3),x_k_pos(4),x_k_pos(1),x_k_pos(2)];
    DR_GNSS_Integration(i,:) = [time(i), x_k_C(1:2)*rad_to_deg, x_k_C(3:4)];
end
DR_GNSS_Integration(:,6)=DR_Pos_Vel(:,6)
writematrix(DR_GNSS_Integration,'zzz_DR_GNSS_output.csv')

%plot(DR_GNSS_Integration(:,2),DR_GNSS_Integration(:,3))
%plot(GNSS_Pos_Vel(:,2),GNSS_Pos_Vel(:,3))
%plot(DR_GNSS_Integration(:,2),DR_GNSS_Integration(:,3))
%plot(DR_Pos_Vel(:,2),DR_Pos_Vel(:,3))
%plot(DR_Pos_Vel(:,2),DR_Pos_Vel(:,3),GNSS_Pos_Vel(:,2),GNSS_Pos_Vel(:,3))
%plot(DR_GNSS_Integration(:,2),DR_GNSS_Integration(:,3),DR_Pos_Vel(:,2),DR_Pos_Vel(:,3),GNSS_Pos_Vel(:,2),GNSS_Pos_Vel(:,3))


