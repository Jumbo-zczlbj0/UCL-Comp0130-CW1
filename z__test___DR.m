clc
clear
format long;
% z_initialize;
Define_Constants;
% Dead reckoning
csv_1 = 'Dead_reckoning.csv';
speed_heading = csvread(csv_1);

time = speed_heading(:,1);
% corrected heading
csv_2 = 'corrected_heading.csv';
corrected_heading = csvread(csv_2);

% initial output
single_DR_output = zeros(850,6);

% latitude (time=0)
L_k_1=deg_to_rad*51.509254461345272;
% longitude (time=0)
lambda_k_1=deg_to_rad*-0.161045484902720;
% Heading when time = 0
Heading_DR_k_1=speed_heading(1,7);

for i = 1:850
    % Rear wheel Velocity
    V_er=(speed_heading(i+1,4)+speed_heading(i+1,5))/2;

    % Speed of the North and East
    % V_Nk=(cos(phi)-0.5*dot_phi*delta_t*sin(phi)*Rear_wheel_Velocity)*V_er*delta_t+
    % (cos(phi)*l_bry+sin(phi)*l_brx)*dot_phi*delta_t
    V_Nk=(cos(corrected_heading(i,2)*deg_to_rad)-0.5*speed_heading(i,6)* ...
        0.5*sin(corrected_heading(i,2)*deg_to_rad))*V_er*0.5 + (cos(corrected_heading(i,2)*deg_to_rad)* ...
        0.2+sin(corrected_heading(i,2)*deg_to_rad)*0.2)*0.5*speed_heading(i,6);

    % V_Ek=(sin(phi)+0.5*dot_phi*delta_t*cos(phi)*Rear_wheel_Velocity)*V_er*delta_t+
    % (sin(phi)*l_bry-cos(phi)*l_brx)*dot_phi*delta_t
    V_Ek=(sin(corrected_heading(i,2)*deg_to_rad)+0.5*speed_heading(i,6)* ...
        0.5*cos(corrected_heading(i,2)*deg_to_rad))*V_er*0.5 + (sin(corrected_heading(i,2)*deg_to_rad)* ...
        0.2-cos(corrected_heading(i,2)*deg_to_rad)*0.2)*0.5*speed_heading(i,6);

    % get R_N R_E
    [R_N,R_E]= Radii_of_curvature(L_k_1);
    
    % height
    h =38.825814539333805;

    % Calculating latitude and longitude
    % L_k_1 = L_k +(V_Nk/(R_N+h))
    L_k=L_k_1+(V_Nk/(R_N+h));
    % lambda_k_1 = lambda_k +(V_Ek/((R_E+h)*cos(L_k_1)))
    lambda_k=lambda_k_1+(V_Ek/((R_E+h)*cos(L_k_1)));

    % update
    L_k_1 = L_k;
    lambda_k_1 = lambda_k;
    
    single_DR_output(i,1)=time(i,1);
    single_DR_output(i,2)=rad_to_deg*L_k;
    single_DR_output(i,3)=rad_to_deg*lambda_k;
    single_DR_output(i,4)=V_Nk;
    single_DR_output(i,5)=V_Ek;
    single_DR_output(i,6)=corrected_heading(i,2);  

end

writematrix(single_DR_output,'single_DR.csv')

csv_3 = 'GNSS_Pos_Vel_NED.csv';
GNSS_output = csvread(csv_3);


%plot(GNSS_output(:,2),GNSS_output(:,3))
plot(single_DR_output(:,2), single_DR_output(:,3))
%plot(output_car(:,2), output_car(:,3),GNSS_output(:,2),GNSS_output(:,3))