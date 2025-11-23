% Created by Alex Kim

clear all       % clears the workspace
clc             % clears the command window
close all       % closes all figures

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Constants

Data = load("flight_data4_7B_Space.mat");
GPS_sat_loc = Data.GPS_sat_loc;

time = Data.time_data;

Calibration0 = load("Calibration_0_7B_Space.mat");
Calibration1 = load("Calibration_1_7B_Space.mat");
Calibration2 = load("Calibration_2_7B_Space.mat");
Calibration3 = load("Calibration_3_7B_Space.mat");
Calibration4 = load("Calibration_4_7B_Space.mat");

pseudo0 = Calibration0.pseudo_range;
pseudo1 = Calibration1.pseudo_range;
pseudo2 = Calibration2.pseudo_range;
pseudo3 = Calibration3.pseudo_range;
pseudo4 = Calibration4.pseudo_range;

PDOP = zeros(1,length(Data.pseudo_range));
TDOP = zeros(1,length(Data.pseudo_range));
NDOP = zeros(1,length(Data.pseudo_range));
EDOP = zeros(1,length(Data.pseudo_range));
DDOP = zeros(1,length(Data.pseudo_range));

% Calibration Points
Point0 = [0 0 0];
Point1 = [0 5 0];
Point2 = [-5 5 0]; % Point 2 was an outlier
Point3 = [-5 0 0];
Point4 = [-10 0 0];

% Takes the std of the calibration pseudo ranges
for i=1:6
    stds0(i) = std(pseudo0(:,i));
    stds1(i) = std(pseudo1(:,i));
    stds2(i) = std(pseudo2(:,i));
    stds3(i) = std(pseudo3(:,i));
    stds4(i) = std(pseudo4(:,i));
end

% Mean standard deviation of the pseudo ranges (sigma_p)
std_total = mean([stds0; stds1; stds3; stds4], 'all');

% Position of the pseudolites (TO BE ENTERED)

r_S1 = [GPS_sat_loc(1,1);GPS_sat_loc(2,1);GPS_sat_loc(3,1)];
r_S2 = [GPS_sat_loc(1,2);GPS_sat_loc(2,2);GPS_sat_loc(3,2)];
r_S3 = [GPS_sat_loc(1,3);GPS_sat_loc(2,3);GPS_sat_loc(3,3)];
r_S4 = [GPS_sat_loc(1,4);GPS_sat_loc(2,4);GPS_sat_loc(3,4)];
r_S5 = [GPS_sat_loc(1,5);GPS_sat_loc(2,5);GPS_sat_loc(3,5)];
r_S6 = [GPS_sat_loc(1,6);GPS_sat_loc(2,6);GPS_sat_loc(3,6)];

% Initial guess of receiver position
r_R_bar = [0;0;0];

for i=1:length(Data.pseudo_range)
    % Ranges between the pseudolites and the receiver (TO BE ENTERED)
    rho1 = Data.pseudo_range(i,1);
    rho2 = Data.pseudo_range(i,2);
    rho3 = Data.pseudo_range(i,3);
    rho4 = Data.pseudo_range(i,4);
    rho5 = Data.pseudo_range(i,5);
    rho6 = Data.pseudo_range(i,6);
    
    % Maximum number of iterations for Newton's method
    max_iters = 10;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve Using Newton's Method
    
    % Initialized loop variables
    iters = 0;
    stop = 0;
    while stop == 0
	    iters = iters+1; % update number of iterations
    
	    % Compute range based on current estimate of receiver position
	    rho1_bar = norm(r_R_bar-r_S1);
	    rho2_bar = norm(r_R_bar-r_S2);
	    rho3_bar = norm(r_R_bar-r_S3);
        rho4_bar = norm(r_R_bar-r_S4);
        rho5_bar = norm(r_R_bar-r_S5);
        rho6_bar = norm(r_R_bar-r_S6);
    
	    % Compute matrices from linearization for Newton's method
	    A = [(r_R_bar(1)-r_S1(1))/rho1_bar (r_R_bar(2)-r_S1(2))/rho1_bar (r_R_bar(3)-r_S1(3))/rho1_bar 1;...
		    (r_R_bar(1)-r_S2(1))/rho2_bar (r_R_bar(2)-r_S2(2))/rho2_bar (r_R_bar(3)-r_S2(3))/rho2_bar 1;...
		    (r_R_bar(1)-r_S3(1))/rho3_bar (r_R_bar(2)-r_S3(2))/rho3_bar (r_R_bar(3)-r_S3(3))/rho3_bar 1; ...
            (r_R_bar(1)-r_S4(1))/rho4_bar (r_R_bar(2)-r_S4(2))/rho4_bar (r_R_bar(3)-r_S4(3))/rho4_bar 1; ...
            (r_R_bar(1)-r_S5(1))/rho5_bar (r_R_bar(2)-r_S5(2))/rho5_bar (r_R_bar(3)-r_S5(3))/rho5_bar 1; ...
            (r_R_bar(1)-r_S6(1))/rho6_bar (r_R_bar(2)-r_S6(2))/rho6_bar (r_R_bar(3)-r_S6(3))/rho6_bar 1];
	    b = [norm(rho1 - rho1_bar); norm(rho2 - rho2_bar); norm(rho3 - rho3_bar); norm(rho4 - rho4_bar); norm(rho5 - rho5_bar); norm(rho6 - rho6_bar)];
    
	    % Solve for updated receiver position estimate
	    delta_x = A\b;
        delta_r = delta_x(1:3);
        c_delta_Tr = delta_x(4);
	    r_R_bar = r_R_bar + delta_r;
    
	    % Stopping criteria based on convergence on solution within 1e-6 m or
	    % exceeding maximum number of iterations
	    if ((norm(delta_x) < 1e-6) || (iters > max_iters)) 
		    stop = 1;
	    end
    end
    
    Q = (A'*A)^-1;
    PDOP(i) = sqrt(Q(1,1)+Q(2,2)+Q(3,3));
    TDOP(i) = sqrt(Q(4,4));
    NDOP(i) = sqrt(Q(1,1));
    EDOP(i) = sqrt(Q(2,2));
    DDOP(i) = sqrt(Q(3,3));
    
    % Estimated position of the receiver
    r_R_sol = r_R_bar;
    r_R_sol_Arr_x(i) = r_R_sol(1);
    r_R_sol_Arr_y(i) = r_R_sol(2);
    r_R_sol_Arr_z(i) = r_R_sol(3);
    
    % Estimated receiver clock offset
    c_delta_Tr_sol = c_delta_Tr;
    c_delta_Tr_sol_Arr{i} = c_delta_Tr_sol;

end

% Calculate the error
sigma_north = NDOP*std_total;
sigma_east = EDOP*std_total;
sigma_down = DDOP*std_total;

upper_x = r_R_sol_Arr_x + 3*sigma_north;
lower_x = r_R_sol_Arr_x - 3*sigma_north;

upper_y = r_R_sol_Arr_y + 3*sigma_east;
lower_y = r_R_sol_Arr_y - 3*sigma_east;

upper_z = r_R_sol_Arr_z + 3*sigma_down;
lower_z = r_R_sol_Arr_z - 3*sigma_dwon;


% 3D trajectory plot
figure;
plot3(r_R_sol_Arr_x, r_R_sol_Arr_y, r_R_sol_Arr_z)
grid on
xlabel('North (m)')
ylabel('East (m)')
zlabel('Down (m)')

% Individual plots
figure;

subplot(3,1,1)
plot(time, r_R_sol_Arr_x)
grid on
ylabel('Est. North Pos. (m)')

subplot(3,1,2)
plot(time, r_R_sol_Arr_y)
grid on
ylabel('Est. East Pos. (m)')

subplot(3,1,3)
plot(time, r_R_sol_Arr_z)
grid on
ylabel('Est. Down Pos. (m)')
xlabel('Time (s)')