clc
clear all
close all

load('results_RNT.mat');

%azimuth = linspace(1,360,360);
%%
b_1_time_hist = reshape(rotor_C_T_time_hist(1,24,:),[360,1]);
% b_2_time_hist = reshape(rotor_C_T_time_hist(2,24,:),[360,1]);
% b_3_time_hist = reshape(rotor_C_T_time_hist(3,24,:),[360,1]);
% b_4_time_hist = reshape(rotor_C_T_time_hist(4,24,:),[360,1]);
%plot(rotor_C_T_time_hist(1,24,:))

%figure
plot(azimuth(2:360), b_1_time_hist(2:360));
% hold on
% plot(azimuth(2:360), b_2_time_hist(2:360));
% hold on
% plot(azimuth(2:360), b_3_time_hist(2:360));
% hold on
% plot(azimuth(2:360), b_4_time_hist(2:360));
xlabel("azimuth (degrees)")
ylabel("blade dC_T")
%legend("blade 1","blade 2","blade 3","blade 4")
%%
C_l_wing_0_per = reshape(wing_C_l_time_hist(1,1,:), [360,1]);
C_l_wing_25_per = reshape(wing_C_l_time_hist(1,2,:), [360,1]);
C_l_wing_50_per = reshape(wing_C_l_time_hist(1,4,:), [360,1]);
C_l_wing_75_per = reshape(wing_C_l_time_hist(1,6,:), [360,1]);
C_l_wing_100_per = reshape(wing_C_l_time_hist(1,8,:), [360,1]);

figure
plot(azimuth(2:360), C_l_wing_0_per(2:360));
xlabel("azimuth (degrees)")
ylabel("wing dC_l y/b = 0.0")
%legend("y/b = 0","y/b = 1.0")

figure
plot(azimuth(2:360), C_l_wing_25_per(2:360));
xlabel("azimuth (degrees)")
ylabel("wing dC_l y/b = 0.25")

figure
plot(azimuth(2:360), C_l_wing_50_per(2:360));
xlabel("azimuth (degrees)")
ylabel("wing dC_l y/b = 0.5")

figure
plot(azimuth(2:360), C_l_wing_75_per(2:360));
xlabel("azimuth (degrees)")
ylabel("wing dC_l y/b = 0.75")

figure
plot(azimuth(2:360), C_l_wing_100_per(2:360));
xlabel("azimuth (degrees)")
ylabel("wing dC_l y/b = 1.0")


l = 0:1:7;
theta = l*pi/8;
y = (1-cos(theta))/4;
y_wing = [y];
%plot(wing_C_l_time_hist(1,:,1));
wing_C_l = wing_part_0_C_l;

figure
plot(y_wing, wing_C_l)

%%

wake_x_rotor_only = reshape(wake_0_trajectory(1,1,:),[1,1024]);
wake_y_rotor_only = reshape(wake_0_trajectory(1,2,:),[1,1024]);
wake_z_rotor_only = reshape(wake_0_trajectory(1,3,:),[1,1024]);


load("results_RNT.mat")

wake_x_wing_plus_rotor = reshape(wake_0_trajectory(1,1,:),[1,1024]);
wake_y_wing_plus_rotor = reshape(wake_0_trajectory(1,2,:),[1,1024]);
wake_z_wing_plus_rotor = reshape(wake_0_trajectory(1,3,:),[1,1024]);

plot3(wake_x_rotor_only,wake_y_rotor_only,wake_z_rotor_only)
hold on
plot3(wake_x_wing_plus_rotor,wake_y_wing_plus_rotor, wake_z_wing_plus_rotor)
xlabel("x/R")
ylabel("y/R")
zlabel("z/R")
legend("rotor only", "rotor + wing")
