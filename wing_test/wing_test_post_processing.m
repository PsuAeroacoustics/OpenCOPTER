clc
clear all
close all

load('results_1.mat');

azimuth = linspace(1,360,360);

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

C_l_wing_0_per = reshape(wing_C_l_time_hist(1,8,:), [360,1]);
C_l_wing_25_per = reshape(wing_C_l_time_hist(1,4,:), [360,1]);
C_l_wing_50_per = reshape(wing_C_l_time_hist(1,1,:), [360,1]);
C_l_wing_75_per = reshape(wing_C_l_time_hist(2,4,:), [360,1]);
C_l_wing_100_per = reshape(wing_C_l_time_hist(2,8,:), [360,1]);

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
y_wing = [y,y+0.5];
%plot(wing_C_l_time_hist(1,:,1));
wing_C_l = [fliplr(wing_part_0_C_l), wing_part_1_C_l];

figure
plot(y_wing, wing_C_l)
