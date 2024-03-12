clear all;
clc;

%% Baseline values
% psi = deg2rad([ ...
%   19.225
%   34.245
%   49.265
%   64.275
%   79.285
%   94.055
%  109.055
%  124.045
%  139.055
%  154.055
%  169.045
%  199.035
%  214.055
%  229.035
%  244.065
%  259.045
%  274.045
%  289.045
%  304.035
%  319.045
%  334.045
%  349.045]);
% 
% el_twist = deg2rad([ ...
% -0.74
% -0.64
% -0.63
% -0.53
% -0.15
% -0.01
% -0.23
% -0.47
% -1.06
% -1.23
% -1.32
% -0.93
% -0.77
% -0.39
% -0.11
%  0.15
%  0.30
%  0.24
% -0.03
% -0.39
% -0.66
% -0.83]);

%% Minimum Noise Values
psi = deg2rad([ ...
 19.055
 34.045
 49.055
 64.055
 79.055
 94.055
109.045
124.055
139.045
154.055
169.055
199.045
214.055
229.055
244.055
259.055
274.055
289.065
304.055
319.065
334.045
349.055]);

el_twist = deg2rad([ ...
 0.35
-0.50
-1.52
-1.76
-1.14
-0.35
-0.01
 0.55
 0.40
-0.62
-1.72
-2.18
-1.11
 0.11
 1.02
 1.16
 0.48
-0.60
-1.34
-1.33
-0.63
 0.08]);

%% Minimum Vibration
% psi = deg2rad([ ...
%  19.055
%  34.055
%  49.055
%  64.045
%  79.055
%  94.045
% 109.055
% 124.045
% 139.045
% 154.035
% 169.035
% 199.045
% 214.045
% 229.045
% 244.055
% 259.045
% 274.055
% 289.045
% 304.055
% 319.045
% 334.035
% 349.045]);
% 
% el_twist = deg2rad([ ...
% -2.24
% -1.92
% -1.02
%  0.21
%  1.53
%  1.72
%  0.73
% -0.83
% -2.14
% -2.70
% -2.09
%  0.18
%  0.68
%  0.44
% -0.35
% -1.03
% -1.15
% -0.35
%  0.63
%  0.86
%  0.48
% -0.52]);

series_size = 6;

upper_bounds = zeros([(2*series_size + 2), 1]);
lower_bounds = zeros([(2*series_size + 2), 1]);
for i=1:(2*series_size + 2)
    upper_bounds(i) = inf;
    lower_bounds(i) = -inf;
end

upper_bounds(end) = 1.0;
lower_bounds(end) = 1.0;

fit_type = sprintf('fourier%d', series_size);

fit_options = fitoptions(fit_type, 'Lower', lower_bounds, 'Upper', upper_bounds);

curve_fit = fit(psi, el_twist, fit_type, fit_options);

fit_psi = linspace(0, 2*pi, 1000);

plot(rad2deg(psi), el_twist, 'r', rad2deg(fit_psi), curve_fit(fit_psi), 'b');

curve_fit
