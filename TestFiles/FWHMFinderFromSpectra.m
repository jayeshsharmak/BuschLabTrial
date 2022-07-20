clc;    % Refresh command window 
close all;  % Close figs etc. 
clear;  % clearvars also works to clear existing variables
workspace;  % Opens workspace in case it's closed

tmp = readOceanOpticsSpectroV2('20220708_835nm.txt')
x = tmp (:,1) %wavelengthsX 
y = tmp (:,2) %PhotonCount
% Initialization steps.

%% 

amp = 50;
mu = 35;
sigma = 20;
fontSize = 12; %ancillary variables for future figure

hFig = figure;
plot(x, y, 'b-', 'LineWidth', 2);
grid on;
xlabel('x', 'FontSize', fontSize);
ylabel('y', 'FontSize', fontSize);
%hFig.WindowState = 'maximized';
%hFig.Position = [10 10 550 400]; 

% Find the half height - midway between min and max y values.
halfHeight = (min(y) + max(y)) / 2;
hold on;
yline(halfHeight, 'Color', 'g', 'LineWidth', 2);

% Find left edge
index1 = find(y >= halfHeight, 1, 'first');
x1 = x(index1)
line([x1, x1], [0, y(index1)], 'Color', 'r', 'LineWidth', 2);
text(x1+1, 5, sprintf('x1 = %.2f', x1), 'FontSize', fontSize, 'Color', 'r');

% Find right edge
index2 = find(y >= halfHeight, 1, 'last');
x2 = x(index2)
line([x2, x2], [0, y(index2)], 'Color', 'r', 'LineWidth', 2);
text(x2+1, 5, sprintf('x2 = %.2f', x2), 'FontSize', fontSize, 'Color', 'r'); 
text(x2+3)

% Compute the full width, half max.
fwhm1 = x2 - x1
text(mu, halfHeight+2, sprintf('width = %.2f', fwhm1), 'FontSize', fontSize, 'Color', 'r', 'HorizontalAlignment', 'center');
caption = sprintf('Full Width, Half Max = %.2f', x2 - x1);
title(caption,  'FontSize', fontSize);
 
