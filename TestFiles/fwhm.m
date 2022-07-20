% Initialization steps.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.

format long g;
format compact;
fontSize = 24;


x = linspace(0, 100, 1024);
%% 

amp = 50;
mu = 35;
sigma = 20;
y = amp * exp(-(x - mu).^2 / (2 * sigma^2));


hFig = figure;
plot(x, y, 'b-', 'LineWidth', 2);
grid on;
xlabel('x', 'FontSize', fontSize);
ylabel('y', 'FontSize', fontSize);
hFig.WindowState = 'maximized';

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

% Compute the full width, half max.
fwhm1 = x2 - x1
text(mu, halfHeight+2, sprintf('width = %.2f', fwhm1), 'FontSize', fontSize, 'Color', 'r', 'HorizontalAlignment', 'center');
caption = sprintf('Full Width, Half Max = %.2f', x2 - x1);
title(caption,  'FontSize', fontSize);