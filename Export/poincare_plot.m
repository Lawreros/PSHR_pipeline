function [SD1, SD2] = poincare_plot(mat)
% Generates a poincare plot from the data
% Inputs:
%   mat: [n-by-1 array], vector containing all RR values

% Returns:
%   SD1: [int] The standard deviation along the y=x axis.
%   SD2: [int] The standard deviation perpindicular to the y=x axis


    cut = [];
    for i=1:length(mat)              % Eliminate NaN's
        if mat(i,1) == 0 || isnan(mat(i,1))
            cut = [cut,i];
        end
    end
    mat(cut) = [];
    mat(:,2) = [mat(2:end,1);0];
    mat(end,:) = [];

    % Calculate SD1 and SD2
    xc = sum(mat(1:end-1,1))/(length(mat)-1);
    yc = sum(mat(2:end,2))/(length(mat)-1);

    SD1 = sqrt((1/length(mat))*nansum(((mat(:,1)-mat(:,2))-nanmean(mat(:,1)-mat(:,2))).^2)/2);
    SD2 = sqrt((1/length(mat))*nansum(((mat(:,1)+mat(:,2))-nanmean(mat(:,1)+mat(:,2))).^2)/2);

    % Making a rotated elipsoid to display SD1 and SD2 https://www.mathworks.com/matlabcentral/answers/342148-how-can-i-rotate-an-ellipse
    x = zeros(1000,1);
    y = zeros(1000,1);
    theta = linspace(0, 2*pi, 1000);
    for k = 1:1000
        x(k) = SD2*cos(theta(k));
        y(k) = SD1*sin(theta(k));
    end
    alpha = pi/4;
    R = [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
    rCoords = R*[x'; y'];
    xr = rCoords(1,:)';
    yr = rCoords(2,:)';

    min_RR = nanmin(mat(:,1));
    max_RR = nanmax(mat(:,1));


    scatter(mat(:,1), mat(:,2), 15)
    axis([min_RR-50 max_RR+50 min_RR-50 max_RR+50])
    xlabel('RR_n (ms)');
    ylabel('RR_n_+_1 Interval (ms)');
    title('Poincare Plot placeholder title')
    hold on;
    plot(xr+xc, yr+yc, 'r');
    plot([0:1600],[0:1600],'r');
    hold off;

end