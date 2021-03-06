% interp - Program to interpolate data using Lagrange 
% polynomial to fit quadratic to three data points
%script and intrpf.m function originally by Frank Toffoletto,
%edited by Benjamin Klimko, PHYS 416, Spring 2018

clear all;   % Clear memory
%* Initialize the data points to be fit by quadratic
disp('Enter data points as x,y pairs (e.g., [1 2])');
disp('Type ''done'' when finished inputting data');
idx = 1;
temp = input('Enter data point: ', 's');
while ~strcmp(temp, 'done')
 temp = str2num(temp);
 x(idx) = temp(1);
 y(idx) = temp(2);
 temp = input('Enter data point: ', 's');
 idx = idx + 1;
end

%Set the intial points to be fit
% x=[1 2 3]
% y=[1 12 1]

%* Establish the range of interpolation (from x_min to x_max)
xr = input('Enter range of x values as [x_min x_max]: ');
% xr=[1 4]
%* Find yi for the desired interpolation values xi using
%  the function intrpf
nplot = 100;     % Number of points for interpolation curve
for i=1:nplot
  xi(i) = xr(1) + (xr(2)-xr(1))*(i-1)/(nplot-1);
  yi(i) = intrpf_2(xi(i),x,y);  % Use intrpf_2 function to interpolate
end

%* Plot the curve given by (xi,yi) and mark original data points
figure(1);clf % open the figure and clear the screen
plot(x,y,'*',xi,yi,'-');
xlabel('x');
ylabel('y');
title('Multi-point interpolation');
legend('Data points','Interpolation  ');