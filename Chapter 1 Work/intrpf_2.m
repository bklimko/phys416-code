function yi = intrpf_2(xi,x,y)
% Function to interpolate between data points
% using Lagrange polynomial (quadratic)
% Inputs
%   x    Vector of x coordinates of data points (3 values)
%   y    Vector of y coordinates of data points (3 values)
%   xi   The x value where interpolation is computed
% Output
%   yi   The interpolation polynomial evaluated at xi
%NOTE: intrpf_2 is identical in function to intrpf except that it is
%generalized to n length input vectors x and y

%* Calculate yi = p(xi) using Lagrange polynomial
yi = 0;
for p = 1:length(x) %loop to calculate summation
    tempval = 1;
    for q = 1:length(x) %loop to calculate product
        if p ~= q
            tempval = tempval * (((xi-x(q))/(x(p)-x(q)))*y(p));
        end
    end
    yi = yi + tempval;
end
return;
