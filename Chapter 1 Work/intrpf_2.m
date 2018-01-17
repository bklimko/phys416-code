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
    num = 1;
    denom = 1;
    x_altered = x([1:p-1, p+1:end]);
    for q = x_altered %loop to calculate product
        num = num*(xi-q);
        denom = denom * (x(p) - q);
    end
    y_0 = (num * y(p))/denom;
    yi = yi + y_0;
    
end
return;
