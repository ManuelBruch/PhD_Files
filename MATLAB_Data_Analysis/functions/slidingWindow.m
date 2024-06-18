function rates=slidingWindow(x,y,width,steps)
% if x is a row vector, assume data in y are also in rows
if size(x,1) == 1
    x   =   transpose(x);
    y   =   transpose(y);
end

% number of steps
[n,m]                   =   size(y);
steps_default           =   n-width+1;              % to include last value


if isnan(steps)
    steps   =   steps_default;
end

% preallocate storage structure
rates.coefficients.m    =   zeros(steps,m);    % slope of fitted line --> Variable of Interest
rates.coefficients.b    =   zeros(steps,m);
rates.error_estimates   =   cell(steps,m);
rates.standardError     =   zeros(steps,m);
rates.steps             =   steps;
rates.dimensions        =   size(rates.coefficients.m);
time_windows            =   cell(steps,m);

% Sliding Window algorithm
for i = 1:m
    for j = 1:steps
        [p,S]                       =   polyfit(x(j:j+width-1),y(j:j+width-1,i),1);
        [yy,delta]                  =   polyval(p,x(j:j+width-1),S);
        rates.standardError(j,i)    =   mean(delta);
        rates.coefficients.m(j,i)   =   p(1);
        rates.coefficients.b(j,i)   =   p(2);
        rates.error_estimates{j,i}  =   S;
        time_windows{j,i}           =   x(j:j+width-1);
    end
end
rates.time_windows      =   time_windows(:,1);

end