function Rsq = calc_rsq(x,y,zeroForced)

% calculate R^2
% https://uk.mathworks.com/help/matlab/data_analysis/linear-regression.html
    
% convert x and y in column vectors
if size(x,2) > 1
    x   =   transpose(x);
end

if size(y,2) > 1
    y   =   transpose(y);
end

% zeroForced should be boolean 1 or 0
if zeroForced
    b1      =   x\y;
    yCalc   =   x*b1;
    Rsq     =   1 - sum((y-yCalc).^2)/sum((y-mean(y)).^2);
else
    X       =   [ones(length(x),1) x];
    b1      =   X\y;
    yCalc   =   X*b1;
    Rsq     =   1 - sum((y-yCalc).^2)/sum((y-mean(y)).^2);
end
end