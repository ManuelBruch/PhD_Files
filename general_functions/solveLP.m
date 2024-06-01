function [x, fval, exitflag, output, lambda] = solveLP(NVArgs)
%SOLVELP manually calculate LP solutions
%   Detailed explanation goes here
arguments
    NVArgs.f        (1,:)   double
    NVArgs.A                mat         =   [];
    NVArgs.b        (1,:)   double      =   [];
    NVArgs.Aeq              mat         =   [];
    NVArgs.beq      (1,:)   double      =   [];
    NVArgs.lb               double      =   [];
    NVArgs.ub               double      =   [];
    NVArgs.options
end

[x, fval, exitflag, output, lambda]     =   linprog(NVArgs.f, NVArgs.A, ...
    NVArgs.b, NVArgs.Aeq, NVArgs.beq, NVArgs.lb, NVArgs.ub, NVArgs.options);

end

