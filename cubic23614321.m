function [ rts, info ] = cubic23614321(C)
%cubic23614321x solves cubic = 0
%   C is coefficient array

%unpacking coefficients
a = C(1);
b = C(2);
c = C(3);
d = C(4);

%Quadratic or Less -- Quadratic Formula
if a == 0
    rts = stableQuadratic(b,c,d);
    info = 'quadratic or less -- easy case';
    
%Actual Cubic
else 
    if evalCubic(a,b,c,d,0) == 0
        r1 = 0;
    else
        
    %p0 = bisection(-2e100, 2e100, a,b,c,d);
    %r1 = newtonMethod(p0, [a,b,c,d]);
    
    r1 = newtonMethod(8, [a,b,c,d]);
    if iscellstr(r1)
        r1 = newtonMethod(1e20, [a,b,c,d]);
        if iscellstr(r1)
            r1 = newtonMethod(-1e20, [a,b,c,d]);
        end
    end
    end
    newQuadratic = deconv(C, [1,-r1]);
    rt = stableQuadratic(newQuadratic(1), newQuadratic(2), newQuadratic(3));
    rts = [r1, rt(1), rt(2)];
    info = 'cubic';
end
end

%Auxiliary Functions
function [eval] = evalCubic(a,b,c,d,x)
    %evaluates poly ax^3 + bx^2 + cx + d at x
    eval = a*x^3 + b*x^2 + c*x + d;
end

function [soln] = newtonMethod(p0, funcArray)
    %begining iterator
    i = 1;
    
    %unpacking function coefficients
    a = funcArray(1);
    b = funcArray(2);
    c = funcArray(3);
    d = funcArray(4);
    
    %setting up derivative
    bPrime = 3*a;
    cPrime = 2*b;
    dPrime = c;
    
    %running newton's method
    while i <= 10000000
        
        %possible double root
        %deriv = evalCubic(0,bPrime,cPrime,dPrime,p0);
        %if abs(deriv) <= 10^-6
          %  ...
        %end
        
        p = p0 - evalCubic(a,b,c,d,p0)/evalCubic(0,bPrime,cPrime,dPrime,p0);
        if abs(p - p0) < 10^-10
            soln = p;
            return
        else
            i = i + 1;
            p0 = p;
        end
    end
    soln = 'failed; more than 10000 iterations';
end

function [rts] = stableQuadratic(a, b, c)
    %evaluates quadratic formula
    sqrtDiscriminant = sqrt(b^2 - 4 * a * c);
    if a == 0
        rts = stableLinear(b,c);
    elseif sqrtDiscriminant == 0
        r1 = - b / (2 * a);
        r2 = r1;
        rts = [r1,r2];
    else
        if b < 0
            quadraticPart = -b + sqrtDiscriminant;
        else
            quadraticPart = -b - sqrtDiscriminant;
        end
        r1 = quadraticPart / (2 * a);
        r2 = (2 * c) / quadraticPart;
        rts = [r1, r2];
    end
end

function [rts] = stableLinear(a,b)
    if a == 0
        rts = [];
    else
        rts = (-b/a);
    end
end

%%function [guess] = bisection(left, right, a,b,c,d)
 %   i = 1;
 %   F = evalCubic(a,b,c,d,left);
 %   while i <= 1000000
 %       guess = left + (right - left)/2;
 %       mid = evalCubic(a,b,c,d,guess);
 %       if mid == 0 || (right-left)/2 < 1e-6
 %           return
 %       end
 %       i = i + 1;
 %       if (F * mid) > 0
 %           left = guess;
 %           F = mid;
 %       else
 %           right = guess;
 %       end
 %   end
%%end