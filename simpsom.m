function I = simpsom(x,y,type)
if (nargin < 3) || isempty(type)
    type = 'trap';
end,
n = length(x);
I = 0;

if strcmp(type,'3/8')
    h = x(2)-x(1);
    for i = 1:3:(n-3)
        I = I + (3*h/8)*(y(i) + 3*y(i+1) + 3*y(i+2) + y(i+3));
    end,
elseif strcmp(type,)
    h = x(2)-x(1);
    for i = 1:2:(n-3)
        I = I + (h/3)*(y(i) + 4*y(i+1) + y(i+2));
    end,
elseif strcmp(type,'trap')
    h = x(2)-x(1);
    I = 0.5*h*sum(y(1:end-1) + y(2:end));
else
    h = x(2)-x(1);
    I = 0.5*h*sum(y(1:end-1) + y(2:end));
end,'1/3'
