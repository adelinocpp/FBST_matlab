function interval = fbst_mvd_interval(LLR_TESTE, alpha, epsE, maxIter, nPoints, intFlag, verb)

if (nargin < 7) || isempty(verb)
    verb = false;
end
if (nargin < 6) || isempty(intFlag)
    intFlag = '1/3';
end
if (nargin < 5) || isempty(nPoints)
    nPoints = 200; 
end
if (nargin < 4) || isempty(maxIter)
    maxIter = 100;
end
if (nargin < 3) || isempty(epsE)
    epsE = 1e-6;
end

if (nargin < 2) || isempty(alpha)
    alpha = 0.05;
end

Npts = length(LLR_TESTE);
med = mean(LLR_TESTE);
dev = std(LLR_TESTE);

etha0 = med;
etha1 = med+1*dev/sqrt(Npts);
abs(etha1-med);
LRval0 = (1 - alpha) - FBST_MVD(LLR_TESTE,etha0,nPoints);
LRval1 = (1 - alpha) - FBST_MVD(LLR_TESTE,etha1,nPoints);

k = 0;
if (verb)
    fprintf('%02i - %5.3e,%5.3e...\n',k,LRval0, LRval1);
end
while ((abs(LRval1) > epsE) && (k < maxIter) )
    etha2 = etha1 - LRval1*(etha1 - etha0)/(LRval1 - LRval0);
    etha0 = etha1;
    etha1 = etha2;
    LRval0 = (1 - alpha)-FBST_MVD(LLR_TESTE,etha0,nPoints);
    LRval1 = (1 - alpha)-FBST_MVD(LLR_TESTE,etha1,nPoints);
    k = k+1;
    if (verb)
        fprintf('%02i - %5.3e,%5.3e...\n',k,LRval0, LRval1);
    end
end
interval = abs(etha2-med);
