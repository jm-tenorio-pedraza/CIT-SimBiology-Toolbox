function [logpdf, loggrad] = getlogPdfGrad(p,fun)
logpdf=fun(p);
loggrad=getGradient(p,fun);
end
