function T=detectPD(x, tspan)
indx = find(x>20,1);
if isempty(indx)
    T = tspan(end);
else
    T = tspan(indx);
end
    

