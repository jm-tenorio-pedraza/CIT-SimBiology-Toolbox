function quantiles =quantiles(x,p)
n_i=size(x,2);
quantiles=nan(1,n_i);
n_j=size(x,1);
for i=1:n_i
    x_sorted = sort(x(:,i));
    indx = ceil(p*n_j);
    quantiles(i)=x_sorted(indx);
end
end



