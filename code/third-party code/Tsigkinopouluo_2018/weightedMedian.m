function wMed = weightedMedian(D,W)
% The function calculates the weighted median from a vector of values based
% on their weights. The input is a vector with the values and a vector with
% their corresponding weights.

if nargin ~= 2
    error('weightedMedian:wrongNumberOfArguments', ...
      'Wrong number of arguments.');
end

if size(D) ~= size(W)
    error('weightedMedian:wrongMatrixDimension', ...
      'The dimensions of the input-matrices must match.');
end

wMed = [];

% (line by line) transformation of the input-matrices to line-vectors
d = reshape(D',1,[]);   
w = reshape(W',1,[]);

% sort the vectors
A = [d' w'];
ASortzeros = sortrows(A,1);
ASort=ASortzeros;
ASort(any(ASortzeros(:,2)<=(10^-14),2),:) = [];

dSort = ASort(:,1)';
wSort = ASort(:,2)';

% If there is only one value, it is the weighted median
if length(dSort)== 1
    wMed=dSort;
end

% If there are only two values the weighted median is their mean if they 
% have equal weights, otherwise the weighted median is the value with the
% largest weight.
if length(dSort)== 2
    if wSort(1)== wSort(2)
        wMed= (dSort(1)+dSort(2))/2; 
    elseif wSort(1)> wSort(2)
        wMed= dSort(1);
    else wMed= dSort(2);
    end
end

% If there are more than two values, find the value which is the 50% weighted percentile:
if ((length(dSort)~=1) && (length(dSort)~=2))
i = 1;
j = length(wSort);
Start=wSort(i);
End=wSort(j);

while i < j-1
    if Start-End > 10^-14
       End= End+wSort(j-1);
       j=j-1;
    else
       Start= Start+wSort(i+1);
       i=i+1;
    end
end

% If the 50% weighted percentile falls between two values, the weighted
% median is the average of the two values if their individual weights are
% equal, otherwise the weighted median is the value with the latgest
% individual weight
if abs(Start-End) < (10^-14)
   wMed=(dSort(i)+dSort(j))/2;
elseif Start-End > (10^-13)
   wMed=dSort(i);
else
   wMed=dSort(j);
end
end 

