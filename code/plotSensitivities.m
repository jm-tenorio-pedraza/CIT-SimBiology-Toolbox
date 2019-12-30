function plotSensitivities(S)
s_ij = log10(diag(S)/max(diag(S)));
plot(s_ij, '-d')
hold on
plot(0:length(s_ij),ones(1,1+length(s_ij))*(-1),'-k')
ylim([-5, 0])
title('Singular values of sensitivity matrix')
ylabel('Log_{10} \sigma_i / max \sigma_i')
xlabel('Principal components')
end