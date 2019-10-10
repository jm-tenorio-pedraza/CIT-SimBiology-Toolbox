VarNames = {'Time1' 'Dose_1mgkg' 'Time2' 'Dose_3mgkg' 'Time3' 'Dose_10mgkg'};
VarType = repelem({'double'}, 1,6);
Tumor_PK = importPKfile('/Users/migueltenorio/Documents/Data/CSV/Nedrow_2017_2_Tumor.csv',VarNames,VarType,[3,Inf]);
Tumor_PK = postProcessPKTable(Tumor_PK);

Blood_PK = importPKfile('/Users/migueltenorio/Documents/Data/CSV/Nedrow_2017_2_Blood.csv',VarNames,VarType,[3,Inf]);
Blood_PK = postProcessPKTable(Blood_PK);


Tumor2Blood_PK = importPKfile('/Users/migueltenorio/Documents/Data/CSV/Nedrow_2017_2_Tumor_to_Blood.csv',VarNames,VarType,[3,Inf]);
Tumor2Blood_PK = postProcessPKTable(Tumor2Blood_PK);


writetable(Tumor_PK, '/Users/migueltenorio/Documents/Data/Nedrow_2017_2.xlsx','Sheet', 'Tumor')
writetable(Blood_PK, '/Users/migueltenorio/Documents/Data/Nedrow_2017_2.xlsx','Sheet', 'Blood')
writetable(Tumor2Blood_PK, '/Users/migueltenorio/Documents/Data/Nedrow_2017_2.xlsx','Sheet', 'Tumor_to_Blood')
