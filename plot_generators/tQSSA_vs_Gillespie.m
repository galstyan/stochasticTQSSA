path = strsplit(pwd, '/');
[a, b] = size (path)
for i = 1 : b-1
    path(i) = strcat(path(i),'/')
end

newpwd = strcat(path{1:b-1})

fid = fopen( strcat(newpwd,'/Data/Interval'))
t = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen( strcat(newpwd,'/data/Exact_Substrate_Mean_PSD'))
gill_psd = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen( strcat(newpwd,'/data/Exact_Substrate_Mean_MSD'))
gill_msd = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen( strcat(newpwd,'/data/Exact_Substrate_Mean'))
gill = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen( strcat(newpwd,'/data/TQSSA_Substrate_Mean_MSD'))
tQSSA_msd = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen( strcat(newpwd,'/data/TQSSA_Substrate_Mean_PSD'))
tQSSA_psd = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen( strcat(newpwd,'/data/TQSSA_Substrate_Mean'))
tQSSA = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

h = figure('Position', [100, 100, 650, 500]);
plot(t,gill_msd,'k--',t,gill,'k', t, gill_psd, 'k--', t, tQSSA_msd, 'k--*', t, tQSSA, 'k-*', t, tQSSA_psd, 'k--*');

xlabel('Time (sec)', 'FontSize',12.5, 'Interpreter','latex');
ylabel('Substrate', 'FontSize',12.5, 'Interpreter','latex');

ylim([-1 inf])

print (h, '-dpsc', strcat(newpwd,'/figures/tQSSA_vs_Gillespie.eps'))