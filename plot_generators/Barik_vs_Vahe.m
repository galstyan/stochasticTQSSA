path = strsplit(pwd, '/');
[a, b] = size (path)
for i = 1 : b-1
    path(i) = strcat(path(i),'/')
end

newpwd = strcat(path{1:b-1})

fid = fopen(strcat(newpwd,'/data/Cond_Exp_Substrate'))
s = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen(strcat(newpwd,'/data/Cond_Exp_Barik'))
barik = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen(strcat(newpwd,'/data/Cond_Exp_TQSSA'))
tQSSA = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen(strcat(newpwd,'/data/Exact_ODE_Cond_1'))
ExactCond1 = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen(strcat(newpwd,'/data/Exact_ODE_Cond_2'))
ExactCond2 = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen(strcat(newpwd,'/data/Exact_ODE_Cond_3'))
ExactCond3 = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

h = figure('Position', [100, 100, 650, 500]);
plot(s,barik,'k--',s,tQSSA,'k', s, ExactCond1, 'k-o', s, ExactCond2, 'k-o', s, ExactCond3, 'k-o');
xlabel('Total Substrate $\bar{s}$', 'FontSize',12.5, 'Interpreter','latex');
ylabel('Conditional Expectation $\langle c |\bar{s} \rangle $', 'FontSize',12.5, 'Interpreter','latex');

print (h, '-dpsc', strcat(newpwd,'/figures/Barik_vs_Vahe.eps'))