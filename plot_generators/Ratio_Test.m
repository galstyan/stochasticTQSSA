path = strsplit(pwd, '/');
[a, b] = size (path)
for i = 1 : b-1
    path(i) = strcat(path(i),'/')
end

newpwd = strcat(path{1:b-1})

fid = fopen(strcat(newpwd,'/data/Exact_ODE_ratio_tSubstrate'))
s = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen(strcat(newpwd,'/data/Exact_ODE_ratio_1'))
ratio1 = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen(strcat(newpwd,'/data/Exact_ODE_ratio_2'))
ratio2 = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

fid = fopen(strcat(newpwd,'/data/Exact_ODE_ratio_3'))
ratio3 = fscanf(fid, '%g %g', [1 inf]);
fclose(fid);

h = figure('Position', [100, 100, 650, 500]);
plot(s,ratio1,'k-o',s,ratio2,'k-s', s, ratio3, 'k-+');
xlabel('Total Substrate $\bar{s}$', 'FontSize',12.5, 'Interpreter','latex');
ylabel('$\sigma^2_{c}(\bar{s})/ e_0 \bar{s}$ ratio', 'FontSize',12.5, 'Interpreter','latex');

print (h, '-dpsc', strcat(newpwd,'/figures/RatioTest.eps'))