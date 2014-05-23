files = dir('*.out');

% load all files
for k = 1:length(files)
    load(files(k).name, '-ascii')
end


hold on
% find and plot electron densities
myvars = whos;

color = hsv(15);
cc = 1
for i = 1:length(myvars)
    %myvars(i).name
    if (strncmp(char(myvars(i).name),'eDens_array_ex_40000_0_8_',25))
        plot(x_q_array_40000,eval(myvars(i).name), 'DisplayName',strrep(myvars(i).name,'_','\_'),'Color',color(cc,:))
        cc = cc +1
    end
end


% plot conduction band
figure()
hold on
cc = 1
for i = 1:length(myvars)
    %myvars(i).name
    if (strncmp(char(myvars(i).name),'pot_tot_array_p_ex_40000_0_8_',29))
        plot(x_array_40000,eval(myvars(i).name), 'DisplayName',strrep(myvars(i).name,'_','\_'),'Color',color(cc,:))
        cc = cc +1
    end
end