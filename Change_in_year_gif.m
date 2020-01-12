%%%%%% change in year %%%%%%
% year=200;    % year input

%%%%%%%%%% clear %%%%%%%%%%%%%%%%
clear Tree_information_I

%%%%%%%%%%% format tranformation %%%%%%%%%%
% ID Type Age dbh xdot ydot
cnames={'ID','Type','Age','dbh','xdot','ydot'};
for year=150:200
clear Tree_information_I
Tree_information_I(:,1)=cell2mat(Tree_information(year,1));
Tree_information_I(:,2)=cell2mat(Tree_information(year,2));
Tree_information_I(:,3)=cell2mat(Tree_information(year,3));
Tree_information_I(:,4)=cell2mat(Tree_information(year,4));
Tree_information_I(:,5)=cell2mat(Tree_information(year,5));
Tree_information_I(:,6)=cell2mat(Tree_information(year,6));

figure(1)
subplot(221)

x_linespace=1:100;
y_linespace=100:-1:1;
pcolor(x_linespace,y_linespace,LF_annual_accumulation_invasive_I);
colorbar
title(['Invasive leave accumulation for the ',num2str(year),'th year']);

subplot(222)
pcolor(x_linespace,y_linespace,LF_annual_accumulation_local_I);
colorbar
title(['local leave accumulation for the ',num2str(year),'th year']);


% figure(2)
subplot(223)
I=find(Tree_information_I(:,2)==1);
scatter(Tree_information_I(I,5),Tree_information_I(I,6),'.','r'); % invasive tree location
hold on;
I=find(Tree_information_I(:,2)==2);
scatter(Tree_information_I(I,5),Tree_information_I(I,6),'.','b');
hold off;
title(['Year ',num2str(year)]);
pause(1)
end
