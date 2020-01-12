%% user interface
%%%%%%%%%%%%%%%%%  input%%%%%%%%%%%%%%
year=42;    % year input

%%%%%%%%%% clear %%%%%%%%%%%%%%%%
clear Tree_information_I

%%%%%%%%%%% format tranformation %%%%%%%%%%
% ID Type Age dbh xdot ydot
cnames={'ID','Type','Age','dbh','xdot','ydot'};
Tree_information_I(:,1)=cell2mat(Tree_information(year,1));
Tree_information_I(:,2)=cell2mat(Tree_information(year,2));
Tree_information_I(:,3)=cell2mat(Tree_information(year,3));
Tree_information_I(:,4)=cell2mat(Tree_information(year,4));
Tree_information_I(:,5)=cell2mat(Tree_information(year,5));
Tree_information_I(:,6)=cell2mat(Tree_information(year,6));

LF_annual_accumulation_invasive_I=cell2mat(LF_annual_accumulation_invasive(year));
LF_annual_accumulation_local_I=cell2mat(LF_annual_accumulation_local(year));

data=dataset({Tree_information_I,'ID','Type','Age','dbh','xdot','ydot'});

for i = 1:Nyear
    temp = cell2mat(Tree_information(i,2));
    IN_of_trees(i)=length(find(temp==1));
    LN_of_trees(i)=length(find(temp==2));
end

%%
%%%%%%%%%%% display %%%%%%%%%%%%%

figure(1)
% subplot(121)
x_linespace=1:100;
y_linespace=100:-1:1;
pcolor(x_linespace,y_linespace,LF_annual_accumulation_invasive_I);
colorbar
title(['invasive leave accumulation for the ',num2str(year),'th year']);
ylabel('grid 1:100');
xlabel('grid 1:100');
figure(2)
% subplot(122)
pcolor(x_linespace,y_linespace,LF_annual_accumulation_local_I);
colorbar
title(['local leave accumulation for the ',num2str(year),'th year']);
ylabel('grid 1:100');
xlabel('grid 1:100');

% I suggest to see the data in excel
% figure(3)
% uitable('Data',Tree_information_I,'ColumnName',cnames);
% clear x_linespace y_linespace

% scatter plot
figure(4)
x=[1:Nyear]';XI=[x,IN_of_trees'];XL=[x,LN_of_trees'];
scatter(x,IN_of_trees);
xlim([0,year]);%ylim([0,4000]);
title(['invasive total tree trending ',num2str(year),'th year']);
ylabel('tree total');
xlabel('Nyear');

figure(5)
scatter(x,LN_of_trees);
xlim([0,year]);%ylim([0,1000]);
title(['local total tree trending ',num2str(year),'th year']);
ylabel('tree total');
xlabel('Nyear');

figure(6)
I=find(Tree_information_I(:,2)==1);
scatter(Tree_information_I(I,5),Tree_information_I(I,6),'.','r'); % invasive tree location
hold on;
I=find(Tree_information_I(:,2)==2);
scatter(Tree_information_I(I,5),Tree_information_I(I,6),'.','b');
title(['total tree location ',num2str(year),'th year']);
ylabel('grid 1:100');
xlabel('grid 1:100');

%% export to excel 
excel_name=['Tree information for the ', num2str(year),'th year'];
xlswrite([excel_name '.xlsx'],cnames,'sheet1','A1');
xlswrite([excel_name '.xlsx'],Tree_information_I,'sheet1','A2');
clear excel_name