%% ManHam model vegetation dynamics part,
%% written on Jan 2010 by Jiang, updated on May 10 2010
%% MelaHam Model vegetation dynamics,
%% adjust on December 2019 by Lu
% field area:100*100 meters, cell for vadose zone 5*5 meters

% This file describes the dynamics of two competing tree species.

clear; clc;

% load deltz
Vs = zeros(20,20,5);   %   Vadose zone salinity, will be connect to another 
                       %   matlab code which simulate salinity. But now
                       %   I set to constant, make sure this code works
                       
LF_accumulation_invasive=zeros(100,100); % initialize invasive leaf accumlation
LF_accumulation_local=zeros(100,100); % initialize local species leaf acculation

%  LF_annualLF_invasive=zeros(100,100,t);
%  LF_annualLF_local=zeros(100,100,t);

deltz=zeros(20,20);
deltz=deltz+1000;

Nspec = 2;                % Number of species;
Nyear =200;              % Time span in years

% Ntree_max = 1000;          % Maximum tree numbers for each species;
phi = 0.2;                % Light extinction coefficient                 
                         
G = [292 267]/12;        % Growth constants for the two species
D_max=[140 100];         % Maximum d.b.h. for the two species(cm)
H_max=[4500 4000];       % Maximum height for the two species (cm)
AGE_max=[300 250];       % Maximum age for the two species (year)
b2=[48.04 77.26];        % Constant in height to d.b.h. relationship
b3=[0.172 0.396];        % Constant in height to d.b.h. relationship
cLW=[38.90 27.55];       % Constant in leaf weight to d.b.h. relationship
bLW=[1.62 1.79];         % Constant in leaf weight to d.b.h. relationship
cLA=6.181;               %       leaf area to leaf weight, both species
bLA=0.9726;              %       leaf area to leaf weight, both species
cRT=0.1193;              %       tree biomass to dbh
bRT=2.393;               %       tree biomass to dbh
eta = 0.3;               %   allocate to root biomass
turnover_invasive = 0.34; % turnover rate for invasive
turnover_local = 0.85;    % turnover rate for local


c1 = 10;                 % coefficient for ZOI
c2 = 0.5;                % coefficient for ZOI
c3 = [0.35,0.5];         % dispersal coefficient, invasive species disperse farther

% Initially generate random location for each tree j.
rand('state',4);
Ntot = 0;                % Total number of trees
KNT = 0;

% initialize tree number
Ntree(1)=30;         % invasive
Ntree(2)=500;        % local

for i = 1: Nspec         % Iterate over the two species
%     Ntree(i)=30;         % Number of trees for species i, set to 10 %local tree
    Ntot = Ntot + Ntree(i);
    NL = KNT + 1;
    NU = KNT + Ntree(i);
    for j = NL: NU
        Ntype(j) = i;        % Identify each tree by species (1 or 2)
        xdot(j) = 100*abs(rand);  % Tree position along x axis,(meter)
        ydot(j) = 100*abs(rand);  % Tree position along y axis
        dbh(j) = 10.27;      % Initial d.b.h. cm
        Iage(j) = 0;         % Initial age
        alive(j) = 1;        % Flag for trees that are alive (=1 when alive)
    end
    KNT = NU;
end

xdot_ori=xdot;
ydot_ori=ydot;


% Find the neighbor of each tree, in order of nearest to furthest. 
% Only record distances lower than 5.0m. Assume trees at distance further than 5.0m do not compete 
for i = 1: Ntot % Iterate over all trees
    Neiber = i; % "Neiber" records ID number, first neighbor is itself
    dis = 0;    % "dis" record distance, nearest neighbor distance is zero
    for j = 1: Ntot  %Iterate over all other trees to find neighbors
        if i ~= j
            distemp = sqrt((xdot(i)-xdot(j))^2 + (ydot(i)-ydot(j))^2); %Distance of neighber trees 
            len = length(dis);  
            for k = len: -1: 1
                if distemp > 5.0 break; end
                if distemp > dis(k)
                    dis(k+1) = distemp;
                    Neiber(k+1) = j;
                    break
                else
                    dis(k+1) = dis(k);   %Store distance to neighbors
                    Neiber(k+1) = Neiber(k);  %Stor ID of neighbors
                end
            end
        end
    end
%Store information on IDs and distances of neighbors in two matrices
    Neib(1:length(dis),i) = Neiber;  % Record all vector info to matrix
    distance(1:length(dis),i) = dis;
end
Neib = sparse(Neib);          % Change to sparse matrix, store only non-zeros
distance = sparse(distance);  % Change to sparse matrix to save space

%Iterate over time in time steps of months
% time_i=1;



for t = 1: Nyear
   
%Next few statements are for debugging purposes    
tic

if t==400
    aaaa=1;
end

% %Matrices for root biomasses of two tree types in 3-D space
%     MB_root=zeros(20,20,5); % Biomass of roots for mangrove
%     HB_root=zeros(20,20,5); % biomass of root for hammock

    
%%%%%%%%%%%%%%%%%% Tree Growth Module  %%%%%%%%%%%%%%%%%%%%%

% This iterates over all trees and computes their growth
LF_grid1 = zeros(100,100);         %litter fall for invasive
LF_grid2 = zeros(100,100);         %litter fall for local
for j = 1: Ntot                    %Iterate over all trees
    Ntp = Ntype(j);                    %Species of the tree                     
    LW = cLW(Ntp)*dbh(j)^bLW(Ntp);     % Leaf mass of tree, grams
    LA = cLA*LW^bLA/8000;              % Leaf area of tree
    AL = exp(-phi*LA);                 % Available light
    rs = 1-exp(-4.64*(AL-0.05));       % Calculate shade-tolerant tree growth reduction
    if rs<0 rs=0; end
    ri=2.24*(1-exp(-1.136*(AL-0.08))); % Calculate shade-intolerant tree growth reduction
    if ri<0 ri=0; end
    
% Note that the model uses two horizontal coordinate systems.  One is 
% continuous,for the locations of trees.  The other coordinate system 
% is a 100 x 100 grid sytem for the vadose zone water and salinity
% dynamics. 
% We need to determine which spatial grid cell(s) each tree is in, so here
% we change continuous coordinates of tree locations to row-column matrix #.
% A tree may overlap more than one grid cell
    col = ceil(xdot(j)/5);       % Column location of tree
    row = 20 - ceil(ydot(j)/5) + 1;   % Row location of tree
    dz = deltz(row,col);   % Soil depth, mm (keep same with "vadose.m")
    
% We compute the amount of root biomass of each tree types in each
% spatial grid cell. This is needed in order to compute the water and
% salinity dynamics of each cell.
    B_root = eta*cRT*dbh(j)^bRT*1000;      % Root biomass, grams
    if Ntp==1
        % Melaleuca
        Max_rl = dbh(j)^2/(41.143+0.11789*dbh(j)^2);  % Maximum horizontal extension of roots, meters
%         Max_rd = Max_rl;                    % Maximum vertical extentsion of roots        
    else
        % Hammock
        Max_rl = dbh(j)^2/(41.143+0.15789*dbh(j)^2); % Maximum horizontal extension of roots
%         Max_rd = Max_rl*1.5; % Maximum vertical extentsion of roots
    end

% Calculate fraction of root biomass occupying each cell 5m*5m, "rootfrc.m"
% return cell identify number locii(row), locjj(column), layer number
% and Biomass of roots for each cell, sum of B_rootfcc should == B_root


%     [locii(time_i,j) locjj(time_i,j) B_rootfcc(time_i,j) LF(time_i,j)] = rootfrc(xdot(j),ydot(j),B_root,Max_rl,Ntp,dbh(j));
%   litter fall according to rootfrc function   
[locii locjj LF] = rootfrc(xdot(j),ydot(j),B_root,Max_rl,Ntp,dbh(j)); 
      for ii = 1:length(LF)
        if Ntp == 1
            LF_grid1(locii(ii),locjj(ii)) = LF_grid1(locii(ii),locjj(ii)) + LF(ii);
        else
            LF_grid2(locii(ii),locjj(ii)) = LF_grid2(locii(ii),locjj(ii)) + LF(ii);
        end
      end
      %from rootfrc.m to get locii locjj LF
    

%     tran = [];
%     for len = 1:length(layerr)
%         % Salinity of each cell and vertical layer, which roots occupy.
%         sa = Vs(locii(len),locjj(len),layerr(len)); 
%         if Ntp==1
%        % Water uptake (evapotranspiration) of this individual tree, 
%        % summed over the different cells it occupies. 
%         tran(len)=mantran(sa)*B_rootfcc(len);
%        % Sum up and store all the root biomass that occupies the same cell.
%         MB_root(locii(len),locjj(len),layerr(len))= MB_root(locii(len),locjj(len),layerr(len))+ B_rootfcc(len); 
%         else
%         tran(len)=hamtran(sa)*B_rootfcc(len);
%         HB_root(locii(len),locjj(len),layerr(len))= HB_root(locii(len),locjj(len),layerr(len))+ B_rootfcc(len); 
%         end
%     end
%     %This is the sum of evapotranspiration for a single tree, summed over
%     %all of the cells it occupies
%     tran = sum(tran);
%     %This is a multiplier, with the negative effect of salinity on growth.
%     r_salt = tran*2000/(2.6*B_root);% when sa=0,tran=2.6/2000
%     Tree height as a function of d.b.h.  
    Ht = 137+b2(Ntp)*dbh(j)-b3(Ntp)*dbh(j)^2;        % Height,cm
    if Ntp==1 r_light=rs; end
    if Ntp==2 r_light=ri; end
  
  % The "field of neighborhood" approach of Berger, Hildenbrandt, and Grimm
  % (2002) is used to simulate competition between neighboring trees.
  % Next, use function "fon.m" calculate neighbor competition, include shading by
  % neighbor and nutrient competition, et cetera
  % FS calculates the shading effects of all neighbors.
    FS_sum = 0;
        r2 = dbh(j)/200; % radius of tree, meters 
        KR = c1*(r2)^c2; % radius of zone of influence (ZOI) for object tree.
        NR_max = c1*(140/200)^c2; % potential neighbor
        di = distance(:,j); % load all neighbor distances
        if length(di)>1
            k = 2;
            while di(k)~=0 & di(k)< KR+NR_max
                r1 = dbh(Neib(k,j))/200; % r1 is neighbor tree
                FS_sum = FS_sum + fon(di(k),r1,KR);
                k = k + 1;
                if k >length(di)
                    break
                end
            end
        end
    FS_ave = FS_sum/(pi*KR^2);

    %FS average is used to calculate a multipier, r_fon, that reduces the
    %growth rate of a tree due to shading.
    r_fon = 1-2*FS_ave; % assume FS_ave bigger 0.5, FON multiplier=0;
    if r_fon<0
        %alive(j)=0;  
        r_fon = 0;
    end        
    
% This calculates the increase in tree diameter at each time step   
   
    DNC(j) = G(Ntp)*dbh(j)*(1-dbh(j)*Ht/(D_max(Ntp)*H_max(Ntp)))/ ...
        (274+3*b2(Ntp)*dbh(j)-4*b3(Ntp)*dbh(j)^2);
    r_salt = 1; % Assume no salt effect
%     r_light = 1; %Assume no light effect 
    DNC(j) = r_salt.*r_light.*r_fon.*DNC(j);
 
    dbh(j) = dbh(j) + DNC(j);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% LF accumulation %%%%%%%%%%%%%%%%%%

%LFaccumulation.m to get invasive and local leaf fall accumulation
%litter fall with turnover rate
LF_accumulation_invasive = (1-turnover_invasive)*LF_accumulation_invasive + LF_grid1;
LF_accumulation_local = (1-turnover_local)*LF_accumulation_local + LF_grid2;

%%%%%%%%%%%% information save %%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%% before Birth and dispersal %%%%%%%%%%%%

LF_annual_accumulation_invasive(t)={LF_accumulation_invasive};  % invasive litter information
LF_annual_accumulation_local(t)={LF_accumulation_local};        % local litter information
Tree_information(t,1)={1:Ntot'};    % ID information
Tree_information(t,2)={Ntype'};     % type information
Tree_information(t,3)={Iage'};      % age  information
Tree_information(t,4)={dbh'};       % dbh information
Tree_information(t,5)={xdot'};      % location information
Tree_information(t,6)={ydot'};      % location information


% loc_annual(t,1)={xdot};
% loc_annual(t,2)={ydot};


%%%%%%%%%%%%%%%%%%% Birth and dispersal %%%%%%%%%%%%%%%%%%%%

% This continues the iteration over all tree to calculate propagule
% production and dispersal to locations

Nnew = 0;           % This keeps track of new propagules
for j = 1 : Ntot    % Iterate over all existing tree to determine if they reproduce.
    Ntp = Ntype(j); % different species, different repreduction
    temp = rand;
    switch Ntp      % melaleuca
        case 1
            [x1,x2]=boundary(xdot(j),[1],[100],[2]);  % boundary prevent exceeding
            [y1,y2]=boundary(ydot(j),[1],[100],[2]);
            if (sum(sum(LF_accumulation_invasive([x1:x2],[y1:y2])))+sum(sum(LF_accumulation_local([x1:x2],[y1:y2])))) < 0.08
            brate = 0.028 * 2; % Reproductive probability each time step
            % depend on how thick the litter, less litter, more likely to
            % survive
            else
            brate = 0.015;
            end
%         end
        case 2      % hammock
            [x1,x2]=boundary(xdot(j),[1],[100],[2]);  % boundary prevent exceeding
            [y1,y2]=boundary(ydot(j),[1],[100],[2]);
            if (sum(sum(LF_accumulation_invasive([x1:x2],[y1:y2])))+sum(sum(LF_accumulation_local([x1:x2],[y1:y2])))) < 0.05
            brate = 0.0084 * 2;  % depend on how thick the litter, less litter, more likely to
            % survive
            else
            brate = 0.00001;
            end
    end
     
    
    if temp < brate
        for brate_i=1:2
        Nnew = Nnew + 1;
        %Ntree(Ntp) = Ntree(Ntp) + 1;
        Ntype(Ntot+Nnew)=Ntp; 
        dbh(Ntot+Nnew)=1.37;   %Initial d.b.h.
        Iage(Ntot+Nnew)=0;     %Initial age
        alive(Ntot+Nnew)=1;    %Assigned a value of 1, meaning alive

    % Compute the random dispersal distance, m 
    % If a propagules is dispersed out of the model area, it is excluded
    switch Ntp
        case 1
        disper = -log(rand)/c3(Ntp)*4;
        case 2
        disper = -log(rand)/c3(Ntp);
    end
%         while disper>10
%             disper = -log(rand)/c3;  % dispersal distance, m
%         end
        alpha = 2*pi*rand;
        xdot(Ntot+Nnew) = xdot(j) + disper*cos(alpha);  %Location on x-axis
        if xdot(Ntot+Nnew)<=0 | xdot(Ntot+Nnew)>100
            %xdot(Ntot+Nnew)=0;
            alive(Ntot+Nnew)=0; %If outside of model area, seedling excluded
        end
        ydot(Ntot+Nnew) = ydot(j) + disper*sin(alpha);  %Location on y-axis
        if ydot(Ntot+Nnew)<=0 | ydot(Ntot+Nnew)>100
            %ydot(Ntot+Nnew)=0;
            alive(Ntot+Nnew)=0;  %If outside of model area, seedling excluded
        end

     % Store new seedling in the Neib matrix; that is, the matrices of
     % IDs of neighbors and distances from neighbors
        disnew = 0;
        Neibernew = Ntot+Nnew;
        for i = 1 : Ntot+Nnew-1
            dis = distance(:,i);
            Neiber = Neib(:,i); 
            
            distemp = sqrt((xdot(i)-xdot(Ntot+Nnew))^2 + (ydot(i)-ydot(Ntot+Nnew))^2);
            len = length(dis);
            for k = len: -1: 1
                if distemp > 5.0 break; end
                if k==1
                    dis(k+1) = distemp;
                    Neiber(k+1) = Ntot+Nnew;
                    break
                end
                if dis(k)~=0
                    if distemp > dis(k)
                        dis(k+1) = distemp;
                        Neiber(k+1) = Ntot+Nnew;
                        break
                    else
                        dis(k+1) = dis(k);
                        Neiber(k+1) = Neiber(k);
                    end
                end
            end
             Neib_temp(1:length(dis),i) = Neiber;
             distance_temp(1:length(dis),i) = dis;
             
             lennew = length(disnew);
             for k = lennew:-1:1
                 if distemp > 5.0 break; end
                 if distemp > disnew(k)
                     disnew(k+1) = distemp;
                     Neibernew(k+1) = i;
                     break
                 else
                     disnew(k+1) = disnew(k);
                     Neibernew(k+1) = Neibernew(k);
                 end
             end
             
        end
        Neib = Neib_temp;
        distance = distance_temp;
        Neib(1:length(disnew),Ntot+Nnew) = Neibernew;
        distance(1:length(disnew),Ntot+Nnew) = disnew;
                Neib_temp =[];
                distance_temp=[];
            
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Tree Mortality %%%%%%%%%%%%%%%%%%

% This continues the iteration over all trees and allows some probability
% of mortality that depends on d.b.h.

for j = 1 : Ntot
    Ntp = Ntype(j);
    temp = rand;
    mc1 = 0.00102;                   % baseline mortality
    mc2 = 0.001*exp(-0.05*dbh(j));   % mortality decreases with size
    mc3 = 0.039*exp(-15.0*DNC(j));   % mortality due to poor growth rate
    mc4 = 0.0061*exp(0.0007*Ntot);        % density depedent death rate

    mc = mc1 + mc2 + mc3 + mc4;
%     mc = mc4
    if temp < mc
        alive(j) = 0;
        dbh(j) = -1;
    end
%     if dbh(j)<10.0 & DNC(j)/dbh(j) < 0.05
%         temp = rand;
%         if temp < 0.5
%             alive(j) = 0;
%             dbh(j) = -1;
%         end        
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Update Tree Information %%%%%%%%%%%%%%%%%%%

% This discards dead trees and updates the indexing of the populations

Ntot = Ntot + Nnew;
subx = find(alive==0); % Index of trees that are not alive
Ntot = Ntot - length(subx);
% delete info of trees that are not alive
Ntype(subx)=[];
xdot(subx)=[];
ydot(subx)=[];
dbh(subx)=[];
alive(subx)=[];
Iage(subx)=[];
Iage = Iage + 1;
DNC = 0;
 
% re-organize Neighbor trees
Neib(:,subx)=[];
distance(:,subx)=[];
%Neib = reshape(Neib,1,[]);
%distance = reshape(distance,1,[]);
for i = 1 : length(subx)
    suby = find(Neib==subx(i));
    y = size(Neib);
    while mod(suby,y(1)) ~= 0
        Neib(suby) = Neib(suby+1);
        distance(suby)=distance(suby+1);
        suby = suby+1;
    end
    Neib(suby)=0;
    distance(suby)=0;
    suby = find(Neib>subx(i));
    Neib(suby) = Neib(suby)-1;
    subx = subx - 1;
end
%Neib = reshape(Neib,[],Ntot);
%distance=reshape(distance,[],Ntot);

% Pass info to vadose zone water movement and return salinity 
%Vs = vadose(t,MB_root,HB_root);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time_i=time_i+1;
t
run_time=toc
end % end of t(year)
