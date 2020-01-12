 function [locii_cell locjj_cell LF_cell] = rootfrc(xdot,ydot,B_root,Max_rl,Ntp,dbh)

% calculate fraction of root biomass occupy current cell 5m*5m
% return cell identify number locii(row), locjj(column), layer number
% and Biomass of roots for each cell, sum of B_rootfcc should == B_root


Max_rd=1000;
dz=1000;
locii = [];
locjj = [];
% layerr = [];
B_rootfcc = [];
LF=[];

if Ntp==1 %invasive
    LFi=0.0055*dbh; %leaf fall equation
elseif Ntp==2 %local
    LFi=0.0005*dbh; %leaf fall equation
end


    B_v0 = B_root*5.776/(1-exp(-5.776*Max_rd)); % Surface root biomass
    B_rooti = [0 0 0 0 0];
    depth=dz;
    layer = 1;
%     % while loop calculate root biomass for each layer
%     while depth < Max_rd 
%         B_rooti(layer) = B_v0*(exp(-5.776*(depth-dz))-exp(-5.776*depth))/5.776;
%         depth = depth + dz;
%         layer = layer + 1;
%         if layer==5
%             break
%         end
%     end
%     if Max_rd>dz*5 Max_rd=dz*5; end
    B_rooti(layer) = B_v0*(exp(-5.776*(depth-dz))-exp(-5.776*Max_rd))/5.776;
%     
%     % calculate root extention range for each layer, surface is Max_rl
    alph = B_rooti(1)/(Max_rl*Max_rl);
    for i = 1:1
        radius(i) = sqrt(B_rooti(i)/alph);
    end
    
    % calculate fractional root biomass in each cell, B_rootfcc
    % assume root extension range is rectangle, percentage of root cut by
    % cell column is 'tempj' from 'c_left' to 'c_right'; percentage of root
    % cut by cell row is 'tempi' from 'r_up' to 'r_down'. 
    
    k=1;
 while B_rooti(k) ~= 0

    % Extention range
       
    if Max_rl<=2.5
    radius(k)=Max_rl; 
    else 
    radius(k)=2.5;
    end
        
    left = xdot - radius(k);% simplified the parameter
    if left<0 left=0; end
    right = xdot + radius(k);
    if right>100 right=100; end
    up = ydot + radius(k);
    if up>100 up=100; end
    down = ydot - radius(k);
    if down<0 down=0; end
    
    % change x-y coordinate to cell identify row-column number
%%%%%%%%%%%%%%%%%%%%%%%%    
    %c_left = ceil(left/5);
    c_left = ceil(left);
    if c_left<1 c_left=1; end
    c_right = ceil(right);    
    r_up = 100 - ceil(up) + 1;
    r_down = 100 - ceil(down) + 1;
    if r_down>100 r_down=100; end
        left_temp = left;        
        for j=c_left:c_right
            %right_temp = j*5;
            right_temp = j;
            if right_temp>right right_temp=right; end
            tempj = (right_temp-left_temp)/(2*radius(k));
            left_temp = right_temp;   
            
            up_temp = up;
            for i=r_up:r_down
                locii = [locii; i];
                locjj = [locjj; j];
%                 layerr = [layerr; k];
                %down_temp = (20-i)*5;
                down_temp = (100-i);
                if down_temp<down down_temp=down; end
                tempi = (up_temp-down_temp)/(2*radius(k));
                up_temp = down_temp;
                LF=[LF;tempj*tempi*LFi(k)];
                B_rootfcc = [B_rootfcc; tempj*tempi*B_rooti(k)];
            end
        end
        k = k + 1;
        if k>5 break; end
        locii_cell=locii; % using cell to save locii locjj
        locjj_cell=locjj;
%         B_rootfcc_cell={B_rootfcc};
        LF_cell=LF;
 end
 