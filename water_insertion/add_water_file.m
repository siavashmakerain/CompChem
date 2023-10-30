function [Coords,databond,dataangles,n_wat,natom]=add_water_file(fid,nlayer,nbonds,nangles,co2_conc,natom,wat_dist,aa2,bb2,ab_thet,vc)
num_wat=208; % For reax use 64 with 14 layers, for comb 47 with 19 layers

mydata1 = textscan(fid,'%s',natom,'headerlines',20,'delimiter', '\n');
%databond= textscan(fid,'%s',nbonds,'headerlines',6+natom,'delimiter', '\n');
%dataangles= textscan(fid,'%s',nangles,'headerlines',3,'delimiter', '\n');
%databond=databond{1};
%databond= cellfun(@(x)str2num(x),databond,'UniformOutput',false);
%databond = cell2mat(databond);
%dataangles=dataangles{1};
%dataangles= cellfun(@(x)str2num(x),dataangles,'UniformOutput',false);
%dataangles = cell2mat(dataangles);
% Cells above are commented because the substrate does not have bonds or
% angles
databond=[];
dataangles=[];
mydata1 = mydata1{1};
mydata2 = cell2mat(cellfun(@(x)str2num(x),mydata1,'UniformOutput',false));
mydata2=sortrows(mydata2,1);
%databond=sortrows(databond,1);
% Refine the condition for every atom on top
% for i=1:length(mydata2)
% %    if i==1377
%     if mydata2(i,7)>= 39
%         if mydata2(i,3) ==14 ||  mydata2(i,3) ==15 || mydata2(i,3) ==1 || mydata2(i,3) ==8 || mydata2(i,3) ==6 || mydata2(i,3) ==9
%         mydata2(i,5:7)=mydata2(i,5:7)-vc';
%         end
%  
%     end
%  %   end
% end
Coords = mydata2(:,4:6);
nmol=1;
% Convert to ReaxFF atom types
type_reax=zeros(length(mydata2),1);
for i=1:length(mydata2)
    if mydata2(i,3)==3 || mydata2(i,3)==2 || mydata2(i,3)==10
        type_reax(i)=8;
    end
    if mydata2(i,3)==1 
        type_reax(i)=6;
    end
    if mydata2(i,3)==4 || mydata2(i,3)==7 || mydata2(i,3)==8 || mydata2(i,3)==9  
        type_reax(i)=3;
    end
    if mydata2(i,3)==5 || mydata2(i,3)==6  
        type_reax(i)=2;
    end
end

%mydata2(all_rem,:)=[];
%mydata2=sortrows(mydata2,1);

natom=size(Coords,1);
%% Water insertion

%aa2 = aa2'; bb2 = bb2';

[z_top,z_top_id] = max(Coords(:,3));
z_top=28;
wat_space = 2.5; % the approximate space assigned for a water molecule
n_adiv = round(norm(aa2)/wat_space); %num of divisions in the a direction
n_bdiv = round(norm(bb2)*sind(ab_thet)/wat_space);

binx = norm(aa2)/n_adiv;
biny = norm(bb2)*sind(ab_thet)/n_bdiv;
bin = [binx,biny,0];
%choose the cells randomly
[nx,ny] = meshgrid(0:n_adiv-1,0:n_bdiv-1);
n_tot = [nx(:),ny(:),zeros(length(nx(:)),1)];
chosen_top = randsample(1:length(nx(:)),num_wat);%the chosen spaces for water

%make the coordinates for the zero'th cell
wat_coord(1,1) = binx/2;
wat_coord(1,2) = biny/2-0.5859/2;
wat_coord(1:3,3) = z_top + 2;
wat_coord(2,1) = wat_coord(1,1) - 0.757;
wat_coord(3,1) = wat_coord(1,1) + 0.757;
wat_coord(2:3,2) = wat_coord(1,2) + 0.5859;

for iwat = 1:length(chosen_top)
    cwat = n_tot(chosen_top(iwat),:);%current water
    wat_coords([3*iwat-2,3*iwat-1,3*iwat],:) = wat_coord + repmat(cwat.*bin,3,1);
%     elem_layer(end+1:end+3) = {'Ow','Hw','Hw'};
end


for i = 1:nlayer
    tmp = wat_coords;
    tmp(:,3) = tmp(:,3) + (i-1)*wat_dist;
    Coords = [Coords;tmp];
    
end
n_wat=1/3*nlayer*length(wat_coords);
%bonds_csh=0;
angles_csh=0;
% Replacing CO2 based on concentration
no_co2=co2_conc/(co2_conc+1)*n_wat;
if co2_conc==0
    every=1000000000;
else
every=round(n_wat/no_co2);
end
natom1=natom;
for i=1:n_wat
    if mod(i,every)==0
        databond(bonds_csh+2*(i-1)+1,:)=[bonds_csh+2*(i-1)+1 3 natom1+3*(i-1)+1 natom1+3*(i-1)+2];
        databond(bonds_csh+2*(i-1)+2,:)=[bonds_csh+2*(i-1)+2 3 natom1+3*(i-1)+1 natom1+3*i];
    else
        databond(bonds_csh+2*(i-1)+1,:)=[bonds_csh+2*(i-1)+1 1 natom1+3*(i-1)+1 natom1+3*(i-1)+2];
        databond(bonds_csh+2*(i-1)+2,:)=[bonds_csh+2*(i-1)+2 1 natom1+3*(i-1)+1 natom1+3*i];
    end
end
for i=1:n_wat
    if mod(i,every)==0
    dataangles(angles_csh+i,:)=[angles_csh+i 2 natom1+3*(i-1)+2 natom1+3*(i-1)+1 natom1+3*i];
    else
      dataangles(angles_csh+i,:)=[angles_csh+i 1 natom1+3*(i-1)+2 natom1+3*(i-1)+1 natom1+3*i];
    end
end
for i=natom+1:3:length(Coords)
    if mod(i-natom,every)==0
    type_reax(i)=1;type_reax(i+1)=3;type_reax(i+2)=3;
    else
        type_reax(i)=3;type_reax(i+1)=2;type_reax(i+2)=2;
    end
end
%file2 = fopen('csh170_pps_ali.txt','w');
file4 = fopen('csh170_water_pps_ali_reax.txt','w');
file4i = fopen('csh170_water_pps_ali_cshff.txt','w');
file5 = fopen('bonds_csh170_water_pps_ali_reax.txt','w');
file6 = fopen('angles_csh170_water_pps_ali_reax.txt','w');
coords_csh=mydata2(:,1:6);
coords_csh_mol=[];
for i=1:size(coords_csh)
    nmol=nmol+1;
   coords_csh_mol(i,:)=[coords_csh(i,1) nmol coords_csh(i,2:6)];
end
coords_wc=[];
%nmol=natom;
for i=1:n_wat
    count=natom1+(i-1)*3+1;
    if mod(i,every)==0
        coords_wc=[coords_wc;count nmol+i 11 0.6512 Coords(natom+3*(i-1)+1,:);...
            count+1 nmol+i 12  -0.3256 Coords(natom+3*(i-1)+1,1)-1.149 Coords(natom+3*(i-1)+1,2) Coords(natom+3*(i-1)+1,3);...
            count+2 nmol+i 12  -0.3256 Coords(natom+3*(i-1)+1,1)+1.149 Coords(natom+3*(i-1)+1,2) Coords(natom+3*(i-1)+1,3)];
    else
       % coords_wc=[coords_wc;count 3 -0.8476 Coords(natom+3*(i-1)+1,:);...
       %     count+1  2   0.4238 Coords(natom+3*(i-1)+2,:);...
       %     count+2  2   0.4238 Coords(natom+3*(i-1)+3,:)];
        coords_wc=[coords_wc;count nmol+i 3 -1.1794 Coords(natom+3*(i-1)+1,:);...
            count+1 nmol+i 4   0.5897 Coords(natom+3*(i-1)+2,:);...
            count+2 nmol+i 4   0.5897 Coords(natom+3*(i-1)+3,:)];
    end
end
%for i=1:length(coords_wc)
%    coords_wc(i,6)=coords_wc(i,6)-1.5;
%end
Coords=[coords_csh_mol;coords_wc];
%Coords=coords_wc;
%Coords=[coords_csh;coords_wc];




%count=0;
% for i = 1:length(Coords)
%     count=count+1;
%    
%     fprintf(file4,'%4.0f\t  %4.0f\t 0.0 %4.8f\t %4.8f\t  %4.8f\t\r\n',count,type_reax(i),Coords(i,:));
% end
% for i=1:natom
%     fprintf(file4i,'%4.0f\t  %4.0f\t %4.0f\t %4.8f\t %4.8f\t %4.8f\t  %4.8f\t\r\n',mydata2(i,1:7));
% end
% for i=1:n_wat
%     count=natom+(i-1)*3+1;
%     if mod(i,every)==0
%         fprintf(file4i,'%4.0f\t  %4.0f\t 11 0.6512 %4.8f\t %4.8f\t  %4.8f\t\r\n',count,nmol+i,Coords(natom+3*(i-1)+1,:));
%         fprintf(file4i,'%4.0f\t  %4.0f\t 12  -0.3256 %4.8f\t %4.8f\t  %4.8f\t\r\n',count+1,nmol+i,Coords(natom+3*(i-1)+2,:));
%         fprintf(file4i,'%4.0f\t  %4.0f\t 12  -0.3256 %4.8f\t %4.8f\t  %4.8f\t\r\n',count+2,nmol+i,Coords(natom+3*(i-1)+3,:));
%     else
%         fprintf(file4i,'%4.0f\t  %4.0f\t 4 -0.82 %4.8f\t %4.8f\t  %4.8f\t\r\n',count,nmol+i,Coords(natom+3*(i-1)+1,:));
%         fprintf(file4i,'%4.0f\t  %4.0f\t 5  0.41 %4.8f\t %4.8f\t  %4.8f\t\r\n',count+1,nmol+i,Coords(natom+3*(i-1)+2,:));
%         fprintf(file4i,'%4.0f\t  %4.0f\t 5  0.41 %4.8f\t %4.8f\t  %4.8f\t\r\n',count+2,nmol+i,Coords(natom+3*(i-1)+3,:));
%     end
% end
% for i=1:length(databond)
%     fprintf(file5,'%4.0f\t %4.0f\t %4.0f\t %4.0f\t\n',databond(i,:));
% end
% for i=1:length(dataangles)
%     fprintf(file6,'%4.0f\t %4.0f\t %4.0f\t %4.0f\t %4.0f\t\n',dataangles(i,:));
% end