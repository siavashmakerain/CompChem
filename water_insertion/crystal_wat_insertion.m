clear all;
%Cell parameters
xlo=0;
xhi=36.792;
ylo=0;
zhi=27;
zlo=0;
yhi=63.72561331;
xy=0;
xz=0;
yz=0;
v_c=[xz;yz;zhi-zlo];
aa2=[xhi-xlo;0];
bb2=[xy;yhi-ylo];
ab_thet = acosd(dot(aa2,bb2)/norm(aa2)/norm(bb2));
wat_dist = 2.5; 

% no. of bonds, angles, atoms,
nbonds=0;
nangles=0;
natom=4096;
co2_conc=[0.0 0.01 0.02 0.05 0.1 0.15];   % CO2 concentration
co2_num=[3 4];   % CO2 concentration
co2_const=[0.33 0 0 0.12];   % CO2 concentration
nlayer = [13]; % MAIN LAYERS
%nlayer = [13 16 27 33]; % EXTENDED POINTS
%nlayer = [7 8]; % EXTENDED POINTS 6-9
cnt2=0;
num_wat=64;
for i=1:length(co2_num)
     num_co2=co2_num(i);
     for j=1:length(nlayer)
fid = fopen('data.agi_slab-35','r');
co2i=0.0;
cnt2=cnt2+1;
nlayerj=nlayer(j);
%[Coords,databond,dataangles,n_co2(cnt2),natom]=add_co2(fid,nlayerj,nbonds,nangles,co2i,natom,wat_dist,aa2,bb2,ab_thet,v_c);
%[Coords,databond,dataangles,n_wat]=add_co2(fid,nlayer,nbonds,nangles,co2_conc,natom,wat_dist,aa2,bb2,ab_thet,num_wat);
[Coords,databond,dataangles,n_wat,natom]=add_water_file(fid,nlayer,nbonds,nangles,co2_conc(1),natom,wat_dist,aa2,bb2,ab_thet,v_c);
        zmax=max(Coords(:,7))+1.5;
        zhip=zmax+5.5;
        cratio=(zhip-zlo)/(zhi-zlo);
        v_cp=v_c*cratio;
        
       
         txtfile=['AgI-0001-tip4p-implicit-' num2str(nlayerj) '-' num2str(num_co2) '.data'];
         txtfile_2=['forst_co2-cell-' num2str(nlayerj) '-' num2str(num_co2) '.data'];
        file=fopen(txtfile,'w');
        file2=fopen(txtfile_2,'w');
        fprintf(file,'\n');
        fprintf(file,'%.0f atoms \r\n',length(Coords));
        fprintf(file,'5 atom types \r\n');
        fprintf(file,'%.0f bonds \r\n',length(databond));
        fprintf(file,'2 bond types \r\n');
        fprintf(file,'%.0f angles \r\n',length(dataangles));
        fprintf(file,'1 angle types \r\n');
        
        fprintf(file,'\n');
        fprintf(file,'%4.8f %4.8f xlo xhi\r\n',xlo,xhi);
        fprintf(file,'%4.8f %4.8f ylo yhi\r\n',ylo,yhi);
        fprintf(file,'%4.8f %4.8f zlo zhi\r\n',zlo,120);
        fprintf(file,'%4.8f %4.8f %4.8f xy xz yz\r\n',xy,v_cp(1),v_cp(2));
        fprintf(file,'\n');
        fprintf(file,'Masses \r\n');
        fprintf(file,'\n');
          fprintf(file,['1 107.868 \r\n',...
            '2 126.904 \r\n',...
            '3 15.9994 \r\n',...
            '4 1.008\r\n',...
            '5 1.0e-100\r\n']);
        fprintf(file,'\n');
        fprintf(file,'Bond Coeffs # harmonic \r\n');
        fprintf(file,'\n');
        fprintf(file,'1 554.135 0.9572 \r\n');
        fprintf(file,'2 554.135 0.125 \r\n');
       % fprintf(file,'3 2017.89000 1.162000 \r\n');
        fprintf(file,'\n');
        fprintf(file,'Angle Coeffs # harmonic \r\n');
        fprintf(file,'\n');
        fprintf(file,'1 45.7696 104.52 \r\n');
      
        fprintf(file,'\n');
        fprintf(file,'Atoms \r\n');
        fprintf(file,'\n');
        
       % cah=find(Coords(:,3)==3);Coords(cah,3)=13;
       % ohp=find(Coords(:,3)==1);Coords(ohp,3)=9;
       % hp=find(Coords(:,3)==2);Coords(hp,3)=6;
        for k=1:length(Coords)
            fprintf(file,'%.0f %.0f %.0f %4.8f %4.8f  %4.8f %4.8f\r\n',Coords(k,1:7));
        end
        fprintf(file,'\n');
        fprintf(file,'Bonds \r\n');
        fprintf(file,'\n');
        for k=1:length(databond)
            fprintf(file,'%4.0f\t %4.0f\t %4.0f\t %4.0f\t\n',databond(k,:));
        end
         fprintf(file,'\n');
         fprintf(file,'Angles \r\n');
         fprintf(file,'\n');
         for k=1:length(dataangles)
            fprintf(file,'%4.0f\t %4.0f\t %4.0f\t %4.0f\t %4.0f\t\n',dataangles(k,:));
        end
        fclose(file);
        
        % Printing cell information to file
        v_ap=[xhi-xlo 0 0];
        v_bp=[xy yhi-ylo 0];
        
            fprintf(file2,'%4.8f %4.8f  %4.8f\r\n',v_ap(:));
            fprintf(file2,'%4.8f %4.8f  %4.8f\r\n',v_bp(:));
            fprintf(file2,'%4.8f  %4.8f %4.8f\r\n',v_cp(:));
 fclose(file2);
        fclose(fid);
        Coords=sortrows(Coords,1);
selems = {};
selems (Coords(:,3)==1) = {'st'};
selems (Coords(:,3)==2) = {'cao'};
selems (Coords(:,3)==3) = {'cao'};
selems (Coords(:,3)==4) = {'o*'};
selems (Coords(:,3)==5) = {'h*'};
selems (Coords(:,3)==6) = {'ho'};
selems (Coords(:,3)==7) = {'ob'};
selems (Coords(:,3)==8) = {'obts'};
selems (Coords(:,3)==9) = {'oh'};
selems (Coords(:,3)==10) = {'cao'};
selems (Coords(:,3)==11) = {'C_EPM2'};
selems (Coords(:,3)==12) = {'O_EPM2'};
selems (Coords(:,3)==13) = {'cah'};
selems (Coords(:,3)==14) = {'mg'};
selems (Coords(:,3)==15) = {'ohs'};
file=fopen('towhee_coords','w');
for k=1:length(Coords)
    fprintf(file,'%4.8f  %4.8f %4.8f\r\n',Coords(k,5:7));
end
fclose(file);
txtfile=['towhee_connectivity_ch-' num2str(nlayerj) '-' num2str(num_co2) '.data'];
file=fopen(txtfile,'w');
cnt=0;
cntmol=0;
for k=1:natom
    cnt=cnt+1;
    %fprintf(file,'\n');
    fprintf(file,'unit ntype \r\n');
    %fprintf(file,'\n');
    fprintf(file,'%4.0f ''%s '' %4.8f\r\n',Coords(k,1),selems{k},Coords(k,4));
    if Coords(k,3)==4 && Coords(k,1) > natom
        cntmol=cntmol+1;
        fprintf(file,'vibration\r\n');
        fprintf(file,'2\r\n');
        bondatoms=find(databond(:,3)==Coords(k,1));
        fprintf(file,'%4.0f %4.0f\r\n',databond(bondatoms(1),4),databond(bondatoms(2),4));
        fprintf(file,'improper torsion\r\n');
        fprintf(file,'0\r\n');
    elseif Coords(k,3)==5 && Coords(k,1) > natom
        fprintf(file,'vibration\r\n');
        fprintf(file,'1\r\n');
        bondatoms=find(databond(:,4)==Coords(k,1));
        fprintf(file,'%4.0f\r\n',databond(bondatoms,3));
        fprintf(file,'improper torsion\r\n');
        fprintf(file,'0\r\n');
    elseif Coords(k,3)==11 && Coords(k,1) > natom
        cntmol=cntmol+1;
        fprintf(file,'vibration\r\n');
        fprintf(file,'2\r\n');
        bondatoms=find(databond(:,3)==Coords(k,1));
        fprintf(file,'%4.0f %4.0f\r\n',databond(bondatoms(1),4),databond(bondatoms(2),4));
        fprintf(file,'improper torsion\r\n');
        fprintf(file,'0\r\n');
    elseif Coords(k,3)==12 && Coords(k,1) > natom
        fprintf(file,'vibration\r\n');
        fprintf(file,'1\r\n');
        bondatoms=find(databond(:,4)==Coords(k,1));
        fprintf(file,'%4.0f\r\n',databond(bondatoms,3));
        fprintf(file,'improper torsion\r\n');
        fprintf(file,'0\r\n');
    else
        fprintf(file,'vibration\r\n');
        fprintf(file,'0\r\n');
        fprintf(file,'improper torsion\r\n');
        fprintf(file,'0\r\n');
    end
    
end
fclose(file);
txtfile=['towhee_connectivity_wat-' num2str(nlayerj) '-' num2str(num_co2) '.data'];
file=fopen(txtfile,'w');
cnt=0;
cntmol=0;

for k=natom+1:length(Coords)
    cnt=cnt+1;
    %fprintf(file,'\n');
    fprintf(file,'unit ntype \r\n');
    %fprintf(file,'\n');
    fprintf(file,'%4.0f ''%s '' %4.8f\r\n',Coords(k,1)-natom,selems{k},Coords(k,4));
    if Coords(k,3)==4 && Coords(k,1) > natom
        cntmol=cntmol+1;
        fprintf(file,'vibration\r\n');
        fprintf(file,'2\r\n');
        bondatoms=find(databond(:,3)==Coords(k,1));
        fprintf(file,'%4.0f %4.0f\r\n',databond(bondatoms(1),4)-natom,databond(bondatoms(2),4)-natom);
        fprintf(file,'improper torsion\r\n');
        fprintf(file,'0\r\n');
    elseif Coords(k,3)==5 && Coords(k,1) > natom
        fprintf(file,'vibration\r\n');
        fprintf(file,'1\r\n');
        bondatoms=find(databond(:,4)==Coords(k,1));
        fprintf(file,'%4.0f\r\n',databond(bondatoms,3)-natom);
        fprintf(file,'improper torsion\r\n');
        fprintf(file,'0\r\n');
    else
        fprintf(file,'vibration\r\n');
        fprintf(file,'0\r\n');
        fprintf(file,'improper torsion\r\n');
        fprintf(file,'0\r\n');
    end
    
end
fclose(file);
     end
end

file=fopen('num_h2o.txt','w');
fprintf(file,'%4.0f',n_co2(:));