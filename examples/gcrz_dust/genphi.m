%fid=fopen('data/g067238.06118','r');
%fid=fopen('data/g1100223012.01149_261','r');
fid=fopen('data/g056129.05550_1','r');
textscan(fid,'%*s',6);
nxefit=fscanf(fid,'%d',[1]);
nyefit=fscanf(fid,'%d',[1]);
gdata=fscanf(fid,'%e',[20]);
fclose(fid);

%%This mesh is generated the same as gfile for magnetic field which vary slowly in RZ domain
xdim=gdata(1);
zdim=gdata(2);
rgrid1=gdata(4);
zmid=gdata(5);
gridr=zeros(1,nxefit);
gridz=zeros(1,nyefit);
for i=1:nxefit
	gridr(i) = (rgrid1 + xdim*(i-1)/(nxefit-1));
end
for j=1:nyefit
	gridz(j) = ((zmid-0.5*zdim)+zdim*(j-1)/(nyefit-1));
end
[meshr,meshz]=meshgrid(gridr,gridz);

%%However, electric potential vary fast in RZ domain, and need fine mesh
rmin=1.3;
rmax=2.3;
zmin=-1.2;
zmax=0.8;
nr_fine=601;
nz_fine=1201;
dr=(rmax-rmin)/(nr_fine-1);
dz=(zmax-zmin)/(nz_fine-1);
for i=1:nr_fine
	gridr_fine(i) = rmin+dr*(i-1);
end
for j=1:nz_fine
	gridz_fine(j) = zmin+dz*(j-1);
end
[meshr_fine,meshz_fine]=meshgrid(gridr_fine,gridz_fine);

fid=fopen('data/prodata.dat','r');
TTN=fscanf(fid,'%e',[6,inf]);
fclose(fid);

TTN1=griddata(TTN(1,:),TTN(2,:),TTN(3,:),meshr_fine,meshz_fine,'v4');
TTN2=griddata(TTN(1,:),TTN(2,:),TTN(4,:),meshr_fine,meshz_fine,'v4');
TTN3=griddata(TTN(1,:),TTN(2,:),TTN(5,:),meshr_fine,meshz_fine,'v4');
TTN4=griddata(TTN(1,:),TTN(2,:),TTN(6,:),meshr_fine,meshz_fine,'v4');

fid=fopen('data/ttn1.dat','w');
fprintf(fid,'%20.10e %20.10e %20.10e %20.10e %20.10e\n',TTN1);
fclose(fid);

fid=fopen('data/ttn2.dat','w');
fprintf(fid,'%20.10e %20.10e %20.10e %20.10e %20.10e\n',TTN2);
fclose(fid);

fid=fopen('data/ttn3.dat','w');
fprintf(fid,'%20.10e %20.10e %20.10e %20.10e %20.10e\n',TTN3);
fclose(fid);

fid=fopen('data/ttn4.dat','w');
fprintf(fid,'%20.10e %20.10e %20.10e %20.10e %20.10e\n',TTN4);
fclose(fid);

figure,surf(meshr_fine,meshz_fine,TTN1);
