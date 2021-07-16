%relax the initial supercell using quasi newton method
%hold transverse direction and optimize longitudinal direction for a fixed
%lattice constant

clear
clc

tic

%build a framework for the superlattice
n1=1;%number of 8-atom cells in InAs region
n2=1;%number of 8-atom cells in AlSb region
ndl=n1+n2;%total number of 8-atom cells
na=ndl*8;%number of atoms in each supercell
 
%fixed lattice constant of the supercell
a=5.80;
xlatt=a;
zlatt=a;
ylatt=ndl*a;%initial value
aa=a*0.05;

%mass in atomic mass unit "g/mol"
mAl=26.98154;
mIn=114.82;
mAs=74.9216;
mSb=121.75;


%load parameters updated 2021
% parameters2021_v4

%look up table for atoms in superlattice
mass=zeros(na,1);
type=zeros(na,1);

%assign atoms to each site
for n=1:n1
    mass(1+(n-1)*8)=mIn;
    mass(2+(n-1)*8)=mIn;
    mass(3+(n-1)*8)=mIn;
    mass(4+(n-1)*8)=mIn;
    
    type(1+(n-1)*8)=1;
    type(2+(n-1)*8)=1;
    type(3+(n-1)*8)=1;
    type(4+(n-1)*8)=1;
    
    mass(5+(n-1)*8)=mAs;
    mass(6+(n-1)*8)=mAs;
    mass(7+(n-1)*8)=mAs;
    mass(8+(n-1)*8)=mAs;
    
    type(5+(n-1)*8)=2;
    type(6+(n-1)*8)=2;
    type(7+(n-1)*8)=2;
    type(8+(n-1)*8)=2;
end

for n=(n1+1):(n1+n2)
    mass(1+(n-1)*8)=mAl;
    mass(2+(n-1)*8)=mAl;
    mass(3+(n-1)*8)=mAl;
    mass(4+(n-1)*8)=mAl;
    
    type(1+(n-1)*8)=3;
    type(2+(n-1)*8)=3;
    type(3+(n-1)*8)=3;
    type(4+(n-1)*8)=3;
    
    mass(5+(n-1)*8)=mSb;
    mass(6+(n-1)*8)=mSb;
    mass(7+(n-1)*8)=mSb;
    mass(8+(n-1)*8)=mSb;
    
    type(5+(n-1)*8)=4;
    type(6+(n-1)*8)=4;
    type(7+(n-1)*8)=4;
    type(8+(n-1)*8)=4;
end

clear("n")

%mass matrix to normalize dynamical matrix
mass3dn=zeros(3*na);
for i=1:na
    mass3dn(3*i-2,3*i-2)=1/sqrt(mass(i));
    mass3dn(3*i-1,3*i-1)=1/sqrt(mass(i));
    mass3dn(3*i,3*i)=1/sqrt(mass(i));
end

%relative coordinates of atoms in a 8-atom cell
%(set to default values, we can include lattice distortion here)
%(length scale on "a" of ffc)
tao=cell(8,1);
tao{1}=[0;0;0];
tao{2}=[1/2;1/2;0];
tao{3}=[1/2;0;1/2];
tao{4}=[0;1/2;1/2];
tao{5}=tao{1}+[1/4;1/4;1/4];
tao{6}=tao{2}+[1/4;1/4;1/4];
tao{7}=tao{3}+[1/4;1/4;1/4];
tao{8}=tao{4}+[1/4;1/4;1/4];

% %initialize neighboring site information
% nncount=zeros(na,1);
% nnR=cell(na,16);%nearest neighor unit cell coordinates
% nntype=zeros(na,16);%nearest neighbor atom type
% nnlatt=cell(na,16);%nearest neighbor lattice coordinates
% nncell=zeros(na,16);%nearest neighbor unit cell index

%initial position coodinates in supercell, applied 'a' scaled on Angstrom
RR=cell(na,1);
for i=1:ndl
    for p=1:8
        nindex=8*(i-1)+p;
        RR{nindex}=[0;a*(i-1);0]+a.*tao{p};        
    end
end