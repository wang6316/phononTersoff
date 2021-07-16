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
a=6.15;
aa=a*0.05;
xlatt=a;
ylatt=2*a;
zlatt=a;

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
tao{1}=[0,0,0];
tao{2}=[1/2,1/2,0];
tao{3}=[1/2,0,1/2];
tao{4}=[0,1/2,1/2];
tao{5}=tao{1}+[1/4,1/4,1/4];
tao{6}=tao{2}+[1/4,1/4,1/4];
tao{7}=tao{3}+[1/4,1/4,1/4];
tao{8}=tao{4}+[1/4,1/4,1/4];

%define intermediate variables for quasi Newton method
nf=3*na+3;%add three degree of freedom for lattice constant
xk=zeros(nf,1);%position vector on the run
Hkinv=eye(nf);%initialize inverse of the Hesian matrix as identity
sk=zeros(nf,1);%virtual displacement of current step
gradxk=zeros(nf,1);
gradxk1=zeros(nf,1);
qk=zeros(nf,1);%difference in gradient

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
        RR{nindex}=[0,a*(i-1),0]+a.*tao{p};
        xk((nindex-1)*3+1)=RR{nindex}(1);
        xk((nindex-1)*3+2)=RR{nindex}(2);
        xk((nindex-1)*3+3)=RR{nindex}(3);%initialize the xk vector, with y component of full 3D position 
    end
end

clear("i","p","nindex")

% energy=Etotal(xk,xlatt,ylatt,zlatt,na,type)

xk(3*na+1)=xlatt;
xk(3*na+2)=ylatt;
xk(3*na+3)=zlatt;%initial value
xk1=xk;
Hk1inv=Hkinv;

%iteration starts here
%using BFGS method to update Hessian
it=0;%counting iterations
go=true;

while go==true    %interation count
    it=it+1;    
   
%     %update corresponding 3D vector RR
%     for i=1:na
%         RR{i}(1)=xk((i-1)*3+1);
%         RR{i}(2)=xk((i-1)*3+2);
%         RR{i}(3)=xk((i-1)*3+3);
%     end
%     
    %update lattice constant in y direction
%     distance1=RR{10}(2)-RR{1}(2);
%     distance2=RR{9}(2)-RR{6}(2);
%     distance3=RR{14}(2)-RR{10}(2);
%     ylatt=distance1+distance2+distance3;
    
%      %update neighboring site information for this iteration
%      nncount=zeros(na,1);
%     for i=1:na%loop i runs over all atoms of the supercell 
%         inn=0;
%         %searching neighbors in 3*3 lattice with the supercell
%         for ix=-1:1:1
%             for iy=-1:1:1
%                 for iz=-1:1:1%lattice index of x,y,z               
%                     for j=1:na%supercell index for neighbor                   
%                         tempR=RR{j}+[xlatt*ix,ylatt*iy,zlatt*iz];
%                         distance=norm(tempR-RR{i});  %scaled on a                 
%                         if (distance>1.5 && distance<4.5)%loosely based on cutoff of Tersoff potential
%                             inn=inn+1;
%                             nnR{i,inn}=tempR;%position vector of this neighbor
%                             nncell(i,inn)=j;%unit cell index of neighbor j
%                             nnlatt{i,inn}=[xlatt*ix,ylatt*iy,zlatt*iz];%lattice coordinates of this neighbor
%                             nntype(i,inn)=type(j);%atom type of this neighbor
%                         end
%                     end
%                 end
%             end
%         end
%         nncount(i)=inn;
%     end

% calculate the derivative of vector d
    xlatt=xk(3*na+1);
    ylatt=xk(3*na+2);
    ylatt=xk(3*na+3);

parfor i=1:nf %runs over na atoms %runs over x,y,z
    xktemp1=xk;
    xktemp1(i)=xk(i)+a*0.01; %right displaced vector for derivative
    xktemp2=xk;
    xktemp2(i)=xk(i)-a*0.01; %left displaced vector of derivative
    %calculate ith elements of the gradient vector
    gradxk(i)=(Etotal(xktemp1,na,type)-Etotal(xktemp2,na,type))/(2*0.01*a);
end

Ek=Etotal(xk,na,type);

Et(it)=Ek;%keep track of total energy for each iteration

%walk using BFGS method 
sk=-Hkinv*gradxk;
xk1=xk+sk*aa/norm(sk);%sk normalized, aa controls the walk radius

if norm(gradxk)==0
    go=false;
    Efinal=Ek;
    break;
end

xlatt=xk1(3*na+1);
ylatt=xk1(3*na+2);
zlatt=xk1(3*na+3);


parfor i=1:nf
    xk1temp1=xk1;   
    xk1temp1(i)=xk1(i)+a*0.01;
    xk1temp2=xk1;
    xk1temp2(i)=xk1(i)-a*0.01;
    %calculate derivaitive at the displaced position
    gradxk1(i)=(Etotal(xk1temp1,na,type)-Etotal(xk1temp2,na,type))/(2*0.01*a);
end

Ek1=Etotal(xk1,na,type);

qk=gradxk1-gradxk;

Hk1inv=Hkinv+(1+(transpose(qk)*Hkinv*qk)/(transpose(sk)*qk))*(sk*transpose(sk))/(transpose(sk)*qk)-(sk*transpose(qk)*Hkinv+Hkinv*qk*transpose(sk))/(transpose(sk)*qk);

   if(abs(Ek-Ek1)<10^-8)        
        Efinal=Ek
        go=false;
        it
        break;
   else
       %update vector xk and inverse of Hessian
        xk=xk1;
        Hkinv=Hk1inv;
        %change walk radius if go across the minimum, direction is
        %automatically altered in sk
        if((Ek1-Ek)>0)
            aa=aa/2;
        end
    
   end
end

%write the relaxed xk vector back in the form of na*3d atom positions
for ia=1:na
    RR{ia}=xk(3*ia-2:3*ia);
end

xlattice=xk(3*na+1);
ylattice=xk(3*na+2);
zlattice=xk(3*na+3);
       
toc