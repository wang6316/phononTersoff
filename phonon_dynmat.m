%following relax_newton_3d
%main routing calculating phonon structrue of the superlattice using
%Tersoff potential

tic

%evaluate interatomic force constant using Tersoff potential 
%fill in informations about 4 nearest neighbors and 12 next nearest neighbors 

KblockAll=cell(na,16);

for i=1:na%loop i runs over all atoms of the supercell 
    inn=0;
    %searching neighbor in 3*3 lattice with supercell
    for ix=-1:1:1
        for iy=-1:1:1
            for iz=-1:1:1%lattice index of x,y,z               
                for j=1:na%supercell index for neighbor                   
                    tempR=RR{j}+[xlatt*ix;ylatt*iy;zlatt*iz];%scaled on 'a'
                    distance=norm(tempR-RR{i});                   
                    if (distance>0.4 && distance<0.75)
                        inn=inn+1;
                        nnR{i,inn}=a.*tempR;%position vector of this neighbor
                        nncell(i,inn)=j;%unit cell index of neighbor j
                        nnlatt{i,inn}=[xlatt*ix,ylatt*iy,zlatt*iz];%lattice coordinates of this neighbor
                        nntype(i,inn)=type(j);%atom type of the neighbor
                    end
                end
            end
        end
    end
end

%load parameters 
% parametersNOV19
parameters2021_v4

%finite step for numerical derivative
dh=0.01*a;

%five-point stencil 
fivepoint=[-2*dh,-dh,0,dh,2*dh];



% with the information above, finding force constant using Tersoff potential
%use five-point stencil to calculate second order derivative
for i=1:na%loop for possible atom i, runs over a supercell
    icoord0=a.*RR{i};%asign realistic value, length in angstrom
    itype=type(i);
    for j=1:16 %loops for 16 possible neighbors j of  atom i
        jcoord0=nnR{i,j};
        jtype=nntype(i,j);
        %above choose i,j in supercell
        %below create 5*5 grid for i_alpha,j_beta
        
        Kblock=zeros(3,3);%3*3 block for di_alpha,dj_beta 
       
        for dim1=1:3%alpha runs over x,y,z
            for dim2=1:3%beta runs over j_x,y,z
                Eij=zeros(5);
%                 der_j=zeros(5,1);%grid for di_alpha
                for ifive1=1:5%five points for alpha
                    if dim1==1
                        icoord=icoord0+[fivepoint(ifive1),0,0];
                    elseif dim1==2
                        icoord=icoord0+[0,fivepoint(ifive1),0];
                    elseif dim1==3
                        icoord=icoord0+[0,0,fivepoint(ifive1)];
                    end
                        
                    
                    for ifive2=1:5%five points for beta
                        if dim2==1
                            jcoord=jcoord0+[fivepoint(ifive2),0,0];
                        elseif dim2==2
                            jcoord=jcoord0+[0,fivepoint(ifive2),0];
                        elseif dim2==3
                            jcoord=jcoord0+[0,0,fivepoint(ifive2)];
                        end
        
                        vecij=jcoord-icoord;
                        rij=norm(vecij);
        

                        %assign parameters for certain combinations of i,j
                        if(itype==1 && jtype==1)
                            ij_InIn
                        elseif(itype==2 && jtype==2)
                            ij_AsAs
                        elseif(itype==3 && jtype==3)
                            ij_AlAl
                        elseif(itype==4 && jtype==4)
                            ij_SbSb
                        elseif((itype==1 && jtype==2)||(itype==2 && jtype==1))
                            ij_InAs
                        elseif((itype==1 && jtype==3)||(itype==3 && jtype==1))
                            ij_InAl
                        elseif((itype==1 && jtype==4)||(itype==4 && jtype==1))
                            ij_InSb
                        elseif((itype==2 && jtype==3)||(itype==3 && jtype==2))
                            ij_AlAs
                        elseif((itype==2 && jtype==4)||(itype==4 && jtype==2))
                            ij_AsSb
                        elseif((itype==3 && jtype==4)||(itype==4 && jtype==3))
                            ij_AlSb                   
                        end
        
                        xi=0;%initialize xi for all 3-body correction
        
                        for k=1:16%loop for 16 possible neighbors k of atom i
                            if(k~=j)               
                                kcoord=nnR{i,k};
                                ktype=nntype(i,k);
                                vecik=kcoord-icoord;%fix ik in derivative
                                rik=norm(vecik);
                
                                %assign parameters for certain combinations of i,k
                                if(itype==1 && ktype==1)
                                    ik_InIn
                                elseif(itype==2 && ktype==2)
                                    ik_AsAs
                                elseif(itype==3 && ktype==3)
                                    ik_AlAl
                                elseif(itype==4 && ktype==4)
                                    ik_SbSb
                                elseif((itype==1 && ktype==2)||(itype==2 && ktype==1))
                                    ik_InAs
                                elseif((itype==1 && ktype==3)||(itype==3 && ktype==1))
                                    ik_InAl
                                elseif((itype==1 && ktype==4)||(itype==4 && ktype==1))
                                    ik_InSb
                                elseif((itype==2 && ktype==3)||(itype==3 && ktype==2))
                                    ik_AlAs
                                elseif((itype==2 && ktype==4)||(itype==4 && ktype==2))
                                    ik_AsSb
                                elseif((itype==3 && ktype==4)||(itype==4 && ktype==3))
                                    ik_AlSb                   
                                end
                
                                costheta=dot(vecij,vecik)/(rij*rik);
                
                                gtheta=gamma_ik*(1+c^2/d^2-c^2/(d^2+(costheta+h)^2));
                                incre=fc(rik,Rik,Dik)*gtheta*exp(lambda3^m*(rij-rik)^m);
                                xi=xi+incre;  
                            end
                        end
        
                        bij=(1+(beta2*xi)^n)^(-1/(2*n));
                        Eij(ifive1,ifive2)=Vij(rij,R,D,bij,De,S,beta,Re);
                    end
                end
                
                %partial i partial j
                Exmid=zeros(5,1);
                for ifive1=1:5
                    Exmid(ifive1)=(-Eij(ifive1,5)+8*Eij(ifive1,4)-8*Eij(ifive1,2)+Eij(ifive1,1))/(12*dh);%five point derivative on j
%                     Exmid(ifive1)=(Eij(ifive1,4)-Eij(ifive1,2))/(2*dh);%two point derivative on j
                end
                dEdij=(-Exmid(5)+8*Exmid(4)-8*Exmid(2)+Exmid(1))/(12*dh);%five point derivative on i
%                 dEdij=(Exmid(4)-Exmid(2))/(2*dh);%two point derivative on i
                Kblock(dim1,dim2)=dEdij;         
            end             
        end        
        KblockAll{i,j}=Kblock(1:3,1:3);

    end
end


%KblockAll generate asymmtric result, trouble shooting


%define k space grid
nk=101;
dk=2*pi/(nk-1);
kvecx=zeros(nk,1);

%phonon frequency
nb=3*na;%number of branches
omega=zeros(nk,nb);

% kvec=[0,pi/3,0];
%build the dynamical matrix using 3*3 matrix blocks
for ik=1
    %initialize dynamical matrix
    Kmatrix=zeros(3*na);
    kvecx(ik)=-pi+(ik-1)*dk;
    kvec=[kvecx(ik),0,0];
    for i=1:na
        for j=1:16
            jj=nncell(i,j);
            Kmatrix(3*i-2:3*i,3*jj-2:3*jj)=Kmatrix(3*i-2:3*i,3*jj-2:3*jj)+KblockAll{i,j}.*exp(1i*dot(kvec,nnlatt{i,j}));
            Kmatrix(3*i-2:3*i,3*i-2:3*i)=Kmatrix(3*i-2:3*i,3*i-2:3*i)-KblockAll{i,j};
        end
    end
    dyn=mass3dn*Kmatrix*mass3dn;%mass normalized
    eigenvalue=sort(abs(eig(dyn)));
    omega(ik,1:nb)=sqrt(eigenvalue(1:nb));
end     

%convert phonon unit to cm^-1
omega(:,:)=521.461464*omega(:,:);

plot(kvecx(1:nk),omega(1:nk,:),'blue'),xlabel('kx(-\pi/a1,\pi/a1)'),ylabel('\omega(cm^-1)')
set(gca,'FontSize',16)
axis([-pi pi 0 450])
title('Transverse direction')
toc

        
     
        
        
        
        
        
        
        
                
            


    








    
