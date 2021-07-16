%plot potential energy surface of 1 atom in the relaxed lattice
tic

parameters2021_v4

%five-point stencil 
dh=0.01*a;
PE=zeros(15);%five point value for derivative

Eij=zeros(1,16);

     %update neighboring site information for this iteration
    nncount=zeros(na,1);
    for i=1:na%loop i runs over all atoms of the supercell 
        inn=0;
        %searching neighbors in 3*3 lattice with the supercell
        for ix=-1:1:1
            for iy=-1:1:1
                for iz=-1:1:1%lattice index of x,y,z               
                    for j=1:na%supercell index for neighbor                   
                        tempR=RR{j}+[xlatt*ix;ylatt*iy;zlatt*iz];
                        distance=norm(tempR-RR{i});  %scaled on Angstrom                
                        if (distance>0.5 && distance<4.5)%loosely based on cutoff of Tersoff potential
                            inn=inn+1;
                            nnR{i,inn}=tempR;%position vector of this neighbor
                            nncell(i,inn)=j;%unit cell index of neighbor j
                            nnlatt{i,inn}=[xlatt*ix,ylatt*iy,zlatt*iz];%lattice coordinates of this neighbor
                            nntype(i,inn)=type(j);%atom type of this neighbor
                        end
                    end
                end
            end
        end
        nncount(i)=inn;
    end

% with the information above, finding force constant using Tersoff potential
%choose xy grid 

for i=5%loop for possible atom i, runs over a supercell
    icoord0=RR{i};%asign realistic value, length in angstrom
    itype=type(i);
%     sumforce=zeros(3,1);
%     singleforce=zeros(3,1);
    netforce=zeros(3,1);
    
        for ix=-7:1:7%five points for alpha
            for iy=-7:1:7
                icoord=icoord0+[dh;0;0]*ix+[0;dh;0]*iy;                
            for j=1:nncount(i) %loops for 16 possible neighbors j of atom i
                jcoord=nnR{i,j};
                jtype=nntype(i,j);
                %above choose i,j in supercell               
                vecij=jcoord-icoord;
                rij=norm(vecij);
                    
                xi=0;%initialize xi for all 3-body correction
        
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
               
                        for k=1:nncount(i)%loop for 16 or fewer possible neighbors k of atom i
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
                        Eij(j)=Vij(rij,R,D,bij,De,S,beta1,Re);%potential energy atom i feels contributed by atom j
%                         ;%total potential energy of atom i at current displaced position
            end
                PE(ix+8,iy+8)=sum(Eij);   
            end
        end
end

surf(PE)
xlabel("perpendicular to interface"),ylabel("parallel to interface"),zlabel("PES")

toc