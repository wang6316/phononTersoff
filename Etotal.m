function [Etotal] = Etotal(xk,na,type)
%calculate total energy of the superlattice given the position vector "xk"
nnR=cell(na,16);
nntype=zeros(na,16);
nncount=zeros(na,1);

xlatt=xk(3*na+1);
ylatt=xk(3*na+2);
zlatt=xk(3*na+3);


%  make the nearest neighbour table
 for ia=1:na
     R=xk(3*ia-2:3*ia);
     inn=0;%count neighbors
     for ix=-1:1:1
         for inj=-1:1:1
             for iz=-1:1:1                 
                for ja=1:na
                    tempR=xk(3*ja-2:3*ja)+[xlatt*ix;ylatt*inj;zlatt*iz];
                    distance=norm(tempR-R);
                    if(distance>0.5 && distance<4.5)%loosely based on cutoff
                        inn=inn+1;
                        nnR{ia,inn}=tempR;%position vector of this neighbor
                        nntype(ia,inn)=type(ja);%atom type of this neighbor
                    end
                end
             end
         end
     end
     nncount(ia)=inn;
 end
 
 %now calculate the total energy
 Etotal=0;%initialize total energy

 parameters2021_v4%load the parameters
 
 for ia=1:na
     icoord=xk(3*ia-2:3*ia);
     itype=type(ia);
     for inj=1:nncount(ia)
         jcoord=nnR{ia,inj};
         jtype=nntype(ia,inj);
         
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
         
         for ink=1:nncount(ia)%loop for 16 possible neighbors k of atom i
             if(ink~=inj)
                 kcoord=nnR{ia,ink};
                 ktype=nntype(ia,ink);
                 
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
         Eij=Vij(rij,R,D,bij,De,S,beta1,Re);%evaluate potential energy between this pair of ij
         Etotal=Etotal+0.5*Eij;         
     end         
 end
                        


