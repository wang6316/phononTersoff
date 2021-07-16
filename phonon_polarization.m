%following phonon_dynmat
%check the polarization vectors 

na=16;

%define k space grid 
nk=101;
dk=pi/(nk-1);
kvecy=zeros(nk,1);

%phonon frequency
nb=3*na;
omega=zeros(nk,nb);
polvec=cell(nk,nb);%define mass normalized polarization vectors

%build the dynamical matrix using 3*3 matrix blocks
for ik=50
    %initialize dynamical matrix
    Kmatrix=zeros(3*na);
    kvecy(ik)=(ik-(nk+1)/2)*dk;
    kvec=[0,kvecy(ik),0];
    for i=1:na
        for j=1:16
            jj=nncell(i,j);
            Kmatrix(3*i-2:3*i,3*jj-2:3*jj)=Kmatrix(3*i-2:3*i,3*jj-2:3*jj)+KblockAll{i,j}.*exp(1i*dot(kvec,nnlatt{i,j}));
            Kmatrix(3*i-2:3*i,3*i-2:3*i)=Kmatrix(3*i-2:3*i,3*i-2:3*i)-KblockAll{i,j};
        end
    end
    dyn=mass3dn*Kmatrix*mass3dn;%mass normalized
    [eigenvector,eigenvalue]=eig(dyn);%no need to sort
    for ib=1:nb
        omega(ik,ib)=sqrt(eigenvalue(ib,ib));
        polvec{ik,ib}=eigenvector(:,ib);
    end
    
end     
