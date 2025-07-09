function [u, b, C]= VGJD(u0,Img, b, Ksigma, timestep,epsilon, iter,alfa,beta)

u=u0;
KB1 = conv2(b,Ksigma,'same');
KB2 = conv2(b.^2,Ksigma,'same');
C =updateC(Img, u, KB1, KB2, epsilon);

u = updateLSF(Img,u, C, b, timestep, epsilon, iter,alfa,beta);
Hu=Heaviside(u,epsilon);
M(:,:,1)=Hu;
M(:,:,2)=1-Hu;
b =updateB(Img, C, M,  Ksigma);


% update level set function
function u = updateLSF(Img, u0, C, b,  timestep, epsilon, iter,alfa,beta)
u=u0;
Hu=Heaviside(u,epsilon);
M(:,:,1)=Hu;
M(:,:,2)=1-Hu;
N_class=size(M,3);
e=zeros(size(M));
u=u0;
for kk=1:N_class
    e(:,:,kk) = (Img-b.*C(kk)).*log(Img./(b.*C(kk))+1e-10);
    %e(:,:,kk) = (Img).*log(Img./(b.*C(kk))+1e-10).^2;
end

for kk=1:iter
    u=NeumannBoundCond(u);
    DiracU=Dirac(u,epsilon);
    e1 = (e(:,:,1)-e(:,:,2));
    dataForce = -alfa*(e1./beta)./(sqrt(1+(e1./beta).^2));
    ImageTerm = DiracU.*dataForce;
    u=u+timestep*(ImageTerm);
end

% update b
function b =updateB(Img, C, M,  Ksigma)

PC1=zeros(size(Img));
PC2=PC1;
N_class=size(M,3);
for kk=1:N_class
    PC1=PC1+C(kk)*M(:,:,kk);
    PC2=PC2+C(kk)^2*M(:,:,kk);
end
KNm1 = conv2(PC1.*Img,Ksigma,'same');
KDn1 = conv2(PC2,Ksigma,'same');

b = KNm1./KDn1+1e-10;

% Update C
function C_new =updateC(Img, u, Kb1, Kb2, epsilon)
Hu=Heaviside(u,epsilon);
M(:,:,1)=Hu;
M(:,:,2)=1-Hu;
N_class=size(M,3);
for kk=1:N_class
    Nm2 = Kb1.*Img.*M(:,:,kk);
    Dn2 = Kb2.*M(:,:,kk);
    C_new(kk) = sum(Nm2(:))/sum(Dn2(:))+1e-10;
end



% Make a function satisfy Neumann boundary condition
function g = NeumannBoundCond(f)
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  


function h = Heaviside(x,epsilon)    
h=0.5*(1+(2/pi)*atan(x./epsilon));

function f = Dirac(x, epsilon)    
f=(epsilon/pi)./(epsilon^2.+x.^2);

