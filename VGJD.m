function [u, b, C]= VGJD(u0,Img, b, Ksigma, timestep,epsilon, iter,alfa,beta)

u=u0;

%% Update constants c inside and outside the contour
convB1 = conv2(b,Ksigma,'same');
convB2 = conv2(b.^2,Ksigma,'same');
C =updateC(Img, u, convB1, convB2, epsilon);

%% Update level set function
u = updateLSF(Img,u, C, b, timestep, epsilon, iter,alfa,beta);
Hu=Heaviside(u,epsilon);
M(:,:,1)=Hu;
M(:,:,2)=1-Hu;

%% Update bias field b
b =updateB(Img, C, M,  Ksigma);

%% Function: Update level set
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
end
for kk=1:iter
    u=NeumannBoundCond(u);
    diracFunc=Dirac(u,epsilon);
    e1 = (e(:,:,1)-e(:,:,2));
    dataDriven = -alfa*(e1./beta)./(sqrt(1+(e1./beta).^2)); 
    imageForce = diracFunc.*dataDriven; % Eq.13
    u=u+timestep*(imageForce); % Eq.14
end

%% Function: Update bias field b
function b =updateB(Img, C, M,  Ksigma)
sumMu=zeros(size(Img));
sumMu2=sumMu;
N_class=size(M,3);
for kk=1:N_class
    sumMu=sumMu+C(kk)*M(:,:,kk);
    sumMu2=sumMu2+C(kk)^2*M(:,:,kk);
end
convMuImg = conv2(sumMu.*Img,Ksigma,'same');
convMu2 = conv2(sumMu2,Ksigma,'same');
b = convMuImg./convMu2+1e-10;

%% Function: Update constants c inside and outside the contour
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

%% Applies Neumann boundary conditions to maintain continuity at matrix edges
function g = NeumannBoundCond(f)
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  

%% Smooth approximation of Heaviside step function
function h = Heaviside(x,epsilon)    
h=0.5*(1+(2/pi)*atan(x./epsilon));

%% Smooth approximation of Dirac delta function
function f = Dirac(x, epsilon)    
f=(epsilon/pi)./(epsilon^2.+x.^2);

