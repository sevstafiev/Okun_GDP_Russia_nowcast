function logf=loglike_clark_okun4(params)
global y 


beta=params(1);
phi1=params(2);
phi2=params(3);
sigma2nu=params(4)^2;
sigma2v=params(5)^2;
sigma2eps=params(6)^2;
sigma2wsy=params(7)^2;
sigma2wsu=params(8)^2;
sigma2ksi=params(9)^2;
sigma2del=params(10)^2;



ARmatrix=[phi1 phi2;
           1    0];
       e=eig(ARmatrix);
if max(abs(e))>=1 
logf=10^7;
else

t=length(y);   

%z(t)=[tao(t) mu(t) us(t) sy(t) sy(t-1) sy(t-2) sy(t-3) su(t) su(t-1) su(t-2) su(t-3) c(t) c(t-1) uc(t)]'    
H=[1 0 0 1 0 0 0 0 0 0 0 1 0 0 
   0 0 1 0 0 0 0 1 0 0 0 0 0 1]';

F=[ 1 1 0  0  0  0  0 0 0 0 0 0 0 0
    0 1 0  0  0  0  0 0 0 0 0 0 0 0 
    0 0 1  0  0  0  0 0 0 0 0 0 0 0 
    0 0 0 -1 -1 -1  0 0 0 0 0 0 0 0 
    0 0 0  1  0  0  0 0 0 0 0 0 0 0 
    0 0 0  0  1  0  0 0 0 0 0 0 0 0 
    0 0 0  0  0  1  0 0 0 0 0 0 0 0 
    0 0 0  0  0  0  0 -1 -1 -1 0 0 0 0 
    0 0 0  0  0  0  0  1  0  0 0 0 0 0
    0 0 0  0  0  0  0 0   1  0 0 0 0 0
    0 0 0  0  0  0  0 0   0  1 0 0 0 0
    0 0 0  0  0  0  0 0   0  0  0 phi1 phi2 0 
    0 0 0  0  0  0  0 0   0  0  0 1 0 0 
    0 0 0  0  0  0  0 0   0  0  0 beta*phi1 beta*phi2 0];

Q=zeros(14,14);
Q(1,1)=sigma2nu+sigma2v;
Q(1,2)=sigma2v;
Q(2,1)=sigma2v;
Q(2,2)=sigma2v;
Q(3,3)=sigma2eps;
Q(4,4)=sigma2wsy;
Q(8,8)=sigma2wsu;
Q(12,12)=sigma2ksi;
Q(12,14)=beta*sigma2ksi;
Q(14,12)=beta*sigma2ksi;
Q(14,14)=beta^2*sigma2ksi+sigma2del;

statet_tm1=zeros(14,t+1);
Pt_tm1=zeros(14,14,t+1);
yt_tm1=zeros(2,t);
loglikle=zeros(1,t);

%%%initialisation
%z(t)=[tao(t) mu(t) us(t) sy(t) sy(t-1) sy(t-2) sy(t-3) su(t) su(t-1) su(t-2) su(t-3) c(t) c(t-1) uc(t)]'    
statet_tm1(:,1)=[y(1,1) y(1,2)-y(1,1) y(2,1)-0.0054 0 0 0 0 0.0054   -0.0009   -0.0049    0.0004 0 0 0]';
vecP10=(eye(size(kron(F(12:14,12:14),F(12:14,12:14))))-kron(F(12:14,12:14),F(12:14,12:14)))^(-1)*vec(Q(12:14,12:14));
P10=reshape(vecP10,3,3);
Pt_tm1(:,:,1)=[1000 0 0 0 0 0 0 0 0 0 0 0 0 0
               0 1000 0 0 0 0 0 0 0 0 0 0 0 0
               0 0 1000 0 0 0 0 0 0 0 0 0 0 0
               0 0 0 1000 0 0 0 0 0 0 0 0 0 0
               0 0 0 0 1000 0 0 0 0 0 0 0 0 0
               0 0 0 0 0 1000 0 0 0 0 0 0 0 0
               0 0 0 0 0 0 1000 0 0 0 0 0 0 0
               0 0 0 0 0 0 0 1000 0 0 0 0 0 0
               0 0 0 0 0 0 0 0 1000 0 0 0 0 0
               0 0 0 0 0 0 0 0 0 1000 0 0 0 0
               0 0 0 0 0 0 0 0 0 0 1000 0 0 0
               0 0 0 0 0 0 0 0 0  0  0 P10(1,1) P10(1,2) P10(1,3)
               0 0 0 0 0 0 0 0 0  0  0 P10(2,1) P10(2,2) P10(2,3)
               0 0 0 0 0 0 0 0 0  0  0 P10(3,1) P10(3,2) P10(3,3)];
           
for i=1:t
yt_tm1(:,i)=H'*statet_tm1(:,i);
statet_tm1(:,i+1)=F*statet_tm1(:,i)+F*Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
Pt_tm1(:,:,i+1)=F*(Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*H'*Pt_tm1(:,:,i))*F'+Q;
loglikle(i)=-0.5*log(2*pi)-0.5*log(det((H'*Pt_tm1(:,:,i)*H)))-0.5*(y(:,i)-H'*statet_tm1(:,i))'*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
end
loglikleGlobe=sum(loglikle);
logf=-loglikleGlobe;
if imag(logf)~=0
logf=10000000000000;
end
end    
end
