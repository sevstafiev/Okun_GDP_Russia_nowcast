function logf=loglike_clark_okun6(params)
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

%z(t)=[tao(t) mu(t) us(t) sy(t) sy(t-1) sy(t-2) sy(t-3) su(t) su(t-1) su(t-2) su(t-3) c(t) c(t-1) uc(t) nu(t) eps(t) wsy(t) wsu(t)]'    
H=[1 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0
   0 0 1 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0]';

F=[ 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 -1 -1 -1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 phi1 phi2 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 beta*phi1 beta*phi2 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

Q=zeros(18,18);
Q(1,1)=sigma2nu+sigma2v;
Q(1,2)=sigma2v;
Q(2,1)=sigma2v;
Q(1,15)=sigma2nu;
Q(15,1)=sigma2nu;
Q(15,15)=sigma2nu;
Q(2,2)=sigma2v;
Q(3,3)=sigma2eps;
Q(3,16)=sigma2eps;
Q(16,3)=sigma2eps;
Q(16,16)=sigma2eps;
Q(4,4)=sigma2wsy;
Q(4,17)=sigma2wsy;
Q(17,4)=sigma2wsy;
Q(17,17)=sigma2wsy;
Q(8,8)=sigma2wsu;
Q(8,18)=sigma2wsu;
Q(18,8)=sigma2wsu;
Q(18,18)=sigma2wsu;
Q(12,12)=sigma2ksi;
Q(12,14)=beta*sigma2ksi;
Q(14,12)=beta*sigma2ksi;
Q(14,14)=beta^2*sigma2ksi+sigma2del;


statet_tm1=zeros(18,t+1);
Pt_tm1=zeros(18,18,t+1);
yt_tm1=zeros(2,t);
loglikle=zeros(1,t);

%%%initialisation
%z(t)=[tao(t) mu(t) us(t) sy(t) sy(t-1) sy(t-2) sy(t-3) su(t) su(t-1) su(t-2) su(t-3) c(t) c(t-1) uc(t) nu(t) eps(t) wsy(t) wsu(t)]'    
statet_tm1(:,1)=[y(1,1) y(1,2)-y(1,1) y(2,1) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
vecP10=(eye(size(kron(F(12:18,12:18),F(12:18,12:18))))-kron(F(12:18,12:18),F(12:18,12:18)))^(-1)*vec(Q(12:18,12:18));
P10=reshape(vecP10,7,7);
Pt_tm1(:,:,1)=[1000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
               0 1000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
               0 0 1000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
               0 0 0 1000 0 0 0 0 0 0 0 0 0 0 0 0 0 0
               0 0 0 0 1000 0 0 0 0 0 0 0 0 0 0 0 0 0
               0 0 0 0 0 1000 0 0 0 0 0 0 0 0 0 0 0 0
               0 0 0 0 0 0 1000 0 0 0 0 0 0 0 0 0 0 0
               0 0 0 0 0 0 0 1000 0 0 0 0 0 0 0 0 0 0
               0 0 0 0 0 0 0 0 1000 0 0 0 0 0 0 0 0 0
               0 0 0 0 0 0 0 0 0 1000 0 0 0 0 0 0 0 0
               0 0 0 0 0 0 0 0 0 0 1000 0 0 0 0 0 0 0
               0 0 0 0 0 0 0 0 0  0  0 P10(1,1) P10(1,2) P10(1,3) P10(1,4) P10(1,5) P10(1,6) P10(1,7)
               0 0 0 0 0 0 0 0 0  0  0 P10(2,1) P10(2,2) P10(2,3) P10(2,4) P10(2,5) P10(2,6) P10(2,7)
               0 0 0 0 0 0 0 0 0  0  0 P10(3,1) P10(3,2) P10(3,3) P10(3,4) P10(3,5) P10(3,6) P10(3,7)
               0 0 0 0 0 0 0 0 0  0  0 P10(4,1) P10(4,2) P10(4,3) P10(4,4) P10(4,5) P10(4,6) P10(4,7)
               0 0 0 0 0 0 0 0 0  0  0 P10(5,1) P10(5,2) P10(5,3) P10(5,4) P10(5,5) P10(5,6) P10(5,7)
               0 0 0 0 0 0 0 0 0  0  0 P10(6,1) P10(6,2) P10(6,3) P10(6,4) P10(6,5) P10(6,6) P10(6,7)
               0 0 0 0 0 0 0 0 0  0  0 P10(7,1) P10(7,2) P10(7,3) P10(7,4) P10(7,5) P10(7,6) P10(7,7)];
           
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
