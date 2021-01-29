% Оценка стандартной модели
% без авторегрессии в цикле безработицы
% c сезонностью
% наукаст при известной безработице за 1, 2, 3 месяца
% без ограничений на тренды выпуска и безработицы

close all
clear all

global y
cd C:\Users\evsta\OneDrive\Desktop\работа\clark_okun

plus_year=2 % 0-1999q1, 1-2000q1
year_start=1999+plus_year;
un_month_known=3; % 2 - два месяца квартала известны, 1 - 1 месяц, 3 - 3
h_forecast_length=28; % сколько последних точек прогнозить

minus_last_quart=1;

data=xlsread('GDP_Russia_wseas_1999q1_2020q3.xlsx');
% на подачу идет месячная безработица с тестового месяца
unempl = xlsread('U_Russia_wseas_1999m1_2020m9.xlsx');
forecast_arima=xlsread('arima_wseas_Russia_2013q4_2020q3.xlsx');

data=data(1:end-minus_last_quart,:);

% на подачу идет месячная безработица с тестового месяца

unempl = unempl(1:end-minus_last_quart*3,:);


forecast_arima=forecast_arima(1:end-minus_last_quart,:);


forecast=zeros(h_forecast_length,1);

unempl = unempl(:,1);
means1 = zeros(1,ceil(length(unempl)/3));
means2 = zeros(1,ceil(length(unempl)/3));
means3 = zeros(1,ceil(length(unempl)/3));


for i = 1:3:length(unempl)
    means1(ceil(i/3)) = mean(unempl(i:(i)));
end

for i = 1:3:length(unempl)
    means2(ceil(i/3)) = mean(unempl(i:(i+1)));
end

for i = 1:3:length(unempl)
    means3(ceil(i/3)) = mean(unempl(i:(i+2)));
end

means1 = means1/100;  %средние по безработице за каждый квартал
means2 = means2/100;
means3 = means3/100;

data(:,2) = means3*100';

y=[log(data(:,1)) data(:,2)/100]';    
tinit=1+4*plus_year; % 1 - 1999q1
tend=length(data) % 86 - 2020q2, 
tt=year_start+(tinit*0.25-0.25):0.25:year_start+(tend*0.25-0.25);

y=y(:,tinit:tend);
t=length(y);


options = optimset('Display','iter',...
    'TolFun',1e-6,... 
    'HessUpdate', 'bfgs',...
            'MaxFunEvals',3000,...
    'TolX',1e-6);


params_start=[-0.2 1.5 -0.8 0.1 0.1 0.1 0.1 0.1 0.1 0.1];

[betaest,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(@loglike_clark_okun4,params_start,options);
betaest;
std=((diag(pinv(1/t*HESSIAN))/t).^0.5)';
coefs=[betaest
    abs(std)
betaest./abs(std)];


params=betaest
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
yt_tm11=zeros(1,t);
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


statet_t=zeros(14,t);
Pt_t=zeros(14,14,t);
Jt=zeros(14,14,t);

statet_T=zeros(14,t);
Pt_T=zeros(14,14,t);

for i=1:t
      yt_tm1(:,i)=H'*statet_tm1(:,i);
      statet_t(:,i)=statet_tm1(:,i)+Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
      statet_tm1(:,i+1)=F*statet_tm1(:,i)+F*Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
      Pt_t(:,:,i)=Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*H'*Pt_tm1(:,:,i);
      Pt_tm1(:,:,i+1)=F*(Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*H'*Pt_tm1(:,:,i))*F'+Q;
      Jt(:,:,i)=Pt_t(:,:,i)*F'*pinv(Pt_tm1(:,:,i+1));
end

statet_T(:,t)=statet_t(:,t);
Pt_T(:,:,t)=Pt_t(:,:,t);

for i=t-1:-1:1
statet_T(:,i)=statet_t(:,i)+Jt(:,:,i)*(statet_T(:,i+1)-statet_tm1(:,i+1)); 
Pt_T(:,:,i)=Pt_t(:,:,i)+Jt(:,:,i)*(Pt_T(:,:,i+1)-Pt_tm1(:,:,i+1))*Jt(:,:,i)';
end

stder=zeros(14,t);
for j=1:14
stder(j,:)=Pt_T(j,j,:).^0.5;
end

%%
%z(t)=[tao(t) mu(t) us(t) sy(t) sy(t-1) sy(t-2) sy(t-3) su(t) su(t-1) su(t-2) su(t-3) c(t) c(t-1) uc(t)]'  
figure
subplot(2,1,1)
plot(tt(1:end),y(1,1:end),'linewidth',2)
grid on
hold on
plot(tt(1:end),y(1,1:end)-statet_T(4,1:end),'linewidth',2,'Color',[0.9290, 0.6940, 0.1250])
plot(tt(1:end),statet_T(1,1:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980])
plot(tt(1:end),statet_T(1,1:end)-stder(1,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
plot(tt(1:end),statet_T(1,1:end)+stder(1,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
hold off
legend('Логарифм ВВП','Логарифм ВВП без сезонности','Трендовая компонента ВВП','68% доверительный интервал','Location','SE')
set(gca,'FontSize',14)
title('Логарифм ВВП')

subplot(2,1,2)
plot(tt(1:end),statet_T(4,1:end),'linewidth',2)
grid on
title('Сезонность ВВП')
set(gca,'FontSize',14)
%%
figure
plot(tt(5:end),400*statet_T(2,(5:end)),'linewidth',2)
grid on
hold on
plot(tt(5:end),400*(statet_T(2,(5:end))-stder(2,(5:end))),'linewidth',1,'Color',[0, 0.4470, 0.7410],'Linestyle','--')
plot(tt(5:end),400*(statet_T(2,(5:end))+stder(2,(5:end))),'linewidth',1,'Color',[0, 0.4470, 0.7410],'Linestyle','--')
hold off
set(gca,'FontSize',14)
title('Темп роста трендовой компоненты (% в год)')

%%
figure
plot(tt(1:end),y(2,1:end),'linewidth',2)
grid on
hold on
plot(tt(1:end),y(2,1:end)-statet_T(8,1:end),'linewidth',2,'Color',[0.9290, 0.6940, 0.1250])
plot(tt(1:end),statet_T(3,1:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980])
legend('Безработица','Безработица без сезонности','Трендовая компонента безработицы','Location','NE')
%%
figure
subplot(2,1,1)
plot(statet_T(8,:),'linewidth',2)
grid on

figure
plot(y(2,:))
%%
%z(t)=[tao(t) mu(t) us(t) sy(t) sy(t-1) sy(t-2) sy(t-3) su(t) su(t-1) su(t-2) su(t-3) c(t) c(t-1) uc(t)]'  
figure
plot(y(2,:))
hold on
plot(statet_T(3,:)+statet_T(8,:))
hold off
%%
plot(statet_T(9,:))