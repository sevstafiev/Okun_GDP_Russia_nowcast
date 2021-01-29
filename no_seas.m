% ќценка стандартной модели
% без авторегрессии в цикле безработицы
% без сезонности
% наукаст при известной безработице за 1, 2, 3 мес€ца
% без ограничений на тренды выпуска и безработицы

close all
clear all

global y
cd C:\Users\evsta\OneDrive\Desktop\работа\clark_okun

plus_year=0 % 0-1999q1, 1-2000q1
year_start=2001+plus_year;
un_month_known=3; % 2 - два мес€ца квартала известны, 1 - 1 мес€ц, 3 - 3
h_forecast_length=28; % сколько последних точек прогнозить

minus_last_quart=1;

data=xlsread('GDP_Russia_noseas_2001q1_2020q3.xlsx');
% на подачу идет мес€чна€ безработица с тестового мес€ца
unempl = xlsread('U_Russia_noseas_2001m1_2020m9.xlsx');
forecast_arima=xlsread('arima_noseas_Russia_2013q4_2020q3.xlsx');

data=data(1:end-minus_last_quart,:);

% на подачу идет мес€чна€ безработица с тестового мес€ца

unempl = unempl(1:end-minus_last_quart*3,:);


forecast_arima=forecast_arima(1:end-minus_last_quart,:);


forecast=zeros(h_forecast_length,1);

unempl = unempl(:,1);
means1 = zeros(1,ceil(length(unempl)/3));
means2 = zeros(1,ceil(length(unempl)/3));
means3 = zeros(1,ceil(length(unempl)/3));

% remainder = length(unempl)/3 - round(length(unempl)/3);
% 
% if remainder ~= 0
%     if remainder > 0.5
%         unempl = [unempl' mean([unempl(length(unempl)-1),unempl(length(unempl))])];
%     else 
%         unempl = [unempl' unempl(length(unempl)) unempl(length(unempl))];
%     end
% end

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

%
y=[log(data(:,1)) data(:,2)/100]';    
tinit=1+4*plus_year; % 1 - 1999q1
tend=length(data); % 86 - 2020q2, 
tt=year_start+(tinit*0.25-0.25):0.25:year_start+(tend*0.25-0.25);

y=y(:,tinit:tend);
t=length(y);

trial=20;
options = optimset('Display','iter',...
    'TolFun',1e-6,... 
    'HessUpdate', 'bfgs',...
            'MaxFunEvals',3000,...
    'TolX',1e-6);

% beta_trial=zeros(trial,8);
% betaest_trial=zeros(trial,8);
% FVAL_trial=zeros(trial);
% 
% 
% for i=1:trial
%     err=1;
%     while err>0
%     %beta_trial(i,1)=2*abs(randn(1));% beta
%     beta_trial(i,1)=-rand(1);
%     beta_trial(i,2)=2*abs(randn(1));%ar(1)
%     beta_trial(i,3)=2*randn(1);%ar(2)
%     beta_trial(i,4)=abs(randn(1));
%     beta_trial(i,5)=abs(randn(1));
%     beta_trial(i,6)=abs(randn(1));
%     beta_trial(i,7)=abs(randn(1));
%     beta_trial(i,8)=abs(randn(1));
%     
%    
% 
%     phi1=beta_trial(i,1);
%     phi2=beta_trial(i,2);
%     ARmatrix=[phi1 phi2;
%               1    0];
%        e=eig(ARmatrix);
%     if max(abs(e))>=1
%          err=1;
%     else
%         err=0;
%     end
%     end
%         
% end
% 
% 
% 
% 
% for i=1:trial
%     i
% [betaest_trial(i,:),FVAL_trial(i),EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(@loglike_clark_okun2,beta_trial(i,:),options);
% end
% [FVAL_trial,I]=sort(FVAL_trial);
% I(1)
% betaest_trial(I(1),:)
% 
% [betaest,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(@loglike_clark_okun2,betaest_trial(I(1),:),options);



params_start=[-0.2 1.5 -0.8 0.1 0.1 0.1  0.1  0.1];

[betaest,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(@loglike_clark_okun2,params_start,options);


std=((diag(pinv(1/t*HESSIAN))/t).^0.5)';
coefs=[betaest
    abs(std)
betaest./abs(std)]


params=betaest
beta=params(1);
phi1=params(2);
phi2=params(3);
sigma2nu=params(4)^2;
sigma2v=params(5)^2;
sigma2eps=params(6)^2;
sigma2ksi=params(7)^2;
sigma2del=params(8)^2;

t=length(y); 


%z(t)=[tao(t) mu(t) us(t) c(t) c(t-1) uc(t)]'    
H=[1 0 0 1 0 0
   0 0 1 0 0 1]';


F=[ 1 1 0 0 0 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 0 0 phi1 phi2 0
    0 0 0 1 0 0
    0 0 0 beta*phi1 beta*phi2 0];

Q=zeros(6,6);
Q(1,1)=sigma2nu+sigma2v;
Q(1,2)=sigma2v;
Q(2,1)=sigma2v;
Q(2,2)=sigma2v;
Q(3,3)=sigma2eps;
Q(4,4)=sigma2ksi;
Q(4,6)=beta*sigma2ksi;
Q(6,4)=beta*sigma2ksi;
Q(6,6)=beta^2*sigma2ksi+sigma2del;

statet_tm1=zeros(6,t+1);
Pt_tm1=zeros(6,6,t+1);
yt_tm1=zeros(2,t);
loglikle=zeros(1,t);

%%%initialisation
%z(t)=[tao(t) mu(t) us(t) c(t) c(t-1) uc(t)]'
statet_tm1(:,1)=[y(1,1) y(1,2)-y(1,1) y(2,1) 0 0 0]';
vecP10=(eye(size(kron(F(4:6,4:6),F(4:6,4:6))))-kron(F(4:6,4:6),F(4:6,4:6)))^(-1)*vec(Q(4:6,4:6));
P10=reshape(vecP10,3,3);
Pt_tm1(:,:,1)=[1000 0 0 0 0 0
               0 1000 0 0 0 0
               0 0 1000 0 0 0
               0  0  0 P10(1,1) P10(1,2) P10(1,3)
               0  0  0 P10(2,1) P10(2,2) P10(2,3)
               0  0  0 P10(3,1) P10(3,2) P10(3,3)];
           
for i=1:t
     yt_tm1(:,i)=H'*statet_tm1(:,i);
     statet_tm1(:,i+1)=F*statet_tm1(:,i)+F*Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
     Pt_tm1(:,:,i+1)=F*(Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*H'*Pt_tm1(:,:,i))*F'+Q;
     loglikle(i)=-0.5*log(2*pi)-0.5*log(det((H'*Pt_tm1(:,:,i)*H)))-0.5*(y(:,i)-H'*statet_tm1(:,i))'*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
end

statet_t=zeros(6,t);
Pt_t=zeros(6,6,t);
Jt=zeros(6,6,t);

statet_T=zeros(6,t);
Pt_T=zeros(6,6,t);

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

stder=zeros(6,t);
for j=1:6
stder(j,:)=Pt_T(j,j,:).^0.5;
end

figure
plot(tt(1:end),y(1,1:end),'linewidth',2)
grid on
hold on
plot(tt(1:end),statet_T(1,1:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980])
plot(tt(1:end),statet_T(1,1:end)-stder(1,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
plot(tt(1:end),statet_T(1,1:end)+stder(1,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
hold off
legend('Ћогарифм ¬¬ѕ','“рендова€ компонента ¬¬ѕ','68% доверительный интервал','Location','SE')
set(gca,'FontSize',14)
%%
figure
plot(tt(1:end),400*statet_T(2,(1:end)),'linewidth',2)
grid on
hold on
plot(tt(1:end),400*(statet_T(2,(1:end))-stder(2,(1:end))),'linewidth',1,'Color',[0, 0.4470, 0.7410],'Linestyle','--')
plot(tt(1:end),400*(statet_T(2,(1:end))+stder(2,(1:end))),'linewidth',1,'Color',[0, 0.4470, 0.7410],'Linestyle','--')
hold off
set(gca,'FontSize',14)
title('“емп роста трендовой компоненты (% в год)')

%%
%z(t)=[tao(t) mu(t) us(t) c(t) c(t-1) uc(t)]'
figure
plot(tt(1:end),y(2,1:end),'linewidth',2)
grid on
hold on
plot(tt(1:end),statet_T(3,1:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980])
plot(tt(1:end),statet_T(3,1:end)-stder(3,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
plot(tt(1:end),statet_T(3,1:end)+stder(3,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
hold off
legend('Ѕезработица','“рендова€ компонента безработицы','68% доверительный интервал','Location','NE')
set(gca,'FontSize',14)

%%
%z(t)=[tao(t) mu(t) us(t) c(t) c(t-1) uc(t)]'
figure
subplot(2,1,1)
plot(tt(1:end),statet_T(4,1:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980])
hold on
grid on
plot(tt(1:end),statet_T(4,1:end)-stder(4,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
plot(tt(1:end),statet_T(4,1:end)+stder(4,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
hold off
legend('÷иклическа€ компонента выпуска','68% доверительный интервал','Location','SE')
set(gca,'FontSize',14)

subplot(2,1,2)
plot(tt(1:end),statet_T(6,1:end),'linewidth',2,'Color',[0.8500, 0.3250, 0.0980])
hold on
grid on
plot(tt(1:end),statet_T(6,1:end)-stder(6,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
plot(tt(1:end),statet_T(6,1:end)+stder(6,1:end),'linewidth',1,'Color',[0.8500, 0.3250, 0.0980],'Linestyle','--')
hold off
legend('÷иклическа€ компонента безработицы','68% доверительный интервал','Location','SE')
set(gca,'FontSize',14)

%%

std=((diag((1/t*HESSIAN)^-1)/t).^0.5)';

t = betaest(1)/std(1)
B = beta;

if t >= -1.64
    x="коэф-т незначим"
elseif t < -2.58
    x=["***",round(std(1),4),round(B,4)]
elseif t < -1.96
    x=["**",round(std(1),4),round(B,4)]
elseif t < -1.64
    x=["*",round(std(1),4),round(B,4)]
end