% ќценка стандартной модели
% без авторегрессии в цикле безработицы
% c сезонностью
% наукаст при известной безработице за 1, 2, 3 мес€ца
% без ограничений на тренды выпуска и безработицы
%%
close all
clear all

global y
cd C:\Users\evsta\OneDrive\Desktop\работа\clark_okun

plus_year=2 % 0-1999q1, 1-2000q1
year_start=1999+plus_year;
un_month_known=3; % 2 - два мес€ца квартала известны, 1 - 1 мес€ц, 3 - 3
h_forecast_length=23; % сколько последних точек прогнозить

data=xlsread('GDP_Russia_wseas_1999q1_2020q3.xlsx');
% на подачу идет мес€чна€ безработица с тестового мес€ца
unempl = xlsread('U_Russia_wseas_1999m1_2020m9.xlsx');
forecast_arima=xlsread('arima_wseas_Russia_2013q4_2020q3.xlsx');


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


for sss=h_forecast_length:-1:1
%sss=1
y=[log(data(:,1)) data(:,2)/100]';    
tinit=1+4*plus_year; % 1 - 1999q1
tend=length(data)-sss; % 86 - 2020q2, 
tt=year_start+(tinit*0.25-0.25):0.25:year_start+(tend*0.25-0.25);

y1=y(:,tinit:tend+1);
y1=y1(2,1:end);

if un_month_known==1
uq_nc=means1;
uq_nc=uq_nc(tinit:tend+1);
y1(1,end)=uq_nc(end);
end
if un_month_known==2
uq_nc=means2;
uq_nc=uq_nc(tinit:tend+1);
y1(1,end)=uq_nc(end);
end

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
[betaest
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

t1=length(y);
t=length(y1);   


%z(t)=[tao(t) mu(t) us(t) sy(t) sy(t-1) sy(t-2) sy(t-3) su(t) su(t-1) su(t-2) su(t-3) c(t) c(t-1) uc(t)]'    
H=[1 0 0 1 0 0 0 0 0 0 0 1 0 0 
   0 0 1 0 0 0 0 1 0 0 0 0 0 1]';

H1=[0 0 1 0 0 0 0 1 0 0 0 0 0 1]';

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
statet_tm1(:,1)=[y(1,1) y(1,2)-y(1,1) y(2,1) 0 0 0 0 0 0 0 0 0 0 0]';
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
   if i<=t1
     yt_tm1(:,i)=H'*statet_tm1(:,i);
     statet_tm1(:,i+1)=F*statet_tm1(:,i)+F*Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
     Pt_tm1(:,:,i+1)=F*(Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*H'*Pt_tm1(:,:,i))*F'+Q;
     loglikle(i)=-0.5*log(2*pi)-0.5*log(det((H'*Pt_tm1(:,:,i)*H)))-0.5*(y(:,i)-H'*statet_tm1(:,i))'*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
   else
     yt_tm11(:,i)=H1'*statet_tm1(:,i);
     statet_tm1(:,i+1)=F*statet_tm1(:,i)+F*Pt_tm1(:,:,i)*H1*(H1'*Pt_tm1(:,:,i)*H1)^(-1)*(y1(:,i)-H1'*statet_tm1(:,i));
     Pt_tm1(:,:,i+1)=F*(Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H1*(H1'*Pt_tm1(:,:,i)*H1)^(-1)*H1'*Pt_tm1(:,:,i))*F'+Q;
     loglikle(i)=-0.5*log(2*pi)-0.5*log(det((H1'*Pt_tm1(:,:,i)*H1)))-0.5*(y1(:,i)-H1'*statet_tm1(:,i))'*(H1'*Pt_tm1(:,:,i)*H1)^(-1)*(y1(:,i)-H1'*statet_tm1(:,i));
end
end

statet_t=zeros(14,t);
Pt_t=zeros(14,14,t);
Jt=zeros(14,14,t);

statet_T=zeros(14,t);
Pt_T=zeros(14,14,t);

for i=1:t
    if i<=t1
      yt_tm1(:,i)=H'*statet_tm1(:,i);
      statet_t(:,i)=statet_tm1(:,i)+Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
      statet_tm1(:,i+1)=F*statet_tm1(:,i)+F*Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*(y(:,i)-H'*statet_tm1(:,i));
      Pt_t(:,:,i)=Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*H'*Pt_tm1(:,:,i);
      Pt_tm1(:,:,i+1)=F*(Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H*(H'*Pt_tm1(:,:,i)*H)^(-1)*H'*Pt_tm1(:,:,i))*F'+Q;
      Jt(:,:,i)=Pt_t(:,:,i)*F'*pinv(Pt_tm1(:,:,i+1));
    else
      yt_tm11(:,i)=H1'*statet_tm1(:,i);
      statet_t(:,i)=statet_tm1(:,i)+Pt_tm1(:,:,i)*H1*(H1'*Pt_tm1(:,:,i)*H1)^(-1)*(y1(:,i)-H1'*statet_tm1(:,i));
      statet_tm1(:,i+1)=F*statet_tm1(:,i)+F*Pt_tm1(:,:,i)*H1*(H1'*Pt_tm1(:,:,i)*H1)^(-1)*(y1(:,i)-H1'*statet_tm1(:,i));
      Pt_t(:,:,i)=Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H1*(H1'*Pt_tm1(:,:,i)*H1)^(-1)*H1'*Pt_tm1(:,:,i);
      Pt_tm1(:,:,i+1)=F*(Pt_tm1(:,:,i)-Pt_tm1(:,:,i)*H1*(H1'*Pt_tm1(:,:,i)*H1)^(-1)*H1'*Pt_tm1(:,:,i))*F'+Q;
      Jt(:,:,i)=Pt_t(:,:,i)*F'*pinv(Pt_tm1(:,:,i+1));
    end
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

forecast(sss)=statet_T(1,end)+statet_T(4,end)+statet_T(12,end);
end

%%
close all
tt=1999+(tinit*0.25-0.25):0.25:2020.25;
forecast_arima=forecast_arima(end-h_forecast_length+1:end);
tt1=tt((end-h_forecast_length+1:end));
y=[log(data(:,1)) data(:,2)/100]';  
test=y(1,end-h_forecast_length+1:end)'
forecast_uc=fliplr(forecast')'

for minus_last=0:2
e_uc=test(1:end-minus_last)-forecast_uc(1:end-minus_last);
e_arima=test(1:end-minus_last)-forecast_arima(1:end-minus_last);
rmse_uc=(1/(h_forecast_length-minus_last)*sum((e_uc).^2))^0.5;
rmse_arima=(1/(h_forecast_length-minus_last)*sum((e_arima).^2))^0.5;


e_uc1=e_uc./test(1:end-minus_last);
e_arima1=e_arima./test(1:end-minus_last);
mape_uc=100/(h_forecast_length-minus_last)*sum(abs(e_uc1));
mape_arima=100/(h_forecast_length-minus_last)*sum(abs(e_arima1));

[rmse_uc, mape_uc]
end

d1=(e_uc).^2-(e_arima).^2;
V1=var(d1);
d2=abs(e_uc)-abs(e_arima);
V2=var(d2);

dm_sqr=mean(d1)/sqrt(V1/length(d1));
dm_abs=mean(d2)/sqrt(V2/length(d2));


[dm_abs dm_sqr 2*tcdf(-abs(dm_abs),length(e_uc)-1) 2*tcdf(-abs(dm_sqr),length(e_uc)-1)]

%%

for minus_last=0:2
e_uc=test(1:end-minus_last)-forecast_uc(1:end-minus_last);
e_arima=test(1:end-minus_last)-forecast_arima(1:end-minus_last);
rmse_uc=(1/(h_forecast_length-minus_last)*sum((e_uc).^2))^0.5;
rmse_arima=(1/(h_forecast_length-minus_last)*sum((e_arima).^2))^0.5;


e_uc1=e_uc./test(1:end-minus_last);
e_arima1=e_arima./test(1:end-minus_last);
mape_uc=100/(h_forecast_length-minus_last)*sum(abs(e_uc1));
mape_arima=100/(h_forecast_length-minus_last)*sum(abs(e_arima1));

[rmse_arima, mape_arima]
end
%%

figure
plot(tt1,test,'Linewidth',2)
grid on
xlim([tt1(1) tt1(end)+0.5])
hold on
plot(tt1,forecast_uc,'Linewidth',2)
plot(tt1,forecast_arima,'Linewidth',2)
legend('Ћогарифм ¬¬ѕ','‘ильтра  алмана при известной безработице за 3 мес€ца','ARIMA','location','NW')
hold off
set(gca,'FontSize',14)
%%
