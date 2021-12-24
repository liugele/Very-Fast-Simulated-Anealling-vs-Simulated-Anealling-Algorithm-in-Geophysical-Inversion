%Program pemodelan inversi Data Transient Electromagnetic (TEM)/TDEM (Time-domain Electromagnetic) Central Loop Configuration
%menggunakan algoritma Very Fast Simulated Annealing (VFSA)
%Mohammad Rheza Zamani
%Reference : Grandis,H.(2009): Pengantar pemodelan inversi geofisika, Jakarta:HAGI

%Reference : W. Srigutomo, M. Heriyanto, and M. Hilmi Aufa. Gravity Inversion of Talwani Model using Very Fast Simulated Annealing. Journal of Mathematical and Fundamental Sciences, Vol. 51, No. 2, 2019, 177-190. doi: 10.5614/j.math.fund.sci.2019.51.2.7.
clear all;
clc;
format long
%Input data sintetik
t1=linspace(log10(10^-6),log(1),30);
t=10.^t1;
a=25;
I=1;
R = [100 1000 10];
thk = [10 100];
%Pehitungan data sintetik
TEM_sin =fwd_TEM(R,thk,t,a,I);
%Data Latihan
%t = [ 0.00004007 0.0001011 0.0001621 0.0002378 0.0003601 0.000527 0.0008145 0.001558];
%TEM_sin = [ 0.0000054811749 0.00000098437 0.000000396178 0.00000018854 0.000000083357 0.000000038893 0.0000000156912 0.0000000034128];
%Definisi ruang model 
nlayer = 3; %Jumlah lapisan 
nitr = 300; %Jumlah iterasi 
T = 1;
dec = 1;
%Batas atas dan bawah (Diatur 5 kali dari nilai parameter model)
%Batas bawah pencarian nilai resistivitas
LBR = [1 1 1];
%Batas atas pencarian nilai resistivitas
UBR = [200 2000 20];
%Batas bawah pencarian nilai ketebalan
LBT = [1 1];
%Batas atas pencarian nilai resistivitas
UBT = [20 200];
%Membuat model awal acak

rho1(1 , :) = LBR + rand*(UBR - LBR);
thick1(1, :) = LBT + rand*(UBT - LBT);

%Menghitung  dB/dt model semu dan misfit dari model
[TEM1] = fwd_TEM(rho1(1,:),thick1(1,:),t,a,I);
dBdt1(1,:) = TEM1;
[misfit1] = misfit_TEM(TEM_sin,dBdt1(1,:));
E1 = misfit1;

v = VideoWriter('TEM VFSA.avi');
open(v);

for itr = 1 : nitr
    rho_int(1 , :) = LBR + rand*(UBR - LBR)
    thick_int(1, :) = LBT + rand*(UBT - LBT)
    ui = rand;
    yi = sign(ui-0.5)*T*((((1 + (1/T)))^abs(2*ui-1))-1)
    rho2(1 , :) = rho_int + yi*(UBR - LBR)
    thick2(1, :) = thick_int + yi*(UBT - LBT)
    if rho2 < LBR
        rho2 = LBR;
    end
    if rho2>UBR
        rho2 = UBR;
    end
    if thick2 < LBT
        thick2 = LBT;
    end
    if thick2 > UBT
        thick2 = UBT;
    end
    [TEM2] = fwd_TEM(rho2(1,:),thick2(1,:),t,a,I);
    dBdt2(1,:) = TEM2;
    [misfit2] = misfit_TEM(TEM_sin,dBdt2(1,:));
    E2 = misfit2;
    delta_E = E2 -E1;
    if delta_E < 0
        rho1 = rho2;
        thick1 = thick2;
         E1 = E2;
    else
        P = exp((-delta_E)/T);
        if  P >= rand
           rho1 = rho2;
           thick1 = thick2;
           E1 = E2;
        end
    end
    [TEM_new] = fwd_TEM(rho1,thick1,t,a,I);
    Egen(itr)=E1;
    T = T*exp(-dec*(itr)^(1/(2*nlayer)-1));
    Temperature(itr) = T;
    %Ploting model
r_plot = [0, R];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,rho1];
Depth_mod = [0,cumsum(thick1),max(thick1)*100];

figure(1)
subplot(1,6,[1 3])
loglog(t,TEM_new,'r',t,TEM_sin,'ob','MarkerSize',6,'LineWidth',2.5);
ylim([10^-20 10^-2])
legend({'Calculated Data','Observed Data'},'Color','none','FontWeight','Bold');
xlabel('AB/2 (m)','FontSize',8,'FontWeight','Bold');
ylabel('dBdt (V/\itm^{2}) (Ohm.m)','FontSize',8,'FontWeight','Bold');
title(['\bf \fontsize{10}\fontname{Times}TEM Respon  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);

grid on
subplot(1,6,[5 6])
stairs(r_plot,t_plot,'--r','Linewidth',2);
hold on
stairs(r_mod,Depth_mod,'-b','Linewidth',2);
hold off
legend({'Synthetic Model','Inversion Model'},'Color','none','FontWeight','Bold','Location','southwest');
%legend({'Inversion Model'},'Color','none','FontWeight','Bold','Location','southwest');
axis([1 10^5 1 1000])
xlabel('Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
ylabel('Depth (m)','FontSize',8,'FontWeight','Bold');
title('\bf \fontsize{10} Model');
set(gca,'YDir','Reverse');
set(gca, 'XScale', 'log');
set(gcf, 'Position', get(0, 'Screensize'));
grid on
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);



%plot misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

%Plot Temperature
figure(3)
plot(1:nitr,Temperature,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Temperature','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Penurunan Temperature ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

saveas(figure(2),'Misfit TEM VFSA.png')
saveas(figure(3),'Temperature TEM VFSA.png')
