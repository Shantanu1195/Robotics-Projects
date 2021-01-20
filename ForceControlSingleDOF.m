clear all
clc
close all
M=50;
k=500000;
Fd=20;
xc=1;
%integral force control
x=xc;
x_dot=0;
f_real=0;
I=0;
t=0;
t_stack=[];
f_stack=[];
F_stack=[];
while t<5
    dt=0.001;
    t=t+dt;
    if x<xc
        F=0;
    else
        F=k*(x-xc);
    end

    %direct force control (PID)
    k_hat=1000;
    if x<xc
        F_dot_hat=0;
    else
        F_dot_hat=k_hat*x_dot;
    end
    kp=100/k_hat;kd=20/k_hat;ki=0/k_hat;
    I=I+(F-Fd)*dt;%update the integrator
    f=F-kp*M*(F-Fd)-kd*M*F_dot_hat-ki*I;
    
%     %indirect force control (impedance control)
%     Mm=50;Dm=400;Km=10;
%     xd=xc+Fd/Km;
%     f=F+M/Mm*(-Dm*(x_dot)-Km*(x-xd)-F);
    
    %state update based on system dynamics
    f_real=f_real-2*pi*100*dt*(f_real-f);
%     f_real=f;
    x_dot=x_dot+1/M*(f_real-F)*dt;
    x=x+x_dot*dt;
    
    t_stack=[t_stack t];
    f_stack=[f_stack f];
    F_stack=[F_stack F];
end
figure(1)
plot(t_stack,F_stack(1,:),'r');
hold on
plot(t_stack,20.*ones(1,length(t_stack)),'k--');
xlim([0 5]);
xlabel('time (sec)');
ylabel('contact force (N)');
xh=get(gca,'xlabel');
yh=get(gca,'ylabel');
set(xh,'FontSize',16,'FontWeight','bold');
set(yh,'FontSize',16,'FontWeight','bold');
set(gca,'FontSize',16,'FontWeight','bold');
grid on    
    