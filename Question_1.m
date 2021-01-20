function [s,s_dot,s_ddot] = S_Curve(t,t0,s0,sf,vmax,amax)
t = 0;
t0 = 0;
s0 = 0;
sf = 15;
vmax =25;
amax=5;
if sf>s0
    if sf-s0>=vmax^2/amax
        t1=t0+vmax/amax;
        t2=t0+vmax/amax+(sf-s0-vmax^2/amax)/vmax;
        tf=t0+2*vmax/amax+(sf-s0-vmax^2/amax)/vmax;
        if t<t0
            s=s0;
            s_dot=0;
            s_ddot=0;
        elseif t>=t0 & t<=t1
            s=s0+0.5*amax*(t-t0)^2;
            s_dot=amax*(t-t0);
            s_ddot=amax;
        elseif t>t1 & t<=t2
            s=s0+0.5*vmax^2/amax+vmax*(t-t1);
            s_dot=vmax;
            s_ddot=0;
        elseif t>t2 & t<=tf
            s=sf-0.5*amax*(t-tf)^2;
            s_dot=amax*(tf-t);
            s_ddot=-amax;
        else
            s=sf;
            s_dot=0;
            s_ddot=0;
        end
    else
        t1=t0+sqrt((sf-s0)/amax);
        tf=t0+2*sqrt((sf-s0)/amax);
        if t<t0
            s=s0;
            s_dot=0;
            s_ddot=0;
        elseif t>=t0 & t<=t1
            s=s0+0.5*amax*(t-t0)^2;
            s_dot=amax*(t-t0);
            s_ddot=amax;
        elseif t>t1 & t<=tf
            s=sf-0.5*amax*(t-tf)^2;
            s_dot=amax*(tf-t);
            s_ddot=-amax;
        else
            s=sf;
            s_dot=0;
            s_ddot=0;
        end
    end
else
    if s0-sf>=vmax^2/amax
        t1=t0+vmax/amax;
        t2=t0+vmax/amax+(s0-sf-vmax^2/amax)/vmax;
        tf=t0+2*vmax/amax+(s0-sf-vmax^2/amax)/vmax;
        if t<t0
            s=s0;
            s_dot=0;
            s_ddot=0;
        elseif t>=t0 & t<=t1
            s=s0-0.5*amax*(t-t0)^2;
            s_dot=-amax*(t-t0);
            s_ddot=-amax;
        elseif t>t1 & t<=t2
            s=s0-0.5*vmax^2/amax-vmax*(t-t1);
            s_dot=-vmax;
            s_ddot=0;
        elseif t>t2 & t<=tf
            s=sf+0.5*amax*(t-tf)^2;
            s_dot=-amax*(tf-t);
            s_ddot=amax;
        else
            s=sf;
            s_dot=0;
            s_ddot=0;
        end
    else
        t1=t0+sqrt((s0-sf)/amax);
        tf=t0+2*sqrt((s0-sf)/amax);
        if t<t0
            s=s0;
            s_dot=0;
            s_ddot=0;
        elseif t>=t0 & t<=t1
            s=s0-0.5*amax*(t-t0)^2;
            s_dot=-amax*(t-t0);
            s_ddot=-amax;
        elseif t>t1 & t<=tf
            s=sf+0.5*amax*(t-tf)^2;
            s_dot=-amax*(tf-t);
            s_ddot=amax;
        else
            s=sf;
            s_dot=0;
            s_ddot=0;
        end
    end
end



t_stack = [];
s_stack = [];
s_dot_stack = [];
s_ddot_stack = [];


while t<= 10
    [s,s_dot,s_ddot] = S_Curve(t,t0,s0,sf,vmax,amax);
    t = t + 0.02;
    t_stack = [t_stack t];
    s_stack = [s_stack s];
    s_dot_stack = [s_dot_stack s_dot];
    s_ddot_stack = [s_ddot_stack s_ddot];
end

plot(t_stack,s_stack,t_stack,s_dot_stack,t_stack,s_ddot_stack)
xlabel('Time (second)')
legend('postion','velocity','acceleration')


end