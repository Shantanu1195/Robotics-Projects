clear all
clc
close all
Jtype=[0 0];% type of joints. 0 means rotational, 1 means tranlational
q=[0;pi/2];% define initial joint variables
q_dot=[0;0];% define initial joint velocities
a1=1;a2=1;
DH=[[0;q(1);a1;0] [0;q(2);a2;0]]; %define DH parameters
n=2;%number of joints
m1=50;m2=50;
m=[m1 m2]; %define mass constants
Il1=10;Il2=10;
I(:,:,1)=[0 0 0;0 Il1 0;0 0 Il1];I(:,:,2)=[0 0 0;0 Il2 0;0 0 Il2];%define inertia tensor
l1=0.5;l2=0.5;
r_C=[[l1-a1;0;0] [l2-a2;0;0]];%define com coordinates of each link
kv=[50;50];%define viscous friction constants
g=9.8;
g0=[0;-g;0];%define gravity vector
lambdafd=20;%desired contact force
J=Jacobian(DH,Jtype);
JA=J(1:2,:);%initialize JA for calculation of \dot{JA}

Joint{1}.length=0.5;
Joint{1}.radius=0.1;
Joint{1}.zdisplacement=0;%z-coordinate of the actual joint in frame {0}
Joint{2}.length=0.5;
Joint{2}.radius=0.1;
Joint{2}.zdisplacement=0;%z-coordinate of the actual joint in frame {1}
Joint{3}.width=0.1;
Link{1}.radius=0.1;
Link{1}.midpoint=[-0.5;0;0];%mid point of link 1 represented in frame {1}
Link{1}.length=1;%lenght of link 1
Link{1}.R=rot([0;1;0],pi/2);% rotation matrix of the link cylinder 1 wrt frame {1}
Link{2}.radius=0.1;
Link{2}.midpoint=[-0.5;0;0];%mid point of link 2 represented in frame {2}
Link{2}.length=1;%lenght of link 2
Link{2}.R=rot([0;1;0],pi/2);% rotation matrix of the link cylinder 2 wrt frame {2}

Wall.param.length=4;
Wall.param.width=0.05;
Wall.param.height=1;
hWall=createCuboid(rot([0;0;1],pi/4), [1;1;0], Wall.param);

A_lump = Forward_Kinematics(DH);
for i=1:n
    if i==1
        Apre=[eye(3) zeros(3,1);...
       zeros(1,3) 1];
    else
        Apre=A_lump(:,:,i-1);
    end
    A=A_lump(:,:,i);
    %%%draw joints
    if Jtype(i)==0
        param.radius=Joint{i}.radius;
        param.height=Joint{i}.length;
        hjoint{i}=createCylinder(Apre(1:3,1:3), Apre(1:3,4)+Apre(1:3,1:3)*[0;0;Joint{i}.zdisplacement], param);
    else
        param.width=Joint{i}.width;
        param.length=Joint{i}.length;
        param.height=Joint{i}.height;
        hjoint{i} = createCuboid(Apre(1:3,1:3), Apre(1:3,4)+Apre(1:3,1:3)*[0;0;Joint{i}.zdisplacement], param);
    end
    param=[];
    %%%draw links
    Acur=A_lump(:,:,i);
    param.radius=Link{i}.radius;
    param.height=Link{i}.length;
    hlink{i}=createCylinder(Acur(1:3,1:3)*Link{i}.R, Acur(1:3,4)+Acur(1:3,1:3)*Link{i}.midpoint, param);
    param=[];
    %%%draw frames
    param.scale=0.2;
    param.width=0.02;
    if i==1
        hglobalframe=create3DFrame(Apre(1:3,1:3), Apre(1:3,4), param);
    end
    hframe{i}=create3DFrame(A(1:3,1:3), A(1:3,4), param);
end
%%%
xlabel('x');ylabel('y');zlabel('z');
axis equal
axis([-1 2 -1 2 -1 2]);%set display range xyz
view([45 90])
hold on

%store data during simulation
t_stack=[];
lambdaf_stack=[];
ev_stack=[];

t=0;
tic
while t<20
    dt=toc;
    dt=0.02;
    t=t+dt;
    tic
    for i=1:n
        q_ddot_temp=[0;0];
        q_ddot_temp(i)=1;
        B(:,i) = Newton_Euler_Iteration(DH, m, I, r_C, Jtype, [0;0], q_ddot_temp, [0;0;0]);
    end
    theta1=DH(2,1);theta2=DH(2,2);
    Coriolis_Gravity=Newton_Euler_Iteration(DH, m, I, r_C, Jtype, q_dot, [0;0], g0);%compute C(q,\dot{q})\dot{q}+G(q)
    
    
    J=Jacobian(DH,Jtype);
    Sf=[cos(pi/4);sin(pi/4)];
    Sv=[-sin(pi/4);cos(pi/4)];
    JA_new=J(1:2,:);
    JA_dot=(JA_new-JA)./dt;
    JA=JA_new;
    Bm=50;Dm=400;Km=10;
    Ae=A_lump(:,:,n);
    temp=[Sf Sv]'*[Ae(1,4);Ae(2,4)];
    sf=temp(1);sv=temp(2);    
    temp=[Sf Sv]'*JA*q_dot;
    sf_dot=temp(1);sv_dot=temp(2);
    sfc=sqrt(2);
    if sf<sfc
        lambdaf=0;
    else
        k=1000;
        lambdaf=k*(sf-sfc);
    end
    lambdav=0;
    
    %desired trajectory in unconstrained direction
    svd=0.5*sin(2*pi/2*t);
    svd_dot=0.5*2*pi/2*cos(2*pi/2*t);
    svd_ddot=-0.5*(2*pi/2)^2*sin(2*pi/2*t);
    
    %position control law in unconstrained direction
    ev=sv-svd;
    ev_dot=sv_dot-svd_dot;
    Kp=25;Kd=5;
    av=svd_ddot-Kp*ev-Kd*ev_dot;
    
    %direct force control in constrained direction
    k_hat=1000;%estimation of contact stiffness k
    %calculate the estimate of velocity in the constrained direction
    %(lambdaf_dot)
    if sf<sfc
        lambdaf_dot_hat=0;
    else
        lambdaf_dot_hat=k_hat*sf_dot;
    end
    kp=25;kd=5;
    af=-kp/k_hat*(lambdaf-lambdafd)-kd/k_hat*lambdaf_dot_hat;
    
    % impedance control in constrained direction
%     Mm=50;Dm=400;Km=10;
%     sfd=sfc+lambdafd/Km;
%     af=1/Mm*(-Dm*(sf_dot)-Km*(sf-sfd)-lambdaf);
    
    %total control law
    tau=Coriolis_Gravity+kv.*q_dot+JA'*(Sf*lambdaf+Sv*lambdav)+B*JA'*inv(JA*JA'+0.01.*eye(2))*(-JA_dot*q_dot+Sf*af+Sv*av);
    
    %dynamics
    q_ddot=inv(B)*(tau-Coriolis_Gravity-kv.*q_dot+JA'*(-Sf*lambdaf-Sv*lambdav));
    q_dot=q_dot+q_ddot.*dt;%update q_dot
    q=q+q_dot.*dt;%update q
    DH(2,1)=q(1);DH(2,2)=q(2);%update DH parameters
    A_lump = Forward_Kinematics(DH);%get homogeneous transformation matrices
    J = Jacobian(DH,Jtype);%get Jacobian
    v_e = J*q_dot;
    %update drawings
    for i=1:n
        if i==1
            Apre=[eye(3) zeros(3,1);...
                zeros(1,3) 1];
        else
            Apre=A_lump(:,:,i-1);
        end
        Acur=A_lump(:,:,i);
        hjoint{i}=updateRigidBody(Apre(1:3,1:3), Apre(1:3,4)+Apre(1:3,1:3)*[0;0;Joint{i}.zdisplacement],hjoint{i});%update joint drawings
        hlink{i}=updateRigidBody(Acur(1:3,1:3)*Link{i}.R, Acur(1:3,4)+Acur(1:3,1:3)*Link{i}.midpoint,hlink{i});%update link drawings
        hframe{i}=updateRigidBody(Acur(1:3,1:3), Acur(1:3,4),hframe{i});%update frame drawings
    end    
    axis equal
    axis([-1 2 -1 2 -1 2]);%set display range xyz
    view([45 90])
    drawnow;
    
    t_stack=[t_stack t];
    lambdaf_stack=[lambdaf_stack lambdaf];
    ev_stack=[ev_stack ev];
end
figure(2)
subplot(2,1,1)
plot(t_stack,lambdaf_stack(1,:),'r');
xlim([0 20]);
ylabel('contact force in S_f (N)');
xh=get(gca,'xlabel');
yh=get(gca,'ylabel');
set(xh,'FontSize',16,'FontWeight','bold');
set(yh,'FontSize',16,'FontWeight','bold');
set(gca,'FontSize',16,'FontWeight','bold');
grid on
subplot(2,1,2)
plot(t_stack,ev_stack(1,:),'r');
xlim([0 20]);
xlabel('time (sec)');
ylabel('tracking error in S_v (m)');
xh=get(gca,'xlabel');
yh=get(gca,'ylabel');
set(xh,'FontSize',16,'FontWeight','bold');
set(yh,'FontSize',16,'FontWeight','bold');
set(gca,'FontSize',16,'FontWeight','bold');
grid on


