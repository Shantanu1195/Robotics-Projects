clear all
clc
close all
Jtype=[0 0];% type of joints. 0 means rotational, 1 means tranlational
q=[0;pi/6];% define initial joint variables
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

%specify desired target for end-effector
Target.radius=0.05;
Target.center=[0.5;1.5;0];

mE=1;%number of obstacles
Obstacle{1}.radius=0.15;
Obstacle{1}.center=[1;1.0;0];

mP=8;%number of control points on the robot
ControlPoint{1}.LinkNum=1;%which link this point is attached to
ControlPoint{1}.p_rel=[-0.75;0;0];%relative position wrt the frame of that link
ControlPoint{2}.LinkNum=1;
ControlPoint{2}.p_rel=[-0.5;0;0];
ControlPoint{3}.LinkNum=1;
ControlPoint{3}.p_rel=[-0.25;0;0];
ControlPoint{4}.LinkNum=1;
ControlPoint{4}.p_rel=[0;0;0];
ControlPoint{5}.LinkNum=2;
ControlPoint{5}.p_rel=[-0.75;0;0];
ControlPoint{6}.LinkNum=2;
ControlPoint{6}.p_rel=[-0.5;0;0];
ControlPoint{7}.LinkNum=2;
ControlPoint{7}.p_rel=[-0.25;0;0];
ControlPoint{8}.LinkNum=2;
ControlPoint{8}.p_rel=[0;0;0];



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
param=[];
param.radius=Obstacle{1}.radius;
hobstacle=createEllipsoid(eye(3),Obstacle{1}.center,param);
param=[];
param.radius=Target.radius;
htarget=createEllipsoid(eye(3),Target.center,param);
%%%
xlabel('x');ylabel('y');zlabel('z');
axis equal
axis([-1 2 -1 2 -1 2]);%set display range xyz
view([0 90])
hold on

%store data during simulation
t_stack=[];
e_stack=[];
vidObj = VideoWriter('PotentialField2.avi');
open(vidObj);
t=0;
tic
while t<10
    dt=toc;
    dt=0.01;
    t=t+dt;
    tic
    for i=1:n
        q_ddot_temp=[0;0];
        q_ddot_temp(i)=1;
        B(:,i) = Newton_Euler_Iteration(DH, m, I, r_C, Jtype, [0;0], q_ddot_temp, [0;0;0]);
    end
    theta1=DH(2,1);theta2=DH(2,2);
    Coriolis_Gravity=Newton_Euler_Iteration(DH, m, I, r_C, Jtype, q_dot, [0;0], g0);%compute C(q,\dot{q})\dot{q}+G(q)
    
    %potential field motion planning
    A_lump = Forward_Kinematics(DH);
    A_cur=A_lump(:,:,n);
    pe=A_cur(1:3,4);%end-effector position
    ka=1;
    Ua=ka*norm(pe(1:2)-Target.center(1:2));
    ped_dot=-ka.*(pe(1:2)-Target.center(1:2))./norm(pe(1:2)-Target.center(1:2));
    J=Jacobian(DH,Jtype);
    qd_dot=J(1:2,:)'*inv(J(1:2,:)*J(1:2,:)'+0.01.*eye(2))*ped_dot;%first, set qd_dot to be corresponding to the attractive potential
    for i=1:mP
        A_cur=A_lump(:,:,ControlPoint{i}.LinkNum);
        pP(:,i)=A_cur(1:3,4)+A_cur(1:3,1:3)*ControlPoint{i}.p_rel;%calculate the coordinate of i-th control point
        Ji=JacobianArbitraryFrame(DH,Jtype,ControlPoint{i}.LinkNum,ControlPoint{i}.p_rel);%calculate the Jacobian of i-th control point
        Ji=Ji(1:3,:);
        for j=1:mE
            eta0=0.2;
            kr=0.5;
            gamma=2;
            eta=norm(pP(:,i)-Obstacle{j}.center)-Obstacle{1}.radius;
            eta_dotp=(pP(:,i)-Obstacle{j}.center)./norm(pP(:,i)-Obstacle{j}.center);
            if eta<=eta0
                Ur=1/gamma*(1/eta-1/eta0)^gamma;
                Ur_doteta=-1/eta^2*(1/eta-1/eta0)^(gamma-1);
            else
                Ur=0;
                Ur_doteta=0;
            end
            vPijd=-kr.*Ur_doteta.*eta_dotp;
            qd_dot=qd_dot+Ji(1:2,:)'*inv(Ji(1:2,:)*Ji(1:2,:)'+0.01.*eye(2))*vPijd(1:2);%for each i,j pair of, add the repulsive potential terms to qd_dot
        end
    end
    %limit joint velocity to [-1 1]rad/sec
    for i=1:n
        if qd_dot(i)>1
            qd_dot(i)=1;
        elseif qd_dot(i)<-1
            qd_dot(i)=-1;
        end
    end
            
    J=Jacobian(DH,Jtype);
    
    e_dot=q_dot-qd_dot;
    
    % centralized velocity control
    Kd=5.*eye(2);
    tau=-B*Kd*e_dot+Coriolis_Gravity+kv.*q_dot;
    
    %dynamics
    Delta=[-5*sin(2*pi/1*t);-5*sin(2*pi/1*t)];
    q_ddot=inv(B)*(tau-Coriolis_Gravity-kv.*q_dot+Delta);
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
    view([0 90])
    drawnow;
    if mod(t/dt,4)>0 & mod(t/dt,4)<=1
        writeVideo(vidObj,getframe(gcf));
    end
    t_stack=[t_stack t];
end
close(vidObj);

