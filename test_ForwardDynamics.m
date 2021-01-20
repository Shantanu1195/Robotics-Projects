clear all
clc
close all
Jtype=[0 0 1];% type of joints. 0 means rotational, 1 means tranlational
q=[0;pi/20;0.5];% define initial joint variables
q_dot=[0;0;0];% define initial joint velocities
DH=[[0;q(1);0;-pi/2] [2;q(2);0;pi/2] [q(3);0;0;0]]; %define DH parameters
n=3;%number of joints
m=[1 1 1]; %define mass constants
I(:,:,1)=eye(3);I(:,:,2)=eye(3);I(:,:,3)=eye(3);%define inertia tensor
r_C=[[0;0;0] [0;-1;0] [0;0;0]];%define com coordinates of each link
kv=[0.5;0.5;0.5];%define viscous friction constants
g0=[0;0;-9.8];%define gravity vector

Joint{1}.length=0.5;
Joint{1}.radius=0.1;
Joint{1}.zdisplacement=0;%z-coordinate of the actual joint in frame {0}
Joint{2}.length=0.5;
Joint{2}.radius=0.1;
Joint{2}.zdisplacement=2;%z-coordinate of the actual joint in frame {1}
Joint{3}.width=0.1;
Joint{3}.length=0.1;
Joint{3}.height=2;
Joint{3}.zdisplacement=0;%z-coordinate of the actual joint in frame {2}
Link{1}.radius=0.1;
Link{1}.midpoint=[0;0;1];%mid point of link 1 represented in frame {1}
Link{1}.length=2;%lenght of link 1
Link{1}.R=eye(3);% rotation matrix of the link cylinder 1 wrt frame {1}
Link{2}.radius=0.1;
Link{2}.midpoint=[0;0;0];%mid point of link 2 represented in frame {2}
Link{2}.length=0;%lenght of link 2
Link{2}.R=eye(3);% rotation matrix of the link cylinder 2 wrt frame {2}
Link{3}.radius=0.1;
Link{3}.midpoint=[0;0;0];%mid point of link 3 represented in frame {3}
Link{3}.length=0.25;%lenght of link 3
Link{3}.R=eye(3);% rotation matrix of the link cylinder 3 wrt frame {3}

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
%     if i==n        
%         %%%draw end-effector velocity vectors
%         hpdot_e=line([A(1,4) A(1,4)],[A(2,4) A(2,4)],[A(3,4) A(3,4)],'Color','k');
%     end
end
%%%
xlabel('x');ylabel('y');zlabel('z');
axis equal
axis([-1 2 -1 2 -1 2]);%set display range xyz
view([45 45])
hold on

% xbox = RobotRaconteur.Connect('tcp://localhost:5437/Xbox_controllerServer/xbox_controller');

t=0;
tic
while 1
    dt=toc;
    t=t+dt;
    tic
    for i=1:n
        q_ddot_temp=[0;0;0];
        q_ddot_temp(i)=1;
        B(:,i) = Newton_Euler_Iteration(DH, m, I, r_C, Jtype, [0;0;0], q_ddot_temp, [0;0;0]);
    end
    Coriolis_Gravity=Newton_Euler_Iteration(DH, m, I, r_C, Jtype, q_dot, [0;0;0], g0);%compute C(q,\dot{q})\dot{q}+G(q)
    tau=[0;0;0];%actual joint torques
    ks=20;%add a virtual spring to prevent the 3rd prismatic link from sliding down to infinity
    if q(3)>=0.5
        Fs=-ks*(q(3)-0.5);
    elseif q(3)<=-0.5
        Fs=-ks*(q(3)+0.5);
    else
        Fs=0;
    end
    q_ddot=inv(B)*(tau-Coriolis_Gravity-kv.*q_dot+[0;0;Fs]);
    q_dot=q_dot+q_ddot.*dt;%update q_dot
    q=q+q_dot.*dt;%update q
    DH(2,1)=q(1);DH(2,2)=q(2);DH(1,3)=q(3);%update DH parameters
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
    view([45 45])
    drawnow;
end



