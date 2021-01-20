clear all
clc
close all
Jtype=[0 0 1];% type of joints. 0 means rotational, 1 means tranlational
q=[0;0;0];% define joint variables
DH=[[0;q(1);0;-pi/2] [2;q(2);0;pi/2] [q(3);0;0;0]]; %define DH parameters
A_lump = Forward_Kinematics(DH);
for i=1:3
    if i==1
        Apre=[eye(3) zeros(3,1);...
       zeros(1,3) 1];
    else
        Apre=A_lump(:,:,i-1);
    end
    A=A_lump(:,:,i);
    %%%draw joints
    if Jtype(i)==0
        param.radius=0.1;
        param.height=0.5;
        hjoint{i}=createCylinder(Apre(1:3,1:3), Apre(1:3,4), param);
    else
        param.width=0.1;
        param.length=0.1;
        param.height=0.5;
        hjoint{i} = createCuboid(Apre(1:3,1:3), Apre(1:3,4), param);
    end
    %%%draw links
    hlink{i}=line([Apre(1,4) A(1,4)],[Apre(2,4) A(2,4)],[Apre(3,4) A(3,4)]);
    %%%draw frames
    param.scale=0.2;
    param.width=0.02;
    if i==1
        hglobalframe=create3DFrame(Apre(1:3,1:3), Apre(1:3,4), param);
    end
    hframe{i}=create3DFrame(A(1:3,1:3), A(1:3,4), param);
    if i==3        
        %%%draw end-effector velocity vectors
        hpdot_e=line([A(1,4) A(1,4)],[A(2,4) A(2,4)],[A(3,4) A(3,4)],'Color','k');
    end
end
%%%
xlabel('x');ylabel('y');zlabel('z');
axis equal
axis([-1 2 -1 2 -1 2]);%set display range xyz
view([45 45])
hold on

xbox = RobotRaconteur.Connect('tcp://localhost:5437/Xbox_controllerServer/xbox_controller');
t=0;
tic
while 1
    dt=toc;
    t=t+dt;
    tic
    xBoxInput = xbox.controller_input;%get Xbox controller input data
    if xBoxInput.B>0.5 %simulation termination
        break;
    end
    %get input signals from Xbox controller and normalize to the range [-1
    %1]
    Input(1)=double(xBoxInput.left_thumbstick_X)/10000;
    Input(2)=double(xBoxInput.left_thumbstick_Y)/10000;
    Input(3)=double(xBoxInput.right_thumbstick_X)/10000;
    Input(4)=double(xBoxInput.right_thumbstick_Y)/10000;
    q_dot=[Input(1);Input(2);Input(3)];%map input signals to joint velocities
    q=q+q_dot.*dt;%update q
    DH(2,1)=q(1);DH(2,2)=q(2);DH(1,3)=q(3);%update DH parameters
    A_lump = Forward_Kinematics(DH);%get homogeneous transformation matrices
    J = Jacobian(DH,Jtype);%get Jacobian
    v_e = J*q_dot;
    %update drawings
    for i=1:3
        if i==1
            Apre=[eye(3) zeros(3,1);...
                zeros(1,3) 1];
        else
            Apre=A_lump(:,:,i-1);
        end
        A=A_lump(:,:,i);
        hjoint{i}=updateRigidBody(Apre(1:3,1:3), Apre(1:3,4),hjoint{i});%update joint drawings
        set(hlink{i},'XData',[Apre(1,4) A(1,4)],'YData',[Apre(2,4) A(2,4)],'ZData',[Apre(3,4) A(3,4)]);%update link drawings
        hframe{i}=updateRigidBody(A(1:3,1:3), A(1:3,4),hframe{i});%update frame drawings
        if i==3
            set(hpdot_e,'XData',[A(1,4) A(1,4)+v_e(1)],'YData',[A(2,4) A(2,4)+v_e(2)],'ZData',[A(3,4) A(3,4)+v_e(3)]);%update end-effector velocity vector drawings
        end
    end    
    axis equal
    axis([-1 2 -1 2 -1 2]);%set display range xyz
    view([45 45])
    drawnow;
end



