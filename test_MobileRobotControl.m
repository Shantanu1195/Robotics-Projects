clear all
clc
close all

x=0;y=0;theta=pi/2;


RobotBase.width=1;
RobotBase.length=0.5;
RobotBase.height=0.3;
RobotBase.gap=0.1;
RobotWheel.radius=0.2;
RobotWheel.thickness=0.02;
b=0.1;


param.width=RobotBase.width;
param.length=RobotBase.length;
param.height=RobotBase.height;
hRobotBase = createCuboid(rot([0;0;1],theta)*eye(3), [x;y;0]+rot([0;0;1],theta)*[0;0;RobotBase.gap+RobotBase.height/2], param);
param=[];
param.radius=RobotWheel.radius;
param.height=RobotWheel.thickness;
hLeftWheel = createCylinder(rot([0;0;1],theta)*rot([0;1;0],pi/2), [x;y;0]+rot([0;0;1],theta)*[0;RobotBase.width/2+RobotWheel.thickness/2;RobotWheel.radius], param);
hRightWheel = createCylinder(rot([0;0;1],theta)*rot([0;1;0],-pi/2), [x;y;0]+rot([0;0;1],theta)*[0;-RobotBase.width/2-RobotWheel.thickness/2;RobotWheel.radius], param);

t=0;
tic
while 1
    dt=toc;
    dt=0.02;
    t=t+dt;
    
    xBd=1-cos(0.5*t);
    xBd_dot=0.5*sin(0.5*t);
    yBd=sin(0.5*t);
    yBd_dot=0.5*cos(0.5*t);
    xB=x+b*cos(theta);
    yB=y+b*sin(theta);
    
    %controller
    k=10;
    ux=xBd_dot-k*(xB-xBd);
    uy=yBd_dot-k*(yB-yBd);
    v=cos(theta)*ux+sin(theta)*uy;
    omega=-sin(theta)/b*ux+cos(theta)/b*uy;
    
    %kinematic level update
    x=x+dt*cos(theta)*v;
    y=y+dt*sin(theta)*v;
    theta=theta+dt*omega;
    
    hRobotBase = updateRigidBody(rot([0;0;1],theta)*eye(3), [x;y;0]+rot([0;0;1],theta)*[0;0;RobotBase.gap+RobotBase.height/2], hRobotBase);
    hLeftWheel = updateRigidBody(rot([0;0;1],theta)*rot([1;0;0],-pi/2), [x;y;0]+rot([0;0;1],theta)*[0;RobotBase.length/2+RobotWheel.thickness/2;RobotWheel.radius], hLeftWheel);
    hRightWheel = updateRigidBody(rot([0;0;1],theta)*rot([1;0;0],pi/2), [x;y;0]+rot([0;0;1],theta)*[0;-RobotBase.length/2-RobotWheel.thickness/2;RobotWheel.radius], hRightWheel);axis equal
    axis([-1 2 -1 2 -1 2]);%set display range xyz
    view([45 45])
    drawnow
end

