clear all
clc
close all
Jtype=[0 0 1];% type of joints. 0 means rotational, 1 means tranlational
q(1)=0;q(2)=0;q(3)=0;% define joint variables
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
        param.height=0.5
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
    %%%
    xlabel('x');ylabel('y');zlabel('z');
    view([45 45])
    axis equal
    hold on
end