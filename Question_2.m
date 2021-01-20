clear all; clc; close all

% Determine initial conditions.
r = [0;0];                      % These are initial joint variables.
r_dot = [0;0];                  % It is defined as initial joint velocities.

% Put Simulation constants over here.
runtime = 17;                   % Runtime in sex
dt = 0.020;                     % Time Step
g = 9.81;                       % Gravity in m/s^2
v_friction = 0;                 % Viscous Friction Constant

% We are defining two link manipulators over here.
L1_length = 1;                  % Define link 1 length
L2_length = L1_length;          % Define link 2 length
n = 2;                          % Number of Joints
Joint_type = [0 0];             % Joint Type (Here, 0=rotational & 1= translational)
j1_torque = 0;                  % it is defined as a joint 1 torque.
j2_torque = 0;                  % it is defined as joint 2 torque.
m_1 = 5;                         % mass of link 1
m_2 = m_1;                        % mass of link 2
l_1 = 0.5;                       % Center of mass on link 1
l_2 = 0.5;                       % Center of mass of link 2

% lets, define DH parameters over here.
a_1 = L1_length;
a_2 = L2_length;
DH = [[0;r(1);a_1;0] [0;r(2);a_2;0]];
g0 = [0;-g;0];                          % Define Gravity Vector
kv = [v_friction;v_friction];           % Define Viscous Friction Constants
r_C=[[-l_1;0;0] [-l_2;0;0]];              % Define CoM Coordinates of Each Link
m = [m_1 m_2];                            % Define Mass Vector
MI_L1 = (1/3) * m_1 * (L1_length^2);    % Calculate Moment of Inertia, Link 1
MI_L2 = (1/3) * m_2 * (L2_length^2);    % Calculate Moment of Inertia, Link 2
I(:,:,1)=[0 0 0;0 MI_L1 0;0 0 MI_L1];       % Define Inertia Tensor for Link 1
I(:,:,2)=[0 0 0;0 MI_L2 0;0 0 MI_L2];       % Define Inertia Tensor for Link 2



% lets, define joint1 over here.
Joint{1}.length= 0.5;
Joint{1}.radius=0.1;
Joint{1}.zdisplacement=0;

% lets, define joint2 over here.
Joint{2}.length=0.5;
Joint{2}.radius=0.1;
Joint{2}.zdisplacement=0;

% lets, define end-effector over here.
Joint{3}.width=0.1;

% lets, define link1 over here.
Link{1}.radius=0.1;
Link{1}.midpoint=[-0.5 * L1_length;0;0];
Link{1}.length = L1_length;
Link{1}.R=rot([0;1;0],pi/2);

% lets, define link2 over here.
Link{2}.radius=0.1;
Link{2}.midpoint=[-0.5 * L2_length;0;0];
Link{2}.length = L2_length;
Link{2}.R=rot([0;1;0],pi/2);

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
    if Joint_type(i)==0
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
    param.width=0.02;
    param.scale=0.2;
    if i==1
        hglobalframe=create3DFrame(Apre(1:3,1:3), Apre(1:3,4), param);
    end
    hframe{i}=create3DFrame(A(1:3,1:3), A(1:3,4), param);
end
%%%

% Here, we are creating link visual
figure(1)
xlabel('x');ylabel('y');zlabel('z');
axis equal
axis([-(L1_length + L2_length) (L1_length + L2_length)...
        -(L1_length + L2_length) (L1_length + L2_length)...
        -(L1_length + L2_length) (L1_length + L2_length)]);
view([45 45])
hold on

% Here, we are setting initial values.
ThetaStack_1 = [];
ThetaStack_2 = [];
T_Stack = [];
t = 0;

while t <= runtime;
    t = t + dt;
    
    % Dynamics Calculations
    for i = 1:n
        q_ddot_temp = [0;0];
        q_ddot_temp(i) = 1;
        B(:,i) = Newton_Euler_Iteration(DH, m, I, r_C, Joint_type, [0;0], q_ddot_temp, [0;0;0]);
    end
    Theta_1 = DH(2,1);
    Theta_2 = DH(2,2);
    Coriolis_Gravity = Newton_Euler_Iteration(DH, m, I, r_C, Joint_type, r_dot, [0;0], g0); % Compute C(q,\dot{q})\dot{q}+G(q)
    Tau = [j1_torque;j2_torque]; % Actual Joint Torques
    r_ddot = inv(B)*(Tau-Coriolis_Gravity-kv.*r_dot); % Update r_ddot
    r_dot = r_dot+r_ddot.*dt; % Update r_dot
    r = r+r_dot.*dt; % Update r
    DH(2,1) = r(1); % Update DH Parameters
    DH(2,2) = r(2); % Update DH Parameters
    A_lump = Forward_Kinematics(DH); % Calculate Homogeneous Transformation Matrices

    % here, Wr are updating drawings
    for i=1:n
        if i==1
            Apre=[eye(3) zeros(3,1);...
                zeros(1,3) 1];
        else
            Apre=A_lump(:,:,i-1);
        end
        Acur=A_lump(:,:,i);
        hjoint{i}=updateRigidBody(Apre(1:3,1:3), Apre(1:3,4)+Apre(1:3,1:3)*[0;0;Joint{i}.zdisplacement],hjoint{i}); % Update Joint Drawings
        hlink{i}=updateRigidBody(Acur(1:3,1:3)*Link{i}.R, Acur(1:3,4)+Acur(1:3,1:3)*Link{i}.midpoint,hlink{i}); % Update Link Drawings
        hframe{i}=updateRigidBody(Acur(1:3,1:3), Acur(1:3,4),hframe{i}); % Update Frame Drawings
    end    

    drawnow;
    
    % Generate Theta 1, Theta 2, and Time Data
    ThetaStack_1 = [ThetaStack_1 Theta_1];
    ThetaStack_2 = [ThetaStack_2 Theta_2];
    T_Stack = [T_Stack t];

end

% Plot Theta 1 and Theta 2
figure(2)
axis equal
axis([0 runtime 0 5])
plot(T_Stack, ThetaStack_1);
hold on
plot(T_Stack, ThetaStack_2);
legend('Theta 1','Theta 2')
title('Theta 1 and Theta 2 vs. Time')
xlabel('Time (s)')
ylabel('Radians')
