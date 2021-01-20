function tau = Newton_Euler_Iteration(DH, m, I, r_C, Jtype, q_dot, q_ddot, g0)
%Given DH parameters, masses, inertia tensor, com positions, joint types, joint positions, velocities,
%accelerations, and gravity vector, calculate the joint torques
n=length(DH(1,:)); % the number of columns of DH equals the number of joints
% forward iteration of velocity and acceleration
for i=1:n
    A_cur = Forward_Kinematics(DH(:,i));
    A_cur=A_cur(:,:,1);
    R_pre_cur = A_cur(1:3,1:3);% get R^{i-1}_i
    r_pre_cur_i = R_pre_cur'*A_cur(1:3,4);% get r_{i-1,i}^i   
    % compute omega_i^i, \dot{omega}_i^i, \ddot{p}_i^i, and \ddot{p}_{C_i}^i
    % recursively
    if i~=1
        omega_pre_i=omega_i(:,i-1);
        omega_dot_pre_i=omega_dot_i(:,i-1);
        p_ddot_pre_i=p_ddot_i(:,i-1);
    else
        omega_pre_i=[0;0;0];
        omega_dot_pre_i=[0;0;0];
        p_ddot_pre_i=[0;0;0];
    end
    if Jtype(i)==0
        omega_i(:,i)=R_pre_cur'*(omega_pre_i+[0;0;1]*q_dot(i));
        omega_dot_i(:,i)=R_pre_cur'*(omega_dot_pre_i+[0;0;1]*q_ddot(i)+cross(omega_pre_i,[0;0;1])*q_dot(i));
        p_ddot_i(:,i)=R_pre_cur'*p_ddot_pre_i+cross(omega_dot_i(:,i),r_pre_cur_i)+cross(omega_i(:,i),cross(omega_i(:,i),r_pre_cur_i));
    else
        omega_i(:,i)=R_pre_cur'*omega_pre_i;
        omega_dot_i(:,i)=R_pre_cur'*omega_dot_pre_i;
        p_ddot_i(:,i)=R_pre_cur'*(p_ddot_pre_i+[0;0;1]*q_ddot(i))+cross(omega_i(:,i),R_pre_cur'*[0;0;1])*2*q_dot(i)+cross(omega_dot_i(:,i),r_pre_cur_i)+cross(omega_i(:,i),cross(omega_i(:,i),r_pre_cur_i));
    end
    pCi_ddot_i(:,i)=p_ddot_i(:,i)+cross(omega_dot_i(:,i),r_C(:,i))+cross(omega_i(:,i),cross(omega_i(:,i),r_C(:,i)));
end
%backward iteration of forces and torques
f_i(:,n+1)=[0;0;0];
mu_i(:,n+1)=[0;0;0];
for i=n:-1:1
    A_pre_cur = Forward_Kinematics(DH(:,i));
    A_pre_cur=A_pre_cur(:,:,1);
    R_pre_cur = A_pre_cur(1:3,1:3);% get R^{i-1}_i
    A_cur = Forward_Kinematics(DH(:,1:i));
    A_cur=A_cur(:,:,i);
    R_cur = A_cur(1:3,1:3);% get R_i
    % get R^{i}_i+1
    if i==n
        R_cur_next = eye(3);
    else
        A_cur_next = Forward_Kinematics(DH(:,i+1));
        A_cur_next = A_cur_next(:,:,1);
        R_cur_next = A_cur_next(1:3,1:3);
    end
    r_pre_cur_i = R_pre_cur'*A_pre_cur(1:3,4); % get r_{i-1,i}^i  
    %compute f_i, mu_i and tau_i recursively
    f_i(:,i)=R_cur_next*f_i(:,i+1)+pCi_ddot_i(:,i)*m(i)-R_cur'*g0*m(i);
    mu_i(:,i)=-cross(f_i(:,i),r_pre_cur_i+r_C(:,i)) + R_cur_next*mu_i(:,i+1) + cross(R_cur_next*f_i(:,i+1),r_C(:,i)) + I(:,:,i)*omega_dot_i(:,i) + cross(omega_i(:,i),I(:,:,i)*omega_i(:,i));
    if Jtype(i)==0
        tau(i,1)=mu_i(:,i)'*R_pre_cur'*[0;0;1];
    else
        tau(i,1)=f_i(:,i)'*R_pre_cur'*[0;0;1];
    end
end

