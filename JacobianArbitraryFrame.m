function J = JacobianArbitraryFrame(DH,Jtype,LinkNum,p_rel)
%compute Jacobian with the "end-effector" frame rigidly attached to link
%LinkNum, with relative position p_rel w.r.t the origin of the frame LinkNum
n=length(DH(1,:)); % the number of columns of DH equals the number of joints
A_lump = Forward_Kinematics(DH); %get the homogeneous transformation matrices for all the frames
A_LinkNum=A_lump(:,:,LinkNum);
p_LinkNum=A_LinkNum(1:3,4); % calculate p(LinkNum)
R_LinkNum=A_LinkNum(1:3,1:3); % calculate R(LinkNum)
p_e=p_LinkNum+R_LinkNum*p_rel;%calculate p_e
J=[]; %initialize J
for i=1:n
    %% get z_{i-1} and p_{i-1} from the homogeneous transformation matrices
    if i==1
        z_pre=[0;0;1];
        p_pre=[0;0;0];
    else
        A_pre=A_lump(:,:,i-1);
        z_pre=A_pre(1:3,1:3)*[0;0;1];
        p_pre=A_pre(1:3,4);
    end
    %% calculate 6*n Jacobian matrix
    if i>LinkNum
        J=zeros(6,1);
    else
        if Jtype(i)==0
            J=[J [cross(z_pre,p_e-p_pre);z_pre]];
        else
            J=[J [z_pre;zeros(3,1)]];
        end
    end
end

