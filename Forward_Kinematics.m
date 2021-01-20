function A_lump = Forward_Kinematics(DH)
%DH is a 4*n matrix, the i-th column represents the DH parameters from link
%i-1 to link i - d_i, theta_i, a_i, alpha_i
n=length(DH(1,:)); % the number of columns of DH equals the number of joints
A_cur=[eye(3) zeros(3,1);...
       zeros(1,3) 1];
for i=1:n
    d(i)=DH(1,i);theta(i)=DH(2,i);a(i)=DH(3,i);alpha(i)=DH(4,i); % get the DH parameters from link i-1 to link i from DH
    A_cur=A_cur*[cos(theta(i)) -sin(theta(i))*cos(alpha(i)) sin(theta(i))*sin(alpha(i)) a(i)*cos(theta(i));...
        sin(theta(i)) cos(theta(i))*cos(alpha(i)) -cos(theta(i))*sin(alpha(i)) a(i)*sin(theta(i));...
        0 sin(alpha(i)) cos(alpha(i)) d(i);...
        0 0 0 1]; % from link i-1 to link i
    [U S V]=svd(A_cur(1:3,1:3));
     A_cur(1:3,1:3) = U*V';%use Singular value decomposition to enforce an orthogonal (or unitary) matrix
    A_lump(:,:,i)=A_cur; % store A^0_i into the i-th element of A_lump
end

