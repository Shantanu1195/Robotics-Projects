function handle = updateRigidBody(R, t, handle)
    %
    % handle = updateRigidBody(R, t, handle)
    %
    % R is desired orientation of the body
    % t is desired position of the body
    % handle is drawing structure for body undergoing transformation
    %
    % returns updated structure
    
    if isempty(handle);  return;  end
    if nargout ~= 1; error('Must have return value for handle'); end;
    
    % Determine local transformation
    R12 = R*handle.R';
    t12 = t - R12*handle.t;
    
    for i=1:length(handle.bodies)
        V = get(handle.bodies(i),'Vertices');
        dim = size(V);
        V = V*R12' + ones(dim(1),1)*t12';
%         V = (V - ones(dim(1),1)*handle.t')*handle.R*R'+ones(dim(1),1)*t';
        set(handle.bodies(i),'Vertices',V);
    end
    handle.R = R;
    handle.t = t;
end