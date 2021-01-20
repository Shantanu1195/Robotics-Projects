function handle = combineRigidBodies(body1, varargin)
% COMBINERIGIDBODIES
%
%   handle = combineRigidBodies(body1,body2) - combines both bodies into single
%       body, retaining body1's pose
%   handle = combineRigidBodies(body1,body2,body3,...) - combines all bodies
%       into single body, retaining body1's pose
%
%   body1 is a struct with the rigid body information, which will be the 
%       main body to combine other bodies onto
%   body2, body3, etc. are handles to other bodies to attach to body1.
%
%   handle is the struct returned with all the bodies combined together
%       and retaining body1's coordinate frame.
%
%   see also ATTACHPREFIX, 

    % returned handle is body1 by default
    handle = body1;
    
    if nargin == 1;  return;  end
    
    % walk through provided bodies and add them to the return handle
    for i=1:(nargin-1)
        body = varargin{i};
        handle.bodies = [handle.bodies body.bodies];
        handle.labels = [handle.labels body.labels];
    end
    
end