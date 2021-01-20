function handle = create3DFrame(R0, t0, param, varargin)
    %
    % handle = create3DFrame(R0, t0, param, ...)
    %
    % R0 is the 3x3 orientation matrix representation of the 3D frame
    % t0 is the 3x1 center of the 3D frame in the world frame
    % param is struct containing fields
    %       scale (length of each axis)
    %       width (width of each cylinder)
    %
    % Additional Parameters include:
    %       'FaceAlpha': default 1
    %       'EdgeAlpha': default 0
    %
    % depends on the following drawing package files:
    %       createCylinder.m
    %       attachPrefix.m
    %
    % returns handle to drawing structure
    
    % Walk through varargin
    for i=1:2:(nargin-3)
        if strcmp(varargin{i},'FaceAlpha')
            fa = varargin{i+1};
        elseif strcmp(varargin{i},'EdgeAlpha')
            ea = varargin{i+1};
        else
            error(['Parameter not recognized: ' varargin{i}]);
        end
    end
    % Set default values
    if ~exist('fa','var'); fa = 1; end
    if ~exist('ea','var'); ea = 0; end
    
    axis_param.radius = param.width/2;
    axis_param.height = param.scale;
    h = param.scale;
    
    Rx = R0*rot([0;1;0],pi/2);
    tx = t0 + Rx*[0;0;h/2];
    Ry = R0*rot([1;0;0],-pi/2);
    ty = t0 + Ry*[0;0;h/2];
    Rz = R0;
    tz = t0 + Rz*[0;0;h/2];
    
    h_x = createCylinder(Rx,tx,axis_param, ...
                                'FaceColor',[1;0;0], ...
                                'FaceAlpha', fa, ...
                                'EdgeAlpha', ea);
	h_y = createCylinder(Ry,ty,axis_param, ...
                                'FaceColor',[0;1;0], ...
                                'FaceAlpha', fa, ...
                                'EdgeAlpha', ea);
	h_z = createCylinder(Rz,tz,axis_param, ...
                                'FaceColor',[0;0;1], ...
                                'FaceAlpha', fa, ...
                                'EdgeAlpha', ea);
    
    h_y.labels = attachPrefix('Y_',h_x.labels);
    h_x.labels = attachPrefix('X_',h_x.labels);
    h_z.labels = attachPrefix('Z_',h_x.labels);
    handle = combineRigidBodies(h_x,h_y,h_z);
    handle.R = R0;
    handle.t = t0;
end
    