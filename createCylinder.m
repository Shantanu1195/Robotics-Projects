function handle = createCylinder(R0, t0, param, varargin)
    %
    % handle = createCylinder(R0, t0, param, ...)
    %
    % R0 is 3 x 3 matrix for orientation of the cylinder
    % t0 is 3 x 1 vector for center of the cylinder
    % param is struct containing fields
    %       radius 
    %       height 
    %       [opt] radius2 
    %       [opt] radiusX
    %       [opt] radiusY
    %       [opt] radiusX2
    %       [opt] radiusY2
    %
    % possible additional properties are:
    %       'FaceColor'  default: [1;1;1]
    %       'FaceAlpha'  default: 1
    %       'LineWidth'  default: 0.5
    %       'EdgeColor'  default: [0;0;0]
    %       'EdgeAlpha'  default: 1
    %
    % returns handle to drawing structure
    
    % Walk through varargin
    for i=1:2:(nargin-3)
        if strcmp(varargin{i},'FaceColor')
            fc = varargin{i+1};
        elseif strcmp(varargin{i},'FaceAlpha')
            fa = varargin{i+1};
        elseif strcmp(varargin{i},'LineWidth')
            lw = varargin{i+1};
        elseif strcmp(varargin{i},'EdgeColor')
            ec = varargin{i+1};
        elseif strcmp(varargin{i},'EdgeAlpha')
            ea = varargin{i+1};
        else
            error(['Property not recognized: ' varargin{i}]);
        end
    end
    % Set defaults if not already established
    if ~exist('fc','var'); fc = [1;1;1]; end
    if ~exist('fa','var'); fa = 1; end
    if ~exist('lw','var'); lw = 0.5; end
    if ~exist('ec','var'); ec = [0;0;0]; end
    if ~exist('ea','var'); ea = 1; end
    
    % Verify parameters are correct
    invalid = false;
    if isfield(param,'height')
        h = param.height;
    else
        invalid = true;
    end
    if isfield(param,'radius')
        rx = param.radius;
        ry = param.radius;
    elseif isfield(param,'radiusX') && isfield(param, 'radiusY')
        rx = param.radiusX;
        ry = param.radiusY;
    else
        invalid = true;
    end
    
    if isfield(param,'radius2')
        rx2 = param.radius2;
        ry2 = param.radius2;
    elseif isfield(param,'radiusX2') && isfield(param, 'radiusY2')
        rx2 = param.radiusX2;
        ry2 = param.radiusY2;
    else
        rx2 = rx;
        ry2 = ry;
    end
    
    if invalid
        disp('Valid parameterizations are:');
        disp('    height, radius - for circular cylinder');
        disp('    height, radius, radius2 -  for cone');
        disp('    height, radiusX, radiusY - for elliptic cylinder');
        disp('    height, radiusX, radiusY, radiusX2, radiusY2 - for elliptic cone');
        error('Invalid Parameterization for cylinder');
    end
    
    n = 20;
    X = cos(0:2*pi/n:2*pi);
    Y = sin(0:2*pi/n:2*pi);
    
    side_faces = zeros(3,4*n);

    for i=1:n
        side_faces(:,(i-1)*4+1:i*4) = [rx*X(i:i+1) rx2*X(i+1:-1:i); ...
                                       ry*Y(i:i+1) ry2*Y(i+1:-1:i); ...
                                       h/2   h/2  -h/2  -h/2];
    end

    side_faces_rot = t0*ones(1,4*n) + R0*side_faces;
    sides_X = reshape(side_faces_rot(1,:),[4 n]);
    sides_Y = reshape(side_faces_rot(2,:),[4 n]);
    sides_Z = reshape(side_faces_rot(3,:),[4 n]);

    cap_faces = [rx*X(1:n) rx2*X(1:n); ...
                 ry*Y(1:n) ry2*Y(1:n); ...
                 h/2*ones(1,n) -h/2*ones(1,n)];
    cap_faces_rot = t0*ones(1,2*n) + R0*cap_faces;
    caps_X = reshape(cap_faces_rot(1,:),[n 2]);
    caps_Y = reshape(cap_faces_rot(2,:),[n 2]);
    caps_Z = reshape(cap_faces_rot(3,:),[n 2]);

    
    handle.bodies(1) = patch(sides_X,sides_Y,sides_Z,1, ...
                                    'FaceColor',fc, ...
                                    'FaceAlpha',fa, ...
                                    'LineWidth',lw, ...
                                    'EdgeColor',ec, ...
                                    'EdgeAlpha',ea);
    handle.bodies(2) = patch(caps_X,caps_Y,caps_Z,1, ...
                                    'FaceColor',fc, ...
                                    'FaceAlpha',fa, ...
                                    'LineWidth',lw, ...
                                    'EdgeColor',ec, ...
                                    'EdgeAlpha',ea);
    handle.R = R0;
    handle.t = t0;
    handle.labels = {'sides','caps'};
end