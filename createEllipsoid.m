function handle = createEllipsoid(R0,t0,param,varargin)
    %
    % handle = createEllipsoid(R0, t0, param,...)
    %
    % R0 is orientation of the ellipsoid relative to body frame
    % t0 is center of the ellipsoid relative to body frame
    % param is struct containing fields
    %       radius (draws a sphere)
    %       [opt] radiusX
    %       [opt] radiusY
    %       [opt] radiusZ
    % possible additional properties are:
    %       'FaceColor'  default: [1;1;1]
    %       'FaceAlpha'  default: 1
    %       'LineWidth'  default: 0.5
    %       'EdgeColor'  default: [0;0;0]
    %       'EdgeAlpha'  default: 1
    %
    % returns handle to ellipsoid drawing structure
    
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
            error(['Parameter not recognized: ' varargin{i}]);
        end
    end
    % Set defaults if not already established
    if ~exist('fc','var'); fc = [1;1;1];end
    if ~exist('fa','var'); fa = 1; end
    if ~exist('lw','var'); lw = 0.5; end
    if ~exist('ec','var'); ec = [0;0;0]; end
    if ~exist('ea','var'); ea = 1; end
    
    % Verify parameters are correct
    if isfield(param,'radius')
        rx = param.radius;
        ry = param.radius;
        rz = param.radius;
    elseif isfield(param,'radiusX') && isfield(param, 'radiusY') && ...
                isfield(param,'radiusZ')
        rx = param.radiusX;
        ry = param.radiusY;
        rz = param.radiusZ;
    else
        disp('Parameterization needs either:');
        disp('    radius - for sphere');
        disp('    radiusX');
        disp('    radiusY');
        disp('    radiusZ - for ellipsoid');
        handle = [];
        return;
    end
    
    N = 5;
    theta = 0:pi/N:2*pi;
    phi = 0:pi/N:pi;
    
    c_th = cos(theta);
    s_th = sin(theta);
    c_phi = cos(phi);
    s_phi = sin(phi);
    
    x1 = rx*c_th(1:2*N);
    x2 = rx*c_th(2:2*N+1);
    y1 = ry*s_th(1:2*N);
    y2 = ry*s_th(2:2*N+1);
    
    faces = zeros(3,8*N*N);
    for i=1:N
        z1 = rz*c_phi(i)*ones(1,2*N);
        z2 = rz*c_phi(i+1)*ones(1,2*N);
        faces(:,(i-1)*8*N+1:4:i*8*N) = [x1*s_phi(i);y1*s_phi(i);z1];
        faces(:,(i-1)*8*N+2:4:i*8*N) = [x2*s_phi(i);y2*s_phi(i);z1];
        faces(:,(i-1)*8*N+3:4:i*8*N) = [x2*s_phi(i+1);y2*s_phi(i+1);z2];
        faces(:,(i-1)*8*N+4:4:i*8*N) = [x1*s_phi(i+1);y1*s_phi(i+1);z2];
    end
    
    faces = t0*ones(1,8*N*N) + R0*faces;
    faces_X = reshape(faces(1,:),[4 2*N*N]);
    faces_Y = reshape(faces(2,:),[4 2*N*N]);
    faces_Z = reshape(faces(3,:),[4 2*N*N]);
    handle.bodies(1) = patch(faces_X,faces_Y,faces_Z,1, ...
                                    'FaceColor',fc, ...
                                    'FaceAlpha',fa, ...
                                    'LineWidth',lw, ...
                                    'EdgeColor',ec, ...
                                    'EdgeAlpha',ea);
    handle.R = R0;
    handle.t = t0;
    handle.labels = {'sides'};
end
    