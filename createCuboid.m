function handle = createCuboid(R0, t0, param, varargin)
    %
    % handle = createCuboid(R0, t0, param,...)
    %
    % R0 is 3 x 3 matrix for orientation of the cuboid
    % t0 is 3 x 1 vector for center of the cuboid
    % param is struct containing fields
    %       width [X axis]
    %       length [Y axis]
    %       height [Z axis]
    %       [opt] width2
    %       [opt] length2
    %       [opt] height2
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
    
    % verify parameters are correct
    if isfield(param,'width') && ...
            isfield(param,'length') && isfield(param,'height')
        w = param.width; l = param.length; h = param.height;
        if isfield(param,'width2') || isfield(param,'length2') ...
                || isfield(param,'height2')
            if isfield(param,'width2') && ~isfield(param,'length2') ...
                    && ~isfield(param,'height2')
                w2 = param.width2;
                l2 = l;
                h2 = h;
            elseif ~isfield(param,'width2') && isfield(param,'length2') ...
                    && ~isfield(param,'height2')
                w2 = w;
                l2 = param.length2;
                h2 = h;
            elseif ~isfield(param,'width2') && ~isfield(param,'length2') ...
                    && isfield(param,'height2')
                w2 = w;
                l2 = l;
                h2 = param.height2;
            else
                disp('Valid parameterizations are:');
                disp('    width, length, height - for cuboid');
                disp('    [optional parameters]:');
                disp('    width2 OR length2 OR height2 for trapezoidal faces');
                error('Invalid Parameterization for cuboid');
            end
        else
            w2 = w;
            l2 = l;
            h2 = h;
        end
    else
        disp('Valid parameterizations are:');
        disp('    width, length, height - for cuboid');
        disp('    [optional parameters]:');
        disp('    width2 OR length2 OR height2 for trapezoidal faces');
        error('Invalid Parameterization for cuboid');
    end
    
    
    
    right_face = 0.5*[w2 w w w2; -l l l2 -l2; h h -h -h];          % +X
    left_face = 0.5*[-w2 -w -w -w2; -l l l2 -l2; h2 h2 -h2 -h2];   % -X
    front_face = 0.5*[-w -w w w; l l2 l2 l; h2 -h2 -h h];          % +Y
    back_face = 0.5*[-w2 -w2 w2 w2; -l -l2 -l2 -l; h2 -h2 -h h];   % -Y
    top_face = 0.5*[-w -w2 w2 w; l -l -l l; h2 h2 h h];            % +Z
    bottom_face = 0.5*[-w -w2 w2 w; l2 -l2 -l2 l2; -h2 -h2 -h -h]; % -Z

    faces_rot = t0*ones(1,24) + ...
        R0*[left_face right_face front_face back_face ...
            top_face bottom_face];

    faces_X = reshape(faces_rot(1,:),[4 6]);
    faces_Y = reshape(faces_rot(2,:),[4 6]);
    faces_Z = reshape(faces_rot(3,:),[4 6]);
    
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