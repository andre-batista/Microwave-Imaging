function [contrastmap] = getmap(varargin)
%GETMAP Get a contrast map.
%   contrastmap = GETMAP(mapname,model,epsilon_robj,sigma_obj,...) creates
%   a contrast map of geometric shapes. The image is represented as two
%   matrices, for relative permittivity ('epsilon' field) and conductivity 
%   [S/m] ('sigma' field). The size of the image and the background 
%   properties are in the model struct. The contrast of the objects are
%   represented by epsilon_robj and sigma_obj, where they may be scalar or
%   array, according to each case. There are seven map model available with
%   the following map names:
%
%   - contrastmap = GETMAP('triangle',model,epsilon_robj,sigma_obj,l) draws
%     a single equilateral triangle at the center of the image. The input l
%     is the lateral size.
%
%   - contrastmap = GETMAP('star',model,epsilon_robj,sigma_obj,l) draws a
%     hexagram at the center, where l is the lateral size of the two
%     equilateral triangles.
% 
%   - contrastmap = GETMAP('ring',model,epsilon_robj,sigma_obj,ra,rb,rc,
%     del) draws a ring with a circle in at the center where: ra is the
%     outer radius of the ring; rb is the inner radius; rc is the circle
%     radius; and del is the shift from center in both axis.
%     
%   - contrastmap = GETMAP('ellipses',model,epsilon_robj,sigma_obj,la,lb,
%     del) draws two ellipses with semi-minor axis la, semi-major axis lb
%     and a distance del between them.
%
%   - contrastmap = GETMAP('2circles',model,epsilon_robj,sigma_obj,ra,del) 
%     draws two circles with the same radius ra and displacement del 
%     between them. Here the object contrast must be an array with three
%     elements: the contrast of the two circles and the intersection 
%     between them.
%
%   - contrastmap = GETMAP('3objects',model,epsilon_robj,sigma_obj,ra,dela,
%     lb,delb,lc,delc) draws a triangle, a lb,delb,lc,delc) circle and a
%     square. The circle radius is ra and dela is an array with the 
%     displacament of the circle from the center in each axis. The square
%     side is lb and delb is an array with the displacement from the center
%     in each axis. The equilateral triangle size is lc and delc is an 
%     array with the displacement from the center in each axis. The object
%     contrast must be an array with the following order: circle, square 
%     and triangle.
%
%   - contrastmap = GETMAP('filledring',model,epsilon_robj,sigma_obj,ra,rb)
%     draws a ring at the center of the image with outer radius ra and 
%     inner radius rb. The object contrast must be an array in which the 
%     first elements is the ring contrast and the second is the contrast of 
%     the area surrounded by the ring.
%
%   All space variables are in meters.
%
%   Example:
%       
%       contrastmap = getmap('3objects',100,100,.1,.1,1.,.0,[2.,3.,4.],[.1,.2,.3],.6,[2.,2.],3.,[-2.,2],2.,[-2.,-2.]);
%       subplot(1,2,1)
%       imagesc(contrastmap.epsilon_r)
%       title('Relative Permittivity')
%       colorbar
%       subplot(1,2,2)
%       i = magesc(contrastmap.sigma)
%       title('Conductivity')
%       colorbar
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    if nargin < 5
        disp('ERROR GETMAP: less inputs than required!')
        return
    end
    
    model = varargin{2};
        
    switch varargin{1}
                
        case 'triangle'
            
            if nargin ~= 5
                disp('ERROR GETMAP: invalid number of inputs for triangle map')
                return
            end
            
            [epsilon_r,sigma] = build_triangle(model.I,model.J,model.dx,...
                model.dy,model.epsrb,model.sigb,varargin{3},varargin{4},...
                varargin{5});
        
        case 'star'
        
            if nargin ~= 5
                disp('ERROR GETMAP: invalid number of inputs for star map')
                return
            end
            
            [epsilon_r,sigma] = build_star(model.I,model.J,model.dx,...
                model.dy,model.epsrb,model.sigb,varargin{3},varargin{4},...
                varargin{5});
        
        case 'ring'
        
            if nargin ~= 8
                disp('ERROR GETMAP: invalid number of inputs for ring map')
                return
            end
            
            [epsilon_r,sigma] = build_ring(model.I,model.J,model.dx,...
                model.dy,model.epsrb,model.sigb,varargin{3},varargin{4},...
                varargin{5},varargin{6},varargin{7},varargin{8});
        
        case 'ellipses'
        
            if nargin ~= 7
                disp('ERROR GETMAP: invalid number of inputs for ellipses map')
                return
            end
            
            [epsilon_r,sigma] = build_ellipses(model.I,model.J,model.dx,...
                model.dy,model.epsrb,model.sigb,varargin{3},varargin{4},...
                varargin{5},varargin{6},varargin{7});
       
        case '2circles'
        
            if nargin ~= 6
                disp('ERROR GETMAP: invalid number of inputs for 2circles map')
                return
            end
        
            [epsilon_r,sigma] = build_2circles(model.I,model.J,model.dx,...
                model.dy,model.epsrb,model.sigb,varargin{3},varargin{4},...
                varargin{5},varargin{6});
        
        case '3objects'
            
            if nargin ~= 10
                disp('ERROR GETMAP: invalid number of inputs for 2circles map')
                return
            end            
        
            [epsilon_r,sigma] = build_3objects(model.I,model.J,model.dx,...
                model.dy,model.epsrb,model.sigb,varargin{3},varargin{4},...
                varargin{5},varargin{6},varargin{7},varargin{8},...
                varargin{9},varargin{10});
            
        case 'filledring'
            
            if nargin ~= 6
                disp('ERROR GETMAP: invalid number of inputs for filledring map')
                return
            end
            
            [epsilon_r,sigma] = build_filledring(model.I,model.J,model.dx,...
                model.dy,model.epsrb,model.sigb,varargin{3},varargin{4},...
                varargin{5},varargin{6});
            
        otherwise
            disp('ERROR GETMAP: incorrect model name!')
            return
    end

    contrastmap.epsilon_r = epsilon_r;
    contrastmap.sigma = sigma;
    
end

function [epsilon_r,sigma] = build_triangle(I,J,dx,dy,epsilon_rb,sigma_b,...
                                        epsilon_robj,sigma_obj,l)

    epsilon_r = epsilon_rb * ones(I,J);
    sigma = sigma_b * ones(I,J);
    x = (1:I) * dx;
    y = (1:J) * dy;
    Lx = I*dx; Ly = J*dy;

    for i = 1:I
        for j = 1:J

            if x(i) >= Lx/2-l/2 && x(i) <= Lx/2+l/2
                a = x(i)-.5*(Lx-l);
                FLAG = false;
                if y(j) < Ly/2
                    b = y(j)-.5*(Ly-l);
                    v = -.5*a+l/2;
                    if b>=v
                        FLAG = true;
                    end
                else
                    b = y(j)-Lx/2;
                    v = .5*a;
                    if b<=v
                        FLAG = true;
                    end
                end
                if FLAG == true
                    epsilon_r(i,j) = epsilon_robj;
                    sigma(i,j) = sigma_obj;
                end
            end
        end
    end
end

function [epsilon_r,sigma] = build_star(I,J,dx,dy,epsilon_rb,sigma_b, ...
                                        epsilon_robj,sigma_obj,l)

    epsilon_r = epsilon_rb * ones(I,J);
    sigma = sigma_b * ones(I,J);
    x = (1:I) * dx;
    y = (1:J) * dy;
    Lx = I*dx; Ly = J*dy;
    xc = l/6;

    for i = 1:I
        for j = 1:J

            if x(i)+xc >= Lx/2-l/2 && x(i)+xc <= Lx/2+l/2
                a = x(i)+xc-.5*(Lx-l);
                FLAG = false;
                if y(j) < Ly/2
                    b = y(j)-.5*(Ly-l);
                    v = -.5*a+l/2;
                    if b>=v
                        FLAG = true;
                    end
                else
                    b = y(j)-Lx/2;
                    v = .5*a;
                    if b<=v
                        FLAG = true;
                    end
                end
                if FLAG == true
                    epsilon_r(i,j) = epsilon_robj;
                    sigma(i,j) = sigma_obj;
                end
            end
        end
    end


    for i = 1:I
        for j = 1:J

            if x(i)-xc >= Lx/2-l/2 && x(i)-xc <= Lx/2+l/2
                a = x(i)-xc-.5*(Lx-l);
                FLAG = false;
                if y(j) < Ly/2
                    b = y(j)-.5*(Ly-l);
                    v = .5*a;
                    if b>=v
                        FLAG = true;
                    end
                else
                    b = y(j)-Lx/2;
                    v = -.5*a+l/2;
                    if b<=v
                        FLAG = true;
                    end
                end
                if FLAG == true
                    epsilon_r(i,j) = epsilon_robj;
                    sigma(i,j) = sigma_obj;
                end
            end
        end
    end
end

function [epsilon_r,sigma] = build_ring(I,J,dx,dy,epsilon_rb,sigma_b,...
                                        epsilon_robj,sigma_obj,ra,rb,rc,del)

    epsilon_r = epsilon_rb * ones(I,J);
    sigma = sigma_b * ones(I,J);
    x = (1:I) * dx;
    y = (1:J) * dy;
    xc = I*dx/2+del;
    yc = J*dy/2+del;

    for i = 1:I
        for j = 1:J

            r = sqrt((x(i)-xc)^2+(y(j)-yc)^2);
            if r <= ra && r >= rb
                epsilon_r(i,j) = epsilon_robj;
                sigma(i,j) = sigma_obj;
            elseif r <= rc
                epsilon_r(i,j) = epsilon_robj;
                sigma(i,j) = sigma_obj;
            end

        end
    end
end

function [epsilon_r,sigma] = build_ellipses(I,J,dx,dy,epsilon_rb,sigma_b,...
                                        epsilon_robj,sigma_obj,la,lb,del)

    epsilon_r = epsilon_rb * ones(I,J);
    sigma = sigma_b * ones(I,J);
    x = (1:I) * dx;
    y = (1:J) * dy;
    xc = I*dx/2;
    yc = J*dy/2;

    for i = 1:I
        for j = 1:J

            if (x(i)-(xc-del))^2/la^2 + (y(j)-yc)^2/lb^2 <= 1
                epsilon_r(i,j) = epsilon_robj;
                sigma(i,j) = sigma_obj;
            elseif (x(i)-(xc+del))^2/la^2 + (y(j)-yc)^2/lb^2 <= 1
                epsilon_r(i,j) = epsilon_robj;
                sigma(i,j) = sigma_obj;
            end
        end
    end
end

function [epsilon_r,sigma] = build_2circles(I,J,dx,dy,epsilon_rb,sigma_b,...
                                        epsilon_robj,sigma_obj,ra,del)

    epsilon_r = epsilon_rb * ones(I,J);
    sigma = sigma_b * ones(I,J);
    x = (1:I) * dx;
    y = (1:J) * dy;
    
    xc1 = I*dx/2+del;
    yc1 = J*dy/2+del;
    
    xc2 = I*dx/2-del;
    yc2 = J*dy/2+del;

    for i = 1:I
        for j = 1:J

            r1 = sqrt((x(i)-xc1)^2+(y(j)-yc1)^2);
            r2 = sqrt((x(i)-xc2)^2+(y(j)-yc2)^2);
            
            if r1 <= ra && r2 <= ra
                epsilon_r(i,j) = epsilon_robj(3);
                sigma(i,j) = sigma_obj(3);
            elseif r1 <= ra
                epsilon_r(i,j) = epsilon_robj(1);
                sigma(i,j) = sigma_obj(1);
            elseif r2 <= ra
                epsilon_r(i,j) = epsilon_robj(2);
                sigma(i,j) = sigma_obj(2);
            end

        end
    end
end

function [epsilon_r,sigma] = build_3objects(I,J,dx,dy,epsilon_rb,sigma_b,...
                                        epsilon_robj,sigma_obj,ra,dela,lb,...
                                        delb,lc,delc)

    epsilon_r = epsilon_rb * ones(I,J);
    sigma = sigma_b * ones(I,J);
    x = (1:I) * dx;
    y = (1:J) * dy;
    Lx = I*dx; Ly = J*dy;
    
    xca = I*dx/2+dela(1);
    yca = J*dy/2+dela(2);
    
    xcb = I*dx/2+delb(1);
    ycb = J*dy/2+delb(2);
    
    xcc = I*dx/2+delc(1);
    ycc = J*dy/2+delc(2);

    for i = 1:I
        for j = 1:J

            r = sqrt((x(i)-xca)^2+(y(j)-yca)^2);
           
            if r <= ra
                epsilon_r(i,j) = epsilon_robj(1);
                sigma(i,j) = sigma_obj(1);
            elseif x(i) >= xcb-lb/2 && x(i) <= xcb+lb/2 && y(j) >= ycb-lb/2 && y(j) <= ycb+lb/2
                epsilon_r(i,j) = epsilon_robj(2);
                sigma(i,j) = sigma_obj(2);
            elseif x(i) >= xcc-lc/2 && x(i) <= xcc+lc/2 && y(j) >= ycc-lc/2 && y(j) <= ycc+lc/2
                FLAG = false;
                v = lc/2*sqrt(3);
                if y(j) >= ycc
                    a = lc/2/v;
                    b = ycc-lc/2/v*(xcc-v/2);
                    if y(j) <= a*x(i)+b
                        FLAG = true;
                    end
                else
                    yp = ycc + ycc-y(j);
                    a = lc/2/v;
                    b = ycc-lc/2/v*(xcc-v/2);
                    if yp <= a*x(i)+b
                        FLAG = true;
                    end
                end
                if FLAG == true
                    epsilon_r(i,j) = epsilon_robj(3);
                    sigma(i,j) = sigma_obj(3);
                end
            end

        end
    end
end

function [epsilon_r,sigma] = build_filledring(I,J,dx,dy,epsilon_rb,sigma_b,...
                                        epsilon_robj,sigma_obj,ra,rb)

    epsilon_r = epsilon_rb * ones(I,J);
    sigma = sigma_b * ones(I,J);
    x = (1:I) * dx;
    y = (1:J) * dy;
    
    xc = I*dx/2;
    yc = J*dy/2;

    for i = 1:I
        for j = 1:J

            r = sqrt((x(i)-xc)^2+(y(j)-yc)^2);
           
            if r <= ra && r >= rb
                epsilon_r(i,j) = epsilon_robj(1);
                sigma(i,j) = sigma_obj(1);
            elseif r <= rb
                epsilon_r(i,j) = epsilon_robj(2);
                sigma(i,j) = sigma_obj(2);
            end

        end
    end
end