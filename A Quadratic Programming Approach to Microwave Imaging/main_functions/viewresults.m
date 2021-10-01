function [fig] = viewresults (varargin)
%VIEWRESULTS Plot images.
%   fig = viewresults(singlemap,model,filename) plots the image of single
%   contrast map saving it as a .fig file. If the relative permittivity or
%   conductivity is background everywhere, then only the non-background
%   property is plotted.
%
%   fig = viewresults(originalmap,originalmodel,recoveredmap,recoveredmodel,
%   filename) plots the images of the original contrast map and the
%   recovered one. If one of the properties is background everwhere, then
%   only the other is plotted.
%
%   Implemented by Andre Costa Batista (Universidade Federal de Minas
%   Gerais, Brazil). Contact: andre-costa@ufmg.br

    if nargin == 3
        
        fig = figure;
        epsilon_r = varargin{1}.epsilon_r;
        sigma = varargin{1}.sigma;
        model = varargin{2};
        filename = varargin{3};
        
        if all(sigma(:)==model.sigb)
            imagesc(model.x,model.y,epsilon_r')
            set(gca,'YDir','normal')
            c = colorbar;
            caxis([min(epsilon_r(:)), max(epsilon_r(:))]);
            c.Label.String = '\epsilon_r';
            title('Relative Permittivity')
            xlabel('x [m]')
            ylabel('y [m]')
            
        elseif all(epsilon_r(:)==model.epsrb)
            imagesc(model.x,model.y,sigma')
            set(gca,'YDir','normal')
            c = colorbar;
            caxis([min(sigma(:)), max(sigma(:))]);
            c.Label.String = '\sigma [S/m]';
            title('Conductivity')
            xlabel('x [m]')
            ylabel('y [m]')
            
        else
            subplot(1,2,1)
            imagesc(model.x,model.y,epsilon_r')
            set(gca,'YDir','normal')
            c = colorbar;
            caxis([min(epsilon_r(:)), max(epsilon_r(:))]);
            c.Label.String = '\epsilon_r';
            title('Relative Permittivity')
            xlabel('x [m]')
            ylabel('y [m]')
            
            subplot(1,2,2)
            imagesc(model.x,model.y,sigma')
            set(gca,'YDir','normal')
            c = colorbar;
            caxis([min(sigma(:)), max(sigma(:))]);
            c.Label.String = '\sigma [S/m]';
            title('Conductivity')
            xlabel('x [m]')
            ylabel('y [m]')
        end
        
        savefig(filename);
        
    elseif nargin == 5
        
        fig = figure;
        ori_epsr = varargin{1}.epsilon_r;
        ori_sig = varargin{1}.sigma;
        ori_model = varargin{2};
        rec_epsr = varargin{3}.epsilon_r;
        rec_sig = varargin{3}.sigma;
        rec_model = varargin{4};
        filename = varargin{5};
        
        if all(ori_sig(:)==ori_model.sigb) ...
                && all(rec_sig(:)==rec_model.sigb)
            subplot(1,2,1)
            imagesc(rec_model.x,rec_model.y,rec_epsr')
            set(gca,'YDir','normal')
            c = colorbar;
            if min(ori_epsr(:)) ~= max(ori_epsr(:))
                caxis([min(ori_epsr(:)), max(ori_epsr(:))]);
            end
            c.Label.String = '\epsilon_r';
            title('Recovered Relative Permittivity')
            xlabel('x [m]')
            ylabel('y [m]')
            subplot(1,2,2)
            imagesc(ori_model.x,ori_model.y,ori_epsr')
            set(gca,'YDir','normal')
            c = colorbar;
            if min(ori_epsr(:)) ~= max(ori_epsr(:))
                caxis([min(ori_epsr(:)), max(ori_epsr(:))]);
            end
            c.Label.String = '\epsilon_r';
            title('Original Relative Permittivity')
            xlabel('x [m]')
            ylabel('y [m]')
        
        elseif all(ori_epsr(:)==ori_model.epsrb) ...
                && all(rec_epsr(:)==rec_model.epsrb)
            subplot(1,2,1)
            imagesc(rec_model.x,rec_model.y,rec_sig')
            set(gca,'YDir','normal')
            c = colorbar;
            if min(ori_sig(:)) ~= max(ori_sig(:))
                caxis([min(ori_sig(:)), max(ori_sig(:))]);
            end
            c.Label.String = '\sigma [S/m]';
            title('Recovered Conductivity')
            xlabel('x [m]')
            ylabel('y [m]')
            subplot(1,2,2)
            imagesc(ori_model.x,ori_model.y,ori_sig')
            set(gca,'YDir','normal')
            c = colorbar;
            if min(ori_sig(:)) ~= max(ori_sig(:))
                caxis([min(ori_sig(:)), max(ori_sig(:))]);
            end
            c.Label.String = '\sigma [S/m]';
            title('Original Conductivity')
            xlabel('x [m]')
            ylabel('y [m]')
            
        else
            subplot(2,2,1)
            imagesc(rec_model.x,rec_model.y,rec_epsr')
            set(gca,'YDir','normal')
            c = colorbar;
            if min(ori_epsr(:)) ~= max(ori_epsr(:))
                caxis([min(ori_epsr(:)), max(ori_epsr(:))]);
            end
            c.Label.String = '\epsilon_r';
            title('Recovered Relative Permittivity')
            xlabel('x [m]')
            ylabel('y [m]')
            
            subplot(2,2,2)
            imagesc(rec_model.x,rec_model.y,rec_sig')
            set(gca,'YDir','normal')
            c = colorbar;
            if min(ori_sig(:)) ~= max(ori_sig(:))
                caxis([min(ori_sig(:)), max(ori_sig(:))]);
            end
            c.Label.String = '\sigma [S/m]';
            title('Recovered Conductivity')
            xlabel('x [m]')
            ylabel('y [m]')
            
            
            subplot(2,2,3)
            imagesc(ori_model.x,ori_model.y,ori_epsr')
            set(gca,'YDir','normal')
            c = colorbar;
            if min(ori_epsr(:)) ~= max(ori_epsr(:))
                caxis([min(ori_epsr(:)), max(ori_epsr(:))]);
            end
            c.Label.String = '\epsilon_r';
            title('Original Relative Permittivity')
            xlabel('x [m]')
            ylabel('y [m]')
            
            subplot(2,2,4)
            imagesc(ori_model.x,ori_model.y,ori_sig')
            set(gca,'YDir','normal')
            c = colorbar;
            if min(ori_sig(:)) ~= max(ori_sig(:))
                caxis([min(ori_sig(:)), max(ori_sig(:))]);
            end
            c.Label.String = '\sigma [S/m]';
            title('Original Conductivity')
            xlabel('x [m]')
            ylabel('y [m]')
        end
        
        savefig(filename);
        
    else
        disp('ERROR VIEWRESULTS: invalid number of inputs!')
        return
    end
end 