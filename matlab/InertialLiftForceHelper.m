classdef InertialLiftForceHelper
    % InertialLiftForceHelper This class helps with parsing/utilising computed
    %   inertial lift force and migration data for a neutrally buoyant 
    %   spherical particle suspended in flow through curved ducts having a  
    %   square or rectangular (with aspect ratio 2:1 or 4:1) cross-section.
    %   
    %   The class will interpolate the raw data for a desired bend radius
    %   of the duct. Interpolation of particle size is not currently 
    %   supported. As such one must pick a particle size matching the 
    %   provided data.
    %       
    %   Note: use of this helper class requires the files 
    %   'curved_duct_lift_data_mat/square_lift_data.mat'  
    %   'curved_duct_lift_data_mat/rect_2x1_lift_data.mat' 
    %   'curved_duct_lift_data_mat/rect_4x1_lift_data.mat' 
    %   to be located relative to the working directory.
    %
    %   Note: this class will not work with GNU Octave because it 
    %   utilises griddedInterpolant which currently raises a warning
    %   as it has not yet been implemented.
    %   
    %   Please ensure you cite our JFM paper if you use this code/data 
    %   (https://doi.org/10.1017/jfm.2019.323). This code is provided under 
    %   an MIT license (see https://opensource.org/licenses/MIT). However, 
    %   I would appreciate it if you contact me and to let me know if you 
    %   use this code/data. Please also don't hesitate to contact me if you 
    %   have any questions/queries.
    %   
    %   Brendan Harding, 2019.

    properties (Access = private)
        a
        R 
        cs
        raw_data
        aR_lookup
        ar 
        R_values
        rs
        zs
        a_data
        aR_data
        kappa 
        Cr_RBS
        Cz_RBS
        Lr_RBS
        Lz_RBS
        Sr_RBS
        Sz_RBS
        Up_RBS
        Wr_RBS
        Wz_RBS
        Fr_RBS
        Fz_RBS
    end
    properties
        % I do not provide any public properties
    end
    methods (Access = private)
        function obj = load_cs_data(obj)
            % Load data associate with the desired cross-section
            if strcmp(obj.cs,'square')
                % If the  data file is not found, a FileNotFoundError exception will be raised
                mat_file = load('curved_duct_lift_data_mat/square_lift_data.mat');
                obj.aR_lookup = mat_file.aR;
                obj.raw_data = mat_file.data;
                obj.ar = 1.0;
            elseif strcmp(obj.cs,'rect_2x1')
                % If the data file is not found, a FileNotFoundError exception will be raised
                mat_file = load('curved_duct_lift_data_mat/rect_2x1_lift_data.mat');
                obj.aR_lookup = mat_file.aR;
                obj.raw_data = mat_file.data;
                obj.ar = 2.0;
            elseif strcmp(obj.cs,'rect_4x1')
                % If the data file is not found, a FileNotFoundError exception will be raised
                mat_file = load('curved_duct_lift_data_mat/rect_4x1_lift_data.mat');
                obj.aR_lookup = mat_file.aR;
                obj.raw_data = mat_file.data;
                obj.ar = 4.0;
            else
                error("The requestied cross-section is not one of 'square' or 'rect_2x1' or 'rect_4x1'");
            end
        end
        function obj = update_a_data(obj)
            % Pre-process the raw data associated with the desired a
            ais = obj.aR_lookup(:,1)==obj.a;
            if sum(ais)==0
                error(strcat('The desired particle size is not present in the raw data.\n',...
                             'Check what is available with the method get_available_a()'));
            end
            obj.R_values = obj.aR_lookup(ais,2);
            a_data_full = obj.raw_data(ais,:,:,:);
            % strip NaN data around the edges
            r = a_data_full(1,1,:,21);
            z = a_data_full(1,2,21,:);
            rf = isfinite(r);
            zf = isfinite(z);
            obj.rs = r(rf);
            obj.zs = z(zf);
            obj.a_data = a_data_full(:,:,rf,zf);
        end
        function obj = update_R_data(obj)
            % Pre-process the raw data associated with the desired R
            if obj.R<obj.ar
                warning(strcat('The requested bend radius R is non-physical. Choose R>W \n',...
                               '\t(and ideally R>>W) where R is interpreted to be relative \n',...
                               '\tto half the duct height, i.e. read R as 2R/H'));
            end
            if obj.R<20.0 || obj.R>1280.0
                warning(strcat('Warning: Results may be innacurate for the given bend radius, \n',...
                               '\t(20<=R<=1280 is best, where R should be read as 2R/H)'));
            end
            if any(obj.R==obj.R_values)
                % use the respective raw data as is
                obj.aR_data = squeeze(obj.a_data(obj.R==obj.R_values,:,:,:));
            else
                % construct an interpolant
                eps = flip(1.0./obj.R_values(:));
                aR_data_full = squeeze(obj.a_data(1,:,:,:));
                F_array = flip(squeeze(obj.a_data(:,3:1:11,:,:)),1);
                aR_data_full(3:1:11,:,:) = reshape(interp1(eps,F_array,1.0/obj.R,'cubic','extrap'),...
                                                   size(aR_data_full(3:1:11,:,:)));
                obj.aR_data = aR_data_full;
            end
        end
        function obj = setup_interpolants(obj)
            % Constructs interpolants of the appropriate fields
            % griddedInterpolant is not implemented in octave, making this difficult to test
            obj.kappa = 4.0/(obj.R*obj.a.^3);
            RS = squeeze(obj.aR_data(1,:,:));
            ZS = squeeze(obj.aR_data(2,:,:));
            method = 'cubic'; % could try 'cubic', 'pchip' or 'spline' for the method
            extrap = 'none'; % using 'none' means NaN is returned outside the data range
            obj.Cr_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(3,:,:)), method,extrap);
            obj.Cz_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(4,:,:)), method,extrap);
            obj.Lr_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(5,:,:)), method,extrap);
            obj.Lz_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(6,:,:)), method,extrap);
            obj.Sr_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(7,:,:)), method,extrap);
            obj.Sz_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(8,:,:)), method,extrap);
            obj.Up_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(9,:,:)), method,extrap);
            obj.Wr_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(10,:,:)),method,extrap);
            obj.Wz_RBS = griddedInterpolant(RS,ZS,squeeze(obj.aR_data(11,:,:)),method,extrap);
            Fr = squeeze(obj.aR_data(5,:,:)+obj.kappa*obj.aR_data(7,:,:));
            Fz = squeeze(obj.aR_data(6,:,:)+obj.kappa*obj.aR_data(8,:,:));
            obj.Fr_RBS = griddedInterpolant(RS,ZS,Fr,method,extrap);
            obj.Fz_RBS = griddedInterpolant(RS,ZS,Fz,method,extrap);
        end
        function plot_aR_component(obj,ind)
            % Useful for testing/debugging purposes
            RS = squeeze(obj.aR_data(1,:,:));
            ZS = squeeze(obj.aR_data(2,:,:));
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0,0,1+3*obj.ar,3];
            hold on
            contourf(RS,ZS,squeeze(obj.aR_data(ind,:,:)),17);
            colorbar();
            nhH = 1.0/obj.a;
            nhW = obj.ar/obj.a;
            plot([-nhW+1,nhW-1,nhW-1,-nhW+1,-nhW+1],...
                 [-nhH+1,-nhH+1,nhH-1,nhH-1,-nhH+1],'r--');
            set(gca,'DataAspectRatio',[1,1,1],...
                    'PlotBoxAspectRatio',[1,1,1],...
                    'XLim',[-nhW,nhW],...
                    'YLim',[-nhH,nhH]);
            hold off
        end
    end
    methods
        function obj = InertialLiftForceHelper(a_,R_,cs_)
            % Initialise the helper class for a given particle radius (a) and bend radius (R).
            % Both a and R should be read as being relative to half of the duct height, 
            % i.e. read a as 2a/H and R as 2R/H (H being the duct height).
            % The cs argument specifies cross-section shape: 'square' for a square 
            % cross-section, or 'rect_2x1' for a rectangular duct (having aspect ratio 2), 
            % or 'rect_4x1' for a rectangular duct (having aspect ratio 4).
            % (Data for additional cross-sections may be added in the future)
            obj.a = a_;
            obj.R = R_;
            obj.cs = cs_;
            obj = obj.load_cs_data();
            obj = obj.update_a_data();
            obj = obj.update_R_data();
            obj = obj.setup_interpolants();
        end
        function cs_ = get_cross_section(obj)
            % Get the current cross-section
            cs_ = obj.cs;
        end
        function bounds_ = get_bounds(obj)
            % Get the bounds of the current cross-section
            % (given in the form [r_min,r_max,z_min,z_max])
            bounds_ = [-obj.ar/obj.a,obj.ar/obj.a,-1.0/obj.a,1.0/obj.a];
        end
        function bounds_ = get_particle_bounds(obj)
            % Get the bounds for the particle centre in the current 
            % cross-section (given in the form [r_min,r_max,z_min,z_max])
            bounds_ = [-obj.ar/obj.a+1.0,obj.ar/obj.a-1.0,-1.0/obj.a+1.0,1.0/obj.a-1.0];
        end
        function obj = set_cross_section(obj,cs_)
            % Change the cross-section (and update the interpolants)
            obj.cs = cs_;
            obj = obj.load_cs_data();
            obj = obj.update_a_data();
            obj = obj.update_R_data();
            obj = obj.setup_interpolants();
        end
        function alist = get_available_a(obj)
            % Get a list of available particle radii
            alist = unique(obj.raw_data.aR(:,1));
        end
        function a_ = get_a(obj)
            % Get the current particle radius
            a_ = obj.a;
        end
        function obj = set_a(obj,a_)
            % Change the particle radius (and update the interpolants)
            obj.a = a_;
            obj = obj.update_a_data();
            obj = obj.update_R_data();
            obj = obj.setup_interpolants();
        end
        function R_ = get_R(obj)
            % Get the current bend radius
            R_ = obj.R;
        end
        function obj = set_R(obj,R_)
            % Change the bend radius of the duct (and update the interpolants)
            obj.R = R_;
            obj = obj.update_R_data();
            obj = obj.setup_interpolants();
        end
        function F = migration_force(obj,r,z)
            % Get the (net) migration force for a neutrally buoyant
            % spherical particle centred at (r,z) within the cross-section
            % (non-dimensionalised via \rho U_m^2 a^4 / H^2 )
            F = [obj.Fr_RBS(r,z),obj.Fz_RBS(r,z)];
        end
        function V = migration_velocity(obj,r,z)
            % Get the migration velocity for a neutrally buoyant
            % spherical particle centred at (r,z) within the cross-section
            % (non-dimensionalised via U_m a / H )
            V = [obj.Fr_RBS(r,z)/obj.Cr_RBS(r,z), ...
                 obj.Fz_RBS(r,z)/obj.Cz_RBS(r,z)];
        end
        function C = drag_coefficient(obj,r,z)
            % Get the drag coefficients in the r,z directions of a neutrally buoyant
            % spherical particle centred at (r,z) within the cross-section
            % (non-dimensionalised via \mu a )
            C = [obj.Cr_RBS(r,z),obj.Cz_RBS(r,z)];
        end
        function S = secondary_flow_drag(obj,r,z)
            % Get the drag coefficients in the r,z directions of a neutrally buoyant
            % spherical particle centred at (r,z) within the cross-section
            % (non-dimensionalised via \rho U_m^2 a^4 / H^2 )
            S = [obj.Sr_RBS(r,z),obj.Sz_RBS(r,z)];
        end
        function Up = axial_velocity(obj,r,z)
            % Get the terminal/steady axial velocity of a neutrally buoyant 
            % spherical particle centred at (r,z) within the cross-section
            % (non-dimensionalised via U_m a / H )
            Up = obj.Up_RBS(r,z);
        end
        function W = spin_components(obj,r,z)
            % Get the terminal/steady r,z spin components of a neutrally buoyant 
            % spherical particle centred at (r,z) within the cross-section
            % (non-dimensionalised via U_m / H )
            W = [obj.Wr_RBS(r,z),obj.Wz_RBS(r,z)];
        end
        function plot_migration_force(obj)
            % Produces a rough sketch of the magnitude of the migration force including
            % the zero contours of the horizontal and vertical components.
            Fr = squeeze(obj.aR_data(5,:,:)+obj.kappa*obj.aR_data(7,:,:));
            Fz = squeeze(obj.aR_data(6,:,:)+obj.kappa*obj.aR_data(8,:,:));
            RS = squeeze(obj.aR_data(1,:,:));
            ZS = squeeze(obj.aR_data(2,:,:));
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0,0,1+3*obj.ar,3];
            hold on;
            contourf(RS,ZS,(Fr.^2+Fz.^2).^0.5,17);
            colorbar();
            contour(RS,ZS,Fr,[0,0],'k');
            contour(RS,ZS,Fz,[0,0],'w');
            nhH = 1.0/obj.a;
            nhW = obj.ar/obj.a;
            plot([-nhW+1,nhW-1,nhW-1,-nhW+1,-nhW+1],...
                 [-nhH+1,-nhH+1,nhH-1,nhH-1,-nhH+1],'r--');
            set(gca,'DataAspectRatio',[1,1,1],...
                    'PlotBoxAspectRatio',[1,1,1],...
                    'XLim',[-nhW,nhW],...
                    'YLim',[-nhH,nhH]);
            hold off
        end
    end
end
