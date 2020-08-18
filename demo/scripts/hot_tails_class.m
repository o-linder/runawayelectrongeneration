classdef hot_tails_class
    %   A class for the hot-tail runaway generation
    %   initialization (constructor) and the function itself (calculate)
    %   Based on H. M. Smith PoP 15 072502 2008, integrating (5) with tau
    %   from (17).
    %   Valid for exponential T(t) and n(t) evolution
    %   (see the function at the end)
    %
    %   ppg@ipp.mpg.de 2015-08
    %
    %   Initialize the calculation by calling
    %
    %     H = hot_tails_class('t_star',t_star,'n0',n0,'nmax',nmax,'T0',T0,'Tf',Tf);
    %
    %   The input parameters can be vectors or scalars, but must have the
    %   same size.
    %   - n0, nmax: initial and *final* density [m^-3]
    %   - T0, Tf  : initial and *final* temperature [eV].
    %   - t_star  : characteristic time in exp(-t/t_star) [s]
    %
    %   Once initialization is done you can calculate the
    %   __total number of hot-tails generated until "time"__
    %   by calling
    %
    %     ht = H.calculate(time,E)
    %
    %   Where
    %   - time is in the same units as t_star (preferably [s]
    %   - E electric field in [V/m]
    %
    %   NOTE: at constant E-fields the formula gives decreasing HT density.
    %   This is something under investigation.
    %   I advise to use a time-dependent E(t) and always store
    %   HT(t) = max( ht(t),HT(t-1) )
    %   HT current may drop if E drops due to changes in v_c, but that
    %   doesn't mean these particles aren't runaways / relativistic anymore.
    
    
    properties
        Burst = struct('type', 'analitic'); % structure for storing the burst runaway calculation parameters
        PhysConst = struct('qe',1.602*10^(-19), 'me',9.109*10^(-31),...
            'epsilon',8.8542*10^(-12), 'vlight',299792458); % Physical constants
    end
    
    methods
        % Constructor
        function o = hot_tails_class(varargin)
            if nargin <= 1
            elseif nargin >= 2
                k = 1;
                while k<= length(varargin)
                    switch varargin{k}
                        case 't_star'
                            o.Burst.t_star = varargin{k+1};
                        case 't_0'
                            o.Burst.t_star = varargin{k+1};
                            
                        case 'n0'
                            o.Burst.n0 = varargin{k+1};
                        case 'nmax'
                            o.Burst.nmax = varargin{k+1};
                        case 'T0'
                            o.Burst.T0 = varargin{k+1};
                        case 'Tf'
                            o.Burst.Tf = varargin{k+1};
                    end
                    k = k+2;
                end
            end
            
            % Call initialization
            o = o.InitHot_tail;
        end
        
        % Initialization - storing normalized parameters
        function o = InitHot_tail(o)
            
            disp('Initializing hot tail runaway generation')
            
            if  ~isfield(o.Burst,'t_star')
                o.Burst.t_star = input('t_star? : ');
            end
            if  ~isfield(o.Burst,'n0')
                o.Burst.n0 = input('n0? : ');
            end
            if  ~isfield(o.Burst,'nmax')
                o.Burst.nmax = input('nmax? : ');
            end
            if  ~isfield(o.Burst,'T0')
                o.Burst.T0 = input('T0? : ');
            end
            if  ~isfield(o.Burst,'Tf')
                o.Burst.Tf = input('Tf? : ');
            end
            
            % Check sizes
            numels = cellfun(@numel,{o.Burst.t_star,o.Burst.n0,o.Burst.nmax,o.Burst.T0,o.Burst.Tf});
            if ~all(numels==numels(1))
                disp(numels);
                error(['Input sizes do not match: ',num2str(numels)])
            else
                o.Burst.size = size(o.Burst.n0);
            end
            
            % Initial thermal velocity
            o.Burst.v_T0 = sqrt(2*o.Burst.T0.*o.PhysConst.qe./o.PhysConst.me);
            
            o.Burst.lnLambda = o.loglambda(o.Burst.n0,o.Burst.T0);
            
            % initial collision frequency
            o.Burst.nu_0 = o.Burst.n0.*o.PhysConst.qe^4.*o.Burst.lnLambda ./ ...
                (4*pi*o.PhysConst.epsilon^2*o.PhysConst.me^2.*o.Burst.v_T0.^3);
            
            o.Burst.p_star=o.Burst.t_star.*o.Burst.nu_0;
            
            % The following array stores when the temperature starts to change
            % Burst calculation are started from this time instant
            % Currently not used but can be implemented.
            o.Burst.starttime = zeros(o.Burst.size); %(s)
            
            o.Burst.nmax_n0 = o.Burst.nmax ./ o.Burst.n0;
            
            %Constant used later for critical speed v_c^2
            o.Burst.xcconst = o.PhysConst.qe^2 / (8*pi*o.PhysConst.epsilon^2);
            o.Burst.tau_A  = zeros(o.Burst.size); % tau at the actual timestep
            
        end
        
        function [nb, o] = calculate(o, time, E)
            % Returns the number of burst runaways at time t
            % Input arguments:
            %   t: time [s]
            %   E: Electric field [V/m]
            % The return value is the total hot tail runaway density after time t.
            
            nb=zeros(o.Burst.size);
            % We could just use the same E-field for all the grid points,
            % but it is safer to throw a warning and request a vector.
            if numel(E) ~= numel(nb)
                error(['You initialized with vectors. ',...
                    'Please provide E-field as a vector of size: ',num2str(o.Burst.size)]);
            end
            
            % Get T(t) and n(t) from exponential decay.
            % Can be implemented as input parameters too.
            Temperature = o.temperature(time);
            Density = o.density(time);
            
            o.Burst.v_T_A = sqrt(2*Temperature*o.PhysConst.qe/o.PhysConst.me);
            vT_vT0_3 = (o.Burst.v_T_A ./ o.Burst.v_T0).^3;
            
            for k=1:length(nb)
                if o.Burst.starttime(k)<=time
                    % Calculate tau
                    % -----------------------------------------------------
                    t_b = time-o.Burst.starttime(k);
                    
                    p=t_b*o.Burst.nu_0(k);
                    if (p < 2*o.Burst.p_star(k))
                        tau = p^2 / (4*o.Burst.p_star(k)) * o.Burst.nmax_n0(k);
                    else
                        tau = (p-o.Burst.p_star(k))*o.Burst.nmax_n0(k);
                    end
                    o.Burst.tau_A(k) = tau;
                    % calculate the distribution function -------------------------------
                    % the distribution function multiplied by 4*pi*v_T^3*Dx
                    Fa = @(x) (4*o.Burst.n0(k) * vT_vT0_3(k) / sqrt(pi)) ...
                        *exp(- (x.^3*vT_vT0_3(k) + 3*o.Burst.tau_A(k)).^(2/3));
                    
                    % Calculate the number of runaways ----------------------------------
                    % square of the critical velocity.
                    % lnL is the initial value
                    %xc_2 = o.Burst.xcconst .* Density(k) .* o.Burst.lnLambda(k) ./ Temperature(k) ./ E(k);
                    % lnL is the actual value
                    xc_2 = o.Burst.xcconst .* Density(k) .* o.loglambda(Density(k),Temperature(k)) ./ Temperature(k) ./ E(k);
                    
                    % burst runaway density
                    func = @(x) Fa(x) .* (x.^2 - xc_2);
                    nb(k) = integral(func,sqrt(xc_2),Inf);
                end
            end
            
        end
        
        % Functions to calculate "exponential decay" temperature and density evolution
        % Time in seconds; or the same units as t0
        function T = temperature(o,t)
            T = o.Burst.Tf + (o.Burst.T0-o.Burst.Tf).*exp(-(t.*ones(o.Burst.size))./o.Burst.t_star);
        end
        
        function n = density(o,t)
            n = o.Burst.nmax + (o.Burst.n0 - o.Burst.nmax).*exp(-(t.*ones(o.Burst.size))./o.Burst.t_star);
        end
    end
    
    methods (Static)
        function lnL = loglambda(n,T)
            % lnLambda = loglambda(ne, te) returns the Coulomb logarithm for
            % electron-electron collisions
            % Arguments:
            %   T: electron temperature (keV)
            %   n: electron density (m^{-3})
            % Return value: lnLambda = 14.9 - 0.5 ln(n/10^20) + ln(T)
            % ref: Wesson Tokamaks section 14.5
            lnL = 14.9 - 0.5*log(n./1e20) + log(0.001.*T);
        end
    end
end


