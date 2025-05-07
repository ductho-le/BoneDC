function [C, N, rho, dof] = f_CB_Aniso(wave)
% Returns anisotropic stiffness and viscosity tensors for cortical bone
% in tensor form, extracted for specific wave propagation types.
%
% INPUT:
%   wave - a string specifying the wave type:
%          'Lamb'            -> in-plane motion (DOFs 1 and 2)
%          'SH(1layer)'      -> out-of-plane motion only (DOF 3)
%          'Lamb+SH(1layer)' -> all motion components (DOFs 1 to 3)
%
% OUTPUT:
%   C   - 2x2 cell array containing submatrices of the stiffness tensor
%         Each C{i,j} corresponds to c(i,dof,dof,j)
%
%   N   - 2x2 cell array containing submatrices of the viscosity tensor
%         Structured the same way as C
%
%   rho - scalar mass density of cortical bone [kg/m^3]
%
%   dof - degrees of freedom used based on the selected wave type
%
% Notes:
% - C0 and N0 are provided in Voigt notation and mapped to full
%   4th-order tensors (3x3x3x3) for stiffness and viscosity.

    % Voigt notation to tensor notation
    C0 = [31.7 13.4 13.4 0 0 0;
          13.4 20.2 10.7 0 0 0;
          13.4 10.7 20.3 0 0 0;
          0 0 0 4.80 0 0;
          0 0 0 0 6.32 0;
          0 0 0 0 0 6.38]*1e9;

    rho = 1948;

    N0 = [157 121 121 0 0 0;
          121 109 121 0 0 0;
          121 121 109 0 0 0;
          0 0 0 18 0 0;
          0 0 0 0 18 0;
          0 0 0 0 0 18];
    
    c = zeros(3,3,3,3); % Stiffness tensor
    n = zeros(3,3,3,3); % Viscosity tensor

    for i = 1:3 
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    if i==j
                        voigti = i;
                    else
                        voigti = 9-i-j;
                    end
    
                    if k==l
                        voigtj = k;
                    else
                        voigtj = 9-k-l;
                    end
    
                    c(i,j,k,l) = C0(voigti,voigtj);
                    n(i,j,k,l) = N0(voigti,voigtj);
                end
            end
        end
    end

    if strcmp(wave, 'Lamb')
        dof = 1:2;
    elseif strcmp(wave, 'SH(1layer)')      % fluid doesn't support Shear waves
        dof = 3;
    elseif strcmp(wave, 'Lamb+SH(1layer)') % fluid doesn't support Shear waves
        dof = 1:3;
    else
        error("Invalid wave input!!!")
    end

    C = cell(2,2);
    C{1,1} = squeeze(c(1,dof,dof,1));
    C{1,2} = squeeze(c(1,dof,dof,2));
    C{2,1} = squeeze(c(2,dof,dof,1));
    C{2,2} = squeeze(c(2,dof,dof,2));

    N = cell(2,2);
    N{1,1} = squeeze(n(1,dof,dof,1)); 
    N{1,2} = squeeze(n(1,dof,dof,2)); 
    N{2,1} = squeeze(n(2,dof,dof,1)); 
    N{2,2} = squeeze(n(2,dof,dof,2)); 
end