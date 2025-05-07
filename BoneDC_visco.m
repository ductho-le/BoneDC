%% %%%%%%%%%%%% Guided Wave Dispersion Curves in Bone by SCM %%%%%%%%%%% %%

% Ductho Le (ductho.le@outlook.com)
% -------------------------------------------------------------------------

clear; clc; close all;
if isempty(gcp("nocreate")), parpool("threads"); end

% 1 = cortical bone (CB)
% 2 = cortical bone (CB) + soft tissue (ST)
% 3 = cortical bone (CB) + soft tissue (ST) + marrow (MR)
layer = 1;

%% Material parameters
% Thickness (m); No. Interpolation Points; Elasticity, Viscousity
h = 4e-3;  N = 20;   [c,n,rho,dof] = f_CB_Aniso('Lamb');            % CB
ht = 1e-3; Nt = 20; rhot = 1e3; clt = 1500; lamt = rhot*clt^2; nt = 1.97; % ST
hb = 1e-3; Nb = 20; rhob = 930; clb = 1480; lamb = rhob*clb^2; nb = 1.97; % MR

%% 3D dispersion curves
Timer = tic;
f0 = linspace(100,1e6,500).';
[AK, SK] = solveWavenumber(f0,layer,c,n,rho,N,h,...
                       lamt,nt,rhot,Nt,ht, ...
                       lamb,nb,rhob,Nb,hb,dof);
fprintf("3D time: %.2f s\n",toc(Timer));
plotDispersion3D(f0,AK,SK,layer);

%% ---------------------- Helper functions --------------------------------
function [AK, SK] = solveWavenumber(f0,layer,c,n,rho,N,h,...
                       lamt,nt,rhot,Nt,ht, ...
                       lamb,nb,rhob,Nb,hb,dof)
    size_mat = 2*N + (layer>1)*Nt + (layer>2)*Nb;
    AK = nan(length(f0), size_mat*2); SK = AK;
    parfor ii = 1:numel(f0)
        w = 2*pi*f0(ii);
        C = cellfun( @(Ci,Ni) Ci + 1i*w*Ni, c, n, 'UniformOutput', false);
        [L2,L1,L0,L] = f_matrix_generation(layer,C,rho,N,h,...
                       lamt+1i*w*nt,sqrt((lamt+1i*w*nt)/rhot),Nt,ht, ...
                       lamb+1i*w*nb,sqrt((lamb+1i*w*nb)/rhob),Nb,hb,dof);

        [U, kii] = polyeig(L0 + w^2*L, L1, L2);

        maskSmode = real(U(N,:)).*real(U(1,:)) > eps ...
                  | imag(U(N,:)).*imag(U(1,:)) > eps;

        sk_loc  = nan(1,size_mat*2); ak_loc  = nan(1,size_mat*2);

        sk_loc(maskSmode) = -1i * kii(maskSmode);
        ak_loc(~maskSmode) = -1i * kii(~maskSmode);

        SK(ii,:) = sk_loc; AK(ii,:) = ak_loc; 
    end
end

function plotDispersion3D(f0,AK,SK,layer)
    dom.kr = [-4 4];   % real(k)  axis limits (rad/mm)
    dom.ki = [-4 4];   % imag(k)  axis limits (Np/mm)
    dom.f  = [ 0 1];   % freq     axis limits (MHz)
    
    % -------- 1. sort the roots and build helper matrices ----------------
    AKsort = sort(AK,2); SKsort = sort(SK,2);
    Fmat  = f0.*ones(size(AKsort));
    
    % Convert to the units used in the plots
    AKrmm  = real(AKsort)/1e3; SKrmm  = real(SKsort)/1e3;         % rad/mm
    AKimm  = imag(AKsort)/1e3; SKimm  = imag(SKsort)/1e3;         % Np/mm
    FMHz  = Fmat/1e6;                                             % MHz
    
    % -------- 2. apply the domain mask *before* classifying roots --------
    maskA =  AKrmm>=dom.kr(1) & AKrmm<=dom.kr(2) & ...
            AKimm>=dom.ki(1) & AKimm<=dom.ki(2);
    maskS = SKrmm>=dom.kr(1) & SKrmm<=dom.kr(2) & ...
            SKimm>=dom.ki(1) & SKimm<=dom.ki(2);
    maskF = FMHz>=dom.f(1)  & FMHz<=dom.f(2);
    
    AKsort(~ maskA & maskF) = NaN;   % cull the out‑of‑window points
    SKsort(~ maskS & maskF) = NaN;
    
    % -------- 3. split into low / high attenuation modes -----------------
    att = 0.1; % threshold attenuation
    maskLS = abs(SKimm) < att; maskLA = abs(AKimm) < att;

    AKL = AKsort(maskLA); SKL = SKsort(maskLS);
    AKH = AKsort(~maskLA); SKH = SKsort(~maskLS);
    
    % -------- 4. draw ----------------------------------------------------
    fig1 = figure("Name","3D Dispersion Curves","Units","normalized", ...
                 "Position",[0.05 0.05 0.6 0.7]);
    tl1  = tiledlayout(fig1,2,3,"Padding","compact","TileSpacing","compact");
    axs = [nexttile(tl1,[1 3]), nexttile(tl1,4), nexttile(tl1,5), nexttile(tl1,6)];
    
    cols = {{'r','b',[.9 .8 0],[0 .6 0]}, {'r','r','b','b'}};
    sel  = 1 + (layer>=2);
    cAKL = cols{sel}{1}; cSKL = cols{sel}{2}; cAKH = cols{sel}{3}; cSKH = cols{sel}{4};

    for ax = axs
        scatter3(ax,real(AKL(:))/1e3,imag(AKL(:))/1e3,FMHz(maskLA),10,'filled', ...
                 'MarkerFaceColor',cAKL); hold(ax,'on');
        scatter3(ax,real(SKL(:))/1e3,imag(SKL(:))/1e3,FMHz(maskLS),10,'filled', ...
                 'MarkerFaceColor',cSKL);
        scatter3(ax,real(AKH(:))/1e3,imag(AKH(:))/1e3,FMHz(~maskLA),10,'filled', ...
                 'MarkerFaceColor',cAKH);
        scatter3(ax,real(SKH(:))/1e3,imag(SKH(:))/1e3,FMHz(~maskLS),10,'filled', ...
                 'MarkerFaceColor',cSKH);
        xlabel(ax,'$k_r$ (rad/mm)','Interpreter','latex');
        ylabel(ax,'$k_i$ (Np/mm)','Interpreter','latex');
        zlabel(ax,'$f$ (MHz)','Interpreter','latex');
        set(ax,'XLim',dom.kr,'YLim',dom.ki,'ZLim',dom.f,'FontSize',12,...
            'XTick',-4:2:4,'YTick',-4:2:4,'ZTick',0:0.25:1);
        grid(ax,'on'); box(ax,'off');
    end
    
    view(axs(1),[35 25]); title(axs(1),'\bf{a)}','Interpreter','latex');
    view(axs(2),[ 0  0]); title(axs(2),'\bf{b)}','Interpreter','latex');
    view(axs(3),[90  0]); title(axs(3),'\bf{c)}','Interpreter','latex');
    view(axs(4),[ 0 90]); title(axs(4),'\bf{d)}','Interpreter','latex');

    % print(gcf, '3D_dispersion_visco.svg', '-dsvg', '-vector', '-r300');
end