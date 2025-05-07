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
% Thickness (m); No. Interpolation Points; Elasticity
h = 4e-3;  N = 20;   [c,~,rho,dof] = f_CB_Aniso('Lamb');            % CB
ht = 1e-3; Nt = 20;  rhot = 1e3; clt = 1500; lamt = rhot*clt^2;     % ST
hb = 1e-3; Nb = 20;  rhob = 930; clb = 1480; lamb = rhob*clb^2;     % MR

%% Build matrices
[L2,L1,L0,L] = f_matrix_generation(layer,c,rho,N,h,...
                                   lamt,clt,Nt,ht, ...
                                   lamb,clb,Nb,hb,dof);

%% 3D dispersion curves (complex k)
Timer3D = tic;
f0 = linspace(100,1e6,500).';
K = solveWavenumber(f0,L2,L1,L0,L);
fprintf("3D time: %.2f s\n",toc(Timer3D));
plotDispersion3D(f0,K);

%% 2D dispersion curves (real k)
Timer2D = tic;
k0 = linspace(1e-2,4000,500);
[SF, SCP, AF, ACP] = solveFrequency(k0,L2,L1,L0,L,N);
fprintf("2D time: %.2f s\n",toc(Timer2D));
plotDispersion2D(k0,SF,SCP,AF,ACP,layer);

%% ---------------------- Helper functions --------------------------------
function K = solveWavenumber(f0,L2,L1,L0,L)
    nK = 2*size(L,2);   K = nan(numel(f0),nK);
    parfor ii = 1:numel(f0)
        w = 2*pi*f0(ii);
        K(ii,:) = -1i*polyeig(L0 + w^2*L, L1, L2).';
    end
end

function [SF, SCP, AF, ACP] = solveFrequency(k0,L2,L1,L0,L,N)
    nM = size(L,2); nK = numel(k0);
    SF  = nan(nM,nK); AF = SF; SCP = SF; ACP = SF;

    parfor ii = 1:nK
        [U,w2] = polyeig((1i*k0(ii))^2*L2 + (1i*k0(ii))*L1 + L0, L);
        w = sqrt(w2);

        maskSmode = real(U(N,:)).*real(U(1,:)) > eps ...
                  | imag(U(N,:)).*imag(U(1,:)) > eps;

        sf_loc  = nan(nM,1); af_loc  = nan(nM,1); scp_loc = nan(nM,1); acp_loc = nan(nM,1);

        sf_loc(maskSmode)   = real(w(maskSmode)) / (2*pi);
        af_loc(~maskSmode)  = real(w(~maskSmode)) / (2*pi);
        scp_loc(maskSmode)  = real(w(maskSmode)) / k0(ii);
        acp_loc(~maskSmode) = real(w(~maskSmode)) / k0(ii);

        sf_loc(sf_loc < 5) = NaN;
        af_loc(af_loc < 5) = NaN;

        SF(:,ii)  = sf_loc; AF(:,ii)  = af_loc; SCP(:,ii) = scp_loc; ACP(:,ii) = acp_loc;
    end
end

function plotDispersion3D(f0,K)
    dom.kr = [-4 4];   % real(k)  axis limits (rad/mm)
    dom.ki = [-4 4];   % imag(k)  axis limits (Np/mm)
    dom.f  = [ 0 1];   % freq     axis limits (MHz)
    
    % -------- 1. sort the roots and build helper matrices ----------------
    Ksort = sort(K,2);
    Fmat  = f0.*ones(size(Ksort));
    
    % Convert to the units used in the plots
    Krmm  = real(Ksort)/1e3;          % rad/mm
    Kimm  = imag(Ksort)/1e3;          % Np/mm
    FMHz  = Fmat/1e6;                 % MHz
    
    % -------- 2. apply the domain mask *before* classifying roots --------
    mask =  Krmm>=dom.kr(1) & Krmm<=dom.kr(2) & ...
            Kimm>=dom.ki(1) & Kimm<=dom.ki(2) & ...
            FMHz>=dom.f(1)  & FMHz<=dom.f(2);
    
    Ksort(~mask) = NaN;   % cull the out‑of‑window points
    FMHz (~mask) = NaN;
    
    % -------- 3. split into real / imaginary / complex root types --------
    KR = Ksort;  KR(abs(imag(Ksort))>=10) = NaN;   % mainly‑real roots
    KI = Ksort;  KI(abs(real(Ksort))>=0.1) = NaN;  % mainly‑imag roots
    KC = Ksort;  KC(abs(imag(Ksort))<10 | abs(real(Ksort))<0.1) = NaN; % complex
    
    % -------- 4. draw ----------------------------------------------------
    fig1 = figure("Name","3D Dispersion Curves","Units","normalized", ...
                 "Position",[0.05 0.05 0.6 0.7]);
    tl1  = tiledlayout(fig1,2,3,"Padding","compact","TileSpacing","compact");
    axs = [nexttile(tl1,[1 3]), nexttile(tl1,4), nexttile(tl1,5), nexttile(tl1,6)];
    
    for ax = axs
        scatter3(ax,real(KC(:))/1e3,imag(KC(:))/1e3,FMHz(:),10,'filled', ...
                 'MarkerFaceColor',[0 .6 0]); hold(ax,'on');
        scatter3(ax,real(KR(:))/1e3,imag(KR(:))/1e3,FMHz(:),10,'filled', ...
                 'MarkerFaceColor','b');
        scatter3(ax,real(KI(:))/1e3,imag(KI(:))/1e3,FMHz(:),10,'filled', ...
                 'MarkerFaceColor','r');
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

    % print(gcf, '3D_dispersion_elas.svg', '-dsvg', '-vector', '-r300');
end

function plotDispersion2D(k0,SF,SCP,AF,ACP,layer)
    fig2 = figure("Name","2D Dispersion Curves","Units","normalized", ...
             "Position",[0.25 0.25 0.6 0.4]);
    tl2  = tiledlayout(fig2,1,2,"Padding","compact","TileSpacing","compact");

    cols = {{'b','r'}, {'k','k'}};
    sel  = 1 + (layer>=2);
    cSF = cols{sel}{1}; cAF = cols{sel}{2};

    nexttile(tl2,1)
    hold on
    scatter(k0/1e3, SF/1e6, 10, cSF, 'filled');
    scatter(k0/1e3, AF/1e6, 10, cAF, 'filled');
    set(gca,"FontSize",12); grid on;
    xlabel('$k$ (rad/mm)','Interpreter','latex');
    ylabel('$f$ (MHz)','Interpreter','latex');
    title("\bf{a)}", Interpreter='latex');
    xlim([0 4]); ylim([0 1]);

    nexttile(tl2,2)
    hold on
    scatter(SF/1e6, SCP/1e3, 10, cSF, 'filled');
    scatter(AF/1e6, ACP/1e3, 10, cAF, 'filled');
    set(gca,"FontSize",12); grid on;
    xlabel('$f$ (MHz)','Interpreter','latex');
    ylabel('$c_p$ (km/s)','Interpreter','latex');
    title("\bf{b)}", Interpreter='latex');
    ylim([0 10]); xlim([0 1]);

    % print(gcf, '2D_dispersion_elas.svg', '-dsvg', '-vector', '-r300');
end
