function [L2, L1, L0, L] = f_matrix_generation(layer, C, rho, N, h, ...
    lamt, clt, Nt, ht, lamb, clb, Nb, hb, dof)
    % INPUTS:
    %   layer – 1 (CB only), 2 (CB + ST), or 3 (ST + CB + MR)
    %   C     – stiffness matrix of cortical bone (CB)
    %   rho   – density of CB
    %   N,h   – n.o. interpolation points (IP) and thickness for CB
    %   lamt,clt,Nt,ht – soft tissue (ST): lambda, speed, n.o. IP, thickness
    %   lamb,clb,Nb,hb – marrow (MR): lambda, speed, n.o. IP, thickness
    %
    % OUTPUTS:
    %   L2, L1, L0, L – coefficient matrices

    % Extract stiffness submatrices for CB
    c11 = C{1,1}; c12 = C{1,2}; c21 = C{2,1}; c22 = C{2,2};
    I = eye(size(c11));

    %% Layer 1: CB
    if layer == 1
        [D1, D2, Id] = getDiffOps(N);
        L2 = h^2 * kron(c11, Id);
        L1 = h   * kron(c12 + c21, D1);
        L0 =       kron(c22, D2);
        L  = h^2 * kron(rho * I, Id);

        % Stress-free boundary conditions (CB)
        B1 = h   * kron(c21, Id([1, N], :));
        B0 =       kron(c22, D1([1, N], :));
        idx = getBoundaryIndices(N, dof);
        L2(idx,:) = 0; L1(idx,:) = B1; L0(idx,:) = B0; L(idx,:) = 0;

    %% Layer 2: ST + CB
    elseif layer == 2
        % CB part
        [D1, D2, Id] = getDiffOps(N);
        M2 = kron(c11, Id); M1 = kron(c12 + c21, D1); M0 = kron(c22, D2);
        M  = kron(rho * I, Id);
        idx = getBoundaryIndices(N, dof);
        B1 = kron(c21, Id([1, N], :)); B0 = kron(c22, D1([1, N], :));
        M2(idx,:) = 0; M1(idx,:) = B1; M0(idx,:) = B0; M(idx,:) = 0;
        normCB = norm(M0, 'fro'); % Normalize CB part
        M2 = M2/normCB; M1 = M1/normCB; M0 = M0/normCB; M = M/normCB;

        % ST part
        [D1t, D2t, Idt] = getDiffOps(Nt);

        % Matrix assembly
        nTotal = Nt + 2*N;
        L2 = zeros(nTotal); L1 = zeros(nTotal); L0 = zeros(nTotal); L = zeros(nTotal);
        L2(1:Nt,1:Nt) = ht^2 * Idt;
        L2(Nt+1:end,Nt+1:end) = h^2 * M2;
        L1(Nt+1:end,Nt+1:end) = h * M1;
        L0(1:Nt,1:Nt) = D2t;
        L0(Nt+1:end,Nt+1:end) = M0;
        L(1:Nt,1:Nt) = ht^2/clt^2 * Idt;
        L(Nt+1:end,Nt+1:end) = h^2 * M;

        % Boundary & interface conditions (ST - CB)
        L2(1,1:Nt) = ht^2 * Idt(1,:); L0(1,1:Nt) = D2t(1,:); L(1,1:Nt) = 0;
        L2(Nt,1:Nt) = 0; L0(Nt,1:Nt) = -D1t(Nt,:); L(Nt,1:Nt) = 0;
        L0(Nt,Nt+N+1:end) = Id(1,:);
        L2(Nt+N+1,1:Nt) = -ht^2 * lamt * Idt(Nt,:) / normCB;
        L0(Nt+N+1,1:Nt) = -lamt * D2t(Nt,:) / normCB;

    %% Layer 3: ST + CB + MR
    elseif layer == 3
        % CB part
        [D1, D2, Id] = getDiffOps(N);
        M2 = kron(c11, Id); M1 = kron(c12 + c21, D1); M0 = kron(c22, D2);
        M  = kron(rho * I, Id);
        idx = getBoundaryIndices(N, dof);
        B1 = kron(c21, Id([1, N], :)); B0 = kron(c22, D1([1, N], :));
        M2(idx,:) = 0; M1(idx,:) = B1; M0(idx,:) = B0; M(idx,:) = 0;
        normCB = norm(M0, 'fro');
        M2 = M2/normCB; M1 = M1/normCB; M0 = M0/normCB; M = M/normCB;

        % ST and MR
        [D1t, D2t, Idt] = getDiffOps(Nt);
        [D1b, D2b, Idb] = getDiffOps(Nb);

        % Total size and initialization
        nTotal = Nt + 2*N + Nb;
        L2 = zeros(nTotal); L1 = zeros(nTotal); L0 = zeros(nTotal); L = zeros(nTotal);

        % ST block
        L2(1:Nt,1:Nt) = ht^2 * Idt;
        L0(1:Nt,1:Nt) = D2t;
        L(1:Nt,1:Nt) = ht^2/clt^2 * Idt;

        % CB block
        L2(Nt+1:Nt+2*N,Nt+1:Nt+2*N) = h^2 * M2;
        L1(Nt+1:Nt+2*N,Nt+1:Nt+2*N) = h * M1;
        L0(Nt+1:Nt+2*N,Nt+1:Nt+2*N) = M0;
        L(Nt+1:Nt+2*N,Nt+1:Nt+2*N) = h^2 * M;

        % MR block
        L2(Nt+2*N+1:end,Nt+2*N+1:end) = hb^2 * Idb;
        L0(Nt+2*N+1:end,Nt+2*N+1:end) = D2b;
        L(Nt+2*N+1:end,Nt+2*N+1:end) = hb^2/clb^2 * Idb;

        % Boundary & interface conditions
        % Top of ST
        L2(1,1:Nt) = ht^2 * Idt(1,:); L0(1,1:Nt) = D2t(1,:); L(1,1:Nt) = 0;

        % ST - CB
        L2(Nt,1:Nt) = 0; L0(Nt,1:Nt) = -D1t(Nt,:); L(Nt,1:Nt) = 0;
        L0(Nt,Nt+N+1:Nt+2*N) = Id(1,:);
        L2(Nt+N+1,1:Nt) = -ht^2 * lamt * Idt(Nt,:) / normCB;
        L0(Nt+N+1,1:Nt) = -lamt * D2t(Nt,:) / normCB;

        % CB - MR
        L2(Nt+2*N+1,Nt+2*N+1:end) = 0;
        L0(Nt+2*N+1,Nt+2*N+1:end) = -D1b(1,:);
        L(Nt+2*N+1,Nt+2*N+1:end) = 0;
        L0(Nt+2*N+1,Nt+N+1:Nt+2*N) = Id(N,:);
        L2(Nt+2*N,Nt+2*N+1:end) = -hb^2 * lamb * Idb(1,:) / normCB;
        L0(Nt+2*N,Nt+2*N+1:end) = -lamb * D2b(1,:) / normCB;

        % Bottom of MR
        L2(end,Nt+2*N+1:end) = hb^2 * Idb(Nb,:);
        L0(end,Nt+2*N+1:end) = D2b(Nb,:);
        L(end,Nt+2*N+1:end) = 0;
    end
end

% ---------- Helper Functions ----------
function [D1, D2, Id] = getDiffOps(N)
    [~, D] = chebdif(N, 2);
    D1 =  2 * D(:,:,1);
    D2 =  4 * D(:,:,2);
    Id = eye(N);
end

function idx = getBoundaryIndices(N, dof)
    idx = [(0:length(dof)-1)*N+1; (1:length(dof))*N];
end
