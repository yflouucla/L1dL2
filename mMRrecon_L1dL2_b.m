function [u,pm] = mMRrecon_L1dL2_b(R, f, pm)
%=============================================================
% ADMM for the L1/L2 constrained model for MRI reconstruction
%
% Solves
%           min  norm(grad(x),1)/norm(grad(x),2)
%           s.t. Ax = b
%
% Reference: "A Scale Invariant Approach for Sparse Signal Recovery" 
%             Yaghoub Rahimi, Chao Wang, Hongbo Dong, Yifei Lou 
% Available at: 
%             https://arxiv.org/abs/1812.08852
% 
% Author: Chao Wang  
% Date: Nov 11. 2018
%============================================================= 

[rows,cols] = size(f);

Num_iter = 3000; rho1 = 50; rho2 = 50;
u0 = zeros(rows,cols); lambda = 50; rho3 = 1;
tol = 1e-5;

if isfield(pm,'rho1'); rho1 = pm.rho1; end
if isfield(pm,'rho2'); rho2 = pm.rho2; end
if isfield(pm,'rho3'); rho3 = pm.rho3; end
if isfield(pm,'tol'); tol = pm.tol; end
if isfield(pm,'lambda'); lambda = pm.lambda; end
if isfield(pm,'Num_iter'); Num_iter = pm.Num_iter; end


u = u0; v = u;


% Build Kernels
scale = sqrt(rows*cols);


uker = zeros(rows,cols);
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker = lambda*(conj(R).*R)+(rho1+rho2)*fft2(uker) + rho3*ones(rows,cols);
% uker(uker<0)=0; uker(uker>1) =1;

ff = f;
f0 = ff; f = f0; murf = ifft2(lambda*R.*f0)*scale;

% Augmented Lagrangian Parameters
dx = zeros(rows,cols);
dy = zeros(rows,cols);
hx = zeros(rows,cols);
hy = zeros(rows,cols);

bx = zeros(rows,cols);
by = zeros(rows,cols);
cx = zeros(rows,cols);
cy = zeros(rows,cols);
h = bx;

List = zeros(Num_iter,4);
for j = 1:Num_iter
        uold = u;

        % u-update
        rhs = murf+rho1*Dxt(dx-bx)+rho1*Dyt(dy-by)+rho2*Dxt(hx-cx)+...
            rho2*Dyt(hy-cy)+ rho3*(v-h);
        u = real(ifft2(fft2(rhs)./uker));
        v = u+h;
        v(v<0)=0; v(v>1) =1;
%         u(u<0)=0; u(u>1) =1;
        % dx,dy-update
        Dxu = Dx(u);
        Dyu = Dy(u);
        hnorm = sqrt(norm(hx(:))^2+norm(hy(:))^2);
        dx = shrink(Dxu+bx, 1/(rho1*hnorm));
        dy = shrink(Dyu+by, 1/(rho1*hnorm));

        % hx,hy-update
        d1 = Dxu + cx;
        d2 = Dyu + cy;
        etha = sqrt(norm(d1(:))^2+norm(d2(:))^2);
        c = norm(dx(:),1)+norm(dy(:),1);
        [hx, hy,tau] = mupdate_h(c,etha,rho2,d1,d2);

        % bx,by,cx,cy-update
        bx = bx+(Dxu-dx);
        by = by+(Dyu-dy);
        cx = cx+(Dxu-hx);
        cy = cy+(Dyu-hy);
        h  = h + u - v;

        f = f+f0-R.*fft2(u)/scale;
        murf = ifft2(lambda*R.*f)*scale;
        List(j,1) = norm(abs(u)-pm.F, 'fro')/norm(pm.F, 'fro');
        List(j,2) = (norm(Dxu(:),1)+norm(Dyu(:),1))/sqrt(norm(Dxu(:)).^2+norm(Dyu(:)).^2);
        List(j,3) = sqrt(norm(Dxu(:)-dx(:))^2+norm(Dyu(:)-dy(:))^2);
        List(j,4) = sqrt(norm(Dxu(:)-hx(:))^2+norm(Dyu(:)-hy(:))^2);
        List(j,5) = norm(uold-u,'fro')/norm(uold,'fro');
        if List(j,5) < tol
            break;
        end
end
pm.i = j;
pm.L = List;
end

function d = Dx(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
end

function d = Dxt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
end

function d = Dy(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
end

function d = Dyt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
d(rows,:) = u(rows,:)-u(1,:);
end


function [tau] = mfindroot(c,eta,rho)

if c == 0 || eta == 0 
    tau = 1;
else
    a = 27*c/(rho*(eta^3)) + 2;
    C = ((a + (a^2 - 4)^0.5)/2)^(1/3);
    tau = (1 + C + 1/C)/3;
end

end


function [hx, hy,tau] = mupdate_h(c,etha,rho,d1,d2)

if etha == 0 
   hx = (c/rho)^(1/3)*ones(size(d1))/sqrt(numel(d1)*2);
   hy = (c/rho)^(1/3)*ones(size(d2))/sqrt(numel(d2)*2);
else
    a = 27*c/(rho*(etha^3)) + 2;
    C = ((a + (a^2 - 4)^0.5)/2)^(1/3);
    tau = (1 + C + 1/C)/3;
    hx = tau*d1;
    hy = tau*d2;
end

end
function z = shrink(x,r)
z = max(0,x - r) - max(0,-x - r);
end

