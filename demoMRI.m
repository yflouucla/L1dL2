clear; close all;
% demo code for MRI reconstruction




%% MRI simulation
N = 256; 
L = 6; 
F = phantom(N);
Mask = fftshift(double(MRImask(N, L)));
data = Mask.*fft2(F)/N;


pm.Num_iter = 3000;
pm.tol = 1e-5;
pm.u_orig = F;

%% L1/L2
tic
pmL1dL2.F = F;
pmL1dL2.rho1 = 1; pmL1dL2.rho2 = 1; pmL1dL2.rho3 = 1; pmL1dL2.lambda = 1000;
[u_L1dL2,pmL1dL2] = mMRrecon_L1dL2_b(Mask,data, pmL1dL2);
timeL1dL2 = toc;
fprintf('Error: %2.12f, runtime: %5.3f, L1/L2 \n',...
    norm(abs(u_L1dL2)-F, 'fro')/norm(F, 'fro'),timeL1dL2);


%% plot
figure;
imshow(abs(u_L1dL2),[]); colormap('gray');
title(['Recovery from ' num2str(L), ' radial lines']);

