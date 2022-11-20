% Kac matrix eigenvalues and eigenvectors
% Jan 17, 2016 -- revision 2 of the original version.
clear all; close all;


NN = 101;                % Rings

ii1 = [1:NN].';
ii2 = [1.5:1:NN-0.5].';    % "Position" of ring-ring couplings

K = sqrt([1:NN-1].' .* [NN-1:-1:1].')/2; % standard Kac couplings expression
%= K = sqrt([NN-1:-1:1] .* [NN-1:-1:1])/2;% more local bandwidth see Sumetsky's paper
%K=ones(1,NN-1)

H = diag(-K,-1) + diag(-K,+1);
[V,D] = eig(H);

figure;
plot(ii2, [2*K -2*K], '-ob', 'LineWidth', 2);
xlabel('Position along coupled-cavity array'); ylabel('Relative frequency, \omega-\omega_o');
hold on;
DD = diag(D);
%KK = ones(1,NN-1).'*2;
plot(ii1, ones(size(ii1)) * DD.', '-r', 'LineWidth', 1.5);
hold off;
ylim([-1 1]*NN);
%legend('Ring-ring coupling','2 \mu','-2 \mu');
title(sprintf('Number of rings = %d', NN));
% Can flip eigenvectors and get a symmetric matrix about the cross-diagional (Uc = U)
% F=[1 1 -1 1 1 1 1 1 -1 1 -1    1 -1 -1 -1 -1 1 1 1 1 1 -1 1 -1];
% Vp=V*diag(F);
% figure; imagesc(Vp); axis image; colorbar; colormap(redbluehilight);
grid on;
xlim([1 NN])
