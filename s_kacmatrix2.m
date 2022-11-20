% Kac matrix eigenvalues and eigenvectors
% Jan 17, 2016 -- revision 2 of the original version.
clear all; close all;


NN =11; 
K = sqrt([1:NN-1] .* [NN-1:-1:1])/2; % standard kadc expression
% K = sqrt([NN-1:-1:1] .* [NN-1:-1:1])/2;% more local bandwidth see Sumetsky's paper
H = diag(-K,-1)+diag(-K,+1);
figure; imagesc(H); axis image; colorbar; colormap(hot);
title('Kac-coupling matrix elements')
[V,D]=eig(H);
figure; imagesc(V); axis image; colorbar; colormap(redbluehilight);
caxis([-1 1]);
title('Kac-coupling matrix eigenvectors')

% Can flip eigenvectors and get a symmetric matrix about the cross-diagional (Uc = U)
% F=[1 1 -1 1 1 1 1 1 -1 1 -1    1 -1 -1 -1 -1 1 1 1 1 1 -1 1 -1];
% Vp=V*diag(F);
% figure; imagesc(Vp); axis image; colorbar; colormap(redbluehilight);


K = (NN)/2*ones(NN-1,1)/2; Ho = diag(-K,-1)+diag(-K,+1);
figure; imagesc(Ho); axis image; colorbar; colormap(hot);
title('Constant-coupling matrix elements')
[Vo,Do]=eig(Ho);
figure; imagesc(Vo); axis image; colorbar; colormap(redbluehilight);
caxis([-0.5 0.5]);
title('Constant-coupling matrix eigenvectors')

figure; plot([diag(Do) diag(D)],'-o');
xlabel('Mode number'); ylabel('Relative resonance frequency, \omega-\omega_o');
legend('Constant coupling (band-edge bunching)','Kac matrix (equispaced)');


% Circular super-resonator
K = sqrt([1:NN-1] .* [NN-1:-1:1])/2; H2 = diag(-K,-1)+diag(-K,+1); H2(1,end) = -NN/128;H2(end,1) = -NN/128;
title('Circular-constant coupling matrix elements')
[V2,D2]=eig(H2);
figure; imagesc(Vo); axis image; colorbar; colormap(redbluehilight);
caxis([-0.5 0.5]);
title('Circular-constant coupling matrix eigenvectors')

figure; plot([diag(Do) diag(D) diag(D2)],'-o');
xlabel('Mode number'); ylabel('Relative resonance frequency, \omega-\omega_o');
legend('Constant coupling (band-edge bunching)','Kac matrix (equispaced)','Circular constant');


% Decompose and propagate pulse

qo = zeros(NN,1); qo(1) = 1;
x = [0.5:NN-0.5].'; qo = exp(-(x/3).^2); qo = qo/sqrt(qo'*qo);     % Gaussian pulse centered at edge ring

VV = V;
DD = D;
ao = VV'*qo;
t = [0:0.05:1.25]*2*pi;
for k = 1:length(t)
    q(:,k) = VV*expm(i*DD*t(k))*ao;
end

%figure; imagesc(abs(q));

figure; w = waterfall(abs(q.'));
view([30 45])
set(w,'EdgeColor',[0 0 0]);
%set(w,'FaceAlpha',0)
title('Pulse propagation');
xlabel('cavity #'); ylabel('time'); zlabel('amplitude');

figure; plot([1:NN],abs(ao));
xlabel('mode #'); ylabel('mode amplitude (a.u.)'); title('Resonator spectrum');

figure; w = waterfall(abs(q));
%view([30-270 45])
view([-45 45])
set(w,'EdgeColor',[0 0 0]);
%set(w,'FaceAlpha',0)
title('Pulse propagation');
ylabel('cavity #'); xlabel('time'); zlabel('amplitude');
