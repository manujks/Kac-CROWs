%% Kac matrix eigenvalues and eigenvectors
% Jan 17, 2016 -- revision 2 of the original version.
clearvars; close all;


%% 1. Setup
%spacing = 100E9;            % IN Hz
%c = 3e8;
%lambda0 = 1.55;                         % in um
%w = (c/lambda0)*1e6;                    % in Hz -- ugh, use f, w is used for rad/s angular frequency

NN = 51; re = 0.001;         % Configuration: NN = # resonators, re = bus coupling decay rate
%product=spacing*(NN-1);
%NN=201;
%spacing=product/(NN-1);

%K = sqrt((1:NN-1) .* (NN-1:-1:1))/2;    % standard Kac matrix coupling expression
K = ones(1,NN-1)*(NN-1)/2;              % uniform coupling distribution
%K = 1:NN-1;                             % not symmetric V matrix
%K = (1:NN-1).*(NN-1:-1:1);                 % seems symmetric
%rng(100); K = rand(NN-1); K = K(1:NN-1).*K(NN-1:-1:1);

% [xx,yy] = meshgrid(0:0.5:2*pi,0:0.5:2*pi);  % Testing how a matrix is plotted/oriented in imagesc plot
% figure; imagesc(xx); title('xx'); colorbar;
% figure; imagesc(yy); title('yy'); colorbar;


%% 2. Compute modes of lossless resonator with no bus waveguides
H = diag(-K,-1) + diag(-K,+1);          % Create system matrix
[V,D] = eig(H);
V = V*diag(sign(V(1,:)));               % Flip all modes to have +ve field on "left" -- without this sign is random
%imagesc(V-V.'); axis image; grid on;    % Attempt at more sophisticated detection of eigenvector sign (the above fails for large NN >= 201)

figure; imagesc(real(V)); axis image; colorbar; colormap(redbluehilight);       % Eigenvectors imagesc colormap
title('Lossless Kac matrix resonator mode fields (eigenvectors)');
xlabel('Mode number'); ylabel('Cavity number');
caxis([-1 1]*max(abs(real(V(:)))));

% figure; plot([1:NN].', real(V), 'LineWidth', 3); grid on;                      % Eigenvectors line plot
% title('Lossless Kac matrix resonator mode fields (eigenvectors)');
% xlabel('Cavity number'); ylabel('Amplitude (a.u.)'); legend(num2str([1:NN].'));

% figure; imagesc(log10(abs(V))); axis image; colorbar; colormap(redbluehilight);       % log10-abs Eigenvectors imagesc colormap
% title('Lossless Kac matrix resonator mode fields (eigenvectors)');
% xlabel('mode number'); ylabel('cavity number');

%return


%% 3. Compute modes with waveguide tap on ith resonator
% for i = 1:NN
%     H = diag(-K,-1) + diag(-K,+1);
%     H(i,i) = 1i*re;     % correction with the coupling sign as Milos suggested
%     [V,D]=eig(H);
% 
%     [d,ind] = sort(diag(real(D)));
%     D = D(ind,ind);
%     V = V(:,ind);
% 
%     D1(:,i) = diag(D)*spacing;      % Store results for ith structure
%     V1(:,:,i) = V;
% %    H(i,i)=0;                      % not needed, start of loop recreates H each time
% end
% 
% % 2b. Step through all supermodes, and add 
% for i = 1:NN
%      D11(i,:) = (D1(i,:));          % ?? What does this do?
%      D1(i,:)  = w + (D1(i,:));      % Just shift frequency by w
% end
%  
% % 2c. Find external Q due to bus at ith ring, and some stats of that.
% for i = 1:NN
%     Q(:,i)    = real(D1(:,i)) ./ (2*abs(imag(D1(:,i))));
%     AVG(i)    = mean(imag(D1(:,i)));
%     STD(i)    = std(imag(D1(:,i)));
%     STDAVG(i) = STD(i)/AVG(i);
% end
% 
% % Plot results
% resonator_field = 1;    % plot field assuming coupling at the given resonator
% 
% figure; imagesc(real(V1(:,:,resonator_field))); axis image; colorbar; colormap(redbluehilight);
% caxis([-1 1]);
% %title('Kac-coupling matrix eigenvectors')
% xlabel('mode number'); ylabel('cavity number');
% title({'Kac-coupling matrix mode fields(eigenvectors) for',[' resonator coupling at cavity #=',num2str(resonator_field)],'(sign normalization with field at ','first resonator in each supermode)'})
% 
% figure; plot([real(D1)*1e-12],'-o');
% %xlabel('Mode number'); ylabel('Resonance frequency, \omega (in THz)');
% xlabel('Mode number'); ylabel('Resonance frequency (THz)');
% legend('Kac matrix (equispaced)');
% 
% nmode = (1:1:NN);
% figure; imagesc(nmode, nmode, log10(abs(Q)));
% axis image; colorbar; colormap(hot); %caxis([0 6]);
% xlabel('Resonator # coupled to bus waveguide'); ylabel('Mode number');
% title({'Plot of log_{10}Q with the given mode number for ',' a given resonator # coupled to bus waveguide',['for re=',num2str(re*1e-9),' GHz and for no. of resonators N=',num2str(NN) ],['for re=',num2str(re*spacing*1e-9),' GHz for ','comb spacing=',num2str(spacing*1e-9),' Ghz (',num2str(re),' for unity comb spacing)']})
% 


%% 4. Time evolution of pulses

% Set spectrum of field in comb
nvec = [-(NN-1)/2:(NN-1)/2].';
nvecp = nvec/(NN-1);
%qvec0 = exp(-(nvecp/0.25).^2) .* exp(j * nvecp * pi*(NN-1)/2);  % Set the original spectrum (complex amplitude vs. mode number)
%qvec0 = diag(eye(NN)).* exp(j * nvecp * pi*(NN-1)/2);
%qvec0 = hamming(NN).* exp(j * nvecp * pi*(NN-1)/2);
qvec0 = tukeywin(NN,0.65).* exp(j * nvecp * pi*(NN-1)/2);
w = 0; spacing = 1;
VV = V; DD = w*eye(NN) + D*spacing;     % [MP] Scale and shift eigenvalues

avec0 = VV' * qvec0;                    % Find the spatial amplitude distribution (vs. cavity number i.e. "z")

%figure; plot(1:NN, abs(qvec0).^2, '-o', 'LineWidth', 2);
figure; stem(1:NN, abs(qvec0).^2, '-o', 'LineWidth', 2);
xlim([ 0 NN+1 ]); ylim([0 1.5*max(abs(qvec0).^2)]);
xlabel('Mode number'); ylabel('Spectral power'); grid on;
figure; plot((1:NN).', [real(avec0) abs(avec0)], '-o', 'LineWidth', 2);
xlabel('Cavity number'); ylabel('Power'); grid on;


dwot = [0:0.05/5:1.5]*2*pi;               % Define time evolution for each supermode
TT = exp(j * nvec*dwot);                % 2pi between adjacent resonators is the full cycle

qt = TT.*(qvec0*ones(size(dwot)));      % Time-evolve the supermodes
at = VV' * qt;                          % Convert to cavity basis

figure;
hw = waterfall((1:NN), dwot/(2*pi), abs(at.')); view([30 45]);
set(hw,'EdgeColor',[0 0 0]); %set(w,'FaceAlpha',0)
title({['Pulse propagation for no. of resonators N = ', num2str(NN)]});
xlabel('Cavity number'); ylabel('Normalized time, \delta\omega_o t/2\pi (round trips)'); zlabel('Cavity amplitude (a.u.)');
xlim([1 NN]); ylim auto;

figure;
plot((1:NN), abs(at(:,1:end)), 'LineWidth', 2); grid on;
title({['Pulse propagation for no. of resonators N = ', num2str(NN)]});
xlabel('Cavity number'); ylabel('Cavity amplitude (a.u.)');
xlim([1 NN]);

figure;
hw = waterfall(dwot/(2*pi), (1:NN), abs(at)); view([30+90 45]);
set(hw,'EdgeColor',[0 0 0]); %set(w,'FaceAlpha',0)
title({['Pulse propagation for no. of resonators N = ', num2str(NN)]});
ylabel('Cavity number'); xlabel('Normalized time, \delta\omega_o t/2\pi (round trips)'); zlabel('Cavity amplitude (a.u.)');
ylim([1 NN]); xlim auto;

figure;
plot(dwot/(2*pi), abs(at(1:2:end,:).'), 'LineWidth', 2); grid on;
title({['Pulse propagation for no. of resonators N = ', num2str(NN)]});
xlabel('Normalized time, \delta\omega_o t/2\pi (round trips)'); ylabel('Cavity amplitude (a.u.)');


n=1:1:NN;
for k = 1:length(dwot)
%    q(:,k) = VV*expm(1i*DD*t(k))*ao;
%    qt(:,k)
    n_effective(k)= sum(n(:).*abs(at(:,k)).^2) ./ sum(abs(at(:,k)).^2);
    std(k)= sqrt(sum(n(:).^2.*abs(at(:,k)).^2) ./ sum(abs(at(:,k)).^2)-(sum(n(:).*abs(at(:,k)).^2) ./ sum(abs(at(:,k)).^2)).^2);
end

figure; imagesc(dwot/(2*pi), n, abs(at));
ylabel('Cavity number'); xlabel('Normalized time, \delta\omega t_o/2\pi (round trips)');
title('Pulse evolution in coupled-cavity system');
colorbar; colormap(hot);

% figure; w = waterfall([1:1:NN],t*1E12/(2*pi),abs(q.'));
% view([30 45])
% set(w,'EdgeColor',[0 0 0]);
% %set(w,'FaceAlpha',0)
% title({['Pulse propagation for no. of resonators N=',num2str(NN)]});
% xlabel('cavity #'); ylabel('time(in picosecond)'); zlabel('amplitude');


% n=1:1:NN;
% for k = 1:size(qt,2)
%     %q(:,k) = VV*expm(1i*DD*t(k))*ao;
%     n_effective(k)= sum(n(:).*abs(qt(:,k)).^2) / sum(abs(qt(:,k)).^2);
%     std(k)= sqrt(sum(n(:).^2.*abs(qt(:,k)).^2) ./ sum(abs(qt(:,k)).^2)-(sum(n(:).*abs(qt(:,k)).^2) ./ sum(abs(qt(:,k)).^2)).^2);
% end
% 
% figure; w = waterfall(dwot/(2*pi),[1:NN],abs(qt));
% view([30-270 45])
%  view(2)
%  set(w,'EdgeColor',[0 0 0]);
%  set(w,'FaceAlpha',0)
%  title({['Pulse propagation(mean position) for no. of resonators N=',num2str(NN)],['(top view) z axis is the view axis'],['(n_effective = sum(n |E_n|^2, n = 1..11) / sum(|E_n|^2).]']} );
%  ylabel('cavity #'); xlabel('time(IN Picosecond)'); zlabel('amplitude');

%hold on;
figure;
plot(dwot.'/(2*pi), n_effective.', 'b', 'LineWidth', 2)
hold on;
plot(dwot.'/(2*pi), [n_effective.'+std.' n_effective.'-std.'], 'r', 'LineWidth', 2)
xlabel('Normalized time, \delta\omega_o t/2\pi (round trips)');
ylabel({['pulse width or standard deviation'],['(in units of discretized cavity lengths)']})
title({['standard deviation(or pulse width) for no. of resonators N=',num2str(NN)]} )
ylim([min(n) max(n)]);
grid on;



return


qo = zeros(NN,1); qo(1) = 1;
%x = [0.5:NN-0.5].'; ao = exp(-(x/3).^2); ao = ao/sqrt(ao'*ao);
x = [0.5:NN-0.5].'; qo = exp(-(x/3).^2); qo = qo/sqrt(qo'*qo);     % Gaussian pulse centered at edge ring
%x = [0.5:NN-0.5].'; qo = exp(-(x-((NN+1)/2)).^2); qo = qo/sqrt(qo'*qo);
%qo=V*ao;
x = [0.5:NN-0.5].'; qo = exp(-(x/3).^2); qo = qo/sqrt(qo'*qo);
% VV = V1(:,:,resonator_field);
% DD = diag(D11(:,resonator_field));
VV = V; DD = w*eye(NN) + D*spacing;     % [MP] Scale and shift eigenvalues

ao = (VV)'*qo;
t = [0:0.05:1.5]*(1/spacing)*2*pi;
n=1:1:NN;
for k = 1:length(t)
    %q(:,k) = VV*expm(1i*DD*t(k))*ao;
    n_effective(k)= sum(n(:).*abs(q(:,k)).^2) ./ sum(abs(q(:,k)).^2);
    std(k)= sqrt(sum(n(:).^2.*abs(q(:,k)).^2) ./ sum(abs(q(:,k)).^2)-(sum(n(:).*abs(q(:,k)).^2) ./ sum(abs(q(:,k)).^2)).^2);
end

%figure; imagesc(abs(q));

figure; w = waterfall([1:1:NN],t*1E12/(2*pi),abs(q.'));
view([30 45])
set(w,'EdgeColor',[0 0 0]);
%set(w,'FaceAlpha',0)
title({['Pulse propagation for no. of resonators N=',num2str(NN)]});
xlabel('cavity #'); ylabel('time(in picosecond)'); zlabel('amplitude');



figure; w = waterfall(t./(2*pi)*1E12,[1:NN],abs(q));
%view([30-270 45])
view([-45 45])
set(w,'EdgeColor',[0 0 0]);
%set(w,'FaceAlpha',0)
title({['Pulse propagation for no. of resonators N=',num2str(NN)]} );
ylabel('cavity #'); xlabel('time(IN Picosecond)'); zlabel('amplitude');


figure; w = waterfall(t./(2*pi)*1E12,[1:NN],abs(q));
%view([30-270 45])
view(2)
set(w,'EdgeColor',[0 0 0]);
%set(w,'FaceAlpha',0)
title({['Pulse propagation(mean position) for no. of resonators N=',num2str(NN)],['(top view) z axis is the view axis'],['(n_effective = sum(n |E_n|^2, n = 1..11) / sum(|E_n|^2).]']} );
ylabel('cavity #'); xlabel('time(IN Picosecond)'); zlabel('amplitude');

hold on;
plot(t./(2*pi)*1E12,n_effective)
figure;
plot(t./(2*pi)*1E12,std)
xlabel('time(IN Picosecond)');
ylabel({['pulse width or standard deviation'],['(in units of discretized cavity lengths)']})
title({['standard deviation(or pulse width) for no. of resonators N=',num2str(NN)]} )
figure; plot([1:NN],abs(ao));
xlabel('mode #'); ylabel('mode amplitude (a.u.)'); title('Resonator mode spectrum for incident mode field');

figure; plot([1:NN],abs(qo));
xlabel('mode #'); ylabel('CAVITY #'); title('Resonator cavities field for incident mode field');
