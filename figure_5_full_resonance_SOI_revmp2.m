% Simulated linear transmission through passive Kac comb
% Mar 27, 2021
clear all; close all;

%% 1. Plotting settings
figurePosition = [325 120 825 550];

greenDark = [0.0 0.4 0.0];  green     = [0.0 0.6 0.2];  % Plot colors
blueDark  = [0.0 0.4 0.6];  blue      = [0.0 0.5 0.7];
redDark   = [0.6 0.0 0.0];  red       = [0.8 0.0 0.0];
yellow    = [1 0 0];


%% 2. Settings
% 2a. General
c = 299792458;

% 2b. Resonator
dfspacing = 100;                        % Comb spacing in GHz
Qo = 1e6;                               % Ring Qo
omega0 = 2*pi*c/1550e-9 / 1e9;          % in Grad/s
ro = omega0/(2*Qo);                     % in Grad/s
domega = [-20:0.0005:20] * 2*pi*dfspacing;

NNset = [8 16 24 32];                   % Set of combs to plot

figure;

for nn = 1:length(NNset)
    % Kac comb parameters
    NN = NNset(nn);

    % Kac comb CMT model
    mu1 = sqrt([1:NN-1].*(NN-[1:NN-1])) / 2 * (2*pi*dfspacing);   % Coupling vector for Kac resonator
    % First find supermode re's for lossless Kac resonator
    Ho = diag(-mu1,-1) + diag(-mu1,+1); % + diag(omega0*ones(1,NN)); % + 1j*diag(ro*ones(1,NN));

    [V,D] = eig(Ho);
    wsupermodes = diag(D); %resvec = imag(wsupermodes); resvecavg = mean(resvec);

    resovre = abs(V).^2/2*2;              % Supermode external Qs (along row) for bus coupled to ith ring (row i)
%    resovreavg = mean(resovre,2);       % Find average supermode external Q for each i choice of ring coupled to bus
    resovremax = max(resovre,[],2);
    resovremin = min(resovre,[],2);
    resovrelogspan = log10(resovremax) - log10(resovremin);
    resovregeomean = 10.^((log10(resovremax) + log10(resovremin))/2);    % Find average supermode external Q for each i choice of ring coupled to bus
    [resovremin resovremax resovrelogspan 10.^resovrelogspan resovregeomean];

    % So, which is the best and worst ring to couple to?  (decide by min or max log span of supermode decay rates)
    idxbestring(nn) = find(resovrelogspan == min(resovrelogspan),1);     % Find only the first min and max
    idxworstring(nn) = find(resovrelogspan == max(resovrelogspan),1);

%     re = 2*pi*10;
%     ro = 2*pi*10*0.01;
%     Q_calc = omega0/(2*ro)            % in Grads
%     pin   = 3;
%     pdrop = NN-pin+1;

    % So, let's set the input and drop to the best ring
    pin = idxbestring(nn)      % Best ring
%    pin = idxworstring(nn);     % Worst ring
%    pin = 1;                    % First ring
 %  pin = bus(nn);     % Best ring
    TsqdB(nn) = 10*log10( ((1 - 10^(resovrelogspan(pin)/2))/(1 + 10^(resovrelogspan(pin)/2)))^2 );

%     % Simulate exact Q to verify
%     retest = 0.001*ro;
%     dHo = zeros(size(Ho)); dHo(pin,pin) = j*retest;
%     [V,D] = eig(Ho+dHo);
%     [wsupermodes,ix] = sort(diag(D),'ComparisonMethod','real');  % Sort the modes by real
%     II = eye(size(D)); P = II(ix,:);    % Permutation matrix that sorts the eigenvalues (inv(P) = P.', so P*P.'=I)
%     %Dp = P*D*inv(P); Dp = P*D*P.';
%     D = P*D*P.'; V = V*P.';             % Permute the eigenvectors accordingly.


    %    pin = idxworstring(nn);
    resovremeantarget = resovregeomean(pin);
%    re = 41.224/2/pi; %7; %ro / resovremeantarget;
    re = ro ./ resovremeantarget; % /2/pi;
%    re = ro / (0.02734375*2); m = 2;

    % Now find lossy resonator
    H = Ho + 1j * diag(ro*ones(1,NN));
    H(pin,pin)     = H(pin,pin)   + 1j * re;
%    H(pdrop,pdrop) = H(pdrop,pdrop) + 1j * re;

%     R = eye(2);                                     % 2x2 identity matrix
    R = 1;

%    Mi = zeros(NN,2);
    Mi = zeros(NN,1);
    Mi(pin,1)   = sqrt(2*re);         % Set input coupling for add port (1) to desired resonator number
%    Mi(pdrop,2) = sqrt(2*re);         % Set input coupling for drop port (2) to desired resonator number

    Mo = Mi.';

    for kk = 1:length(domega)
    %    S(:,:,kk) = R + 1j*Mo * inv( omega(kk)*eye(NN) - H ) * Mi;  % Compute S matrix at each frequency
        S(:,:,kk) = R + 1j*Mo * (( domega(kk)*eye(NN) - H ) \ Mi);  % Compute S matrix at each frequency
    %    trans(kk)=1-(abs(S(1,1,kk))).^2;
    end

    thru = abs( squeeze(S(1,1,:)) ).^2;
%    drop = abs( squeeze(S(2,1,:)) ).^2;

    subplot(4,1,nn)
%    plot(domega/(2*pi), 10*log10(thru), '-o', 'Color', redDark, 'LineWidth', 0.75);
    plot(domega/(2*pi), 10*log10(thru), '-', 'Color', redDark, 'LineWidth', 0.75);
    ylabel('Transmission (dB)');
    xlabel('Frequency detuning from isolated ring resonance, (\omega-\omega_{0})/2\pi (GHz)');
    xlim([-2000 2000]); ylim([-50 2]); grid on;
end

figure;

for nn = 1:length(NNset)
    % Kac comb parameters
    NN = NNset(nn);

    % Kac comb CMT model
    mu1 = sqrt([1:NN-1].*(NN-[1:NN-1])) / 2 * (2*pi*dfspacing);   % Coupling vector for Kac resonator
    % First find supermode re's for lossless Kac resonator
    Ho = diag(-mu1,-1) + diag(-mu1,+1); % + diag(omega0*ones(1,NN)); % + 1j*diag(ro*ones(1,NN));

    [V,D] = eig(Ho);
    wsupermodes = diag(D); %resvec = imag(wsupermodes); resvecavg = mean(resvec);

    resovre = abs(V).^2/2*2;              % Supermode external Qs (along row) for bus coupled to ith ring (row i)
%    resovreavg = mean(resovre,2);       % Find average supermode external Q for each i choice of ring coupled to bus
    resovremax = max(resovre,[],2);
    resovremin = min(resovre,[],2);
    resovrelogspan = log10(resovremax) - log10(resovremin);
    resovregeomean = 10.^((log10(resovremax) + log10(resovremin))/2);    % Find average supermode external Q for each i choice of ring coupled to bus
    [resovremin resovremax resovrelogspan 10.^resovrelogspan resovregeomean];

    % So, which is the best and worst ring to couple to?  (decide by min or max log span of supermode decay rates)
    idxbestring(nn) = find(resovrelogspan == min(resovrelogspan),1);     % Find only the first min and max
    idxworstring(nn) = find(resovrelogspan == max(resovrelogspan),1);

%     re = 2*pi*10;
%     ro = 2*pi*10*0.01;
%     Q_calc = omega0/(2*ro)            % in Grads
%     pin   = 3;
%     pdrop = NN-pin+1;

    % So, let's set the input and drop to the best ring
    pin = idxbestring(nn)      % Best ring
%    pin = idxworstring(nn);     % Worst ring
%    pin = 1;                    % First ring
 %  pin = bus(nn);     % Best ring
    TsqdB(nn) = 10*log10( ((1 - 10^(resovrelogspan(pin)/2))/(1 + 10^(resovrelogspan(pin)/2)))^2 );

%     % Simulate exact Q to verify
%     retest = 0.001*ro;
%     dHo = zeros(size(Ho)); dHo(pin,pin) = j*retest;
%     [V,D] = eig(Ho+dHo);
%     [wsupermodes,ix] = sort(diag(D),'ComparisonMethod','real');  % Sort the modes by real
%     II = eye(size(D)); P = II(ix,:);    % Permutation matrix that sorts the eigenvalues (inv(P) = P.', so P*P.'=I)
%     %Dp = P*D*inv(P); Dp = P*D*P.';
%     D = P*D*P.'; V = V*P.';             % Permute the eigenvectors accordingly.


    %    pin = idxworstring(nn);
    resovremeantarget = resovregeomean(pin);
%    re = 41.224/2/pi; %7; %ro / resovremeantarget;
    re = ro ./ resovremeantarget; % /2/pi;
%    re = ro / (0.02734375*2); m = 2;

    % Now find lossy resonator
    H = Ho + 1j * diag(ro*ones(1,NN));
    H(pin,pin)     = H(pin,pin)   + 1j * re;
%    H(pdrop,pdrop) = H(pdrop,pdrop) + 1j * re;

%     R = eye(2);                                     % 2x2 identity matrix
    R = 1;

%    Mi = zeros(NN,2);
    Mi = zeros(NN,1);
    Mi(pin,1)   = sqrt(2*re);         % Set input coupling for add port (1) to desired resonator number
%    Mi(pdrop,2) = sqrt(2*re);         % Set input coupling for drop port (2) to desired resonator number

    Mo = Mi.';

    for kk = 1:length(domega)
    %    S(:,:,kk) = R + 1j*Mo * inv( omega(kk)*eye(NN) - H ) * Mi;  % Compute S matrix at each frequency
        S(:,:,kk) = R + 1j*Mo * (( domega(kk)*eye(NN) - H ) \ Mi);  % Compute S matrix at each frequency
    %    trans(kk)=1-(abs(S(1,1,kk))).^2;
    end

    
%    drop = abs( squeeze(S(2,1,:)) ).^2;
   res_ovr_ro=resovre(:,pin)./ resovremeantarget;
    subplot(4,1,nn)
    
    T=[(1-10.^log10(res_ovr_ro))./(1+10.^log10(res_ovr_ro))].^2;
%    plot(domega/(2*pi), 10*log10(thru), '-o', 'Color', redDark, 'LineWidth', 0.75);
    plot(log10(res_ovr_ro), 10*log10(T), '*', 'Color', redDark, 'LineWidth', 0.75);
    ylabel('Transmission (dB)');
    xlabel('Frequency detuning from isolated ring resonance, (\omega-\omega_{0})/2\pi (GHz)');
    xlim([-3 3]); ylim([-50 2]); grid on;
end