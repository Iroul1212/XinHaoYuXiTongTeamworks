% =========================================================================
% Transformer winding fault diagnosis via impulsive sound (paper-oriented)
% Paper: Diagnosis on Winding Failure Through Impulsive Sound of Power Transformer
% Method summary:
%   1) Complex wavelet time-frequency analysis (for visualization)
%   2) 3-level wavelet packet decomposition (db5)
%   3) Fault features: Sound Entropy (SE) and Odd-Even Ratio (OER)
% =========================================================================

clc;
clear;
close all;

%% 0) Configuration
fs = 1024;           % Paper setting: 1024 S/s
T = 1.0;             % Paper setting: 1 second record
waveletName = 'db5'; % Paper setting
wpLevel = 3;         % Paper setting: 3-layer WPD

% Data mode:
%   'synthetic' : run demo with generated impulsive sounds
%   'matfile'   : load user data from MAT file
dataMode = getenv('IMPULSIVE_SOUND_DATA_MODE');
if isempty(dataMode)
    dataMode = 'synthetic';
end

% If dataMode = 'matfile', MAT must include one of:
%   A) signalMatrix [N x M], each column is one test sample
%   B) singleSignal [N x 1], one sample only
matFilePath = '';

t = (0:1/fs:T-1/fs).';

%% 1) Load or create signals
if strcmpi(dataMode, 'matfile')
    if isempty(matFilePath) || ~isfile(matFilePath)
        error('Please set matFilePath to a valid MAT file when dataMode = matfile.');
    end

    S = load(matFilePath);
    if isfield(S, 'signalMatrix')
        signalMatrix = S.signalMatrix;
    elseif isfield(S, 'singleSignal')
        signalMatrix = S.singleSignal(:);
    else
        error('MAT file must contain signalMatrix or singleSignal.');
    end

    if size(signalMatrix, 1) < numel(t)
        error('Signal length is shorter than 1 second at fs=1024.');
    end
    signalMatrix = signalMatrix(1:numel(t), :);
else
    % Synthetic demo to mimic paper trend across repeated short-circuit tests:
    % early stage: looseness -> SE and OER increase
    % later stage: deformation -> possible SE decrease and OER variation
    rng(2026, 'twister');
    nTests = 30;
    signalMatrix = zeros(numel(t), nTests);

    for k = 1:nTests
        loosenessFactor = min(1, (k-1)/14);   % increases in first half
        deformFactor = max(0, (k-16)/14);     % increases in later stage

        % Baseline components around 50 Hz harmonics.
        a50  = 1.00 + 0.65*loosenessFactor - 0.20*deformFactor;
        a100 = 0.70 + 0.30*loosenessFactor - 0.15*deformFactor;
        a150 = 0.30 + 0.40*loosenessFactor + 0.05*deformFactor;
        a200 = 0.20 + 0.20*loosenessFactor - 0.08*deformFactor;
        a250 = 0.12 + 0.22*loosenessFactor + 0.03*deformFactor;

        % Impulsive envelope: rapid onset + decaying free vibration.
        env = (1 - exp(-50*t)) .* exp(-4.5*t);

        sig = env .* ( ...
            a50 *sin(2*pi*50*t)  + ...
            a100*sin(2*pi*100*t) + ...
            a150*sin(2*pi*150*t) + ...
            a200*sin(2*pi*200*t) + ...
            a250*sin(2*pi*250*t) );

        % Add weak wideband disturbance.
        sig = sig + 0.02*randn(size(sig));

        signalMatrix(:, k) = sig;
    end
end

nTests = size(signalMatrix, 2);

%% 2) Time-frequency analysis with complex wavelet (visualization)
% Show the first sample as an example (paper Fig.4/Fig.5 style).
exampleSig = signalMatrix(:, 1);

figure('Name', 'Impulsive sound signal (example)', 'Position', [100 100 900 280]);
plot(t, exampleSig, 'k');
xlabel('Time (s)');
ylabel('Amplitude');
title('Impulsive sound (example sample)');
xlim([0 T]);

figure('Name', 'Time-frequency map by complex wavelet', 'Position', [100 420 920 420]);
[cfs, f] = cwt(exampleSig, fs, 'amor');
imagesc(t, f, abs(cfs));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Complex wavelet time-frequency map');
ylim([0 400]);
colormap(jet);
colorbar;

%% 3) Wavelet packet features: SE and OER
SE = zeros(nTests, 1);
OER = zeros(nTests, 1);

% Paper node order in frequency sequence at level-3:
% 7, 8, 10, 9, 13, 14, 12, 11
freqOrderedNodes = [7 8 10 9 13 14 12 11];

% Odd/Even groups in paper:
oddNodes = [7 10 13 12];
evenNodes = [8 9 14 11];

for k = 1:nTests
    x = signalMatrix(:, k);

    % Build wavelet packet tree.
    Tnode = wpdec(x, wpLevel, waveletName);

    % Node energies in the exact order used by paper.
    E = zeros(size(freqOrderedNodes));
    for i = 1:numel(freqOrderedNodes)
        nodeIdx = freqOrderedNodes(i);
        coef = wpcoef(Tnode, nodeIdx);
        E(i) = sum(coef.^2);
    end

    Eall = sum(E) + eps;
    p = E / Eall;

    % Sound Entropy (SE): -sum(p_i * log2(p_i)).
    SE(k) = -sum(p .* log2(p + eps));

    % Odd-Even Ratio (OER): odd band energy ratio over even band ratio.
    Eodd = 0;
    Eeven = 0;
    for i = 1:numel(freqOrderedNodes)
        nodeIdx = freqOrderedNodes(i);
        if any(nodeIdx == oddNodes)
            Eodd = Eodd + E(i);
        elseif any(nodeIdx == evenNodes)
            Eeven = Eeven + E(i);
        end
    end
    OER(k) = (Eodd / Eall) / max(Eeven / Eall, eps);
end

%% 4) Feature trend and diagnosis logic
figure('Name', 'SE and OER across tests', 'Position', [1050 100 750 460]);
subplot(2,1,1);
plot(1:nTests, SE, 'k-o', 'LineWidth', 1.2, 'MarkerSize', 4);
grid on;
xlabel('Test index');
ylabel('SE');
title('Sound Entropy (SE) trend');

subplot(2,1,2);
plot(1:nTests, OER, 'k-s', 'LineWidth', 1.2, 'MarkerSize', 4);
grid on;
xlabel('Test index');
ylabel('OER');
title('Odd-Even Ratio (OER) trend');

% Relative change against baseline (first sample)
SE_base = SE(1);
OER_base = OER(1);
SE_dev = (SE - SE_base) / max(abs(SE_base), eps) * 100;
OER_dev = (OER - OER_base) / max(abs(OER_base), eps) * 100;

fprintf('\n=== Feature deviation relative to sample #1 ===\n');
fprintf('Index\tSE_dev(%%)\tOER_dev(%%)\tDiagnosis\n');
for k = 1:nTests
    diagnosis = localDiagnosis(SE_dev(k), OER_dev(k));
    fprintf('%d\t%8.2f\t%9.2f\t%s\n', k, SE_dev(k), OER_dev(k), diagnosis);
end

%% 5) Optional: inspect level-3 band energies of one sample
inspectIdx = min(1, nTests);
x = signalMatrix(:, inspectIdx);
Tnode = wpdec(x, wpLevel, waveletName);
E = zeros(size(freqOrderedNodes));
for i = 1:numel(freqOrderedNodes)
    E(i) = sum(wpcoef(Tnode, freqOrderedNodes(i)).^2);
end

figure('Name', 'Level-3 WPD band energies (example)', 'Position', [1050 620 750 300]);
bar(E / sum(E + eps));
grid on;
xticklabels(cellstr(num2str(freqOrderedNodes(:))));
xlabel('Terminal node index');
ylabel('Normalized band energy');
title('3-level WPD band energy distribution (db5)');


function label = localDiagnosis(seDevPct, oerDevPct)
% Paper-inspired rule:
%   SE and OER increase -> looseness tendency
%   SE and OER decrease -> deformation tendency
%   opposite trend -> mixed or transition state
th = 5; % simple tolerance for weak fluctuation

seUp = seDevPct > th;
seDown = seDevPct < -th;
oerUp = oerDevPct > th;
oerDown = oerDevPct < -th;

if seUp && oerUp
    label = 'Looseness tendency';
elseif seDown && oerDown
    label = 'Deformation tendency';
else
    label = 'Mixed/transition state';
end
end
