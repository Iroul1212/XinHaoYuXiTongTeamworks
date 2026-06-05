% =========================================================================
% 变压器铁芯直流偏磁振动的时频分析 (基于 HHT)
% 参考文献复刻代码
% =========================================================================

clc; clear; close all;

%% 0. 运行模式设置
% strict-paper: 更贴近论文复现（默认）
% robust-analysis: 更稳健的工程分析
run_mode = getenv('HHT_RUN_MODE');
if isempty(run_mode)
    run_mode = 'strict-paper';
end
rng(2026, 'twister'); % 固定随机种子，保证结果可复现

%% 1. 参数设置与模拟信号构造
fs = 10000;           % 采样频率 (Hz)
t = 0:1/fs:0.35;      % 时间向量 (0 到 0.35s，与文献一致)

% 模拟正常工作下的振动信号 (主要为 100Hz 和 200Hz)
f1 = 100; amp1_n = 1.0;
f2 = 200; amp2_n = 0.5;
sig_normal = amp1_n*sin(2*pi*f1*t) + amp2_n*sin(2*pi*f2*t);

% 模拟直流偏磁下的振动信号 (引入奇次谐波 50Hz, 150Hz, 250Hz，且振幅增大)
f0 = 50;  amp0_dc = 0.8;
f1 = 100; amp1_dc = 2.0; % 100Hz 振幅增大
f3 = 150; amp3_dc = 0.6;
f2 = 200; amp2_dc = 1.2; % 200Hz 振幅增大
f4 = 250; amp4_dc = 0.4;
sig_dc = amp0_dc*sin(2*pi*f0*t) + amp1_dc*sin(2*pi*f1*t) + ...
         amp3_dc*sin(2*pi*f3*t) + amp2_dc*sin(2*pi*f2*t) + amp4_dc*sin(2*pi*f4*t);

% 组合信号：0.18s 为偏磁注入点
t0 = 0.18;
if strcmpi(run_mode, 'strict-paper')
    % 严格复现模式：采用阶跃注入
    vibration_sig = sig_normal;
    dc_index = t >= t0;
    vibration_sig(dc_index) = sig_dc(dc_index);
else
    % 稳健分析模式：采用短时平滑过渡，降低伪高频
    transition_dur = 0.01; % 10 ms 过渡窗
    w = zeros(size(t));
    w(t >= (t0 + transition_dur/2)) = 1;
    in_trans = (t > (t0 - transition_dur/2)) & (t < (t0 + transition_dur/2));
    w(in_trans) = 0.5 * (1 + sin(pi * (t(in_trans) - t0) / transition_dur));
    vibration_sig = (1 - w).*sig_normal + w.*sig_dc;
end

% 严格复现模式默认不额外加噪声；稳健模式加入少量噪声
if strcmpi(run_mode, 'robust-analysis')
    vibration_sig = vibration_sig + 0.1*randn(size(t));
end

%% 2. 绘制原始时域信号 (对应文献 Figure 6)
figure('Name', 'Vibration Signal', 'Position', [100, 500, 800, 300]);
plot(t, vibration_sig, 'k');
xline(0.18, 'r--', 'DC bias injected', 'LabelVerticalAlignment', 'bottom');
xlabel('time / s');
ylabel('Vibration signal (V)');
title('Simulated Transformer Core Vibration Signal (z-direction equivalent)');
xlim([0, 0.35]);

%% 3. 经验模态分解 (EMD) (对应文献 Figure 7 左侧)
% 使用 MATLAB 内置的 emd 函数
[imf, residual] = emd(vibration_sig, 'Interpolation','pchip');

% 绘制前 6 个 IMFs 和残余分量
num_imf_plot = min(6, size(imf, 2));
figure('Name', 'EMD Decomposition', 'Position', [100, 100, 800, 800]);
for k = 1:num_imf_plot
    subplot(num_imf_plot+1, 1, k);
    plot(t, imf(:, k), 'k');
    ylabel(['IMF ', num2str(k)]);
    xlim([0, 0.35]);
end
subplot(num_imf_plot+1, 1, num_imf_plot+1);
plot(t, residual, 'k');
ylabel('Residue');
xlabel('time / s');
xlim([0, 0.35]);

%% 4. Hilbert-Huang 频谱 (HHT Spectrum) (对应文献 Figure 10)
figure('Name', 'Hilbert-Huang Spectrum', 'Position', [500, 100, 800, 600]);
if strcmpi(run_mode, 'strict-paper')
    % 论文复现优先：按常见做法舍弃 IMF1 高频噪声，取 IMF2~IMF6 进行 HHT
    num_imf = size(imf, 2);
    i_start = min(2, num_imf);
    i_end = min(6, num_imf);
    imf_analyzed = imf(:, i_start:i_end);
else
    % 稳健分析：保留在目标频带中能量占比高的 IMF
    num_imf = size(imf, 2);
    band = [30, 400];
    band_ratio = zeros(1, num_imf);
    for k = 1:num_imf
        [pxx, f_psd] = periodogram(imf(:, k), [], [], fs);
        idx_band = (f_psd >= band(1)) & (f_psd <= band(2));
        e_band = trapz(f_psd(idx_band), pxx(idx_band));
        e_all = trapz(f_psd, pxx) + eps;
        band_ratio(k) = e_band / e_all;
    end
    sel = band_ratio >= 0.2;
    if ~any(sel)
        [~, idx_best] = max(band_ratio);
        sel(idx_best) = true;
    end
    imf_analyzed = imf(:, sel);
end
[hs, f_hht, t_hht] = hht(imf_analyzed, fs);

% 显式绘制 HHT 谱（当 hht 有输出参数时，MATLAB 不会自动绘图）
imagesc(t_hht, f_hht, hs);
axis xy;
xlabel('time / s');
ylabel('frequency / Hz');

% 限制频率显示范围以便于观察 0-1000Hz
ylim([0, 1000]);
title('Hilbert-Huang Spectrum');
hold on;
xline(0.18, 'w--', '0.18s', 'LineWidth', 2, 'FontSize', 14, 'Color', 'white');
hold off;
colorbar;
colormap('jet'); 

%% 5. 边际谱计算与绘制 (Marginal Spectrum) (对应文献 Figure 11)
% 边际谱是 HHT 谱在时间尺度上的积分
marginal_spec = trapz(t_hht, hs, 2);

figure('Name', 'Marginal Spectrum', 'Position', [900, 500, 600, 400]);
plot(f_hht, marginal_spec, 'k', 'LineWidth', 1.5);
xlabel('frequency / Hz');
ylabel('amplitude / V');
title('Marginal Spectrum');
xlim([0, 1000]); % 限制在 1000Hz 以内观察特征频率
grid on;

% 标出重要的频率点
hold on;
important_freqs = [50, 100, 150, 200, 250];
for i = 1:length(important_freqs)
    xline(important_freqs(i), 'r:', [num2str(important_freqs(i)) 'Hz']);
end
hold off;