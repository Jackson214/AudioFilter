%Jackson Spray

clearvars;

[x,Fs] = audioread("noisyaudio.wav");
samples = 129423; % obtained by size(x)

X = fft(x);
wx = linspace(-Fs/2,Fs/2,samples);

figure(1),subplot(2,2,1),plot(wx,(2/samples)*abs(fftshift(X))),title('DFT vs Frequency');
normDFT = X/max(X);
figure(1),subplot(2,2,2),semilogy(wx,abs(fftshift(normDFT))),title('Normalized DFT');

wp = 0.25*pi;
ws = 0.65*pi;
sp = 0.3;
ss = 0.1;

kp = 1/((1-sp)^2) - 1;
ks = 1/((ss^2)) - 1;

Omegap = wp;
Omegas = ws;

N = (1/2)*((log10(ks/kp))/(log10(Omegas/Omegap)));
N = ceil(N);

OmegaC = Omegap/(kp^(1/2*N));

B = OmegaC^3;

A = poly([OmegaC*-1 OmegaC*exp(1i*(pi+(pi/3))) OmegaC*exp(1i*(pi-(pi/3)))]);

[H,W] = freqs(B,A,100);

[R, P, K] = residue(B, A);

Hc = tf(B,A);
Hd = c2d(Hc,1);

[Num,Den] = tfdata(Hd, 'v');
syms z;
sys_syms=poly2sym(Num,z)/poly2sym(Den,z);

z = exp(1i*2*pi*(1:101)/101);
Hz = freqresp(Hd,z);

t = linspace(0,Fs/2,samples);
HaS = 20*log10(sqrt(1./(1+(t/OmegaC).^(2*N))));
figure(1),subplot(2,2,3),semilogy(t,HaS),semilogx(t,HaS),title('Logarithmic Gain of Freq. Response'),ylabel('Gain (dB)'),xlabel('Frequency (Hz)');

wc = atan(OmegaC/2);
Wn = wc/pi;

[b,a] = butter(N,Wn);
audio = filter(b,a,x);
newDFT = fft(audio);
figure(1),subplot(2,2,4),plot(wx,(2/N)*abs(fftshift(newDFT))),title('Filtered Audio DFT'),ylabel('DFT'),xlabel('Frequency (Hz)');

sound(audio, Fs);
audiowrite("filteredAudio.wav", audio, Fs);
