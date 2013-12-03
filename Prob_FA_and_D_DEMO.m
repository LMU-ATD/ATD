% Probability of False Alarm and Detection Demo
%   MATLAB m-file
%   Author: Terry Kong
%   Variables: 
%       A = constant signal, 
%       vt = threshold voltage, 
%       ymin = plot lower bound, 
%       ymax = plot upper bound, 
%       orange = colormap for orange, 
%       n = plot resolution, 

% Pd and Pfa Determined Numerically by running a time domain simulation 
% Generate noise and signals for SNR=13dB 
clear all; clc; orange = [1,0.5,0.2];
myf = figure;
set(myf,'NumberTitle','off','Name','Signal(w/o noise): I(t) = 5cos(t), Q(t) = 5sin(t), s(t) = 5') 

A = 5;
ymin = -0.2;
ymax = 0.4;
% Determine the probabilities for different thresholds 
vt=3; vty = linspace(ymin,ymax,100); vtx = zeros(size(vty));
n = 100000;

a = (1:n); 
x = randn(size(a)); 
sigi = A*sin(a/1000); 
y = randn(size(a)); 
sigq = A*cos(a/1000);
% Determine the envelope 
c = sqrt(x.^2+y.^2); 
csig = sqrt((x+sigi).*(x+sigi)+(y+sigq).*(y+sigq)); 
%Plot the distributions 
[f1,x1] = ksdensity(c); 
[f2,x2] = ksdensity(csig);

subplot(3,1,1);hold on;grid on;title('PDF of Noise only')
h1 = plot(x1,f1,'r');
set(h1,'LineWidth',1);
I = find(x1 > vt);
h11 = area(x1(I),f1(I),'FaceColor','r'); legend(h11,sprintf('Probability of FA: P(Amplitude > %0.1f = Vt)',vt));
xlim([-2,12]); ylim([ymin,0.8]); xlabel('Amplitude(V)'); ylabel('Probability')
line([vt,vt],[ymin,ymax],'LineStyle',':','LineWidth',2,'Color',orange)
text(2.05,0.5,'Detection Threshold','FontWeight','bold','FontSize',12,'Color',orange)

subplot(3,1,2);hold on;grid on;title(sprintf('PDF of Noise + Signal(s(t) = %0.1f)',A))
h2 = plot(x2,f2,'b');
set(h2,'LineWidth',1);
I = find(x2 > vt);
h22 = area(x2(I),f2(I),'FaceColor','b'); legend(h22,sprintf('Probability of Detection: P(Amplitude > %0.1f = Vt)',vt));
xlim([-2,12]); ylim([ymin,0.8]); xlabel('Amplitude(V)'); ylabel('Probability')
line([vt,vt],[ymin,ymax],'LineStyle',':','LineWidth',2,'Color',orange)
text(2.05,0.5,'Detection Threshold','FontWeight','bold','FontSize',12,'Color',orange)

subplot(3,1,3);hold on;grid on;
h3 = plot(x1,f1,'r',x2,f2,'b');
set(h3,'LineWidth',1)
legend('Noise only',sprintf('Noise + Signal(s(t) = %0.1f)',A))
xlim([-2,12]); ylim([ymin,0.8]); xlabel('Amplitude(V)'); ylabel('Probability')
line([vt,vt],[ymin,ymax],'LineStyle',':','LineWidth',2,'Color',orange)
text(2.05,0.5,'Detection Threshold','FontWeight','bold','FontSize',12,'Color',orange)


I = find(x2 > vt);
h1 = area(x2(I),f2(I),'FaceColor','b');
I = find(x1 > vt);
h2 = area(x1(I),f1(I),'FaceColor','r');

% Look for noise peaks above the threshold 
nfa = find(c>vt); 
pfa = length(nfa)/n
% Look for S+N peaks above the threshold 
nd = find(csig>vt); 
pd = length(nd)/n
