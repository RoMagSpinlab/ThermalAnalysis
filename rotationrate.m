%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation rate calulation 2.0
% 
% Yufan Xu
% 
% This routine calculates the rotation rate of romag from the background 
% magnetic field data. The sinFit funtion is by Peter Seibold
% (https://www.mathworks.com/matlabcentral/fileexchange/66793-sine-fitting)
% 
% Matlab Version: R2019b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mise en place
clear all
close all
clc

%% Data acquisition
filename = uigetfile('*.txt');
fileID = fopen(filename,'r');
startRow = 4;
formatSpec = '%11s%11s%11s%11s%s%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', ...
    'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, ...
    'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Assign arrays
time = cell2mat(raw(:, 1));
%Bz = cell2mat(raw(:, 2));
Bx = cell2mat(raw(:, 2));
By = cell2mat(raw(:, 3));
Bz = cell2mat(raw(:, 4));
Bmod = cell2mat(raw(:, 5));

%% Plot
hFig1 = figure(1);
    set(hFig1, 'Position', [100 10 700 400])
    set(gca,'OuterPosition',[0 0 1 0.95]);
    % timelim = [min(time),max(time)*1.1];
    % time= time(1:10000);
    % Bz = Bz(1:10000);
    smoothB = movmean(By,20);
    plot(time,smoothB,'LineWidth',1); 
    hold on
    grid on
    legend('By')
    
%     % Guesses
%     g_amp = 0.03; % mT
%     g_per = 60/20.154; % 60/RPM 
%     g_offset = -0.0336; % offset
%     g_phase = 5.8754; % phase
%     
%    fit = @(b,time)  g_amp.*(sin(2*pi*time./b(1) + b(2))) + b(3); 
%         % Function to fit
%     fcn = @(b) sum((fit(b,time) - movmean(By,20)).^2);                              
%         % Least-Squares cost function
%     s = fminsearch(fcn, [g_per;  g_phase;  g_offset]);                       
%         % Minimize Least-Squares
%     rpm=60/s(1)
        
    xp = linspace(min(time),max(time),1000); 
    SineParams=sineFit(time,By,0);
    fit = @(SineParams,time)  SineParams(1)+SineParams(2).*...
        (sin(2*pi*time.*SineParams(3) + SineParams(4))) + SineParams(5); 
    plot(xp,fit(SineParams,xp), 'r'); % plot fitted line
    rmp = 60*SineParams(3) 
    legend('smoothed B', 'Best fit');
    set(gca,'fontsize',18, 'FontName','Times');
    
    
%% Funtions    
function SineParams=sineFit(x,y,varargin)
%Purpose: Estimation of noisy sine curve parameters by FFT and non linear fitting.
%No toolbox required.
%
% Syntax:
%       [SineParams]=sineFit(x,y,optional)
%       Input: x and y values, y=offs+amp+sin(2*pi*f*x+phi)+noise
%              optional: plot graphics if ommited. Do not plot if 0.
%       Output: SineParams(1): offset (offs)
%               SineParams(2): amplitude (amp)
%               SineParams(3): frequency (f)
%               SineParams(4): phaseshift (phi)
%               SineParams(5): MSE , if negative then SineParams are from FFT 
%       yOut=offs+amp*sin(2*pi*f*x+phi)
%
% Example:
% % generate y(x)
% x=-4:5;
% y=1+2*(sin(2*pi*0.1*x+2)+0.3*randn(size(x)));%Sine + noise
% [SineP]=sineFit(x,y)
% figure;
% xx=x(1):(x(end)-x(1))/222:x(end);%better resolution
% plot(x,y,xx,SineP(1)+SineP(2)*sin(2*pi*SineP(3)*xx+SineP(4)));
% %uncomment following lines if you want to save y=f(x) and run it sineFitDemo
% %paramsClean=[1,2,0.1,2];
% % save('xy.mat','x','y','paramsClean');
%
%Author: Peter Seibold

if nargin>2 && varargin{1}==0
  boolGraphic=false;
else
  boolGraphic=true;
end
%% FFT
pi2=2*pi;
NumSamples=length(x);
T=x(2)-x(1);
fNy=1/(2*T);%Nyquist frequency
offs=mean(y);%DC value, do not take (max(y)+min(y))/2!
y_m=y-offs;%FFT much better without offset
n = 128*2^nextpow2(NumSamples);%heavy zero padding
Y = fft(y_m,n);%Y(f)
n2=floor(n/2);
P2 = abs(Y/NumSamples);
P1 = P2(1:n2+1);
P1(2:end-1) = 2*P1(2:end-1);
fs = (0:n2)/n/T;% frequency scale
% %FFT parameters at peak
[maxFFT,maxFFTindx]=max(P1);%Peak magnitude and location
fPeak=fs(maxFFTindx);% f at peak
Phip=angle(Y(maxFFTindx))+pi/2;%Phi-Peak is for cos, sin(90°+alpha)=cos(betta), alpha=-betta
Phip=Phip-x(1)*fPeak*pi2;%shift for phi at x=0
if Phip<0;Phip=2*pi+Phip;end
%% Fitting
if numel(x)<12 || x(end)-x(1)<5/fPeak
  %Low sample number or low period number
  offs=(max(y)+min(y))/2;%Better results with offset by peak
elseif abs(offs)<0.1 && maxFFT>0.9
  offs=0;%Priority to 0
end
paramsFFTp=[offs,maxFFT,fPeak,Phip];
P1P=P1(1:maxFFTindx);%FFT only from f=0 to peak
if maxFFTindx>0.99*n2 && ~(maxFFTindx+2<n2 && P1(n2-2)<0.7*maxFFT)
  %FFT peak close to nyquist frequency and not a sharp peak
  fIndxExtra=[round(maxFFTindx*0.995);find(P1P<0.7*maxFFT,1,'last');find(P1P<0.5*maxFFT,1,'last')];
  AExtra=(max(y)-min(y))/2;
  Numf=3;%number of frequencies
  PeakInds = find(diff(sign(diff(P1))))+1;%Indices of FFT peaks
  PeakVals=sortrows([PeakInds(:),fs(PeakInds)',P1(PeakInds)'],3);%Values at peaks sorted
  if numel(PeakInds)>1 && PeakVals(end,3)*.95<PeakVals(end-1,3)
    %Second FFT peak nearly as large as max. peak
    %Evaluate extra peak
    fIndxExtra=[fIndxExtra;PeakVals(end-1,1)];
    Numf=4;%number of frequencies
  end
elseif x(end)-x(1)<1/fPeak
  %Period propably < 1
  fIndxExtra=[round(maxFFTindx*0.9);find(P1P<0.6*maxFFT,1,'last')];
  AExtra=(max(y)-min(y));
  Numf=2;%number of frequencies
else
  %FFT peak not at f-Nyquist and propably more than 1 sine period
  %Evaluate only at FFT peak
  Numf=1;%number of frequencies
  paramsFFT=paramsFFTp;
end
if Numf>1
  fExtra=fs(fIndxExtra);
  PhiExtra=(angle(Y(fIndxExtra))+pi/2-x(1)*fExtra*pi2);
  AExtra=repelem(AExtra,Numf);
  offExtra=repelem(offs,Numf);
  paramsFFT=[offExtra',AExtra',fExtra',PhiExtra'];
end
paramsOut=zeros(Numf,6);%for regression outputs
%fminsearch ======================================================
for i=1:Numf
  fun = @(SineParams)sseval(SineParams,x,y);
  [SineParams,SE] =fminsearch(fun,paramsFFT(i,:),...%SE: squared error
    optimset('MaxFunEvals',200000,'MaxIter',200000));
  %make frequency positive
  if SineParams(3)<0
    SineParams(3)=-SineParams(3);
    SineParams(4)=pi-SineParams(4);%sin(2*pi*-f-phi)=sin(2*pi*f+phi+pi)
  end
  %make amplitude positive
  if SineParams(2)<0
    SineParams(2)=-SineParams(2);
    SineParams(4)=SineParams(4)+pi;
  end
  MSE=SE/numel(x);
  paramsOut(i,:)=[SineParams,MSE,MSE];
  if SineParams(3)>fNy
    %f larger than nyquist limit
    paramsOut(i,5)=Inf;%set MSE to terrible
  end
end
%% take best manipulated score
[MSEmin,MSEminIndx]=min(paramsOut(:,5));
SineParams=paramsOut(MSEminIndx,1:4);
%  Determine max allowed amplitude by MSEmin
if MSEmin<=0.00001 || ...%extremly good MSE
    NumSamples<5 || ... %no MSE with nlinfit for less than 5 samples
    (NumSamples==5 && SineParams(3)<0.8*paramsFFT(1,3)) ||... %num period propably <1
    (MSEmin<1 && x(end)-x(1)<0.5/SineParams(3))%propably less than 0.5 periods
  maxAmp=66*maxFFT;%max allowed amplitude
elseif MSEmin>0.3
  maxAmp=4*maxFFT;
elseif MSEmin>0.01
  maxAmp=6*maxFFT;
elseif MSEmin>0.001
  maxAmp=18*maxFFT;
else
  %very good MSE, 0.00001 < MSE <0.001
  maxAmp=33*maxFFT;
end
if SineParams(2)>maxAmp || SineParams(3)>fNy
  %Best regression has too big amplitude or is over Nyquist limit,
  %take original FFT result
  SineParams=paramsFFTp;
  MSE=(sum((y - (SineParams(1)+SineParams(2)*sin(2*pi*SineParams(3)*x+SineParams(4)))).^2))/numel(x);
  SineParams(5)=-MSE;%Negative indicates SineParams are from FFT
else
  SineParams(5)=paramsOut(MSEminIndx,6);%for PlotResults
end
%make phase between 0 and 2 pi
SineParams(4)=rem(SineParams(4),pi2);
if SineParams(4)<0
  SineParams(4)=SineParams(4)+pi2;
end
if boolGraphic
  PlotResults(x,y,SineParams,paramsFFTp,fs,P1,maxFFTindx,maxFFT);
end
end

function sse = sseval(SineParams,x,y)
offs = SineParams(1);
A =SineParams(2);
f=SineParams(3);
Phi=SineParams(4);
sse = sum((y - (offs+A*sin(2*pi*f*x+Phi))).^2);
end

%% Plot results (optional)
function PlotResults(x,y,SineParams,paramsFFTp,fs,P1,maxFFTindx,maxFFT)
xstart=x(1);
xend=x(end);
xLength=xend-xstart;
xSstep=min(xLength/100,1/SineParams(3)*0.1);
xS=xstart:xSstep:xend;
xFFT=xstart-xLength*.02:xLength/100:xend+xLength*.02;%FFTcurve longer
y5=SineParams(1)+SineParams(2)*sin(2*pi*SineParams(3)*xS+SineParams(4));%result
yFFT=paramsFFTp(1)+paramsFFTp(2)*sin(2*pi*paramsFFTp(3)*xFFT+paramsFFTp(4));
%time plot:
hFigPlotSin = findobj( 'Type', 'Figure', 'Tag', 'Fig$PlotSin' );
if isempty(hFigPlotSin)
  screensize=get(0, 'MonitorPositions');
  hFigPlotSin=figure('Tag','Fig$PlotSin','Name','Sinus',...
    'OuterPosition',[960,screensize(1,4)/2,screensize(1,3)-960,screensize(1,4)/2]);
  drawnow
end
figure(hFigPlotSin(1));
cla reset;
plot(x,y,'k.');%time series as dots
xlabel('Time [s]');
hold on;
pIn=plot(x,y,'r-');%time series as line
pFFT=plot(xFFT,yFFT,'color',[0.7 0.7 0.7]);
pResult=plot(xS,y5,'b-');%result
legend([pIn,pResult,pFFT],'Input','Result', 'FFT peak');
hold off;
grid on;
%FFT plot:
hFigPlotFFT = findobj( 'Type', 'Figure', 'Tag', 'Fig$PlotFFT' );
if isempty(hFigPlotFFT)
  hFigPlotFFT=figure('Tag','Fig$PlotFFT','Name','FFT',...
    'OuterPosition',[960,40,screensize(1,3)-960,screensize(1,4)/2-45]);
  drawnow
end
figure(hFigPlotFFT(1));
cla reset;
pFFTin=plot(fs,P1,'r-');
xlabel('Frequency [Hz]');
ylabel('Amplitude')
hold on;
pFFTmax=plot(fs(maxFFTindx),maxFFT,'r+','MarkerSize',12);%max FFT
pFFTresult=plot(SineParams(3),SineParams(2),'b+','LineWidth',2);
plot([SineParams(3),SineParams(3)],[0,max(max(P1)*1.01,SineParams(2))],'b-');
hLeg=legend([pFFTin,pFFTresult,pFFTmax],'Input',...
  ['Result:     ' num2str(SineParams(2),3) ', ' num2str(SineParams(3),3) ' Hz'],...
  ['max FFT:  ' num2str(maxFFT,3) ', ' num2str(fs(maxFFTindx),3) ' Hz'],...
  'Location','best');
title(hLeg,'        amplitude | frequency','FontSize',8);
hold off;
grid on;
disp(['Result:        y= ' num2str(SineParams(1)) ' + ' num2str(SineParams(2)) ...
  ' * sin(2*pi*' num2str(SineParams(3)) '+' num2str(SineParams(4)) ')   MSE: ' num2str(SineParams(5))]);
disp(['FFT:           y= ' num2str(paramsFFTp(1)) ' + ' num2str(paramsFFTp(2)) ...
  ' * sin(2*pi*' num2str(paramsFFTp(3)) '+' num2str(paramsFFTp(4)) ')' ]);
end