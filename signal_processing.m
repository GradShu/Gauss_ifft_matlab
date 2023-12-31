clear;close all
%The simulation use for signal processing in wavenumber scanning domains
% of white-light interferometers. In the simulation, The IFT(Inverse Fourie
% -r Transform) of I can be get.
%Creat guass function:
%      iSigma =  a1*exp(-((sigma-b1)/c1)^2)
%Coefficients (with 95% confidence bounds):
%        a1 =       568.8
%        b1 =        1.45  
%        c1 =      0.3129
%Creat a phase functoion:
%      pSigma = -4*pi*z0*sigma;

%Related parameters
%intensityMark,phaseMark,ifftIntenMark,ifftPhaseMark:image identification
intensityMark = 221;
phaseMark = 223;
ifftIntenMark = 222;
ifftPhaseMark = 224;
%npoint:sampling point 
nPoint = 2^16;
z0 = 10;
sigma = 0:0.0083:(nPoint-1)*0.0083;
fSigma = 568.8*exp(-((sigma-1.45)/0.3129).^2);
pSigma = -4*pi*z0*sigma;
iSigma = fSigma.*exp(1i*pSigma);
%calculate ifft of iSigma
[intensity,phase,zDataInter]=ifftOfIsigma(iSigma,sigma(1),sigma(end),nPoint);
zeroPhase = findZeroPhase(zDataInter,phase,0);

figure(1);
subplot(2,2,1);
setImage(221,sigma,fSigma,'Wavenumber \sigma(\mum^{-1})','Intensity(a.u.)');
subplot(2,2,3);
setImage(223,sigma,pSigma,'Wavenumber \sigma(\mum^{-1})','Phase(rad)');
subplot(2,2,2);
setImage(222,zDataInter,intensity,'Position z(\mum)','Intensity(a.u.)');
subplot(2,2,4);
setImage(224,zDataInter,phase,'Position z(\mum)','Phase(rad)',zeroPhase);

function [intensity,phase,zDataInter]=ifftOfIsigma(input,sigmaStart,sigmaEnd,nPoint)
%This function use for calculating ifft of iSigma.
%The output is intensity and phase after ifft.
deltaSigma = (sigmaEnd-sigmaStart)/(nPoint-1);
deltaZ = 1/(2*nPoint*deltaSigma);
zData = (1:nPoint)*deltaZ;
intensity = abs(ifft(input));
phase = angle(ifft(input));
%interpolation
zDataInter = linspace(zData(1),zData(end),nPoint*10);
intensity = interpn(zData,intensity,zDataInter,'linear');
phase = interpn(zData,phase,zDataInter,'linear');
end

function zeroPhase = findZeroPhase(zData,phaseData,zeroPoint)
diff = 0.01;
%limit the range of zData,about (zp-lambda/8,zp+lambda/8)
LimZData = zData(zData>(10-(1/1.45)/8)&zData<(10+(1/1.45)/8));
phaseData = phaseData(zData>(10-(1/1.45)/8)&zData<(10+(1/1.45)/8));
% phaseData = phaseData()
positivePhase = phaseData(phaseData>0);
negativePhase = phaseData(phaseData<0);
while true 
    if length(find(abs(positivePhase-zeroPoint)<diff))>1
       lastDiff = diff;
       diff = diff/2;
    elseif (isempty(find(abs(positivePhase-zeroPoint)<diff, 1)))
        diff = (diff + lastDiff)/2;
    else
        zeroPhase1= positivePhase(abs(positivePhase-zeroPoint)<diff);
        break
    end
end
while true 
    if length(find(abs(negativePhase-zeroPoint)<diff))>1
       lastDiff = diff;
       diff = diff/2;
    elseif (isempty(find(abs(negativePhase-zeroPoint)<diff, 1)))
        diff = (diff + lastDiff)/2;
    else
        zeroPhase2= negativePhase(abs(negativePhase-zeroPoint)<diff);
        break
    end
end
zeroData1 = LimZData(phaseData==zeroPhase1);
zeroData2 = LimZData(phaseData==zeroPhase2);
zeroPhase = ((0-zeroPhase2)*(zeroData1-zeroData2))/(zeroPhase1-zeroPhase2)...
            +zeroData2;
end

function setImage(imageMark,xData,yData,xlabel,ylabel,zeroPhase)
%This function setImage use for setting image property, such as axes, line,
%colour and so on. 
%ampMaxZ:The abscissa where the maximum amplitude value is located.
global ampMaxZ
yDataLine = plot(xData,yData);
%Find the max of Intensity.
[yDataMax,yDataMaxindex]= max(yData);
%Set the location of the tick marks along the axis
if imageMark == 221
    hold on
    ampMaxZ = xData(yDataMaxindex);
    yMaxLine = plot([xData(yDataMaxindex),xData(yDataMaxindex)],...
               [0,yDataMax]);
    set(gca,'Xlim',   [xData(yDataMaxindex)-0.7 xData(yDataMaxindex)+0.7]);
    set(yDataLine,'color',      'blue',...
                  'LineStyle',  '-',...
                  'LineWidth',   2);
    set(yMaxLine, 'color',       'k',...
                  'LineStyle',  '--',...
                  'LineWidth',   2);
elseif imageMark == 223
    set(gca,'Xlim',   [ampMaxZ-0.7 ampMaxZ+0.7]);
    set(yDataLine,'color',      'blue',...
                  'LineStyle',  '-',...
                  'LineWidth',   2);
elseif imageMark == 222
    hold on
    ampMaxZ = xData(yDataMaxindex);
    yMaxLine = plot([xData(yDataMaxindex),xData(yDataMaxindex)],...
               [0,yDataMax]);
    set(yDataLine,'color',      'blue',...
                  'LineStyle',  '-',...
                  'LineWidth',   2);
    set(yMaxLine, 'color',       'k',...
                  'LineStyle',   '--',...
                  'LineWidth',    2);
    set(gca,'Xlim',   [xData(yDataMaxindex)-2 xData(yDataMaxindex)+2]);
else
    hold on
    zeroLine = plot([ampMaxZ-2,ampMaxZ+2],[0,0]);
    zeroPhaseLine = plot([zeroPhase,zeroPhase],[-3.14,0]);
    set(gca,'Xlim', [ampMaxZ-2 ampMaxZ+2]);
    set(gca,'Ylim', [-3.14 3.14]);
    set(yDataLine,'color',      'blue',...
                  'LineStyle',  '-',...
                  'LineWidth',   2);
    set(zeroLine, 'color',       'k',...
                  'LineStyle',   '--',...
                  'LineWidth',    2);
    set(zeroPhaseLine, 'color',  'k',...
                  'LineStyle',   '--',...
                  'LineWidth',    2);
end

set(gca,'Xcolor',     [0 0 0],...
        'Ycolor',     [0 0 0],...
        'Color' ,     [1 1 1],...
        'FontName',   'Arial',...
        'FontAngle',  'normal',...
        'FontSize',    13);
%Set label names for X coordinates and Y coordinates
set(get(gca,'XLabel'), 'String',xlabel,...
                       'FontName','Arial',...
                       'FontSize',14);
set(get(gca,'YLabel'), 'String',ylabel,...
                       'FontName','Arial',...
                       'FontSize',14);
end