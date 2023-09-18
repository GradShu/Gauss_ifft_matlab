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
nPoint = 8192;
z0 = 10;
sigma = 0:0.0083:8191*0.0083;
fSigma = 568.8*exp(-((sigma-1.45)/0.3129).^2);
pSigma = -4*pi*z0*sigma;
iSigma = fSigma.*exp(1i*pSigma);
%calculate ifft of iSigma
[intensity,phase,zData]=ifftOfIsigma(iSigma,sigma(1),sigma(end),nPoint);

figure(1);
subplot(2,2,1);
setIntensityImage(sigma,fSigma,'Wavenumber \sigma(\mum^{-1})','Intensity(a.u.)');
subplot(2,2,3);
setPhaseImage(sigma,pSigma,'Wavenumber \sigma(\mum^{-1})','Phase(rad)');
subplot(2,2,2);
setIntensityImage(zData,intensity,'Position z(\mum)','Intensity(a.u.)');
subplot(2,2,4);
setPhaseImage(zData,phase,'Position z(\mum)','Phase(rad)');

function setIntensityImage(xData,yData,xlabel,ylabel)
%This function setImage use for setting image property, such as axes, line,
%colour and so on. 
yDataLine = plot(xData,yData);
%Find the max of Intensity.
[yDataMax,yDataMaxindex]= max(yData);
hold on
yMaxLine = plot([xData(yDataMaxindex),xData(yDataMaxindex)],[0,yDataMax]);
%Set the location of the tick marks along the axis
set(gca,'Xcolor',     [0 0 0],...
        'Ycolor',     [0 0 0],...
        'Color' ,     [1 1 1],...
        'FontName',   'Arial',...
        'FontAngle',  'normal',...
        'FontSize',    13,...
        'Xlim',       [xData(yDataMaxindex)-0.7 xData(yDataMaxindex)+0.7]);
%Set label names for X coordinates and Y coordinates
set(get(gca,'XLabel'),      'String',xlabel,...
                            'FontName','Arial',...
                            'FontSize',14);
set(get(gca,'YLabel'),      'String',ylabel,...
                            'FontName','Arial',...
                            'FontSize',14);
set(yDataLine,'color',      'blue',...
              'LineStyle',  '-',...
              'LineWidth',   2);
set(yMaxLine,'color',       'k',...
              'LineStyle',  '-.',...
              'LineWidth',   2);
end

function setPhaseImage(xData,yData,xlabel,ylabel)
%This function setImage use for setting image property, such as axes, line,
%colour and so on. 
yDataLine = plot(xData,yData);
% zeroPhase = ;
%Set the location of the tick marks along the axis
set(gca,'Xcolor',     [0 0 0],...
        'Ycolor',     [0 0 0],...
        'Color' ,     [1 1 1],...
        'FontName',   'Arial',...
        'FontAngle',  'normal',...
        'FontSize',   13);
%Set label names for X coordinates and Y coordinates
set(get(gca,'XLabel'),'String',xlabel,...
                      'FontName','Arial',...
                      'FontSize',14);
set(get(gca,'YLabel'),'String',ylabel,...
                      'FontName','Arial',...
                      'FontSize',14);
set(yDataLine,'color',      'blue',...
              'LineStyle',  '-',...
              'LineWidth',   2);
end

function [intensity,phase,zData]=ifftOfIsigma(input,sigmaStart,sigmaEnd,nPoint)
%This function use for calculating ifft of iSigma.
%The output is intensity and phase after ifft.
deltaSigma = (sigmaEnd-sigmaStart)/(nPoint-1);
deltaZ = 1/(2*nPoint*deltaSigma);
zData = (1:nPoint)*deltaZ;
intensity = abs(ifft(input));
phase = angle(ifft(input));
end

