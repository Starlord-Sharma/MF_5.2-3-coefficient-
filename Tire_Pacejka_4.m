
tiledlayout(4,1) %To create 4 axes in one figure
nexttile %1st axes
kPa2psi = 0.145038;
P = P*kPa2psi;
plot(ET,P) %Plot P (Inflation Pressure) vs ET
ylabel('Inf. Pressure [psi]')
xlabel('Elapsed Time [s]')
nexttile %2nd axes
plot(ET,IA) %Plot IA vs ET
ylabel('Inc. Angle [deg]')
xlabel('Elapsed Time [s]')
nexttile %3rd axes
plot(ET,FZ) %Plot FZ vs ET
ylabel('Norm. Load [N]')
xlabel('Elapsed Time [s]')
nexttile %4th axes
plot(ET,SA) %Plot SA vs ET
ylabel('Slip Angle [deg]')
xlabel('Elapsed Time [s]')
%% 
%First we generate the logical array which contain 1 only in the
%position in which the condition is satisfied
idx = ET > 480 & ET < 630
%Then, we use this array as index to extract
%the correct data section from SA, FZ and FY
SAext = SA(idx)
FZext = FZ(idx)
FYext = FY(idx)
tiledlayout(2,2) %To create 4 axes in one figure
nexttile %1st axes
plot(ET(idx),P(idx)) %Plot P (Inflation Pressure) vs ET
ylabel('Inf. Pressure [psi]')
xlabel('Elapsed Time [s]')
nexttile %2nd axes
plot(ET(idx),IA(idx)) %Plot IA vs ET
ylabel('Inc. Angle [deg]')
xlabel('Elapsed Time [s]')
nexttile %3rd axes
plot(ET(idx),FZ(idx)) %Plot FZ vs ET
ylabel('Norm. Load [N]')
xlabel('Elapsed Time [s]')
nexttile %4th axes
plot(ET(idx),SA(idx)) %Plot SA vs ET
ylabel('Slip Angle [deg]')
xlabel('Elapsed Time [s]')
%% 

%%function fy = Pacejka4_Model(P,X)

%%D = (P(1) + P(2)/1000*X(:,2)).*X(:,2);
%%fy = D.*sin(P(4)*atan(P(3).*X(:,1)));

%%end
%% 
figure
plot3(SAext,FZext,FYext,'k-')
xlabel('Slip Angle [deg]')
ylabel('Normal Load [N]')
zlabel('Lateral Force [N]')
%% 
xdata = [SAext FZext]; %Input variables array
ydata = FYext; %Independent variable
x0 = [1 1 1 1] %Initial guesses
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',10e9,'StepTolerance',1e-10);
%% 
[FY4,f_resnorm,f_residual] = lsqcurvefit('Pacejka4_Model',x0,xdata,ydata,[],[],options)
fz = linspace(-200,-1600,100)' %Range for FZ
sa = linspace(-13,13,100)' %Range for SA
fy = zeros(length(sa),length(fz)); %Array to store FY
for i=1:length(fz) %For each value of FZ
fy(i,:) = Pacejka4_Model(FY4,[sa fz(i)*ones(100,1)]);
end
fy
%Raw data plot
figure
plot3(SAext,FZext,FYext,'k-')
xlabel('Slip Angle [deg]')
ylabel('Normal Load [N]')
zlabel('Lateral Force [N]')
hold on
%Fitted data plot
surf(sa,fz,fy)
hold off
x0 = [0.0817 -0.5734 -0.5681 -0.1447];
options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',10e9,'StepTolerance',1e-10);
[FY4(2,:),f_resnorm(2),f_residual(:,2)] = lsqcurvefit('Pacejka4_Model',x0,xdata,ydata,[],[],options)
for i=1:length(fz) %For each value of FZ
fy(i,:,2) = Pacejka4_Model(FY4(2,:),[sa fz(i)*ones(100,1)]);
end
plot(f_residual(:,1))
hold on
plot(f_residual(:,2))
hold off
legend("Fit 1","Fit 2")
%Raw data plot
figure
plot3(SAext,FZext,FYext,'k-')
xlabel('Slip Angle [deg]')
ylabel('Normal Load [N]')
zlabel('Lateral Force [N]')
hold on
%Fitted data plot
surf(sa,fz,fy(:,:,2))
hold off
%The guesses were generated randomly using the following command
%x0 = -1*rand(1,4).^randi([1 2]);
%The fitting was then run until a reasonable solution was achieved
