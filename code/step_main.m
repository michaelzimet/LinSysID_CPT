%Make predictions for Impulse and Ramp Response
%after running physio_sysID.m and eye_sysID.m with iOrder=1 and iCond=1

A = thisA;
B = thisB;
C = thisC;
l=1;
m=1/90;
tau = 90;
tmax=150;

z=figure;
t = tiledlayout(2,2);
xlabel(t, 'Time (s)','Interpreter', 'latex')

% Step Function Response
for t=1:89
    %Step On
    Xsoln(:,t) = expm(A*t)*ic_cont + inv(A)*(expm(A*t)-eye(nx))*B;
end
for t=90:tmax
    %Step Off
    Xsoln(:,t) = expm(A*t)*ic_cont + expm(A*(t-tau))*inv(A)*(expm(A*tau) - eye(nx))*B;
end
ysoln =  C*Xsoln;

% Impulse Response

for t=1:tmax
    ysoln_imp(t) = C * expm(A*t) * ic_cont + l * C * expm(A*t) * B;
end

for t=1:89
    %Ramp On
    Xsoln_ramp(:,t) = expm(A*t) *ic_cont + m*((inv(A))^2*(expm(A*t)-eye(nx))-inv(A)*t)*B;
end
for t=90:tmax
    %Ramp Off
    Xsoln_ramp(:,t) = expm(A*t)*ic_cont + expm(A*(t-tau))*m*((inv(A))^2*(expm(A*t)-eye(nx))-inv(A)*t)*B;
end
ysoln_ramp =  C*Xsoln_ramp;

%openfig('figs/imp_ramp.fig');
figure(z)
nexttile;
plot(1:tmax, ysoln, 'b', 1:tmax, ysoln_imp, 'r', 1:tmax, ysoln_ramp, 'g', 'LineWidth',3);
ylabel([meas ', ' thisYLabel], 'Interpreter','latex');
ylim([-0.5 1.5])
xlim([0 150])
xticks(0:30:150)
yticks(-0.5:0.5:1.5)
saveas(gcf, 'figs/imp_ramp.fig','fig')

