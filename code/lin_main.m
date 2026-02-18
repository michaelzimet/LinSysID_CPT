% second order model of pupil 
% after running eye_sysID.m with nx=2, iOrder=1, and iCond=1

global A B;

A = thisA;
B = thisB;
C = thisC;
u = 1;
xstar = -inv(A)*B*u;

r0 = 1;
th0 = -pi/2;

%initial condition
x10 = xstar(1);
x20 = xstar(2);

[tsoln,Xsoln] = ode45('func_lin',[0 100],[x10,x20]);

x1soln = Xsoln(:,1);
x2soln = Xsoln(:,2);

ysoln =  C(1)*x1soln +C(2)*x2soln;

a=figure;
b=figure;

figure(a)
subplot(2,1,1)
plot(tsoln,-1*x1soln,'r',tsoln,-1*x2soln,'g',tsoln,ysoln,'b','LineWidth',3);
ylim([-0.2 0.5])
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 12)
ylabel([meas ', ' thisYLabel], 'Interpreter','latex', 'FontSize', 12);
saveas(a, ['figs/pupil_traj.eps'],'epsc')


figure(b)
plot(-1*x1soln,-1*x2soln, 'b','LineWidth',3)
x1_plot = 0:0.001:0.25;    
hold on
for y_plot = -0.5:0.1:0.2
    x2_plot = -(C(1)/C(2))*x1_plot + y_plot/C(2);
    plot(x1_plot,x2_plot, 'r')
end
axis([0,0.25,-0.1,0.5])
xlabel('$x_1$', 'Interpreter','latex')
ylabel('$x_2$', 'Interpreter','latex')
saveas(b, ['figs/pupil_cont.eps'],'epsc')

