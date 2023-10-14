syms u v

f = 0.5*(0.25 *sin(4*u)* sin(v) *(cos(-4*v) + cos(-v))-0.25*sin(2*u)*sin(2*u)*(cos(-5*v) -cos(-3*v)) -sin(4*u)*sin(-4*v)*(cos(v)+1) -0.5 * sin(4*u)*sin(-4*v) *(cos(v)+1) -0.5* sin(4*u) *sin(-v))

f2 = matlabFunction(f);
H = hessian(f,[u,v]);
g = gradient(f, [u,v]);

beta = 1./(1 + norm(g).^2);

S = beta.*(trace(H).^2 - trace(H*H)) + 2.*beta.^2.*(transpose(g)*(H*H - trace(H).*H)*g);
colormap parula;
S2 = matlabFunction(S);
u2 = linspace(-1.7, 1.7, 100);
v2 = linspace(-1.7, 1.7, 100);

[U,V] = meshgrid(u2,v2);

surf(U,V, f2(U,V), S2(U,V))
c = colorbar('TickLabelInterpreter', 'latex');

c.Limits = [-1, 1];
caxis([-1 1]);
xticks([-pi/2, -pi/4, 0, pi/4, pi/2]);
yticks([-pi/2, -pi/4, 0, pi/4, pi/2]);
zticks([-1.5, 0, 1.5]);
xticklabels({'$-0.5\pi$', '$-\pi/4$', '0', '$\pi/4$', '$\pi/2$'});
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
c.Label.Interpreter = 'latex';
zlim([-3, 3]); % Set the z-axis limits directly
xlabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 31); % Use LaTeX for labels
ylabel('$\gamma$', 'Interpreter', 'latex', 'FontSize', 31);
zlabel('$C(\gamma, \beta)$', 'Interpreter', 'latex', 'FontSize', 31);
set(gca, 'XTickLabel', {'$-0.5\pi$', '$-0.25\pi$', '$0$','$0.25\pi$' '$0.5\pi$'}, 'FontSize', 31); % Change font size for x-axis tick labels
set(gca, 'YTickLabel', {'$-0.5\pi$', '$-0.25\pi$', '$0$','$0.25\pi$' '$0.5\pi$'}, 'FontSize', 31); % Change font size for y-axis tick labels
set(gca, 'ZTickLabel', {'$-1$', '$0$', '$1$'}, 'FontSize', 32); % Change font size for y-axis tick labels
xtickangle(0);
ytickangle(0);
saveas(gcf, 'my_plot.png');





