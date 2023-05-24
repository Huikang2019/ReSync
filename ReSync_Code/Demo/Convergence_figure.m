
% parameters with uniform topology
n =200; p = (log(n) / n)^(1 / 3); q = (log(n) / n)^(1 / 3); sigma=0; model='uniform';

% generate data with uniform topology
model_out = Uniform_Topology(n,q,1-p,sigma,model);

Ind = model_out.Ind; % matrix of edge indices (m by 2)
RijMat = model_out.RijMat; % given corrupted and noisy relative rotations
ErrVec = model_out.ErrVec; % ground truth corruption levels
R_orig = model_out.R_orig; % ground truth rotations

% set ReSync defult parameters
ReSync_parameters.max_iter = 300;
ReSync_parameters.stepsize = 1 / (n*p*q);
ReSync_parameters.stop_threshold = 1e-30;

% run ReSync
R_SP_co = SpectrIn(Ind, RijMat);

ReSync_parameters.decay = 0.98;
[~, Dist1] = ReSync(Ind , RijMat, R_SP_co, R_orig, ReSync_parameters);
ReSync_parameters.decay = 0.95;
[~, Dist2] = ReSync(Ind , RijMat, R_SP_co, R_orig, ReSync_parameters);
ReSync_parameters.decay = 0.90;
[~, Dist3] = ReSync(Ind , RijMat, R_SP_co, R_orig, ReSync_parameters);
ReSync_parameters.decay = 0.85;
[~, Dist4] = ReSync(Ind , RijMat, R_SP_co, R_orig, ReSync_parameters);
ReSync_parameters.decay = 0.80;
[~, Dist5] = ReSync(Ind , RijMat, R_SP_co, R_orig, ReSync_parameters);
ReSync_parameters.decay = 0.7;
[~, Dist6] = ReSync(Ind , RijMat, R_SP_co, R_orig, ReSync_parameters);

fig = figure;

semilogy(Dist1,'-s','LineWidth',2,'MarkerIndices', 1:30:ReSync_parameters.max_iter,'MarkerSize',8);
hold on
box on
semilogy(Dist2,'-o','LineWidth',2,'MarkerIndices', 1:30:ReSync_parameters.max_iter,'MarkerSize',8);
semilogy(Dist3,'-hexagram','LineWidth',2,'MarkerIndices', 1:30:ReSync_parameters.max_iter,'MarkerSize',8);
semilogy(Dist4,'-^','LineWidth',2,'MarkerIndices', 1:30:ReSync_parameters.max_iter,'MarkerSize',8);
semilogy(Dist5,'-d','LineWidth',2,'MarkerIndices', 1:30:ReSync_parameters.max_iter,'MarkerSize',8);
semilogy(Dist6,'-*','LineWidth',2,'MarkerIndices', 1:20:ReSync_parameters.max_iter,'MarkerSize',8);
hold off
set(gcf, 'Color', 'white');
ylim([5e-8 3])
set(gca, 'LineWidth' , 1.7, 'FontName', 'Times New Roman','FontSize',18);
legend('$\gamma = 0.98$','$\gamma = 0.95$','$\gamma = 0.9$','$\gamma = 0.85$','$\gamma = 0.8$','$\gamma = 0.7$','Interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','NorthEast')
xlabel('Iteration Number','Interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('dist$(\textbf{X} - \textbf{X}^\star) / \sqrt{n}$','Interpreter','latex','FontName','Times New Roman','FontSize',20)


