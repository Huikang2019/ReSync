
% parameters with uniform topology
n = 200; q=0.2; sigma=1; model='uniform';

num_sim = 10;
num_iter = 9;

final_beta = 5; % temp

% set CEMP defult parameters
CEMP_parameters.max_iter = 6;
CEMP_parameters.reweighting = 2.^((1:6)-1);
CEMP_parameters.nsample = 50;
CEMP_parameters.gcw_beta = final_beta;

% set DESC default parameters
lr = 0.01;
DESC_parameters.iters = 100;
DESC_parameters.learning_rate = lr;
DESC_parameters.make_plots = false;
DESC_parameters.Gradient = ConstantStepSize(lr);

% set MPLS default parameters
MPLS_parameters.stop_threshold = 1e-3;
MPLS_parameters.max_iter = 100;
MPLS_parameters.reweighting = CEMP_parameters.reweighting(end);
MPLS_parameters.thresholding = [0.95,0.9,0.85,0.8];
MPLS_parameters.cycle_info_ratio = 1./((1:MPLS_parameters.max_iter)+1);

% set ReSync defult parameters
ReSync_parameters.max_iter = 400;
ReSync_parameters.decay = 0.95;
ReSync_parameters.stop_threshold = 1e-8;

dist_MPLS = zeros(num_sim, num_iter);
dist_CEMP = zeros(num_sim, num_iter);
dist_IRLS_L12 = zeros(num_sim, num_iter);
dist_DESC = zeros(num_sim, num_iter);
dist_ReSync = zeros(num_sim, num_iter);

for sim = 1 : num_sim
    for iter = 1 : num_iter
        p = (iter + 1) / 10;

        % generate data with uniform topology
        model_out = Uniform_Topology(n,q,1-p,sigma,model);

        Ind = model_out.Ind; % matrix of edge indices (m by 2)
        RijMat = model_out.RijMat; % given corrupted and noisy relative rotations
        ErrVec = model_out.ErrVec; % ground truth corruption levels
        R_orig = model_out.R_orig; % ground truth rotations

        % run MPLS and CEMP+MST
        [R_MPLS, R_CEMP_MST] = MPLS(Ind,RijMat,CEMP_parameters, MPLS_parameters);
        dist_MPLS(sim, iter) = Dist2(R_MPLS, R_orig);

        % run ReSync
        R_SP_co = SpectrIn(Ind, RijMat);

        ReSync_parameters.stepsize = 1 / (n*p*q);
        [R_ReSync, ~] = ReSync(Ind , RijMat, R_SP_co, R_orig, ReSync_parameters );

        dist_ReSync(sim, iter) = Dist2(R_ReSync, R_orig);

        % run CEMP+GCW
        R_CEMP_GCW = CEMP_GCW(Ind, RijMat, CEMP_parameters);

        dist_CEMP(sim, iter) = Dist2(R_CEMP_GCW, R_orig);

        % run IRLS
        % R_IRLS_GM = IRLS_GM(RijMat, Ind);
        R_IRLS_L12 = IRLS_L12(RijMat,Ind);
        dist_IRLS_L12(sim, iter) = Dist2(R_IRLS_L12, R_orig);        

        % run DESC
        [R_DESC, R_DESC_init, ~] = DESC(Ind, RijMat, DESC_parameters, R_orig, ErrVec);
        dist_DESC(sim, iter) = Dist2(R_DESC, R_orig);
    end
end

fig = figure;

x_label = (2 : 1 : num_iter + 1) / 10;
plot(x_label, sum(dist_MPLS, 1) / num_sim,'-s','LineWidth',2,'MarkerSize',8);
hold on
box on
plot(x_label, sum(dist_CEMP, 1) / num_sim,'-o','LineWidth',2,'MarkerSize',8);
plot(x_label, sum(dist_IRLS_L12, 1) / num_sim,'-hexagram','LineWidth',2,'MarkerSize',8);
plot(x_label, sum(dist_DESC, 1) / num_sim,'-^','LineWidth',2,'MarkerSize',8);
plot(x_label, sum(dist_ReSync, 1) / num_sim,'-d','LineWidth',2,'MarkerSize',8);
hold off
%title('$\sigma = 0.0$','Interpreter','latex')
set(gcf, 'Color', 'white');
set(gca, 'LineWidth' , 1.7, 'FontName', 'Times New Roman','FontSize',18);
legend('MPLS','CEMP\_GCW','IRLS\_L12','DESC','ReSync','FontName','Times New Roman','FontSize',20,'Location','NorthEast')
xlabel('$p$','Interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('dist$(\textbf{X} - \textbf{X}^\star) / \sqrt{n}$','Interpreter','latex','FontName','Times New Roman','FontSize',20)



