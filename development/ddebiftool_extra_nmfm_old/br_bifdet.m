function [newbranch, success] = br_bifdet(funcs, branch, varargin)
%% Compute stability, detect and correct bifurcations along branch
% function [newbranch, success] = br_bifdet(branch)
% Purpose:
%   This routine computes stability, detects possible bifurcations,
%   corrects bifurcation candidates to genuine bifurcation points,
%   and adds them to the branch as a "native" point of the branch
% INPUT:
%  funcs: defining functions
%	branch: continued branch (contains method, initial points and all new points)
% OUTPUT:
%	newbranch: extended branch
%  success: 1 if a point was detected, 0 otherwise
%
% $Id: br_bifdet.m 109 2015-08-31 23:45:11Z jansieber $
%

%% Process method parameters
default={'plotaxis',[]};
options=dde_set_options(default,varargin);
tolcorr = branch.method.bifurcation.correction_tolerance;
imagthresh=branch.method.bifurcation.imagthreshold; % threshold for treating a root as complex
isimag=@(x)x>imagthresh;
monitor_eigenvalues=branch.method.bifurcation.monitor_eigenvalues;
plot_testfunctions=branch.method.bifurcation.plot_testfunctions;
pause_on_bifurcation=branch.method.bifurcation.pause_on_bifurcation;
if monitor_eigenvalues && isempty(options.plotaxis)
    options.plotaxis=gca;
end
%% Set variables
success = 0;
kind = branch.point(1).kind;
branch = br_bifinit(branch); % Make sure all points get flag, nmfm, nvec (hopf)
newbranch = branch;
% method = branch.method.point;
stmethod = branch.method.stability;
free_par = branch.parameter.free;

for i=1:length(branch.point)
    p=branch.point(i);
    
    if monitor_eigenvalues
        p_splot(p,'plotaxis',options.plotaxis);
        pause(.1);
    end
    
    
    if strcmp(kind,'stst')
        %% stst branch
        par1(i) = branch.point(i).parameter(free_par); % needed for plotting later
        phi_f(i) = det(ch_matrix(funcs,p.x,p.parameter,0));
        phi_H(i) = nmfm_smrp(funcs, p,stmethod,'remove_omega',false,'threshold',isimag);
        
        if(i==1)
            continue;
        end
        
        if abs(sign(phi_f(i))+sign(phi_f(i-1)))==0
            %% fold
            fprintf('BR_BIFDET: Fold detected near par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)));
            % correction
            method=df_mthod(funcs,'fold');
            method.stability=branch.method.stability;
            method.minimal_accuracy = tolcorr;
            fold=p_tofold(funcs,p);
            [fold, located] = p_correc(funcs,fold,free_par,[],method.point);
            if located
                fold.stability=p_stabil(funcs,fold,method.stability);
                fold = nmfm_cusp(funcs, fold); % normal form computation for fold
                [newbranch, added] = add_to_branch(funcs,newbranch,branch,i,fold); % add to new branch
                if added
                    fprintf('BR_BIFDET: Fold located at  par(%d) = %.10f.\n',...
                        free_par(1), fold.parameter(free_par(1)));
                    fprintf('BR_BIFDET: Normal form coefficient: b = %.10f\n\n', fold.nmfm.b);
                end
            else
                disp('BR_BIFDET: Failed to correct fold.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
        
        if abs(sign(phi_H(i))+sign(phi_H(i-1)))==0 && phi_H(i) ~= 0
            %% Hopf
            % correction
            fprintf('BR_BIFDET: Hopf detected near par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)));
            method=df_mthod(funcs,'hopf');
            method.minimal_accuracy = tolcorr;
            hopf=p_tohopf(funcs,p);
            [hopf, located] = p_correc(funcs,hopf,free_par,[],method.point);
            if located
                hopf.stability=p_stabil(funcs,hopf,method.stability);
                hopf = nmfm_hopf(funcs, hopf,  branch.point(i-1)); % normal form computation for hopf
                [newbranch,added] = add_to_branch(funcs,newbranch,branch,i,hopf); % add to new branch
                if added                   
                    fprintf('BR_BIFDET: Hopf located at  par(%d) = %.10f.\n',...
                        free_par(1), hopf.parameter(free_par(1)));
                    fprintf('BR_BIFDET: Normal form coefficient: L1 = %.10f\n\n', hopf.nmfm.L1);
                end
            else
                fprintf('BR_BIFDET: Failed to correct Hopf.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
    elseif strcmp(kind,'fold')
        %% fold branch
        p=branch.point(i);
        fold = nmfm_fold(funcs, p);
        
        par1(i) = p.parameter(free_par(1)); % needed for plotting
        par2(i) = p.parameter(free_par(2));
        
        % test functions
        phi_cusp(i) = fold.nmfm.b;
        phi_bt(i) = phi_f_bt(funcs,p);
        phi_zh(i) = nmfm_smrp(funcs, p,stmethod,'remove_omega',false,'threshold',isimag);
        
        if(i==1)
            continue;
        end
        
        if(i>1)
            phi_bt_dir(i) = phi_bt(i)-phi_bt(i-1);
            if(i==3)
                phi_bt_dir(1)=phi_bt_dir(2);
            end
        end
        
        if abs(sign(phi_cusp(i))+sign(phi_cusp(i-1)))==0
            %% CUSP
            fprintf('BR_BIFDET: Cusp detected near par(%d) = %.10f, par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)),...
                free_par(2), p.parameter(free_par(2)));
            [cusp,located] = locate_codim2(funcs,branch,fold,i,'cusp'); % correction
            if located
                cusp = nmfm_cusp(funcs,cusp);
                [newbranch,added] = add_to_branch(funcs,newbranch,branch,i,cusp);
                if added
                    fprintf('BR_BIFDET: Cusp located at  par(%d) = %.10f, par(%d) = %.10f.\n',...
                        free_par(1), cusp.parameter(free_par(1)),...
                        free_par(2), cusp.parameter(free_par(2)));
                    fprintf('BR_BIFDET: Normal form coefficients: b = %.10f, c = %.10f.\n\n', cusp.nmfm.b, cusp.nmfm.c);
                end
            else
                fprintf('BR_BIFDET: Failed to correct cusp.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
        
        if abs(sign(phi_bt(i))+sign(phi_bt(i-1)))==0 || (i>2 && ...
                abs(sign(phi_bt(i-1)-phi_bt(i-2))+sign(phi_bt(i)-phi_bt(i-1)))==0 && abs(phi_bt(i-1)-phi_bt(i))~=0)
            %% BT
            fprintf('BR_BIFDET: Bogdanov-Takens detected near par(%d) = %.10f, par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)),...
                free_par(2), p.parameter(free_par(2)));
            [bt,located] = locate_codim2(funcs,branch,p,i,'BT'); % correction
            if located
                [newbranch,added] = add_to_branch(funcs,newbranch,branch,i,bt);
                bt = nmfm_bt(funcs, bt);
                if added  
                    fprintf('BR_BIFDET: a = %.10f, b = %.10f, par(%d) = %.10f, par(%d) = %.10f.\n\n', ...
                        bt.nmfm.a2,bt.nmfm.b2,free_par(1), bt.parameter(free_par(1)),...
                        free_par(2),bt.parameter(free_par(2)));
                end
            else
                fprintf('BR_BIFDET: Failed to correct Bogdanov-Takens point.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
        
        if abs(sign(phi_zh(i))+sign(phi_zh(i-1)))==0 && phi_zh(i) ~= 0
            %% ZH
            fprintf('BR_BIFDET: Zero-Hopf detected near par(%d) = %.10f, par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)),...
                free_par(2), p.parameter(free_par(2)));
            [zh,located] = locate_codim2(funcs,branch,p,i,'ZH_f'); % correction
            if located
                zh = nmfm_zeho(funcs, zh, branch.point(i-1)); % normal form computation
                [newbranch,added] = add_to_branch(funcs,newbranch,branch,i,zh);
                if added
                    fprintf('BR_BIFDET: Zero Hopf located at par(%d) = %.10f, par(%d) = %.10f.\n',...
                        free_par(1), zh.parameter(free_par(1)),...
                        free_par(2), zh.parameter(free_par(2)));
                    fprintf('BR_BIFDET: s = %.10f, theta = %.10f.\n\n', zh.nmfm.s, zh.nmfm.theta);
                end
            else
                fprintf('BR_BIFDET: Failed to correct zero-Hopf.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
    elseif strcmp(kind,'hopf')
        %% Hopf branch
        p=branch.point(i);
        p=nmfm_hopf(funcs,p);
        
        par1(i) = p.parameter(free_par(1)); % needed for plotting
        par2(i) = p.parameter(free_par(2));
        
        % test functions
        phi_h_bt(i) = p.omega;
        phi_h_zh(i) = det(ch_matrix(funcs,p.x,p.parameter,0));
        phi_h_hh(i) = nmfm_smrp(funcs, p,stmethod,'remove_omega',true,'threshold',isimag);
        phi_h_gh(i) = p.nmfm.L1;
        
        if(i==1)
            continue;
        end
        
        if abs(sign(phi_h_bt(i))+sign(phi_h_bt(i-1)))==0
            %% BT
            fprintf('BR_BIFDET: Bogdanov-Takens detected near par(%d) = %.10f, par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)),...
                free_par(2), p.parameter(free_par(2)));
            [bt,located] = locate_codim2(funcs,branch,p,i,'BT'); % correction
            if located
                bt = nmfm_bt(funcs, bt);
                [newbranch,added] = add_to_branch(funcs,newbranch,branch,i,bt);
                if added
                fprintf('BR_BIFDET: a = %.10f, b = %.10f, par(%d) = %.10f, par(%d) = %.10f.\n\n', ...
                    bt.nmfm.a2,bt.nmfm.b2,free_par(1), bt.parameter(free_par(1)),...
                    free_par(2),bt.parameter(free_par(2)));
                end
            else
                fprintf('BR_BIFDET: Failed to correct Bogdanov-Takens point.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
        
        if abs(sign(phi_h_zh(i))+sign(phi_h_zh(i-1)))==0
            %% ZH
            fprintf('BR_BIFDET: Zero-Hopf detected near par(%d) = %.10f, par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)),...
                free_par(2), p.parameter(free_par(2)));
            [zh,located] = locate_codim2(funcs,branch,p,i,'ZH'); % correction
            if located
                zh = nmfm_zeho(funcs, zh, branch.point(i-1)); % normal form computation
                [newbranch,added] = add_to_branch(funcs,newbranch,branch,i,zh);
                if added
                	fprintf('BR_BIFDET: omega = %.10f, par(%d) = %.10f, par(%d) = %.10f.\n',...
                    zh.omega, free_par(1), zh.parameter(free_par(1)),...
                    free_par(2), zh.parameter(free_par(2)));
                fprintf('BR_BIFDET: s = %.10f, theta = %.10f.\n\n', zh.nmfm.s, zh.nmfm.theta);
                end
            else
                fprintf('BR_BIFDET: Failed to correct zero-Hopf.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
        
        if abs(sign(phi_h_hh(i))+sign(phi_h_hh(i-1)))==0 && phi_h_hh(i) ~=0
            %% HH
            fprintf('BR_BIFDET: Double-Hopf detected near par(%d) = %.10f, par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)),...
                free_par(2), p.parameter(free_par(2)));
            [hh,located] = locate_codim2(funcs,branch,p,i,'HH'); % correction
            if located
                hh = nmfm_hoho(funcs, hh, branch.point(i-1)); % normal form computation
                [newbranch,added] = add_to_branch(funcs,newbranch,branch,i,hh);
                if added 
                    fprintf('BR_BIFDET: omega1 = %.10f, omega2 = %.10f, par(%d) = %.10f, par(%d) = %.10f.\n',...
                    hh.omega1, hh.omega2, free_par(1),...
                    hh.parameter(free_par(1)), free_par(2),...
                    hh.parameter(free_par(2)));
                
                fprintf('BR_BIFDET: theta = %.10f, delta = %.10f.\n', hh.nmfm.theta, hh.nmfm.delta);
                fprintf('BR_BIFDET: The eigenvalues are\n');
                disp(hh.stability.l1);
                fprintf('\n');
                end
            else
                fprintf('BR_BIFDET: Failed to correct double-Hopf.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
        
        if abs(sign(phi_h_gh(i))+sign(phi_h_gh(i-1)))==0
            %% GH
            fprintf('BR_BIFDET: Generalized Hopf detected near par(%d) = %.10f, par(%d) = %.10f.\n',...
                free_par(1), p.parameter(free_par(1)),...
                free_par(2), p.parameter(free_par(2)));
            [gh,located] = locate_codim2(funcs,branch,p,i,'GH'); % correction
            if located
                gh = nmfm_genh(funcs, gh,  branch.point(i-1)); % normal form computation
                [newbranch,added] = add_to_branch(funcs,newbranch,branch,i,gh);
                if added
                                    fprintf('BR_BIFDET: L2 = %.10f, omega = %.10f, par(%d) = %.10f, par(%d) = %.10f.\n\n',...
                    gh.nmfm.L2, gh.omega, free_par(1),...
                    gh.parameter(free_par(1)), free_par(2),...
                    gh.parameter(free_par(2)));
                end
                
            else
                fprintf('BR_BIFDET: Failed to correct generalized Hopf.\n\n');
            end
            if pause_on_bifurcation
                pause;
            end
        end
        
    end
end
%% plot test functions
if plot_testfunctions
    if (strcmp(kind,'stst'))
        figure;
        plot(par1,[phi_f; phi_H]);
        hold;
        plot(par1,zeros(length(branch.point)),'r');
        h=legend('$\phi_f$','$\phi_H$');
        title('Test function plot','FontSize',15,'Interpreter','latex');
        xstr = ['free parameter ',num2str(free_par(1))];
        xlabel(xstr,'FontSize',13,'Interpreter','latex');
        set(h,'Interpreter','latex');
        set(h,'FontSize',20);
    elseif strcmp(kind,'fold')
        figure;
        plot(par1, [phi_cusp; phi_bt; phi_bt_dir; phi_zh]);
        h=legend('$\phi^f_{cusp}$','$\phi^f_{bt}$','$\phi^f_{bt\_dir}$','$\phi^f_{zh}$');
        title('Test function plot','FontSize',15,'Interpreter','latex');
        xstr = ['free parameter ',num2str(free_par(1))];
        xlabel(xstr,'FontSize',13,'Interpreter','latex');
        set(h,'Interpreter','latex');
        set(h,'FontSize',20);        
        
        figure;
        plot(par2, [phi_cusp; phi_bt; phi_bt_dir; phi_zh]);
        h=legend('$\phi^f_{cusp}$','$\phi^f_{bt}$','$\phi^f_{bt\_dir}$','$\phi^f_{zh}$');
        title('Test function plot','FontSize',15,'Interpreter','latex');
        xstr = ['free parameter ',num2str(free_par(2))];
        xlabel(xstr,'FontSize',13,'Interpreter','latex');
        set(h,'Interpreter','latex');
        set(h,'FontSize',20);        
    elseif strcmp(kind,'hopf')
        figure;
        plot(par1, [phi_h_gh; phi_h_bt; phi_h_zh; phi_h_hh]);
        h=legend('$\phi^h_{gh}$','$\phi^h_{bt}$','$\phi^h_{zh}$','$\phi^h_{hh}$');
        title('Test function plot','FontSize',15,'Interpreter','latex');
        xstr = ['free parameter ',num2str(free_par(1))];
        xlabel(xstr,'FontSize',13,'Interpreter','latex');
        set(h,'Interpreter','latex');
        set(h,'FontSize',20);        
        
        figure;
        plot(par2, [phi_h_gh; phi_h_bt; phi_h_zh; phi_h_hh]);
        h=legend('$\phi^h_{gh}$','$\phi^h_{bt}$','$\phi^h_{zh}$','$\phi^h_{hh}$');
        title('Test function plot','FontSize',15,'Interpreter','latex');
        xstr = ['free parameter ',num2str(free_par(2))];
        xlabel(xstr,'FontSize',13,'Interpreter','latex');
        set(h,'Interpreter','latex');
        set(h,'FontSize',20);
    end
end

end
