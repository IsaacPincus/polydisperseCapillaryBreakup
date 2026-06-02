function [tau_Simon, longest_relax, longest_relax_EC,...
    t_c, tau_EC, longest_relax_min_ep,...
    tau_EC_fit,tau_EC_gt_v,t_breakup] = ...
    getTimes(a,Az,Ar,f,stress,e,Wi,t,te,ye,ie,...
    chi, mu, N, L, lambda, g, a0,...
    graphs_on)

    eta_s = mu;
    
    ordered_elastic_viscous_stress_ratio = 2/(3*eta_s)*cumsum(g.*lambda.*L.^2);
    tau_Simon = lambda(find(ordered_elastic_viscous_stress_ratio>1,1));
    
    local_e_maxes = find(islocalmax(e)==1);
    t_max_e_all = t(islocalmax(e)==1);
    t_max_e = t_max_e_all(1);
    [~, max_stress_mode] = max(sum(stress(:,t>t_max_e), 2));
    longest_relax = lambda(max_stress_mode);
    if length(local_e_maxes)>=2
        [~, max_stress_mode_EC] = max(sum(stress(:,local_e_maxes(1):local_e_maxes(2)), 2));
        longest_relax_EC = lambda(max_stress_mode_EC);
    else
        longest_relax_EC = nan;
    end
    % % t_c = 
    % t_c = t(find((e*longest_relax<2).*(t>t_max_e)',1));
    indices_EC = (e*longest_relax<5).*(t>t_max_e)'==1;
    if ~(sum(indices_EC)==0)
        times_EC = t(indices_EC);
        t_c = times_EC(1);
    
        index_min_epsilon = find((islocalmin(e)==1).*(t>t_c)',1);
        ep_min = e(index_min_epsilon);
        tau_EC = (2/3)*1/ep_min;
    
    
        [~,highest_stress_mode_min_ep] = max(stress(:,index_min_epsilon));
        longest_relax_min_ep = lambda(highest_stress_mode_min_ep);
        
        indices_fit = find(indices_EC);
        indices_fit = indices_fit(1:round(length(indices_fit)/3));
        [xxx,yyy] = prepareCurveData(t(indices_fit)-t_c,a(indices_fit)/a0);
        fitobject = fit(xxx, yyy, 'exp1');
        tau_EC_fit = -1./fitobject.b;
        
        % sum stresses during EC regime
        stresses_EC = stress(:,indices_EC);
        summed_stress_EC = sum(stresses_EC, 2);
        ordered_stresses_EC = cumsum(summed_stress_EC);
        summed_viscous_stress_EC = sum(3*mu*e(indices_EC));
        E_gt_v_index = find(ordered_stresses_EC>summed_viscous_stress_EC,1);
        if ~isempty(E_gt_v_index)
            tau_EC_gt_v = lambda(E_gt_v_index);
        else
            tau_EC_gt_v = nan;
        end
    
    else
        t_c = t_max_e;
        tau_EC = nan;
        tau_EC_fit = nan;
        tau_EC_gt_v = nan;
        longest_relax_min_ep = nan;
        fitobject = nan;
    end
    
    t_breakup = te(1);

    % fit to EC regime
    if ~(sum(indices_EC)==0)&&graphs_on
        figure();
        hold on
        plot(fitobject,xxx,yyy);
        axes1 = gca;
        axes1.YScale = 'log';
    end

end