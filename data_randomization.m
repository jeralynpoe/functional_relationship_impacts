function [f_samples,q10_samples,cil_samples,ciu_samples,final_q10_counts,random_temps,random_recos] = data_randomization(site_list_daily,site_list_monthly,sites_cil,sites_ciu,site_name,site_num)

site_months = site_list_monthly.(site_name).Month;
site_years = site_list_monthly.(site_name).Year;

npoints_month = (1:30);
current_plot = site_list_daily.(site_name);
for i = 1:500
    for j = 1:size(npoints_month,2)
        for k = 1:size(site_months,1)
        current_month = site_months(k);
        current_year = site_years(k);
        idx = find(current_plot.Year == current_year & current_plot.Month == current_month);
        eml_month_year = current_plot(idx,:);

        %randomize data from current month and year
        rndIDX = randperm(size(eml_month_year,1));

        %if the number of days available is greater than the number of days
        %being evaluated then extract those randomized data points
        %otherwise if the number of days available is less than then number
        %of days being evaluated, then extract all the data points
        %available
        if size(eml_month_year,1) >= npoints_month(j)
            points_per_month = eml_month_year(rndIDX(1:npoints_month(j)), :); 
        else
            points_per_month = eml_month_year(rndIDX(1:size(eml_month_year,1)), :); 
        end

        points_per_month_all{k} = points_per_month;
        points_per_month_all_table = vertcat(points_per_month_all{:});
        data_sample_months{i,j} = points_per_month_all_table;

        end
    end 
end

%construct functional relationships for each random sample
for n = 1:size(data_sample_months,1)
    current_iteration = data_sample_months(n,:);
    for m = 1:size(data_sample_months,2)
        current_random_sample = current_iteration{m};
        current_random_reco = current_random_sample.RECO;
        current_random_temp = current_random_sample.TS_5cm;

        [current_f,current_gof1]=fit(current_random_temp,current_random_reco,'exp1');
        
        f_samples{n,m} = current_f;

        current_q10 = exp(10.*current_f.b);
        q10_samples{n,m} = current_q10;

        current_ci = confint(current_f);
        current_cil = exp(10.*current_ci(1,2));
        current_ciu = exp(10.*current_ci(2,2));

        cil_samples{n,m} = current_cil;
        ciu_samples{n,m} = current_ciu;

        random_temp1{m} = current_random_temp;
        random_reco1{m} = current_random_reco;

        random_temps{n} = random_temp1;
        random_recos{n} = random_reco1;
        
    end
end

cil_daily = sites_cil{site_num};
ciu_daily = sites_ciu{site_num};

%Calculate percentage of Q10s that fall within uncertainty of original data
for p = 1:size(q10_samples,2)
    current_set = q10_samples(:,p);
    current_set = cell2mat(current_set);

    current_q10_range = (current_set < cil_daily) | (current_set > ciu_daily);
    count_true = sum(current_q10_range(:) == 0);
    final_q10_counts(p) = (count_true./500).*100;
end

end