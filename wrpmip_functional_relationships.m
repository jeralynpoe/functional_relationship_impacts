%Benchmark functional relationships from WrPMIP models based on choices
%made during benchmarking process (temporal averaging, temporal extent,
%number of daily observations, and choice of dataset)

%Load in lat/lon from flux tower sites
site_info = readtable('/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Flux Towers/flux_tower_lat_lon.csv','TreatAsMissing','NA');
site_abbrev = site_info.Site;
site_lat = site_info.Latitude;
site_lon = site_info.Longitude;
site_coords = horzcat(site_lat,site_lon);

%Subset WrPMIP model from 2001-2020
start_time = datetime(2001,01,01);
end_time = datetime(2020,12,01);
time_range = (start_time:calmonths(1):end_time)';
isNotLeapDay=~(month(time_range)==2 & day(time_range)==29);
time_range=time_range(isNotLeapDay);

%Load in WrPMIP models
wrpmip_dir = '/Users/jmp838/Desktop/Research/Chapter 3 - WrPMIP Model Evaluation/WrPMIP_Outputs';
wrpmip_info = dir(fullfile(wrpmip_dir,'*.nc'));
wrpmip_nfiles = length(wrpmip_info);
wrpmip_filenames = fullfile(wrpmip_dir,{wrpmip_info.name});
wrpmip_recos_ctl = cell(wrpmip_nfiles,1);
wrpmip_temps_ctl = cell(wrpmip_nfiles,1);
wrpmip_lats = cell(wrpmip_nfiles,1);
wrpmip_lons = cell(wrpmip_nfiles,1);
for i = 1:wrpmip_nfiles
    current_wrpmip_file = wrpmip_filenames{i};

    current_timevar=ncread(current_wrpmip_file,'time');
    current_monthssince = 0:(size(current_timevar,1)-1);

    current_timeunits = ncreadatt(current_wrpmip_file,"time","units");
    current_timeunits2 = extractAfter(current_timeunits,"since ");

    time_contents = count(current_timeunits2,'-');
    if time_contents == 1
        current_timeunits2 = append(current_timeunits2,'-01');
    end

    current_timeunits_dt = datetime(current_timeunits2,'Format','dd-MMM-yyyy HH:mm:ss'); %convert time units to datetime
    time_model = datetime(current_timeunits2,'Format','01-MMM-yyyy 00:00:00');

    if i == 7
        current_ra_ctl = ncread(current_wrpmip_file,'AutoResp',[1 1 1 1],[inf inf inf 1]);
        current_rh_ctl = ncread(current_wrpmip_file,'HeteroResp',[1 1 1 1],[inf inf inf 1]);
        current_reco_ctl = current_ra_ctl + current_rh_ctl;
    else
        current_reco_ctl = ncread(current_wrpmip_file,'TotalResp',[1 1 1 1],[inf inf inf 1]);
    end
    current_temp25layer_ctl = ncread(current_wrpmip_file,'SoilTemp');  %lon x lat x layer x time x sim
    current_temp_5cm_ctl = squeeze(current_temp25layer_ctl(:,:,2,:,1));

    wrpmip_lats{i} = ncread(current_wrpmip_file,'lat');
    wrpmip_lons{i} = ncread(current_wrpmip_file,'lon');

    current_reco_ctl = squeeze(current_reco_ctl);
    current_temp_5cm_ctl = squeeze(current_temp_5cm_ctl);

    if min(current_temp_5cm_ctl,[],'all') > 200
        current_temp_5cm_ctl = current_temp_5cm_ctl - 273.15;
    end

    wrpmip_recos_ctl{i} = current_reco_ctl(:,:,:);
    wrpmip_temps_5cm_ctl{i} = current_temp_5cm_ctl(:,:,:);
end

%make sure model names here are in the same order that were loaded in above
model_names = ["CLASSIC","CLM5-ExIce","CLM5","ecosys","ELM-ECA","ELM-RD",...
    "JSBACH","JULES","LPJ-GUESS-ML","LPJ-GUESS","ORCHIDEE-MICT-teb","ORCHIDEE-MICT","UVic-ESCM"];

%for all WrPMIP models pull out lat/lon of flux tower sites and subset to
%the time ranges at each site
for i = 1:size(site_names,1)
    for j = 1:wrpmip_nfiles
        current_site_years = site_list_monthly2.(site_names{i}).RECO;
        current_ec_lat = site_lat(i);
        current_ec_lon = site_lon(i);
        wrpmip_lat = wrpmip_lats{j};
        wrpmip_lon = wrpmip_lons{j};
        current_wrpmip_reco = wrpmip_recos_ctl{j};
        current_wrpmip_temp_5cm = wrpmip_temps_5cm_ctl{j};
        if all(wrpmip_lon(:)>=0)
            mask = wrpmip_lon >= 180;
            wrpmip_lon(mask) = wrpmip_lon(mask) - 360;
        end
        wrpmip_lat_idx = wrpmip_lat - current_ec_lat;
        wrpmip_lon_idx = wrpmip_lon - current_ec_lon;
        [wrpmip_lat_val,wrpmip_lat_ix]=min(abs(wrpmip_lat_idx)); %find the location closest to zero
        [wrpmip_lon_val,wrpmip_lon_ix]=min(abs(wrpmip_lon_idx)); %find the location closest to zero
        wrpmip_reco_sites = current_wrpmip_reco(wrpmip_lon_ix,wrpmip_lat_ix,:);
        wrpmip_reco_sites = wrpmip_reco_sites(13:252); %subset from 2001-2020
        wrpmip_reco = reshape(wrpmip_reco_sites,240,1);
        wrpmip_reco(isnan(current_site_years)) = NaN;
        wrpmip_reco2 = wrpmip_reco;

        wrpmip_temp_sites = current_wrpmip_temp_5cm(wrpmip_lon_ix,wrpmip_lat_ix,:);
        wrpmip_temp_sites = wrpmip_temp_sites(13:252); %subset from 2001-2020
        wrpmip_temp = reshape(wrpmip_temp_sites,240,1);
        wrpmip_temp(isnan(current_site_years)) = NaN;

        wrpmip_reco(isnan(wrpmip_temp)) = [];
        wrpmip_temp(isnan(wrpmip_reco2)) = [];
        wrpmip_reco(isnan(wrpmip_reco)) = [];
        wrpmip_temp(isnan(wrpmip_temp)) = [];

        wrpmip_sites_temp{j} = wrpmip_temp;
        wrpmip_models_sites_temp{i} = wrpmip_sites_temp;
        wrpmip_sites_reco{j} = wrpmip_reco;
        wrpmip_models_sites_reco{i} = wrpmip_sites_reco;
    end
end

%Create functional relationships for monthly model output
for i=1:16
    for j=1:13

    current_model_sites_reco = wrpmip_models_sites_reco{i};
    current_model_site_reco = current_model_sites_reco{j};
    current_model_sites_temp = wrpmip_models_sites_temp{i};
    current_model_site_temp = current_model_sites_temp{j};

    %monthly
    if size(current_model_site_reco,1) > 2
        [wrpmip_f_monthly,wrpmip_gof1_monthly]=fit(current_model_site_temp,current_model_site_reco,'exp1');
        wrpmip_ci_monthly = confint(wrpmip_f_monthly);
        wrpmip_a_monthly = wrpmip_f_monthly.a;
        wrpmip_b_monthly = wrpmip_f_monthly.b;
        wrpmip_temp_sort_monthly = sort(current_model_site_temp);
        wrpmip_y_fit_monthly = wrpmip_a_monthly.*exp(wrpmip_b_monthly.*wrpmip_temp_sort_monthly);
        wrpmip_y_cbl_monthly = wrpmip_ci_monthly(1,1).*exp(wrpmip_ci_monthly(1,2).*wrpmip_temp_sort_monthly);
        wrpmip_y_cbh_monthly = wrpmip_ci_monthly(2,1).*exp(wrpmip_ci_monthly(2,2).*wrpmip_temp_sort_monthly);
        wrpmip_q10_monthly = exp(10.*wrpmip_f_monthly.b);
        wrpmip_cil_monthly = exp(10.*wrpmip_ci_monthly(1,2));
        wrpmip_ciu_monthly = exp(10.*wrpmip_ci_monthly(2,2));
    else
        wrpmip_f_monthly = NaN;
        wrpmip_temp_sort_monthly = NaN;
        wrpmip_temp_monthly = NaN;
        wrpmip_reco_monthly = NaN;
        wrpmip_y_cbl_monthly = NaN;
        wrpmip_y_cbh_monthly = NaN;
        wrpmip_q10_monthly = NaN;
        wrpmip_cil_monthly = NaN;
        wrpmip_ciu_monthly = NaN;
    end

    wrpmip_fs_sites{j} = wrpmip_f_monthly;
    wrpmip_temp_sorts_sites{j} = wrpmip_temp_sort_monthly;
    wrpmip_temps_sites{j} = current_model_site_temp;
    wrpmip_recos_sites{j} = current_model_site_reco;
    wrpmip_cbls_sites{j} = wrpmip_y_cbl_monthly;
    wrpmip_cbhs_sites{j} = wrpmip_y_cbh_monthly;
    wrpmip_q10s_sites{j} = wrpmip_q10_monthly;
    wrpmip_cils_sites{j} = wrpmip_cil_monthly;
    wrpmip_cius_sites{j} = wrpmip_ciu_monthly;

    wrpmip_fs_monthly{i} = wrpmip_fs_sites;
    wrpmip_temp_sorts_monthly{i} = wrpmip_temp_sorts_sites;
    wrpmip_temps_monthly{i} = wrpmip_temps_sites;
    wrpmip_recos_monthly{i} = wrpmip_recos_sites;
    wrpmip_cbls_monthly{i} = wrpmip_cbls_sites;
    wrpmip_cbhs_monthly{i} = wrpmip_cbhs_sites;
    wrpmip_q10s_monthly{i} = wrpmip_q10s_sites;
    wrpmip_cils_monthly{i} = wrpmip_cils_sites;
    wrpmip_cius_monthly{i} = wrpmip_cius_sites;
    end
end

%Create supplemental figure of model q10s compared to eml (fig S8)
% wrpmip_q10s_data = horzcat(wrpmip_q10s,q10s_daily{10});
% wrpmip_cils_data = horzcat(wrpmip_cils,cils_daily{10});
% wrpmip_cius_data = horzcat(wrpmip_cius,cius_daily{10});
% model_names = {"CLASSIC","CLM5-ExIce","CLM5","ecosys","ELM-ECA","ELM-RD",...
%     "JSBACH","JULES","LPJ-GUESS-ML","LPJ-GUESS","ORCHIDEE-MICT-teb","ORCHIDEE-MICT","UVic-ESCM","Data"};
% [m,ix] = sort(wrpmip_q10s_data); %sort based on median
% models_sort = model_names(ix);
% cils_sort=wrpmip_cils_data(ix); cius_sort=wrpmip_cius_data(ix);
% 
% hold on
% errorbar(1:14,m,m-cils_sort,cius_sort-m,'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
%         'MarkerSize',10);
% errorbar(5,m(5),m(5)-cils_sort(5),cius_sort(5)-m(5),'o','color','r','linewidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r',...
%         'MarkerSize',10);
% xlim([0.25 14.75])
% ylim([0.5 6.5])
% set(gca,'xtick',1:1:14);
% set(gca,'TickLabelInterpreter','none')
% set(gca,'xticklabel',[models_sort],'fontsize',14)
% xtickangle(30)
% ylabel('Intercept','fontsize',14);
% box on;
% ylabel('Q_{10} (unitless)','fontsize',14);

%Create supplemental figure of model q10s compared to 16 flux tower sites
model_names = {"CLASSIC","CLM5-ExIce","CLM5","ecosys","ELM-ECA","ELM-RD",...
    "JSBACH","JULES","LPJ-GUESS-ML","LPJ-GUESS","ORCHIDEE-MICT-teb","ORCHIDEE-MICT","UVic-ESCM","Data"};
tiledlayout(4,4,"TileSpacing","compact")
for i = 1:16
    %for j = 1:13
        nexttile

        wrpmip_q10s_data = horzcat(wrpmip_q10s_monthly{i}{:},q10s_daily{i});
        wrpmip_cils_data = horzcat(wrpmip_cils_monthly{i}{:},cils_daily{i});
        wrpmip_cius_data = horzcat(wrpmip_cius_monthly{i}{:},cius_daily{i});
        [m,ix] = sort(wrpmip_q10s_data); %sort based on median
        models_sort = model_names(ix);
        cils_sort=wrpmip_cils_data(ix); cius_sort=wrpmip_cius_data(ix);

        hold on
        errorbar(1:14,m,m-cils_sort,cius_sort-m,'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',6);
        %find location of data
        data_str = find(cellfun(@(x) x == "Data", models_sort));
        errorbar(data_str,m(data_str),m(data_str)-cils_sort(data_str),cius_sort(data_str)-m(data_str),'o','color','r','linewidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'MarkerSize',6);
        xlim([0.25 14.75])
        set(gca,'xtick',1:1:14);
        set(gca,'TickLabelInterpreter','none')
        set(gca,'xticklabel',[models_sort],'fontsize',9)
        xtickangle(30)
        ylabel('Intercept','fontsize',9);
        box on;
        title(site_names_full{i},'FontSize',9)
        ylabel('Q_{10} (unitless)','fontsize',9);
end

%find min and max q10 from all models across all sites
min_q10 = min(cellfun(@(x) min(cell2mat(x)), wrpmip_q10s_monthly)); %1.1575
max_q10 = max(cellfun(@(x) max(cell2mat(x)), wrpmip_q10s_monthly)); %5.4902
wrpmip_q10s_monthly2 = cellfun(@(x) cell2mat(x), wrpmip_q10s_monthly, 'UniformOutput', false);
wrpmip_q10s_monthly2 = cell2mat(wrpmip_q10s_monthly2);
mean_q10 = nanmean(wrpmip_q10s_monthly2(:)); %2.4376

%difference between model q10 and daily q10 at each location
for i = 1:16
    models_site_diff = cell2mat(wrpmip_q10s_monthly{i}) - q10s_daily{i};
    models_data_diff{i} = models_site_diff;
end

%extract sites from temporal averaging
for i = 1:16
    time_avg_site = [q10s_og{i},q10s_daily{i},q10s_weekly{i},q10s_monthly{i},q10s_yearly{i}];
    time_avg_sites{i} = time_avg_site;
end

%calculate change in temporal extent
for i = 1:16
    time_ext_site = q10_yrs_sites{i};
    time_ext_site2 = cell2mat(vertcat(time_ext_site{:}));
    time_ext_site3 = time_ext_site2(:);
    min_time_ext_q10s = min(time_ext_site2,[],1);
    max_time_ext_q10s = max(time_ext_site2,[],1);
    min_sites_time_ext_q10s{i} = min_time_ext_q10s;
    max_sites_time_ext_q10s{i} = max_time_ext_q10s;
    time_ext_sites{i} = time_ext_site3;
end

load("/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Code/site_data_randomization_5cm_0703.mat",...
    "eml_fs","eml_q10s","eml_cils","eml_cius","eml_q10_counts",...
    "man_q10s","oas_q10s","obs_q10s","ojp_q10s","scb_q10s","scc_q10s","sj2_q10s",...
    "bzf_q10s","bzs_q10s","fcr_q10s","ich_q10s","ics_q10s","prr_q10s","rpf_q10s","uaf_q10s")

sites_rdm_q10s = {man_q10s,oas_q10s,obs_q10s,ojp_q10s,scb_q10s,scc_q10s,sj2_q10s,...
    bzf_q10s,bzs_q10s,eml_q10s,fcr_q10s,ich_q10s,ics_q10s,prr_q10s,rpf_q10s,uaf_q10s};

for i = 1:16
    current_site_rdm_q10 = sites_rdm_q10s{i};
    current_site_q10s2 = cell2mat(current_site_rdm_q10);
    current_site_q10s3 = current_site_q10s2(:);
    min_site_rdm_q10s = min(current_site_q10s2,[],1);
    max_site_rdm_q10s = max(current_site_q10s2,[],1);
    min_sites_rdm_q10s{i} = min_site_rdm_q10s; %calculate 1 min for each day
    max_sites_rdm_q10s{i} = max_site_rdm_q10s; %calculate 1 max for each day
    rdm_q10_sites{i} = current_site_q10s3;
end

%calculate model scores based on level of agreement in q10 between models
%and observations
for h = 1:16
    wrpmip_site_q10_monthly = wrpmip_q10s_monthly{h};
    time_avg_site = time_avg_sites{h};
    current_time_ext_site = time_ext_sites{h};
    current_site_rdm_q10s = rdm_q10_sites{h};
    current_site_tower_q10 = sites_tower_q10s{h};
    current_site_gridded_q10 = sites_gridded_q10s{h};
for i = 1:wrpmip_nfiles
    models_data_diff_score(i) = exp(-norm(wrpmip_site_q10_monthly{i} - q10s_daily{10})./norm(wrpmip_site_q10_monthly{i}));
    wrpmip_time_avg_hh(i) = exp(-norm(wrpmip_site_q10_monthly{i} - time_avg_site(1))./norm(wrpmip_site_q10_monthly{i}));
    wrpmip_time_avg_daily(i) = exp(-norm(wrpmip_site_q10_monthly{i} - time_avg_site(2))./norm(wrpmip_site_q10_monthly{i}));
    wrpmip_time_avg_weekly(i) = exp(-norm(wrpmip_site_q10_monthly{i} - time_avg_site(3))./norm(wrpmip_site_q10_monthly{i}));
    wrpmip_time_avg_monthly(i) = exp(-norm(wrpmip_site_q10_monthly{i} - time_avg_site(4))./norm(wrpmip_site_q10_monthly{i}));
    wrpmip_time_avg_annual(i) = exp(-norm(wrpmip_site_q10_monthly{i} - time_avg_site(5))./norm(wrpmip_site_q10_monthly{i}));

    time_ext_score = [];
    for j = 1:size(current_time_ext_site,1)
        current_time_ext_score = exp(-norm(wrpmip_site_q10_monthly{i} - current_time_ext_site(j))./norm(wrpmip_site_q10_monthly{i}));
        time_ext_score{j} = current_time_ext_score;
        wrpmip_time_ext{i} = time_ext_score;
        wrpmip_sites_time_ext{h} = wrpmip_time_ext;
    end

    rdm_q10_score = [];
    for j = 1:15000
        current_rdm_q10_score = exp(-norm(wrpmip_site_q10_monthly{i} - current_site_rdm_q10s(j))./norm(wrpmip_site_q10_monthly{i}));
        rdm_q10_score{j} = current_rdm_q10_score;
        wrpmip_rdm_q10{i} = rdm_q10_score;
        wrpmip_sites_rdm_q10s{h} = wrpmip_rdm_q10;
    end

    for m = 1:6
        if h == 10
            model_tower_scores(m) = exp(-norm(wrpmip_site_q10_monthly{i} - current_site_tower_q10(m))./norm(wrpmip_site_q10_monthly{i}));
            model_chamber_scores(m) = exp(-norm(wrpmip_site_q10_monthly{i} - chamber_q10s(m))./norm(wrpmip_site_q10_monthly{i}));
            model_gridded_scores(m) = exp(-norm(wrpmip_site_q10_monthly{i} - current_site_gridded_q10(m))./norm(wrpmip_site_q10_monthly{i}));
        else
            model_tower_scores(m) = exp(-norm(wrpmip_site_q10_monthly{i} - current_site_tower_q10(m))./norm(wrpmip_site_q10_monthly{i}));
            model_chamber_scores(m) = nan;
            model_gridded_scores(m) = exp(-norm(wrpmip_site_q10_monthly{i} - current_site_gridded_q10(m))./norm(wrpmip_site_q10_monthly{i}));
        end
    end

    wrpmip_tower_scores{i} = model_tower_scores;
    wrpmip_tower_scores1{i} = model_tower_scores(1);
    wrpmip_chamber_scores{i} = model_chamber_scores;
    wrpmip_gridded_scores{i} = model_gridded_scores;

    models_sites_data_diff_score{h} = models_data_diff_score;
    wrpmip_sites_time_avg_hh{h} = wrpmip_time_avg_hh;
    wrpmip_sites_time_avg_daily{h} = wrpmip_time_avg_daily;
    wrpmip_sites_time_avg_weekly{h} = wrpmip_time_avg_weekly;
    wrpmip_sites_time_avg_monthly{h} = wrpmip_time_avg_monthly;
    wrpmip_sites_time_avg_annual{h} = wrpmip_time_avg_annual;
    wrpmip_sites_tower_scores{h} = wrpmip_tower_scores;
    wrpmip_sites_tower_scores1{h} = wrpmip_tower_scores1;
    wrpmip_sites_chamber_scores{h} = wrpmip_chamber_scores;
    wrpmip_sites_gridded_scores{h} = wrpmip_gridded_scores;
end
end

%find min and max model scores across all models and sites
min_score = min(cell2mat(models_sites_data_diff_score)); %0.3128
max_score = max(cell2mat(models_sites_data_diff_score)); %0.9954
mean_score = nanmean(cell2mat(models_sites_data_diff_score)); %0.7794

models_sites_scores = cell2mat(models_sites_data_diff_score');

for i = 1:16
    model_time_avg = [wrpmip_sites_time_avg_hh{i};wrpmip_sites_time_avg_daily{i};wrpmip_sites_time_avg_weekly{i};wrpmip_sites_time_avg_monthly{i}];
    wrpmip_time_avgs{i} = model_time_avg;
end

%calculate range in model scores
for h = 1:16
for i = 1:wrpmip_nfiles

    current_wrpmip_time_avgs = wrpmip_time_avgs{h};
    current_wrpmip_time_ext = wrpmip_sites_time_ext{h};
    current_wrpmip_rdm_q10 = wrpmip_sites_rdm_q10s{h};
    current_wrpmip_chamber_scores = wrpmip_sites_chamber_scores{h};
    current_wrpmip_tower_scores = wrpmip_sites_tower_scores{h};
    current_wrpmip_gridded_scores = wrpmip_sites_gridded_scores{h};

    current_time_avg_max = max(current_wrpmip_time_avgs(:,i));
    current_time_avg_min = min(current_wrpmip_time_avgs(:,i));

    current_site_time_ext = current_wrpmip_time_ext{i};
    current_time_ext_min = min([current_site_time_ext{:}]);
    current_time_ext_max = max([current_site_time_ext{:}]);

    current_site_rdm_q10 = current_wrpmip_rdm_q10{i};
    current_site_rdm_q102 = cell2mat(current_site_rdm_q10);
    current_rdm_q10_min = min(current_site_rdm_q102);
    current_rdm_q10_max = max(current_site_rdm_q102);

    current_dset_choices = horzcat(current_wrpmip_chamber_scores{i},current_wrpmip_tower_scores{i},current_wrpmip_gridded_scores{i});
    current_dset_choice_min = min(current_dset_choices);
    current_dset_choice_max = max(current_dset_choices);

    %calculate range for each model and category
    range_time_avg = current_time_avg_max-current_time_avg_min;
    range_time_ext = current_time_ext_max-current_time_ext_min;
    range_rdm_q10 = current_rdm_q10_max-current_rdm_q10_min;
    range_dset_choice = current_dset_choice_max-current_dset_choice_min;

    ranges_time_avg{i} = range_time_avg;
    ranges_time_ext{i} = range_time_ext;
    ranges_rdm_q10{i} = range_rdm_q10;
    ranges_dset_choice{i} = range_dset_choice;

    ranges_sites_time_avg{h} = ranges_time_avg;
    ranges_sites_time_ext{h} = ranges_time_ext;
    ranges_sites_rdm_q10{h} = ranges_rdm_q10;
    ranges_sites_dset_choice{h} = ranges_dset_choice;

end
end

%plot range in model skill for each category
model_color=[5,5,5,10,10,10,15,15,15,22,22,22,25];
markers={'diamond','square','^','diamond','square','^','diamond','square','^',...
    'diamond','square','^','diamond'};

for i = 1:16
avg_line = horzcat(nanmean(cell2mat(ranges_sites_time_avg{i})),nanmean(cell2mat(ranges_sites_dset_choice{i})),...
    nanmean(cell2mat(ranges_sites_time_ext{i})),nanmean(cell2mat(ranges_sites_rdm_q10{i})));
avg_lines_sites{i} = avg_line;
end

avg_line_all_sites = cell2mat(avg_lines_sites');
avg_line_all_sites = mean(avg_line_all_sites,1);

%range in scores plot for all sites
model_names = {"CLASSIC","CLM5-ExIce","CLM5","ecosys","ELM-ECA","ELM-RD",...
    "JSBACH","JULES","LPJ-GUESS-ML","LPJ-GUESS","ORCHIDEE-MICT-teb","ORCHIDEE-MICT","UVic-ESCM"};
tiledlayout(4,4,"TileSpacing","compact")
jitter_amount = 0.6;
for h = 1:16
    nexttile
    current_wrpmip_time_avg_range = ranges_sites_time_avg{h};
    current_wrpmip_time_ext_range = ranges_sites_time_ext{h};
    current_wrpmip_rdm_q10_range = ranges_sites_rdm_q10{h};
    current_wrpmip_dset_choice_range = ranges_sites_dset_choice{h};
hold on
for i = 1:13
    x1 = 1 + (rand - 0.5) * jitter_amount;
    x2 = 2 + (rand - 0.5) * jitter_amount;
    x3 = 3 + (rand - 0.5) * jitter_amount;
    x4 = 4 + (rand - 0.5) * jitter_amount;

    swarmchart(x1,cell2mat(current_wrpmip_time_avg_range(i)), 70, model_color(i,:), ...
        'filled','MarkerFaceAlpha', 1,'Marker', markers{i}, 'MarkerEdgeColor','k');
    swarmchart(x2,cell2mat(current_wrpmip_dset_choice_range(i)), 70, model_color(i,:), ...
        'filled','MarkerFaceAlpha', 1,'Marker', markers{i}, 'MarkerEdgeColor','k');
    swarmchart(x3,cell2mat(current_wrpmip_time_ext_range(i)), 70, model_color(i,:), ...
        'filled','MarkerFaceAlpha', 1,'Marker', markers{i}, 'MarkerEdgeColor','k');
    swarmchart(x4,cell2mat(current_wrpmip_rdm_q10_range(i)), 70, model_color(i,:), ...
        'filled','MarkerFaceAlpha', 1,'Marker', markers{i}, 'MarkerEdgeColor','k');
end
plot(avg_lines_sites{h},'k.','MarkerSize',32)
% legend("CLASSIC","","","","CLM5-ExIce","","","","CLM5","","","","ecosys","","","",...
%     "ELM-ECA","","","","ELM-RD","","","","JSBACH","","","","JULES","","","",...
%     "LPJ-GUESS-ML","","","","LPJ-GUESS","","","","ORCHIDEE-MICT-teb","","","","ORCHIDEE-MICT",...
%     "","","","UVic-ESCM","",'location','northwest')
grid on; box on;
xticks(1:4)
if h < 13
    set(gca,'XTickLabel',[],'TickLength',[0 0]);
else
    xticklabels({'Temporal Averaging','Choice of Dataset','Temporal Extent','Number of Observations'})
end
xtickangle(20)
xlim([0.3 4.7])
ylim([0 1])
title(site_names_full{h})
fontsize(12,"points")
end
for i = 1:numel(model_names)
    hLegend(i) = plot(nan, nan, markers{i}, ...
        'MarkerFaceColor', model_color(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, ...
        'LineStyle', 'none', ...
        'DisplayName', model_names{i});
end
lg = legend(hLegend, model_names, 'NumColumns', numel(model_names), ...
            'Location', 'southoutside');
lg.Layout.Tile = 'south'; 

%avg plot of all plots and sites
mat_time_avg = cell2mat(vertcat(ranges_sites_time_avg{:}));
mat_time_ext = cell2mat(vertcat(ranges_sites_time_ext{:}));
mat_rdm_q10 = cell2mat(vertcat(ranges_sites_rdm_q10{:}));
mat_dset_choice = cell2mat(vertcat(ranges_sites_dset_choice{:}));
jitter_amount = 0.6;

model_color=[5,5,5,10,10,10,15,15,15,22,22,22,25];
markers={'diamond','square','^','diamond','square','^','diamond','square','^',...
    'diamond','square','^','diamond'};
model_color = [110/255 92/255 190/255;
    110/255 92/255 190/255;
    110/255 92/255 190/255;
    102/255 155/255.6 253/255;
    102/255 155/255.6 253/255;
    102/255 155/255.6 253/255;
    78/255 206/255 202/255;
    78/255 206/255 202/255;
    78/255 206/255 202/255;
    254/255 209/255 107/255;
    254/255 209/255 107/255;
    254/255 209/255 107/255;
    251/255 251/255 80/255];
for h = 1:16
    current_wrpmip_time_avg_range = ranges_sites_time_avg{h};
    current_wrpmip_time_ext_range = ranges_sites_time_ext{h};
    current_wrpmip_rdm_q10_range = ranges_sites_rdm_q10{h};
    current_wrpmip_dset_choice_range = ranges_sites_dset_choice{h};
hold on
for i = 1:13
    x1 = 1 + (rand - 0.5) * jitter_amount;
    x2 = 2 + (rand - 0.5) * jitter_amount;
    x3 = 3 + (rand - 0.5) * jitter_amount;
    x4 = 4 + (rand - 0.5) * jitter_amount;

    swarmchart(x1,cell2mat(current_wrpmip_time_avg_range(i)), 70, model_color(i,:), ...
        'filled','MarkerFaceAlpha', 1,'Marker', markers{i},'MarkerEdgeColor','k');
    swarmchart(x2,cell2mat(current_wrpmip_dset_choice_range(i)), 70, model_color(i,:), ...
        'filled','MarkerFaceAlpha', 1,'Marker', markers{i}, 'MarkerEdgeColor','k');
    swarmchart(x3,cell2mat(current_wrpmip_time_ext_range(i)), 70, model_color(i,:), ...
        'filled','MarkerFaceAlpha', 1,'Marker', markers{i}, 'MarkerEdgeColor','k');
    swarmchart(x4,cell2mat(current_wrpmip_rdm_q10_range(i)), 70, model_color(i,:), ...
        'filled','MarkerFaceAlpha', 1,'Marker', markers{i}, 'MarkerEdgeColor','k');
end
grid on; box on;
xticks(1:4)
if h < 13
    set(gca,'XTickLabel',[],'TickLength',[0 0]);
else
    xticklabels({'Temporal Averaging','Choice of Dataset','Temporal Extent','Number of Observations'})
end
ylabel("Range in Inferred Model Skill (unitless)")
xtickangle(20)
xlim([0.3 4.7])
ylim([0 1])
fontsize(12,"points")
end
boxchart(ones(208,1)*1,mat_time_avg(:),'BoxFaceColor','k','MarkerStyle','o','MarkerColor','k')
boxchart(ones(208,1)*2,mat_dset_choice(:),'BoxFaceColor','k','MarkerStyle','o','MarkerColor','k')
boxchart(ones(208,1)*3,mat_time_ext(:),'BoxFaceColor','k','MarkerStyle','o','MarkerColor','k')
boxchart(ones(208,1)*4,mat_rdm_q10(:),'BoxFaceColor','k','MarkerStyle','o','MarkerColor','k')
for i = 1:numel(model_names)
    hLegend(i) = plot(nan, nan, markers{i}, ...
        'MarkerFaceColor', model_color(i,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, ...
        'LineStyle', 'none', ...
        'DisplayName', model_names{i});
end
legend(hLegend, model_names, 'Location', 'northwest');

rdm_q10_median = median(mat_rdm_q10(:),'omitmissing');
rdm_q10_mean = mean(mat_rdm_q10(:),'omitmissing');
time_ext_median = median(mat_time_ext(:),'omitmissing');
dset_choice_median = median(mat_dset_choice(:),'omitmissing');
time_avg_median = median(mat_time_avg(:),'omitmissing');

%plot range in model skill from benchmarking choices for each model
model_names2 = ["CLASSIC","CLM5-ExIce","CLM5","ecosys","ELM-ECA","ELM-RD",...
    "JSBACH","JULES","LPJ-GUESS-ML","LPJ-GUESS","ORCHIDEE-MICT-teb","ORCHIDEE-MICT","UVic-ESCM"];

% Convert color code to 1-by-3 RGB array (0~1 each)
blue_str = '#1AADDB';
blue = sscanf(blue_str(2:end),'%2x%2x%2x',[1 3])/255;
yellow_str = '#F6DF28';
yellow = sscanf(yellow_str(2:end),'%2x%2x%2x',[1 3])/255;
green_str = '#66a182';
green = sscanf(green_str(2:end),'%2x%2x%2x',[1 3])/255;

tiledlayout(5,3,"TileSpacing","compact");
for i = 1:wrpmip_nfiles
    nexttile
    b=[1,8,16,28];
    c=[1,7,15,21,28];
    d=flip(1:30);
    hold on
    p1=swarmchart(ones(1,4),wrpmip_time_avgs(:,i),80,b,'filled','MarkerEdgeColor','k');
    p2=swarmchart(ones(1,5)*2,wrpmip_min_time_ext{i},80,c,'filled','MarkerEdgeColor','k');
    p3=swarmchart(ones(1,5)*2,wrpmip_max_time_ext{i},80,c,'filled','MarkerEdgeColor','k');
    p4=swarmchart(ones(1,30)*3,wrpmip_min_rdm_q10s{i},80,d,'filled','MarkerEdgeColor','k');
    p5=swarmchart(ones(1,30)*3,wrpmip_max_rdm_q10s{i},80,d,'filled','MarkerEdgeColor','k');
    p6=swarmchart(ones(1,6)*4,wrpmip_chamber_scores{i}',80,navy,'filled','MarkerEdgeColor','k');
    p7=swarmchart(ones(1,6)*4,wrpmip_tower_scores{i}',80,blue,'filled','MarkerEdgeColor','k');
    p8=swarmchart(ones(1,6)*4,wrpmip_gridded_scores{i}',80,yellow,'filled','MarkerEdgeColor','k');
    xtickangle(20)
    xlim([0.5 4.5])
    ylim([0.0 1])
    grid on; box on;
    p1.XJitterWidth = 0.3;p2.XJitterWidth = 0.3;p3.XJitterWidth = 0.3;p4.XJitterWidth = 0.3;
    p5.XJitterWidth = 0.3;p6.XJitterWidth = 0.3;p7.XJitterWidth = 0.3;p8.XJitterWidth = 0.3;
    title(model_names2{i},'FontSize',14)
    ax2 = gca; 
    ax2.FontSize = 13;
    if i == 1
        ylabel("Q_{10} Score",'FontSize',14)
        set(gca,'XTickLabel',[],'TickLength',[0 0]);
    elseif (i>1) && (i<4)
        set(gca,'XTickLabel',[],'YTickLabel',[],'TickLength',[0 0]);
    elseif i == 4
        ylabel("Q_{10} Score",'FontSize',14)
        set(gca,'XTickLabel',[],'TickLength',[0 0]);
    elseif (i>4) && (i<7)
        set(gca,'XTickLabel',[],'YTickLabel',[],'TickLength',[0 0]);
    elseif i == 7
        ylabel("Q_{10} Score",'FontSize',14)
        set(gca,'XTickLabel',[],'TickLength',[0 0]);
    elseif (i>7) && (i<10)
        set(gca,'XTickLabel',[],'YTickLabel',[],'TickLength',[0 0]);
    elseif i == 10
        ylabel("Q_{10} Score",'FontSize',14)
        set(gca,'XTickLabel',[],'TickLength',[0 0]);
    elseif (i>10) && (i<13)
        set(gca,'XTickLabel',[],'YTickLabel',[],'TickLength',[0 0]);
    elseif i == 13
        ylabel("Q_{10} Score",'FontSize',14)
        xticks(1:1:5)
        xticklabels({'Temporal Averaging','Temporal Extent','Data Amount','Choice of Dataset'})
        cb = colorbar('southoutside','YDir','reverse');
        cb.Ticks = linspace(5,25,2);
        cb.TickLabels = {"Finer Benchmark Resolution","Coarser Benchmark Resolution"};
    end
end
