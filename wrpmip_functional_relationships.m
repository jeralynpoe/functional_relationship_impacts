%Benchmark functional relationships from WrPMIP models based on choices
%made during benchmarking process (temporal averaging, temporal extent,
%number of daily observations, and choice of dataset)

%Load in lat/lon from flux tower sites
site_info = readtable('/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Flux Towers/flux_tower_lat_lon.csv','TreatAsMissing','NA');
site_abbrev = site_info.Site;
site_lat = site_info.Latitude;
site_lon = site_info.Longitude;

%Subset WrPMIP model from 2015-2019
start_time = datetime(2015,01,01);
end_time = datetime(2019,12,01);
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
model_names = ["CLASSIC","CLM5ExIce","CLM5","ecosys","ELMECA","ELMRD",...
    "JSBACH","JULES","LPJGUESSML","LPJGUESS","ORCHIDEEMICTteb","ORCHIDEEMICT","UVicESCM"];

%for all WrPMIP models pull out lat/lon of US-EML and subset from 2015-2019
for i = 1:wrpmip_nfiles
        current_ec_lat = site_lat(10);
        current_ec_lon = site_lon(10);
        wrpmip_lat = wrpmip_lats{i};
        wrpmip_lon = wrpmip_lons{i};
        current_wrpmip_reco = wrpmip_recos_ctl{i};
        current_wrpmip_temp_5cm = wrpmip_temps_5cm_ctl{i};
        if all(wrpmip_lon(:)>=0)
            mask = wrpmip_lon >= 180;
            wrpmip_lon(mask) = wrpmip_lon(mask) - 360;
        end
        wrpmip_lat_idx = wrpmip_lat - current_ec_lat;
        wrpmip_lon_idx = wrpmip_lon - current_ec_lon;
        [wrpmip_lat_val,wrpmip_lat_ix]=min(abs(wrpmip_lat_idx)); %find the location closest to zero
        [wrpmip_lon_val,wrpmip_lon_ix]=min(abs(wrpmip_lon_idx)); %find the location closest to zero
        wrpmip_reco_sites = current_wrpmip_reco(wrpmip_lon_ix,wrpmip_lat_ix,:);
        wrpmip_reco_sites = wrpmip_reco_sites(181:240); %subset from 2015-2019
        wrpmip_reco = reshape(wrpmip_reco_sites,60,1);

        wrpmip_temp_sites = current_wrpmip_temp_5cm(wrpmip_lon_ix,wrpmip_lat_ix,:);
        wrpmip_temp_sites = wrpmip_temp_sites(181:240); %subset from 2015-2019
        wrpmip_temp = reshape(wrpmip_temp_sites,60,1);

        wrpmip_table = table(time_range,wrpmip_reco,wrpmip_temp);
        wrpmip_table_sites{i} = wrpmip_table;
end
wrpmip_struct = cell2struct(wrpmip_table_sites,model_names,2);

%convert each model table to timetable
for i=1:wrpmip_nfiles
    [y,m,d] = ymd(wrpmip_struct.(model_names{i}).time_range);
    wrpmip_struct.(model_names{i}).Year = y;
    wrpmip_struct.(model_names{i}).Month = m;
    wrpmip_struct.(model_names{i}).Day = d;
    wrpmip_struct.(model_names{i}).TIMESTAMP = datetime(wrpmip_struct.(model_names{i}).Year,wrpmip_struct.(model_names{i}).Month,wrpmip_struct.(model_names{i}).Day,'Format','dd-MMM-yyyy');    
    wrpmip_struct.(model_names{i}) = table2timetable(wrpmip_struct.(model_names{i}));
    wrpmip_struct.(model_names{i}) = retime(wrpmip_struct.(model_names{i}),time_range); %insert nans for missing data
end

%Create curves for monthly output
wrpmip_monthly = wrpmip_struct;
for i=1:wrpmip_nfiles
    wrpmip_struct.(model_names{i}) = wrpmip_struct.(model_names{i})((wrpmip_struct.(model_names{i}).Month == 5 | wrpmip_struct.(model_names{i}).Month == 6 | wrpmip_struct.(model_names{i}).Month == 7 | wrpmip_struct.(model_names{i}).Month == 8 | wrpmip_struct.(model_names{i}).Month == 9),:);

    %monthly
    wrpmip_struct.(model_names{i})(any(isnan(wrpmip_struct.(model_names{i}).wrpmip_reco), 2), :) = [];
    wrpmip_struct.(model_names{i})(any(isnan(wrpmip_struct.(model_names{i}).wrpmip_temp), 2), :) = [];
    wrpmip_reco_monthly = wrpmip_struct.(model_names{i}).wrpmip_reco;
    wrpmip_temp_monthly = wrpmip_struct.(model_names{i}).wrpmip_temp;
    [wrpmip_f_monthly,wrpmip_gof1_monthly]=fit(wrpmip_temp_monthly,wrpmip_reco_monthly,'exp1');

    wrpmip_ci_monthly = confint(wrpmip_f_monthly);
    wrpmip_a_monthly = wrpmip_f_monthly.a;
    wrpmip_b_monthly = wrpmip_f_monthly.b;

    wrpmip_temp_sort_monthly = sort(wrpmip_temp_monthly);
    wrpmip_y_fit_monthly = wrpmip_a_monthly.*exp(wrpmip_b_monthly.*wrpmip_temp_sort_monthly);
    wrpmip_y_cbl_monthly = wrpmip_ci_monthly(1,1).*exp(wrpmip_ci_monthly(1,2).*wrpmip_temp_sort_monthly);
    wrpmip_y_cbh_monthly = wrpmip_ci_monthly(2,1).*exp(wrpmip_ci_monthly(2,2).*wrpmip_temp_sort_monthly);

    wrpmip_q10_monthly = exp(10.*wrpmip_f_monthly.b);
    wrpmip_cil_monthly = exp(10.*wrpmip_ci_monthly(1,2));
    wrpmip_ciu_monthly = exp(10.*wrpmip_ci_monthly(2,2));

    wrpmip_fs_monthly{i} = wrpmip_f_monthly;
    wrpmip_temp_sorts_monthly{i} = wrpmip_temp_sort_monthly;
    wrpmip_temps_monthly{i} = wrpmip_temp_monthly;
    wrpmip_recos_monthly{i} = wrpmip_reco_monthly;
    wrpmip_cbls_monthly{i} = wrpmip_y_cbl_monthly;
    wrpmip_cbhs_monthly{i} = wrpmip_y_cbh_monthly;
    wrpmip_q10s_monthly{i} = wrpmip_q10_monthly;
    wrpmip_cils_monthly{i} = wrpmip_cil_monthly;
    wrpmip_cius_monthly{i} = wrpmip_ciu_monthly;
end

wrpmip_q10s_95ci=string(round(cell2mat(wrpmip_q10s_monthly),2))+' ['+string(round(cell2mat(wrpmip_cils_monthly),2))+', '+string(round(cell2mat(wrpmip_cius_monthly),2))+']';

min(cell2mat(wrpmip_q10s_monthly))
max(cell2mat(wrpmip_q10s_monthly))

wrpmip_q10s = cell2mat(wrpmip_q10s_monthly);
wrpmip_cils = cell2mat(wrpmip_cils_monthly);
wrpmip_cius = cell2mat(wrpmip_cius_monthly);

%Create supplemental figure of model q10s compared to eml (fig S8)
wrpmip_q10s_data = horzcat(wrpmip_q10s,q10s_daily{10});
wrpmip_cils_data = horzcat(wrpmip_cils,cils_daily{10});
wrpmip_cius_data = horzcat(wrpmip_cius,cius_daily{10});
model_names = {"CLASSIC","CLM5-ExIce","CLM5","ecosys","ELM-ECA","ELM-RD",...
    "JSBACH","JULES","LPJ-GUESS-ML","LPJ-GUESS","ORCHIDEE-MICT-teb","ORCHIDEE-MICT","UVic-ESCM","Data"};
[m,ix] = sort(wrpmip_q10s_data); %sort based on median
models_sort = model_names(ix);
cils_sort=wrpmip_cils_data(ix); cius_sort=wrpmip_cius_data(ix);

hold on
errorbar(1:14,m,m-cils_sort,cius_sort-m,'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(5,m(5),m(5)-cils_sort(5),cius_sort(5)-m(5),'o','color','r','linewidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r',...
        'MarkerSize',10);
xlim([0.25 14.75])
ylim([0.5 6.5])
set(gca,'xtick',1:1:14);
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[models_sort],'fontsize',14)
xtickangle(30)
ylabel('Intercept','fontsize',14);
box on;
ylabel('Q_{10} (unitless)','fontsize',14);


%difference between model q10 at EML and daily EML q10
models_data_diff = wrpmip_q10s - q10s_daily{10};

%extract EML from temporal averaging
time_avg_eml = [q10s_og{10},q10s_daily{10},q10s_weekly{10},q10s_monthly{10},q10s_yearly{10}];

%calculate change in temporal extent
time_ext_eml = q10_yrs_sites{10};
time_ext_eml = cell2mat(vertcat(time_ext_eml{:}));
min_time_ext_eml = min(time_ext_eml);
max_time_ext_eml = max(time_ext_eml);

load("/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Code/site_data_randomization_5cm.mat",...
    "eml_fs","eml_q10s","eml_cils","eml_cius","eml_q10_counts")

eml_q10s2 = cell2mat(eml_q10s);
min_eml_rdm_q10s = min(eml_q10s2,[],1);
max_eml_rdm_q10s = max(eml_q10s2,[],1);

%calculate model scores based on level of agreement in q10 between models
%and observations
for i = 1:wrpmip_nfiles
    models_data_diff_score(i) = exp(-norm(wrpmip_q10s(i) - q10s_daily{10})./norm(wrpmip_q10s(i)));

    wrpmip_time_avg_hh(i) = exp(-norm(wrpmip_q10s(i) - time_avg_eml(1))./norm(wrpmip_q10s(i)));
    wrpmip_time_avg_daily(i) = exp(-norm(wrpmip_q10s(i) - time_avg_eml(2))./norm(wrpmip_q10s(i)));
    wrpmip_time_avg_weekly(i) = exp(-norm(wrpmip_q10s(i) - time_avg_eml(3))./norm(wrpmip_q10s(i)));
    wrpmip_time_avg_monthly(i) = exp(-norm(wrpmip_q10s(i) - time_avg_eml(4))./norm(wrpmip_q10s(i)));
    wrpmip_time_avg_annual(i) = exp(-norm(wrpmip_q10s(i) - time_avg_eml(5))./norm(wrpmip_q10s(i)));

    for j = 1:5
        min_time_ext(j) = exp(-norm(wrpmip_q10s(i) - min_time_ext_eml(j))./norm(wrpmip_q10s(i)));
        max_time_ext(j) = exp(-norm(wrpmip_q10s(i) - max_time_ext_eml(j))./norm(wrpmip_q10s(i)));
    end
    wrpmip_min_time_ext{i} = min_time_ext;
    wrpmip_max_time_ext{i} = max_time_ext;

    for k = 1:30
        min_rdm_q10s(k) = exp(-norm(wrpmip_q10s(i) - min_eml_rdm_q10s(k))./norm(wrpmip_q10s(i)));
        max_rdm_q10s(k) = exp(-norm(wrpmip_q10s(i) - max_eml_rdm_q10s(k))./norm(wrpmip_q10s(i)));
    end
    wrpmip_min_rdm_q10s{i} = min_rdm_q10s;
    wrpmip_max_rdm_q10s{i} = max_rdm_q10s;

    for m = 1:6
        model_tower_scores(m) = exp(-norm(wrpmip_q10s(i) - eml_tower_q10s(m))./norm(wrpmip_q10s(i)));
        model_chamber_scores(m) = exp(-norm(wrpmip_q10s(i) - chamber_q10s(m))./norm(wrpmip_q10s(i)));
        model_gridded_scores(m) = exp(-norm(wrpmip_q10s(i) - gridded_q10s_st(m))./norm(wrpmip_q10s(i)));
    end
    wrpmip_tower_scores{i} = model_tower_scores;
    wrpmip_tower_scores1{i} = model_tower_scores(1);
    wrpmip_chamber_scores{i} = model_chamber_scores;
    wrpmip_gridded_scores{i} = model_gridded_scores;
end

wrpmip_time_avgs = [wrpmip_time_avg_hh;wrpmip_time_avg_daily;wrpmip_time_avg_weekly;wrpmip_time_avg_monthly];

%calculate range in model scores
for i = 1:wrpmip_nfiles
    current_time_avg_max = max(wrpmip_time_avgs(:,i));
    current_time_avg_min = min(wrpmip_time_avgs(:,i));

    current_time_exts = horzcat(wrpmip_min_time_ext{i},wrpmip_max_time_ext{i});
    current_time_ext_min = min(current_time_exts);
    current_time_ext_max = max(current_time_exts);

    current_rdm_q10s = horzcat(wrpmip_min_rdm_q10s{i},wrpmip_max_rdm_q10s{i});
    current_rdm_q10_min = min(current_rdm_q10s);
    current_rdm_q10_max = max(current_rdm_q10s);

    current_dset_choices = horzcat(wrpmip_chamber_scores{i},wrpmip_tower_scores{i},wrpmip_gridded_scores{i});
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
end

%plot range in model skill for each category
model_color=[5,5,5,10,10,10,15,15,15,22,22,22,25];
markers={'diamond','square','^','diamond','square','^','diamond','square','^',...
    'diamond','square','^','diamond'};

avg_line = horzcat(nanmean(cell2mat(ranges_time_avg)),nanmean(cell2mat(ranges_time_ext)),...
    nanmean(cell2mat(ranges_dset_choice)),nanmean(cell2mat(ranges_rdm_q10)));

jitter_amount = 0.6;
hold on
for i = 1:13
    x1 = 1+(rand-0.5)*jitter_amount;
    x2 = 2+(rand-0.5)*jitter_amount;
    x3 = 3+(rand-0.5)*jitter_amount;
    x4 = 4+(rand-0.5)*jitter_amount;

    swarmchart(x1,cell2mat(ranges_time_avg(i)),80,model_color(i),...
        'filled','MarkerFaceAlpha',0.75,'Marker',markers{i},'MarkerEdgeColor','k');
    swarmchart(x2,cell2mat(ranges_time_ext(i)),80,model_color(i),...
        'filled','MarkerFaceAlpha',0.75,'Marker',markers{i},'MarkerEdgeColor','k');
    swarmchart(x3,cell2mat(ranges_dset_choice(i)),80,model_color(i),...
        'filled','MarkerFaceAlpha',0.75,'Marker',markers{i},'MarkerEdgeColor','k');
    swarmchart(x4,cell2mat(ranges_rdm_q10(i)),80,model_color(i),...
        'filled','MarkerFaceAlpha',0.75,'Marker',markers{i},'MarkerEdgeColor','k');
end
plot(avg_line,'k.','MarkerSize',32)
legend("CLASSIC","","","","CLM5-ExIce","","","","CLM5","","","","ecosys","","","",...
    "ELM-ECA","","","","ELM-RD","","","","JSBACH","","","","JULES","","","",...
    "LPJ-GUESS-ML","","","","LPJ-GUESS","","","","ORCHIDEE-MICT-teb","","","","ORCHIDEE-MICT",...
    "","","","UVic-ESCM","",'location','northwest')
grid on; box on;
xticks(1:4)
xticklabels({'Temporal Averaging','Temporal Extent','Choice of Dataset','Number of Observations'})
ylabel("Range in Inferred Model Skill (unitless)")
xtickangle(20)
xlim([0.3 4.7])
ylim([0 0.75])
fontsize(16,"points")

%plot range in model skill from benchmarking choices for each model (fig S9)
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
