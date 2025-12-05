%How does the choice of dataset impact the functional relationship between soil temp and reco?

%Load in chamber data from CiPEHR
fluxdb_2015 = readtable("/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/CiPEHR/794_EML_AK_CiPEHR_Compiled_Dataset_HH_2015.csv","TreatAsMissing","NA");
fluxdb_2016 = readtable("/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/CiPEHR/794_EML_AK_CiPEHR_Compiled_Dataset_HH_2016.csv","TreatAsMissing","NA");
fluxdb_2017 = readtable("/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/CiPEHR/794_EML_AK_CiPEHR_Compiled_Dataset_HH_2017.csv","TreatAsMissing","NA");
fluxdb_2018 = readtable("/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/CiPEHR/794_EML_AK_CiPEHR_Compiled_Dataset_HH_2018.csv","TreatAsMissing","NA");
fluxdb_2019 = readtable("/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/CiPEHR/794_EML_AK_CiPEHR_Compiled_Dataset_HH_2019.csv","TreatAsMissing","NA");

fluxdb_2015 = fluxdb_2015((fluxdb_2015.treatment == "Control"), :);
fluxdb_2015 = fluxdb_2015((fluxdb_2015.season == 1), :);
fluxdb_2015 = fluxdb_2015(:,{'flux_year','fence','ts','month','week','doy','hour','hourmin','reco','par','t5','t10'});
fluxdb_2016 = fluxdb_2016((fluxdb_2016.treatment == "Control"), :);
fluxdb_2016 = fluxdb_2016((fluxdb_2016.season == 1), :);
fluxdb_2016 = fluxdb_2016(:,{'flux_year','fence','ts','month','week','doy','hour','hourmin','reco','par','t5','t10'});
fluxdb_2017 = fluxdb_2017((fluxdb_2017.treatment == "Control"), :);
fluxdb_2017 = fluxdb_2017((fluxdb_2017.season == 1), :);
fluxdb_2017 = fluxdb_2017(:,{'flux_year','fence','ts','month','week','doy','hour','hourmin','reco','par','t5','t10'});
fluxdb_2018 = fluxdb_2018((fluxdb_2018.treatment == "Control"), :);
fluxdb_2018 = fluxdb_2018((fluxdb_2018.season == 1), :);
fluxdb_2018 = fluxdb_2018(:,{'flux_year','fence','ts','month','week','doy','hour','hourmin','reco','par','t5','t10'});
fluxdb_2019 = fluxdb_2019((fluxdb_2019.treatment == "Control"), :);
fluxdb_2019 = fluxdb_2019((fluxdb_2019.season == 1), :);
fluxdb_2019 = fluxdb_2019(:,{'flux_year','fence','ts','month','week','doy','hour','hourmin','reco','par','t5','t10'});

fluxdb = vertcat(fluxdb_2015,fluxdb_2016,fluxdb_2017,fluxdb_2018,fluxdb_2019);

%Remove rows where co2 or soil temp have nans
fluxdb(any(isnan(fluxdb.reco), 2), :) = [];
fluxdb(any(isnan(fluxdb.t5), 2), :) = [];

fluxdb=fluxdb((fluxdb.par < 5), :); %look at data when par < 5

%only look at data between hours of 18:00 and 6:00
fluxdb_nt_mask = fluxdb.hour >= 6 & fluxdb.hour < 18;
fluxdb(fluxdb_nt_mask,:) = [];

%average data from fences together, group data by site, and store in a struct
cipehr_fence = unique(fluxdb.fence);
cipehr_fence = cellstr(string(cipehr_fence));
fluxdb.fence = grp2idx(fluxdb.fence);
ala1_plots2 = unique(fluxdb.fence);
fence_names = {"Fence_1","Fence_2","Fence_3","Fence_4","Fence_5","Fence_6"}';
fence_names = cellstr(string(fence_names));

cipehr_fences = cell(1,length(cipehr_fence));
cipehr_fence_list = cell(length(cipehr_fence),1);
for i = 1:length(cipehr_fence)
    current_fence = cipehr_fence{i};
    current_plot = ala1_plots2(i);
    cipehr_fence_avg = fluxdb(fluxdb.fence == current_plot,:);
    cipehr_fence_avg2 = groupsummary(cipehr_fence_avg,'ts','mean');
    cipehr_fences{i} = cipehr_fence_avg2;
end
cipehr_fence_list = cell2struct(cipehr_fences,fence_names,2);

%Make all chamber plots have the same time range from 2015-2019 to correspond with time period from US-EML
start_time = datetime(2015,01,01);
end_time = datetime(2019,12,31,23,30,00);
time_range = (start_time:minutes(30):end_time)';
for i=1:numel(fence_names)
    cipehr_fence_list.(fence_names{i}) = table2timetable(cipehr_fence_list.(fence_names{i}));
    cipehr_fence_list.(fence_names{i}) = retime(cipehr_fence_list.(fence_names{i}),time_range); %insert nans for missing data
end

%Calculate functional relationships and q10s for control plots
cipehr_fence_list_daily = cipehr_fence_list;
fs_exp_ctl_sites = cell(size(fence_names,1),1);
recos_ctl = cell(size(fence_names,1),1);
temps_ctl = cell(size(fence_names,1),1);
ctl_q10s = cell(size(fence_names,1),1);
ciu_ctls = cell(size(fence_names,1),1);
cil_ctls = cell(size(fence_names,1),1);
for i = 1:size(fence_names,1)
    cipehr_fence_list_daily.(fence_names{i}) = retime(cipehr_fence_list.(fence_names{i}), 'daily', 'mean');
    current_ctl_plot = cipehr_fence_list_daily.(fence_names{i});
    if isempty(current_ctl_plot)
        continue
    else
        cipehr_fence_list_daily.(fence_names{i})(any(isnan(cipehr_fence_list_daily.(fence_names{i}).mean_reco), 2), :) = [];
        cipehr_fence_list_daily.(fence_names{i})(any(isnan(cipehr_fence_list_daily.(fence_names{i}).mean_t5), 2), :) = [];

        reco_ctl = cipehr_fence_list_daily.(fence_names{i}).mean_reco;
        temp_ctl = cipehr_fence_list_daily.(fence_names{i}).mean_t5;

        [f_exp_ctl,gof1_exp_ctl]=fit(temp_ctl,reco_ctl,'exp1');
        fs_exp_ctl_sites{i} = f_exp_ctl;

        ctl_q10 = exp(10.*f_exp_ctl.b);
        ctl_q10s{i} = ctl_q10;

        recos_ctl{i} = reco_ctl;
        temps_ctl{i} = temp_ctl;

        ci_ctl = confint(f_exp_ctl);
        cil_ctl = exp(10.*ci_ctl(1,2));
        ciu_ctl = exp(10.*ci_ctl(2,2));
        cil_ctls{i} = cil_ctl;
        ciu_ctls{i} = ciu_ctl;

        ctl_a = f_exp_ctl.a;
        ctl_b = f_exp_ctl.b;
        temp_sort_ctl = sort(temp_ctl);
        fit_ctl = ctl_a.*exp(ctl_b.*temp_sort_ctl);
        cbl_ctl = ci_ctl(1,1).*exp(ci_ctl(1,2).*temp_sort_ctl);
        cbh_ctl = ci_ctl(2,1).*exp(ci_ctl(2,2).*temp_sort_ctl);
        fits_ctl{i} = fit_ctl;
        cbls_ctl{i} = cbl_ctl;
        cbhs_ctl{i} = cbh_ctl;
    end
end

%Load in FLUXCOMX yearly files and combine
fluxcomdir_nee = '/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/FLUXCOM-X/NEE';
fluxcominfo_nee = dir(fullfile(fluxcomdir_nee,'*.nc'));
fluxcom_nee_nfiles = length(fluxcominfo_nee);
fluxcom_nee_filenames = fullfile(fluxcomdir_nee,{fluxcominfo_nee.name});
fluxcom_nees = cell(fluxcom_nee_nfiles,1);
fluxcom_nee_times = cell(fluxcom_nee_nfiles,1);
for i = 1:fluxcom_nee_nfiles
  current_file = fluxcom_nee_filenames{i};
  fluxcom_nee_times{i} = ncread(current_file,'time');
  fluxcom_nees{i} = ncread(current_file,'NEE',[1 1 1],[inf inf inf]); %3D - lon, lat, time
end
fluxcom_time = cat(1,fluxcom_nee_times{:}); %concatenate vars to one matrix
fluxcom_nee = cat(3,fluxcom_nees{:}); %concatenate vars to one matrix

fluxcomdir_gpp = '/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/FLUXCOM-X/GPP';
fluxcominfo_gpp = dir(fullfile(fluxcomdir_gpp,'*.nc'));
fluxcom_gpp_nfiles = length(fluxcominfo_gpp);
fluxcom_gpp_filenames = fullfile(fluxcomdir_gpp,{fluxcominfo_gpp.name});
fluxcom_gpps = cell(fluxcom_gpp_nfiles,1);
fluxcom_gpp_times = cell(fluxcom_gpp_nfiles,1);
for i = 1:fluxcom_gpp_nfiles
  current_file = fluxcom_gpp_filenames{i};
  %fluxcom_gpp_times{i} = ncread(current_file,'time'); %,[5],[5]
  fluxcom_gpps{i} = ncread(current_file,'GPP',[1 1 1],[inf inf inf]); %3D - lon, lat, time [1 1 5],[inf inf 5]
end
fluxcom_gpp = cat(3,fluxcom_gpps{:}); %concatenate vars to one matrix
fluxcom_lat = ncread(current_file,'lat');
fluxcom_lon = ncread(current_file,'lon');

fluxcom_reco = (fluxcom_gpp + fluxcom_nee);

%Load in Virkkala yearly files and combine
virkkaladir = '/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/Virkkala_Reco_2015_2019';
virkkalainfo = dir(fullfile(virkkaladir,'*.nc'));
virkkala_nfiles = length(virkkalainfo);
virkkala_filenames = fullfile(virkkaladir,{virkkalainfo.name});
virkkala_recos = cell(virkkala_nfiles,1);
virkkala_times = cell(virkkala_nfiles,1);
for i = 1:virkkala_nfiles
  current_file = virkkala_filenames{i};
  virkkala_recos{i} = ncread(current_file,'reco'); %3D - lon, lat, time
end
virkkala_reco = cat(3,virkkala_recos{:}); %concatenate vars to one matrix
virkkala_reco = virkkala_reco./365;
virkkala_lat = ncread(current_file,'lat');
virkkala_lon = ncread(current_file,'lon');

%Load in MERRA-2
merradir = '/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/raw/MERRA-2';
merrainfo = dir(fullfile(merradir,'*.nc4'));
merra_nfiles = length(merrainfo);
merra_filenames = fullfile(merradir,{merrainfo.name});
merra_filenames = natsortfiles(merra_filenames);
merra_temps = cell(merra_nfiles,1);
for i = 1:merra_nfiles
  current_file = merra_filenames{i};
  merra_temps{i} = ncread(current_file,'TSOIL1'); %layer 1 is 10cm
end
merra_lat = ncread(current_file,'lat');
merra_lon = ncread(current_file,'lon');
merra_temp = cat(3,merra_temps{:}); %concatenate vars to one matrix. units in k
merra_temp = merra_temp - 273.15;

%Load in ERA-5
eradir = '/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Data/raw/ERA-5';
erainfo = dir(fullfile(eradir,'*.nc'));
era_nfiles = length(erainfo);
era_filenames = fullfile(eradir,{erainfo.name});
era_filenames = natsortfiles(era_filenames);
era_temps = cell(era_nfiles,1);
for i = 1:era_nfiles
  current_file = era_filenames{i};
  era_temps{i} = ncread(current_file,'STL1'); %layer 1 is 10cm
end
era_lat = ncread(current_file,'latitude');
era_lon = ncread(current_file,'longitude');
era_temp = cat(3,era_temps{:}); %concatenate vars to one matrix. units in k
era_temp = era_temp - 273.15;
era_lon_mask = era_lon >= 180;
era_lon(era_lon_mask) = era_lon(era_lon_mask) - 360;

%Read in flux tower coordinates
site_info = readtable('/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Flux Towers/flux_tower_lat_lon.csv','TreatAsMissing','NA');
site_abbrev = site_info.Site;
ec_lats = site_info.Latitude;
ec_lons = site_info.Longitude;
site_coords = horzcat(ec_lats,ec_lons);

start_time = datetime(2001,01,01);
end_time = datetime(2020,12,01);
time_range = (start_time:calmonths(1):end_time)';

site_list_monthly2 = site_list_monthly;
for i=1:numel(site_names)
    site_list_monthly2.(site_names{i}) = retime(site_list_monthly2.(site_names{i}),time_range); %insert nan for missing data
end

%Extract flux tower locations from upscaled datasets
for i = 1:size(site_abbrev,1)
    current_site_years = site_list_monthly2.(site_names{i}).RECO;
    lat_idx = fluxcom_lat - ec_lats(i);
    lon_idx = fluxcom_lon - ec_lons(i);
    [lat_val,lat_ix]=min(abs(lat_idx)); %find the location closest to zero
    [lon_val,lon_ix]=min(abs(lon_idx)); %find the location closest to zero
    fluxcom_sites = fluxcom_reco(lon_ix,lat_ix,:);
    fluxcom_sites = reshape(fluxcom_sites,240,1);
    fluxcom_sites(isnan(current_site_years)) = NaN;
    fluxcom_model_sites{i} = fluxcom_sites;
end

%concatenate virkkala et al. growing season outputs with nans for non-growing seaon
nan4=nan(720,113,4);
nan7=nan(720,113,7);
nan3=nan(720,113,3);
virkkala_reco_20yr = cat(3,nan4,virkkala_reco(:,:,1:5),nan7,virkkala_reco(:,:,6:10),nan7,virkkala_reco(:,:,11:15),...
    nan7,virkkala_reco(:,:,16:20),nan7,virkkala_reco(:,:,21:25),nan7,virkkala_reco(:,:,26:30),nan7,virkkala_reco(:,:,31:35),...
    nan7,virkkala_reco(:,:,36:40),nan7,virkkala_reco(:,:,41:45),nan7,virkkala_reco(:,:,46:50),nan7,virkkala_reco(:,:,51:55),...
    nan7,virkkala_reco(:,:,56:60),nan7,virkkala_reco(:,:,61:65),nan7,virkkala_reco(:,:,66:70),nan7,virkkala_reco(:,:,71:75),...
    nan7,virkkala_reco(:,:,76:80),nan7,virkkala_reco(:,:,81:85),nan7,virkkala_reco(:,:,86:90),nan7,...
    virkkala_reco(:,:,91:95),nan7,virkkala_reco(:,:,96:100),nan3); 

%extract flux tower locations from virkkala
for i = 1:size(site_abbrev,1)
    current_site_years = site_list_monthly2.(site_names{i}).RECO;
    lat_idx = virkkala_lat - ec_lats(i);
    lon_idx = virkkala_lon - ec_lons(i);
    [lat_val,lat_ix]=min(abs(lat_idx));
    [lon_val,lon_ix]=min(abs(lon_idx));
    virkkala_sites = virkkala_reco_20yr(lon_ix,lat_ix,:);
    virkkala_sites = reshape(virkkala_sites,240,1);
    virkkala_sites(isnan(current_site_years)) = NaN;
    virkkala_model_sites{i} = virkkala_sites;
end

%extract flux tower locations from merra
for i = 1:size(site_abbrev,1)
    current_site_years = site_list_monthly2.(site_names{i}).RECO;
    lat_idx = merra_lat - ec_lats(i);
    lon_idx = merra_lon - ec_lons(i);
    [lat_val,lat_ix]=min(abs(lat_idx));
    [lon_val,lon_ix]=min(abs(lon_idx));
    merra_sites = merra_temp(lon_ix,lat_ix,:);
    merra_sites = reshape(merra_sites,240,1);
    merra_sites(isnan(current_site_years)) = NaN;
    merra_model_sites{i} = merra_sites;
end

%extract flux tower locations from era
for i = 1:size(site_abbrev,1)
    current_site_years = site_list_monthly2.(site_names{i}).RECO;
    lat_idx = era_lat - ec_lats(i);
    lon_idx = era_lon - ec_lons(i);
    [lat_val,lat_ix]=min(abs(lat_idx));
    [lon_val,lon_ix]=min(abs(lon_idx));
    era_sites = era_temp(lon_ix,lat_ix,:);
    era_sites = reshape(era_sites,240,1);
    era_sites(isnan(current_site_years)) = NaN;
    era_model_sites{i} = era_sites;
end

%calculate q10s using fluxcom and merra
for i = 1:size(site_abbrev,1)
    current_fluxcom = fluxcom_model_sites{i};
    current_merra = merra_model_sites{i};

    current_merra(isnan(current_fluxcom)) = [];
    current_fluxcom(isnan(current_fluxcom)) = [];

    [current_f,current_gof1]=fit(current_merra,current_fluxcom,'exp1');
    current_q10 = exp(10.*current_f.b);
    current_ci = confint(current_f);
    current_cil = exp(10.*current_ci(1,2));
    current_ciu = exp(10.*current_ci(2,2));

    q10_merra_fluxcom{i} = current_q10;
    ci_merra_fluxcom{i} = current_ci;
    cil_merra_fluxcom{i} = current_cil;
    ciu_merra_fluxcom{i} = current_ciu;
end

%calculate q10s using fluxcom and era
for i = 1:size(site_abbrev,1)
    current_fluxcom = fluxcom_model_sites{i};
    current_era = era_model_sites{i};

    current_era(isnan(current_fluxcom)) = [];
    current_fluxcom(isnan(current_fluxcom)) = [];

    [current_f,current_gof1]=fit(current_era,current_fluxcom,'exp1');
    current_q10 = exp(10.*current_f.b);
    current_ci = confint(current_f);
    current_cil = exp(10.*current_ci(1,2));
    current_ciu = exp(10.*current_ci(2,2));

    q10_era_fluxcom{i} = current_q10;
    ci_era_fluxcom{i} = current_ci;
    cil_era_fluxcom{i} = current_cil;
    ciu_era_fluxcom{i} = current_ciu;
end

%calculate q10s using virkkala and merra
for i = 1:size(site_abbrev,1)
    current_virkkala = virkkala_model_sites{i};
    current_merra = merra_model_sites{i};

    current_merra(isnan(current_virkkala)) = [];
    current_virkkala(isnan(current_virkkala)) = [];

    [current_f,current_gof1]=fit(current_merra,current_virkkala,'exp1');
    current_q10 = exp(10.*current_f.b); %2.5659
    current_ci = confint(current_f);
    current_cil = exp(10.*current_ci(1,2));
    current_ciu = exp(10.*current_ci(2,2));

    q10_merra_virkkala{i} = current_q10;
    ci_merra_virkkala{i} = current_ci;
    cil_merra_virkkala{i} = current_cil;
    ciu_merra_virkkala{i} = current_ciu;
end

%calculate q10s using fluxcom and era
for i = 1:size(site_abbrev,1)
    current_virkkala = virkkala_model_sites{i};
    current_era = era_model_sites{i};

    current_era(isnan(current_virkkala)) = [];
    current_virkkala(isnan(current_virkkala)) = [];

    [current_f,current_gof1]=fit(current_era,current_virkkala,'exp1');
    current_q10 = exp(10.*current_f.b); %2.5659
    current_ci = confint(current_f);
    current_cil = exp(10.*current_ci(1,2));
    current_ciu = exp(10.*current_ci(2,2));

    q10_era_virkkala{i} = current_q10;
    ci_era_virkkala{i} = current_ci;
    cil_era_virkkala{i} = current_cil;
    ciu_era_virkkala{i} = current_ciu;
end

plot_letters2={'(a)','(b)'};

%plot dotplot
figure(1)
hold on
errorbar(1:16,cell2mat(m2),cell2mat(m2)-cils_sort2,cius_sort2-cell2mat(m2),'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
xlim([0.25 16.75])
ylim([1.25 4.75])
set(gca,'xtick',1:1:16);
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[sites_sort2],'fontsize',14)
xtickangle(25)
ax = gca;
ax.Units = 'normalized';
axPos = ax.Position;
xLimits = xlim;
yArm = 0.015;
%deciduous broaleaf forest
x1 = 0.5; x2 = 2;
xNorm1 = axPos(1) + axPos(3) * (x1 - xLimits(1)) / diff(xLimits);
xNorm2 = axPos(1) + axPos(3) * (x2 - xLimits(1)) / diff(xLimits);
yBracket1 = axPos(2) - 0.07;
% Clamp positions to [0, 1]
xNorm1 = max(min(xNorm1, 1), 0);
xNorm2 = max(min(xNorm2, 1), 0);
yBracket1 = max(min(yBracket1, 1), 0);
yText1 = max(min(yBracket1 - 0.04, 1), 0);
annotation('line', [xNorm1 xNorm2], [yBracket1 yBracket1], 'Color', 'k', 'LineWidth', 2);
annotation('line', [xNorm1 xNorm1], [yBracket1 yBracket1 + yArm], 'Color', 'k', 'LineWidth', 2);
annotation('line', [xNorm2 xNorm2], [yBracket1 yBracket1 + yArm], 'Color', 'k', 'LineWidth', 2);
annotation('textbox', [0.08,0, 0.2, 0.03], ...
    'String', 'Deciduous Broadleaf Forest', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'Interpreter', 'none', 'FontSize', 14, 'FontName', 'Arial');
%evergreen needleleaf forest
x3 = 2.5; x4 = 10;
xNorm3 = axPos(1) + axPos(3) * (x3 - xLimits(1)) / diff(xLimits);
xNorm4 = axPos(1) + axPos(3) * (x4 - xLimits(1)) / diff(xLimits);
yBracket2 = axPos(2) - 0.07;
xNorm3 = max(min(xNorm3, 1), 0);
xNorm4 = max(min(xNorm4, 1), 0);
yBracket2 = max(min(yBracket2, 1), 0);
yText2 = max(min(yBracket2 - 0.04, 1), 0);
annotation('line', [xNorm3 xNorm4], [yBracket2 yBracket2], 'Color', 'k', 'LineWidth', 2);
annotation('line', [xNorm3 xNorm3], [yBracket2 yBracket2 + yArm], 'Color', 'k', 'LineWidth', 2);
annotation('line', [xNorm4 xNorm4], [yBracket2 yBracket2 + yArm], 'Color', 'k', 'LineWidth', 2);
annotation('textbox', [0.32, yText2, 0.2, 0.03], ...
    'String', 'Evergreen Needleleaf Forest', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'Interpreter', 'none', 'FontSize', 14, 'FontName', 'Arial');
%shrublands
x5 = 10.5; x6 = 13;
xNorm5 = axPos(1) + axPos(3) * (x5 - xLimits(1)) / diff(xLimits);
xNorm6 = axPos(1) + axPos(3) * (x6 - xLimits(1)) / diff(xLimits);
yBracket2 = axPos(2) - 0.07;
xNorm5 = max(min(xNorm5, 1), 0);
xNorm6 = max(min(xNorm6, 1), 0);
yBracket2 = max(min(yBracket2, 1), 0);
yText2 = max(min(yBracket2 - 0.04, 1), 0);
annotation('line', [xNorm5 xNorm6], [yBracket2 yBracket2], 'Color', 'k', 'LineWidth', 2);
annotation('line', [xNorm5 xNorm5], [yBracket2 yBracket2 + yArm], 'Color', 'k', 'LineWidth', 2);
annotation('line', [xNorm6 xNorm6], [yBracket2 yBracket2 + yArm], 'Color', 'k', 'LineWidth', 2);
annotation('textbox', [0.57, yText2, 0.2, 0.03],'String','Open Shrublands','EdgeColor','none',...
    'HorizontalAlignment','center','Interpreter','none','FontSize',14,'FontName','Arial');
%wetlands
x7 = 13.5; x8 = 16;
xNorm7 = axPos(1) + axPos(3) * (x7 - xLimits(1)) / diff(xLimits);
xNorm8 = axPos(1) + axPos(3) * (x8 - xLimits(1)) / diff(xLimits);
yBracket2 = axPos(2) - 0.07;
xNorm7 = max(min(xNorm7, 1), 0);
xNorm8 = max(min(xNorm8, 1), 0);
yBracket2 = max(min(yBracket2, 1), 0);
yText2 = max(min(yBracket2 - 0.04, 1), 0);
annotation('line', [xNorm7 xNorm8], [yBracket2 yBracket2], 'Color', 'k', 'LineWidth', 2);
annotation('line', [xNorm7 xNorm7], [yBracket2 yBracket2 + yArm], 'Color', 'k', 'LineWidth', 2);
annotation('line', [xNorm8 xNorm8], [yBracket2 yBracket2 + yArm], 'Color', 'k', 'LineWidth', 2);
annotation('textbox', [0.71, yText2, 0.2, 0.03],'String','Permanent Wetlands','EdgeColor','none',...
    'HorizontalAlignment','center','Interpreter','none','FontSize',14,'FontName','Arial');
box on;
ylabel('Q_{10} (unitless)','fontsize',14);
text(0.4,4.65, plot_letters2{1}, 'FontWeight', 'bold', 'FontSize', 14)
axes('Position',[0.19 0.63 0.3 0.25])
box on
hold on
errorbar(1:6,cell2mat(ctl_q10s),cell2mat(ctl_q10s)-cell2mat(cil_ctls),cell2mat(ciu_ctls)-cell2mat(ctl_q10s),'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(7,q10s_daily{10},q10s_daily{10}-cils_daily{10},cius_daily{10}-q10s_daily{10},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(8,q10_merra_fluxcom{10},q10_merra_fluxcom{10}-cil_merra_fluxcom{10},ciu_merra_fluxcom{10}-q10_merra_fluxcom{10},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(9,q10_era_fluxcom{10},q10_era_fluxcom{10}-cil_era_fluxcom{10},ciu_era_fluxcom{10}-q10_era_fluxcom{10},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(10,q10_merra_virkkala{10},q10_merra_virkkala{10}-cil_merra_virkkala{10},ciu_merra_virkkala{10}-q10_merra_virkkala{10},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(11,q10_era_virkkala{10},q10_era_virkkala{10}-cil_era_virkkala{10},ciu_era_virkkala{10}-q10_era_virkkala{10},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
xlim([0.25 11.75])
ylim([1.5 5])
set(gca,'xtick',[3,7,10]);
xline(6.5); xline(7.5);
text(0.5,4.7, plot_letters2{2}, 'FontWeight', 'bold', 'FontSize', 14)
set(gca,'xticklabel',["Chamber","Flux Tower","Upscaled"],'fontsize',14)
%xtickangle(45)
title("US-EML")
box on;
ylabel('Q_{10} (unitless)','fontsize',14);

%Plot Q10 based on choice of dataset for all sites
tiledlayout(4,4,"TileSpacing","compact")
for i = 1:size(site_abbrev,1)
    nexttile
    hold on
    if i == 10
        errorbar(1:6,cell2mat(ctl_q10s),cell2mat(ctl_q10s)-cell2mat(cil_ctls),cell2mat(ciu_ctls)-cell2mat(ctl_q10s),'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(7,q10s_daily{i},q10s_daily{i}-cils_daily{i},cius_daily{i}-q10s_daily{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(8,q10_merra_fluxcom{i},q10_merra_fluxcom{i}-cil_merra_fluxcom{i},ciu_merra_fluxcom{i}-q10_merra_fluxcom{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(9,q10_era_fluxcom{i},q10_era_fluxcom{i}-cil_era_fluxcom{i},ciu_era_fluxcom{i}-q10_era_fluxcom{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(10,q10_merra_virkkala{i},q10_merra_virkkala{i}-cil_merra_virkkala{i},ciu_merra_virkkala{i}-q10_merra_virkkala{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(11,q10_era_virkkala{i},q10_era_virkkala{i}-cil_era_virkkala{i},ciu_era_virkkala{i}-q10_era_virkkala{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        xlim([0.25 11.75])
        set(gca,'xtick',[3,7,10.1]);
        xline(6.5); xline(7.5);
        set(gca,'xticklabel',["Chamber","Flux Tower","Upscaled"],'fontsize',14)
    else
        errorbar(1,q10s_daily{i},q10s_daily{i}-cils_daily{i},cius_daily{i}-q10s_daily{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(2,q10_merra_fluxcom{i},q10_merra_fluxcom{i}-cil_merra_fluxcom{i},ciu_merra_fluxcom{i}-q10_merra_fluxcom{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(3,q10_era_fluxcom{i},q10_era_fluxcom{i}-cil_era_fluxcom{i},ciu_era_fluxcom{i}-q10_era_fluxcom{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(4,q10_merra_virkkala{i},q10_merra_virkkala{i}-cil_merra_virkkala{i},ciu_merra_virkkala{i}-q10_merra_virkkala{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        errorbar(5,q10_era_virkkala{i},q10_era_virkkala{i}-cil_era_virkkala{i},ciu_era_virkkala{i}-q10_era_virkkala{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerSize',9);
        xlim([0.25 5.75])
        xline(1.5);
        set(gca,'xtick',[1,3.5]);
        set(gca,'xticklabel',["Flux Tower","Upscaled"],'fontsize',14)
    end
    title(site_abbrev{i})
    box on;
    ylabel('Q_{10} (unitless)','fontsize',14);
    fontsize(12,"points")
end

%see range in q10 for each dataset category
eml_tower_q10s = [q10s_daily{10},nan,nan,nan,nan,nan]';
chamber_q10s = cell2mat(ctl_q10s);
gridded_q10s_st = [q10_merra_fluxcom{10},q10_era_fluxcom{10},q10_merra_virkkala{10},q10_era_virkkala{10},nan,nan]';

for i = 1:16
    site_tower_q10s = [q10s_daily{i},nan,nan,nan,nan,nan]';
    site_gridded_q10s = [q10_merra_fluxcom{i},q10_era_fluxcom{i},q10_merra_virkkala{i},q10_era_virkkala{i},nan,nan]';
    sites_tower_q10s{i} = site_tower_q10s;
    sites_gridded_q10s{i} = site_gridded_q10s;
end

tower_median = nanmedian(eml_tower_q10s); %2.5026
chamber_median = nanmedian(chamber_q10s); %3.6312
gridded_median = nanmedian(gridded_q10s_st); %2.7066
min_chamber_q10 = min(chamber_q10s); %2.8777
max_chamber_q10 = max(chamber_q10s); %4.5835
min_gridded_q10 = min(gridded_q10s_st); %2.3410
max_gridded_q10 = max(gridded_q10s_st); %3.0923
