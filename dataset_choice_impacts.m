%How does the choice of dataset imapct the functional relationship?

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

%Remove rows where reco or soil temp have nans
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

%Make all sites have the same time range from 2015-2019 to correspond 
%with time period from US-EML
start_time = datetime(2015,01,01);
end_time = datetime(2019,12,31,23,30,00);
time_range = (start_time:minutes(30):end_time)';
for i=1:numel(fence_names)
    cipehr_fence_list.(fence_names{i}) = table2timetable(cipehr_fence_list.(fence_names{i}));
    cipehr_fence_list.(fence_names{i}) = retime(cipehr_fence_list.(fence_names{i}),time_range); %insert nans for missing data
end

%Calculate functional relationships and q10s for control sites
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
  fluxcom_nee_times{i} = ncread(current_file,'time',[5],[5]);
  fluxcom_nees{i} = ncread(current_file,'NEE',[1 1 5],[inf inf 5]); %3D - lon, lat, time
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
  %fluxcom_gpp_times{i} = ncread(current_file,'time',[5],[5]);
  fluxcom_gpps{i} = ncread(current_file,'GPP',[1 1 5],[inf inf 5]); %3D - lon, lat, time
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

%Extract EML location (63.87, -149.25) and summer months
fluxcom_summer = repmat([1:12],1,5)';
fluxcom_eml = squeeze(fluxcom_reco(62,53,:)); %2015-2019
virkkala_eml = squeeze(virkkala_reco(62,53,:)); %2015-2019
merra_fluxcom_eml = squeeze(merra_temp(50,308,181:240)); %subset from 2015-2019
merra_virkkala_eml = squeeze(merra_temp(50,308,181:240)); %subset from 2015-2019
era_fluxcom_eml = squeeze(era_temp(844,106,169:228)); %subset from 2015-2019
era_virkkala_eml = squeeze(era_temp(844,106,169:228)); %subset from 2015-2019

merra_fluxcom_eml((fluxcom_summer < 5),:) = NaN; merra_fluxcom_eml((fluxcom_summer > 9),:) = NaN; 
merra_virkkala_eml((fluxcom_summer < 5),:) = NaN; merra_virkkala_eml((fluxcom_summer > 9),:) = NaN; 
era_fluxcom_eml((fluxcom_summer < 5),:) = NaN; era_fluxcom_eml((fluxcom_summer > 9),:) = NaN; 
era_virkkala_eml((fluxcom_summer < 5),:) = NaN; era_virkkala_eml((fluxcom_summer > 9),:) = NaN; 

merra_fluxcom_eml(isnan(merra_fluxcom_eml)) = [];
merra_virkkala_eml(isnan(merra_virkkala_eml)) = [];
era_fluxcom_eml(isnan(era_fluxcom_eml)) = [];
era_virkkala_eml(isnan(era_virkkala_eml)) = [];

%calculate q10 from each combination of upscaled data product
[f_merra_fluxcom,gof1_merra_fluxcom]=fit(merra_fluxcom_eml,fluxcom_eml,'exp1');
q10_merra_fluxcom = exp(10.*f_merra_fluxcom.b);
ci_merra_fluxcom = confint(f_merra_fluxcom);
cil_merra_fluxcom = exp(10.*ci_merra_fluxcom(1,2));
ciu_merra_fluxcom = exp(10.*ci_merra_fluxcom(2,2));

[f_era_fluxcom,gof1_era_fluxcom]=fit(era_fluxcom_eml,fluxcom_eml,'exp1');
q10_era_fluxcom = exp(10.*f_era_fluxcom.b);
ci_era_fluxcom = confint(f_era_fluxcom);
cil_era_fluxcom = exp(10.*ci_era_fluxcom(1,2));
ciu_era_fluxcom = exp(10.*ci_era_fluxcom(2,2));

[f_merra_virkkala,gof1_merra_virkkala]=fit(merra_virkkala_eml,virkkala_eml,'exp1');
q10_merra_virkkala = exp(10.*f_merra_virkkala.b);
ci_merra_virkkala = confint(f_merra_virkkala);
cil_merra_virkkala = exp(10.*ci_merra_virkkala(1,2));
ciu_merra_virkkala = exp(10.*ci_merra_virkkala(2,2));

[f_era_virkkala,gof1_era_virkkala]=fit(era_virkkala_eml,virkkala_eml,'exp1');
q10_era_virkkala = exp(10.*f_era_virkkala.b);
ci_era_virkkala = confint(f_era_virkkala);
cil_era_virkkala = exp(10.*ci_era_virkkala(1,2));
ciu_era_virkkala = exp(10.*ci_era_virkkala(2,2));


%plot q10 values for each flux tower site, cipehr, and upscaled datasets
plot_letters2={'(a)','(b)'};
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
ax=gca;
ax.Units='normalized';
axPos=ax.Position;
xLimits=xlim;
yArm=0.015;
%bracket for deciduous broaleaf forest
x1=0.5; x2=2;
xNorm1 = axPos(1)+axPos(3)*(x1-xLimits(1))/diff(xLimits);
xNorm2 = axPos(1)+axPos(3)*(x2-xLimits(1))/diff(xLimits);
yBracket1 = axPos(2)-0.07;
xNorm1 = max(min(xNorm1,1),0);
xNorm2 = max(min(xNorm2,1),0);
yBracket1 = max(min(yBracket1,1),0);
yText1 = max(min(yBracket1-0.04,1),0);
annotation('line',[xNorm1 xNorm2],[yBracket1 yBracket1],'Color','k','LineWidth',2);
annotation('line',[xNorm1 xNorm1],[yBracket1 yBracket1 + yArm],'Color','k','LineWidth',2);
annotation('line',[xNorm2 xNorm2],[yBracket1 yBracket1 + yArm],'Color','k','LineWidth',2);
annotation('textbox',[0.08,0, 0.2, 0.03],'String','Deciduous Broadleaf Forest',...
    'EdgeColor','none','HorizontalAlignment','center','Interpreter','none','FontSize',14,'FontName', 'Arial');
%bracket for evergreen needleleaf forest
x3=2.5; x4=10;
xNorm3 = axPos(1)+axPos(3)*(x3-xLimits(1))/diff(xLimits);
xNorm4 = axPos(1)+axPos(3)*(x4-xLimits(1))/diff(xLimits);
yBracket2 = axPos(2)-0.07;
xNorm3 = max(min(xNorm3,1),0);
xNorm4 = max(min(xNorm4,1),0);
yBracket2 = max(min(yBracket2,1),0);
yText2 = max(min(yBracket2-0.04,1),0);
annotation('line',[xNorm3 xNorm4],[yBracket2 yBracket2],'Color','k','LineWidth',2);
annotation('line',[xNorm3 xNorm3],[yBracket2 yBracket2 + yArm],'Color','k','LineWidth',2);
annotation('line',[xNorm4 xNorm4],[yBracket2 yBracket2 + yArm],'Color','k','LineWidth',2);
annotation('textbox',[0.32, yText2, 0.2, 0.03],'String','Evergreen Needleleaf Forest',...
    'EdgeColor','none','HorizontalAlignment','center','Interpreter','none','FontSize',14,'FontName','Arial');
%bracket for shrublands
x5=10.5; x6=13;
xNorm5 = axPos(1)+axPos(3)*(x5-xLimits(1))/diff(xLimits);
xNorm6 = axPos(1)+axPos(3)*(x6-xLimits(1))/diff(xLimits);
yBracket2 = axPos(2)-0.07;
xNorm5 = max(min(xNorm5,1),0);
xNorm6 = max(min(xNorm6,1),0);
yBracket2 = max(min(yBracket2,1),0);
yText2 = max(min(yBracket2-0.04,1),0);
annotation('line',[xNorm5 xNorm6],[yBracket2 yBracket2],'Color','k','LineWidth',2);
annotation('line',[xNorm5 xNorm5],[yBracket2 yBracket2 + yArm],'Color','k','LineWidth',2);
annotation('line',[xNorm6 xNorm6],[yBracket2 yBracket2 + yArm],'Color','k','LineWidth',2);
annotation('textbox', [0.57, yText2, 0.2, 0.03],'String','Open Shrublands','EdgeColor','none',...
    'HorizontalAlignment','center','Interpreter','none','FontSize',14,'FontName','Arial');
%bracket for wetlands
x7=13.5; x8=16;
xNorm7 = axPos(1)+axPos(3)*(x7-xLimits(1))/diff(xLimits);
xNorm8 = axPos(1)+axPos(3)*(x8-xLimits(1))/diff(xLimits);
yBracket2 = axPos(2)-0.07;
xNorm7 = max(min(xNorm7,1),0);
xNorm8 = max(min(xNorm8,1),0);
yBracket2 = max(min(yBracket2,1),0);
yText2 = max(min(yBracket2-0.04,1),0);
annotation('line',[xNorm7 xNorm8],[yBracket2 yBracket2],'Color','k','LineWidth',2);
annotation('line',[xNorm7 xNorm7],[yBracket2 yBracket2 + yArm],'Color','k','LineWidth',2);
annotation('line',[xNorm8 xNorm8],[yBracket2 yBracket2 + yArm],'Color','k','LineWidth',2);
annotation('textbox',[0.71, yText2, 0.2, 0.03],'String','Permanent Wetlands','EdgeColor','none',...
    'HorizontalAlignment','center','Interpreter','none','FontSize',14,'FontName','Arial');
box on;
ylabel('Q_{10} (unitless)','fontsize',14);
text(0.4,4.65, plot_letters2{1}, 'FontWeight', 'bold', 'FontSize', 14)
axes('Position',[0.19 0.63 0.3 0.25])
hold on
errorbar(1:6,cell2mat(ctl_q10s),cell2mat(ctl_q10s)-cell2mat(cil_ctls),cell2mat(ciu_ctls)-cell2mat(ctl_q10s),'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(7,q10s_daily{10},q10s_daily{10}-cils_daily{10},cius_daily{10}-q10s_daily{10},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(8,q10_merra_fluxcom,q10_merra_fluxcom-cil_merra_fluxcom,ciu_merra_fluxcom-q10_merra_fluxcom,'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(9,q10_era_fluxcom,q10_era_fluxcom-cil_era_fluxcom,ciu_era_fluxcom-q10_era_fluxcom,'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(10,q10_merra_virkkala,q10_merra_virkkala-cil_merra_virkkala,ciu_merra_virkkala-q10_merra_virkkala,'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(11,q10_era_virkkala,q10_era_virkkala-cil_era_virkkala,ciu_era_virkkala-q10_era_virkkala,'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
xlim([0.25 11.75])
ylim([1.5 5])
set(gca,'xtick',[3,7,10]);
xline(6.5); xline(7.5);
text(0.5,4.7, plot_letters2{2}, 'FontWeight', 'bold', 'FontSize', 14)
set(gca,'xticklabel',["Chamber","Flux Tower","Upscaled"],'fontsize',14)
title("US-EML")
box on;
ylabel('Q_{10} (unitless)','fontsize',14);

%see range in q10 for each dataset category
eml_tower_q10s = [q10s_daily{10},nan,nan,nan,nan,nan]';
chamber_q10s = cell2mat(ctl_q10s);
gridded_q10s_st = [q10_merra_fluxcom,q10_era_fluxcom,q10_merra_virkkala,q10_era_virkkala,nan,nan]';

tower_median = nanmedian(eml_tower_q10s);
chamber_median = nanmedian(chamber_q10s);
gridded_median = nanmedian(gridded_q10s_st);
min_chamber_q10 = min(chamber_q10s);
max_chamber_q10 = max(chamber_q10s);
min_gridded_q10 = min(gridded_q10s_st);
max_gridded_q10 = max(gridded_q10s_st);
