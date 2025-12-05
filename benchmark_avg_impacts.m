%Script to evaluate how temporal averaging, temporal extent, and number of
%daily observations impact the functional relationship from multiple flux
%tower sites

%Load in flux tower data and store in a struct
pth = '/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Flux Towers/Ameriflux/Current_Sites/';
file_list = dir(strcat(pth,'*.csv'));
files = {file_list.name}';
filenames = fullfile(pth,{file_list.name})';
site_list = cell(size(filenames,1),1);
for i = 1:size(files,1)
  %get site abbreviations from file
  underline_location1 = strfind(files, '_');
  underline_location2 = underline_location1{i};
  underline_location3(i) = underline_location2(2);
  current_file = files{i};
  first_part = current_file(8:underline_location3(i)-1);
  first_part_full = current_file(5:underline_location3(i)-1);
  site_names(i,1) = cellstr(first_part);
  site_names_full(i,1) = cellstr(first_part_full);

  %put data into struct
  current_filename = filenames{i};
  sites{i} = readtable(current_filename,'TreatAsMissing','NA');
  site_list{i} = struct(site_names{i},sites{i});  
end
site_list = cell2struct(sites,site_names,2);

%Put all sites in the same time range from 2001-2020 and convert to timetable. 
%Insert -9999 for missing data
start_time = datetime(2001,01,01);
end_time = datetime(2020,12,31,23,30,00);
time_range = (start_time:minutes(30):end_time)';
for i=1:numel(site_names)
    site_list.(site_names{i}).TIMESTAMP_START = num2str(site_list.(site_names{i}).TIMESTAMP_START);
    site_datevec = datevec(site_list.(site_names{i}).TIMESTAMP_START,'yyyymmddHHMM');
    site_list.(site_names{i}).Year = site_datevec(:,1);
    site_list.(site_names{i}).Month = site_datevec(:,2);
    site_list.(site_names{i}).Day = site_datevec(:,3);
    site_list.(site_names{i}).Hour = site_datevec(:,4);
    site_list.(site_names{i}).Minute = site_datevec(:,5);
    site_list.(site_names{i}).Second = site_datevec(:,6);
    site_list.(site_names{i}).TIMESTAMP = datetime(site_list.(site_names{i}).Year,site_list.(site_names{i}).Month,site_list.(site_names{i}).Day,site_list.(site_names{i}).Hour,site_list.(site_names{i}).Minute,site_list.(site_names{i}).Second,'Format','dd-MMM-yyyy hh:mm');
    site_list.(site_names{i}).TIMESTAMP_START = [];
    site_list.(site_names{i}) = table2timetable(site_list.(site_names{i}));
    site_list.(site_names{i}) = retime(site_list.(site_names{i}),time_range); %insert -9999 for missing data
end

%Remove data from year 2008 from Rpf because it's only a few days from Aug.
site_list.Rpf((site_list.Rpf.Year == 2008),:) = [];

%Only look at growing season, insert NaNs for missing values, and rename to common column names
for i=1:numel(site_names)
    site_list.(site_names{i}) = site_list.(site_names{i})((site_list.(site_names{i}).Month == 5 | site_list.(site_names{i}).Month == 6 | site_list.(site_names{i}).Month == 7 | site_list.(site_names{i}).Month == 8 | site_list.(site_names{i}).Month == 9),:);

    %run commented lines only when calculating q10 at various depths for overlapping time periods
    if any(strcmp('TA', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TA((site_list.(site_names{i}).TA == -9999),:) = NaN;
        %site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TA),2), :) = []; %
        site_list.(site_names{i}).TA = site_list.(site_names{i}).TA;
    else any(strcmp('TA_1_1_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TA_1_1_1((site_list.(site_names{i}).TA_1_1_1 == -9999),:) = NaN;
        %site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TA_1_1_1),2), :) = []; %
        site_list.(site_names{i}).TA = site_list.(site_names{i}).TA_1_1_1;
    end

    %soil temp at 2cm depth
    %run commented lines only when calculating q10 at various depths for overlapping time periods
    if any(strcmp('TS_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_1((site_list.(site_names{i}).TS_1 == -9999),:) = NaN;
        %site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_1),2), :) = []; %
        site_list.(site_names{i}).TS_2cm = site_list.(site_names{i}).TS_1;
    elseif any(strcmp('TS_1_1_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_1_1_1((site_list.(site_names{i}).TS_1_1_1 == -9999),:) = NaN;
        %site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_1_1_1), 2), :) = []; %
        site_list.(site_names{i}).TS_2cm = site_list.(site_names{i}).TS_1_1_1;
    elseif any(strcmp('TS_PI_F_1_1_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_PI_F_1_1_1((site_list.(site_names{i}).TS_PI_F_1_1_1 == -9999),:) = NaN;
        %site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_PI_F_1_1_1), 2), :) = []; %
        site_list.(site_names{i}).TS_2cm = site_list.(site_names{i}).TS_PI_F_1_1_1;
    else any(strcmp('TS_PI_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_PI_1((site_list.(site_names{i}).TS_PI_1 == -9999),:) = NaN;
        %site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_PI_1), 2), :) = []; %
        site_list.(site_names{i}).TS_2cm = site_list.(site_names{i}).TS_PI_1;
    end

    %soil temp at 5cm depth
    if any(strcmp('TS_2', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_2((site_list.(site_names{i}).TS_2 == -9999),:) = NaN;
        site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_2), 2), :) = [];
        site_list.(site_names{i}).TS_5cm = site_list.(site_names{i}).TS_2;
    elseif any(strcmp('TS_1_2_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_1_2_1((site_list.(site_names{i}).TS_1_2_1 == -9999),:) = NaN;
        site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_1_2_1), 2), :) = [];
        site_list.(site_names{i}).TS_5cm = site_list.(site_names{i}).TS_1_2_1;
    elseif any(strcmp('TS_PI_F_1_2_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_PI_F_1_2_1((site_list.(site_names{i}).TS_PI_F_1_2_1 == -9999),:) = NaN;
        site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_PI_F_1_2_1), 2), :) = [];
        site_list.(site_names{i}).TS_5cm = site_list.(site_names{i}).TS_PI_F_1_2_1;
    elseif any(strcmp('TS_2_1_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_2_1_1((site_list.(site_names{i}).TS_2_1_1 == -9999),:) = NaN;
        site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_2_1_1), 2), :) = [];
        site_list.(site_names{i}).TS_5cm = site_list.(site_names{i}).TS_2_1_1;
    else any(strcmp('TS_PI_2', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_PI_2((site_list.(site_names{i}).TS_PI_2 == -9999),:) = NaN;
        site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).TS_PI_2), 2), :) = [];
        site_list.(site_names{i}).TS_5cm = site_list.(site_names{i}).TS_PI_2;
    end

    %soil temp at 10cm depth
    if any(strcmp('TS_PI_F_1_3_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_PI_F_1_3_1((site_list.(site_names{i}).TS_PI_F_1_3_1 == -9999),:) = NaN;
        site_list.(site_names{i}).TS_10cm = site_list.(site_names{i}).TS_PI_F_1_3_1;
    elseif any(strcmp('TS_1_3_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).TS_1_3_1((site_list.(site_names{i}).TS_1_3_1 == -9999),:) = NaN;
        site_list.(site_names{i}).TS_10cm = site_list.(site_names{i}).TS_1_3_1;
    else 
        site_list.(site_names{i}).TS_10cm = nan(size(site_list.(site_names{i}),1),1);
    end

    if any(strcmp('RECO_PI', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).RECO_PI((site_list.(site_names{i}).RECO_PI == -9999),:) = NaN;
        site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).RECO_PI), 2), :) = [];
        site_list.(site_names{i}).RECO = site_list.(site_names{i}).RECO_PI;
    elseif any(strcmp('RECO_PI_F', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).RECO_PI_F((site_list.(site_names{i}).RECO_PI_F == -9999),:) = NaN;
        site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).RECO_PI_F), 2), :) = [];
        site_list.(site_names{i}).RECO = site_list.(site_names{i}).RECO_PI_F;
    elseif any(strcmp('RECO_PI_1_1_1', site_list.(site_names{i}).Properties.VariableNames))
        site_list.(site_names{i}).RECO_PI_1_1_1((site_list.(site_names{i}).RECO_PI_1_1_1 == -9999),:) = NaN;
        site_list.(site_names{i})(any(isnan(site_list.(site_names{i}).RECO_PI_1_1_1), 2), :) = [];
        site_list.(site_names{i}).RECO = site_list.(site_names{i}).RECO_PI_1_1_1;
    else
        site_list.(site_names{i}).RECO = nan(size(site_list.(site_names{i}),1),1);
    end
    site_list.(site_names{i}).RECO = site_list.(site_names{i}).RECO.*((1/1000000).*12.0107.*86400); %micromol CO2/m2/s --> gC/m2/day
end

%Only look at nighttime and when PPFD < 10, insert NaNs for missing values, and rename to common column names
site_list_nt = site_list;
for i=1:numel(site_names)
    if any(strcmp('PPFD_IN', site_list_nt.(site_names{i}).Properties.VariableNames))
        site_list_nt.(site_names{i}).PPFD_IN((site_list_nt.(site_names{i}).PPFD_IN == -9999),:) = NaN;
        site_list_nt.(site_names{i}).PPFD_IN((site_list_nt.(site_names{i}).PPFD_IN >= 10),:) = NaN;
        site_list_nt.(site_names{i})(any(isnan(site_list_nt.(site_names{i}).PPFD_IN), 2), :) = [];
        site_list_nt.(site_names{i}).PAR = site_list_nt.(site_names{i}).PPFD_IN;
    elseif any(strcmp('PPFD_IN_1_1_1', site_list_nt.(site_names{i}).Properties.VariableNames))
        site_list_nt.(site_names{i}).PPFD_IN_1_1_1((site_list_nt.(site_names{i}).PPFD_IN_1_1_1 == -9999),:) = NaN;
        site_list_nt.(site_names{i}).PPFD_IN_1_1_1((site_list_nt.(site_names{i}).PPFD_IN_1_1_1 >= 10),:) = NaN;
        site_list_nt.(site_names{i})(any(isnan(site_list_nt.(site_names{i}).PPFD_IN_1_1_1), 2), :) = [];
        site_list_nt.(site_names{i}).PAR = site_list_nt.(site_names{i}).PPFD_IN_1_1_1;
    elseif any(strcmp('PPFD_IN_PI_F', site_list_nt.(site_names{i}).Properties.VariableNames))
        site_list_nt.(site_names{i}).PPFD_IN_PI_F((site_list_nt.(site_names{i}).PPFD_IN_PI_F == -9999),:) = NaN;
        site_list_nt.(site_names{i}).PPFD_IN_PI_F((site_list_nt.(site_names{i}).PPFD_IN_PI_F >= 10),:) = NaN;
        site_list_nt.(site_names{i})(any(isnan(site_list_nt.(site_names{i}).PPFD_IN_PI_F), 2), :) = [];
        site_list_nt.(site_names{i}).PAR = site_list_nt.(site_names{i}).PPFD_IN_PI_F;
    else 
        site_list_nt.(site_names{i}).PAR = nan(size(site_list_nt.(site_names{i}),1),1);
    end

    nt_mask = site_list_nt.(site_names{i}).Hour >= 6 & site_list_nt.(site_names{i}).Hour < 18; %only look at nighttime measurements between 18:00-6:00
    site_list_nt.(site_names{i})(nt_mask,:) = [];
end

%Calculate means over time, construct functional relationships, and calculate q10 values for each site
site_list_daily = site_list_nt;
site_list_weekly = site_list_nt;
site_list_monthly = site_list_nt;
site_list_yearly = site_list_nt;
site_list_ta = site_list_nt;
site_list_2cm = site_list_nt;
site_list_10cm = site_list_nt;
for i=1:numel(site_names)
    site_list_daily.(site_names{i}) = retime(site_list_daily.(site_names{i}), 'daily', 'mean'); 
    site_list_weekly.(site_names{i}) = retime(site_list_weekly.(site_names{i}), 'weekly', 'mean'); 
    site_list_monthly.(site_names{i}) = retime(site_list_monthly.(site_names{i}), 'monthly', 'mean'); 
    site_list_yearly.(site_names{i}) = retime(site_list_yearly.(site_names{i}), 'yearly', 'mean'); 

    site_list_daily.(site_names{i})(any(isnan(site_list_daily.(site_names{i}).RECO), 2), :) = [];
    %site_list_daily.(site_names{i})(any(isnan(site_list_daily.(site_names{i}).TS_2cm), 2), :) = [];
    site_list_daily.(site_names{i})(any(isnan(site_list_daily.(site_names{i}).TS_5cm), 2), :) = [];
    %site_list_daily.(site_names{i})(any(isnan(site_list_daily.(site_names{i}).TA), 2), :) = [];

    %create daily relationship using air temp
    site_list_ta_daily.(site_names{i}) = retime(site_list_ta.(site_names{i}), 'daily', 'mean'); 
    site_list_ta_daily.(site_names{i})(any(isnan(site_list_ta_daily.(site_names{i}).RECO), 2), :) = [];
    site_list_ta_daily.(site_names{i})(any(isnan(site_list_ta_daily.(site_names{i}).TA), 2), :) = [];
    reco_ta = site_list_ta_daily.(site_names{i}).RECO;
    temp_ta = site_list_ta_daily.(site_names{i}).TA;
    [f_ta,gof1_ta]=fit(temp_ta,reco_ta,'exp1');

    ci_ta = confint(f_ta);
    a_ta = f_ta.a;
    b_ta = f_ta.b;
    temp_sort_ta = sort(temp_ta);
    y_fit_ta = a_ta.*exp(b_ta.*temp_sort_ta);
    y_cbl_ta = ci_ta(1,1).*exp(ci_ta(1,2).*temp_sort_ta);
    y_cbh_ta = ci_ta(2,1).*exp(ci_ta(2,2).*temp_sort_ta);
    q10_ta = exp(10.*f_ta.b);
    cil_ta = exp(10.*ci_ta(1,2));
    ciu_ta = exp(10.*ci_ta(2,2));

    fs_ta{i} = f_ta;
    temp_sorts_ta{i} = temp_sort_ta;
    temps_ta{i} = temp_ta;
    recos_ta{i} = reco_ta;
    cbls_ta{i} = y_cbl_ta;
    cbhs_ta{i} = y_cbh_ta;
    q10s_ta{i} = q10_ta;
    cils_ta{i} = cil_ta;
    cius_ta{i} = ciu_ta;

    %construct half-hourly relationship using 5cm soil temp
    reco_og = site_list_nt.(site_names{i}).RECO;
    temp_og = site_list_nt.(site_names{i}).TS_5cm;
    [f_og,gof1_og,opt_og]=fit(temp_og,reco_og,'exp1');

    ci_og = confint(f_og);
    a_og = f_og.a;
    b_og = f_og.b;
    r2_og = gof1_og.rsquare;
    temp_sort_og = sort(temp_og);
    y_fit_og = a_og.*exp(b_og.*temp_sort_og);
    y_cbl_og = ci_og(1,1).*exp(ci_og(1,2).*temp_sort_og);
    y_cbh_og = ci_og(2,1).*exp(ci_og(2,2).*temp_sort_og);
    q10_og = exp(10.*f_og.b);
    cil_og = exp(10.*ci_og(1,2));
    ciu_og = exp(10.*ci_og(2,2));

    fs_og{i} = f_og;
    temp_sorts_og{i} = temp_sort_og;
    temps_og{i} = temp_og;
    recos_og{i} = reco_og;
    y_fits_og{i} = y_fit_og;
    cbls_og{i} = y_cbl_og;
    cbhs_og{i} = y_cbh_og;
    q10s_og{i} = q10_og;
    cils_og{i} = cil_og;
    cius_og{i} = ciu_og;
    r2s_og{i} = r2_og;

    %construct half-hourly relationships using 5cm soil temp during daytime
    reco_dt = site_list.(site_names{i}).RECO;
    temp_dt = site_list.(site_names{i}).TS_5cm;
    [f_dt,gof1_dt]=fit(temp_dt,reco_dt,'exp1');

    ci_dt = confint(f_dt);
    a_dt = f_dt.a;
    b_dt = f_dt.b;
    temp_sort_dt = sort(temp_dt);
    y_fit_dt = a_dt.*exp(b_dt.*temp_sort_dt);
    y_cbl_dt = ci_dt(1,1).*exp(ci_dt(1,2).*temp_sort_dt);
    y_cbh_dt = ci_dt(2,1).*exp(ci_dt(2,2).*temp_sort_dt);
    q10_dt = exp(10.*f_dt.b);
    cil_dt = exp(10.*ci_dt(1,2));
    ciu_dt = exp(10.*ci_dt(2,2));

    fs_dt{i} = f_dt;
    temp_sorts_dt{i} = temp_sort_dt;
    temps_dt{i} = temp_dt;
    recos_dt{i} = reco_dt;
    cbls_dt{i} = y_cbl_dt;
    cbhs_dt{i} = y_cbh_dt;
    q10s_dt{i} = q10_dt;
    cils_dt{i} = cil_dt;
    cius_dt{i} = ciu_dt;

    %construct daily relationships using 2cm soil temp
    site_list_2cm_daily.(site_names{i}) = retime(site_list_2cm.(site_names{i}), 'daily', 'mean'); 
    site_list_2cm_daily.(site_names{i})(any(isnan(site_list_2cm_daily.(site_names{i}).RECO), 2), :) = [];
    site_list_2cm_daily.(site_names{i})(any(isnan(site_list_2cm_daily.(site_names{i}).TS_2cm), 2), :) = [];
    reco_2cm = site_list_2cm_daily.(site_names{i}).RECO;
    temp_2cm = site_list_2cm_daily.(site_names{i}).TS_2cm;
    [f_2cm,gof1_2cm]=fit(temp_2cm,reco_2cm,'exp1');
    ci_2cm = confint(f_2cm);
    a_2cm = f_2cm.a;
    b_2cm = f_2cm.b;
    temp_sort_2cm = sort(temp_2cm);
    y_fit_2cm = a_2cm.*exp(b_2cm.*temp_sort_2cm);
    y_cbl_2cm = ci_2cm(1,1).*exp(ci_2cm(1,2).*temp_sort_2cm);
    y_cbh_2cm = ci_2cm(2,1).*exp(ci_2cm(2,2).*temp_sort_2cm);
    q10_2cm = exp(10.*f_2cm.b);
    cil_2cm = exp(10.*ci_2cm(1,2));
    ciu_2cm = exp(10.*ci_2cm(2,2));

    fs_2cm{i} = f_2cm;
    temp_sorts_2cm{i} = temp_sort_2cm;
    temps_2cm{i} = temp_2cm;
    recos_2cm{i} = reco_2cm;
    cbls_2cm{i} = y_cbl_2cm;
    cbhs_2cm{i} = y_cbh_2cm;
    q10s_2cm{i} = q10_2cm;
    cils_2cm{i} = cil_2cm;
    cius_2cm{i} = ciu_2cm;

    %construct daily relationships using 10cm soil temp
    site_list_10cm.(site_names{i}) = retime(site_list_10cm.(site_names{i}), 'daily', 'mean');
    site_list_10cm.(site_names{i})(any(isnan(site_list_10cm.(site_names{i}).TS_10cm), 2), :) = [];
    site_list_10cm.(site_names{i})(any(isnan(site_list_10cm.(site_names{i}).RECO), 2), :) = [];
    current_ts_10cm = site_list_10cm.(site_names{i}).TS_10cm;
    current_ts_10cm(any(isnan(current_ts_10cm), 2), :) = [];
    if isempty(current_ts_10cm)
        f_10cm = [];
        temp_sort_10cm = NaN;
        temp_10cm = NaN;
        reco_10cm = NaN;
        y_cbl_10cm = NaN;
        y_cbh_10cm = NaN;
        q10_10cm = NaN;
        cil_10cm = NaN;
        ciu_10cm = NaN;
    else
        reco_10cm = site_list_10cm.(site_names{i}).RECO;
        temp_10cm = site_list_10cm.(site_names{i}).TS_10cm;
        [f_10cm,gof1_10cm]=fit(temp_10cm,reco_10cm,'exp1');
        ci_10cm = confint(f_10cm);
        a_10cm = f_10cm.a;
        b_10cm = f_10cm.b;
        temp_sort_10cm = sort(temp_10cm);
        y_fit_10cm = a_10cm.*exp(b_10cm.*temp_sort_10cm);
        y_cbl_10cm = ci_10cm(1,1).*exp(ci_10cm(1,2).*temp_sort_10cm);
        y_cbh_10cm = ci_10cm(2,1).*exp(ci_10cm(2,2).*temp_sort_10cm);
        q10_10cm = exp(10.*f_10cm.b);
        cil_10cm = exp(10.*ci_10cm(1,2));
        ciu_10cm = exp(10.*ci_10cm(2,2));
    end
    fs_10cm{i} = f_10cm;
    temp_sorts_10cm{i} = temp_sort_10cm;
    temps_10cm{i} = temp_10cm;
    recos_10cm{i} = reco_10cm;
    cbls_10cm{i} = y_cbl_10cm;
    cbhs_10cm{i} = y_cbh_10cm;
    q10s_10cm{i} = q10_10cm;
    cils_10cm{i} = cil_10cm;
    cius_10cm{i} = ciu_10cm;

    %construct daily relationships using 5cm soil temp
    reco_daily = site_list_daily.(site_names{i}).RECO;
    temp_daily = site_list_daily.(site_names{i}).TS_5cm;
    [f_daily,gof1_daily]=fit(temp_daily,reco_daily,'exp1');

    ci_daily = confint(f_daily);
    a_daily = f_daily.a;
    b_daily = f_daily.b;
    r2_daily = gof1_daily.rsquare;
    temp_sort_daily = sort(temp_daily);
    y_fit_daily = a_daily.*exp(b_daily.*temp_sort_daily);
    y_cbl_daily = ci_daily(1,1).*exp(ci_daily(1,2).*temp_sort_daily);
    y_cbh_daily = ci_daily(2,1).*exp(ci_daily(2,2).*temp_sort_daily);
    q10_daily = exp(10.*f_daily.b);
    cil_daily = exp(10.*ci_daily(1,2));
    ciu_daily = exp(10.*ci_daily(2,2));

    fs_daily{i} = f_daily;
    temp_sorts_daily{i} = temp_sort_daily;
    temps_daily{i} = temp_daily;
    recos_daily{i} = reco_daily;
    cbls_daily{i} = y_cbl_daily;
    cbhs_daily{i} = y_cbh_daily;
    q10s_daily{i} = q10_daily;
    cils_daily{i} = cil_daily;
    cius_daily{i} = ciu_daily;
    r2s_daily{i} = r2_daily;

    %construct weekly relationships using 5cm soil temp
    site_list_weekly.(site_names{i})(any(isnan(site_list_weekly.(site_names{i}).RECO), 2), :) = [];
    site_list_weekly.(site_names{i})(any(isnan(site_list_weekly.(site_names{i}).TS_5cm), 2), :) = [];
    reco_weekly = site_list_weekly.(site_names{i}).RECO;
    temp_weekly = site_list_weekly.(site_names{i}).TS_5cm;
    [f_weekly,gof1_weekly]=fit(temp_weekly,reco_weekly,'exp1');

    ci_weekly = confint(f_weekly);
    a_weekly = f_weekly.a;
    b_weekly = f_weekly.b;
    r2_weekly = gof1_weekly.rsquare;
    temp_sort_weekly = sort(temp_weekly);
    y_fit_weekly = a_weekly.*exp(b_weekly.*temp_sort_weekly);
    y_cbl_weekly = ci_weekly(1,1).*exp(ci_weekly(1,2).*temp_sort_weekly);
    y_cbh_weekly = ci_weekly(2,1).*exp(ci_weekly(2,2).*temp_sort_weekly);
    q10_weekly = exp(10.*f_weekly.b);
    cil_weekly = exp(10.*ci_weekly(1,2));
    ciu_weekly = exp(10.*ci_weekly(2,2));

    fs_weekly{i} = f_weekly;
    temp_sorts_weekly{i} = temp_sort_weekly;
    temps_weekly{i} = temp_weekly;
    recos_weekly{i} = reco_weekly;
    cbls_weekly{i} = y_cbl_weekly;
    cbhs_weekly{i} = y_cbh_weekly;
    q10s_weekly{i} = q10_weekly;
    cils_weekly{i} = cil_weekly;
    cius_weekly{i} = ciu_weekly;
    r2s_weekly{i} = r2_weekly;

    %construct monthly relationships using 5cm soil temp
    site_list_monthly.(site_names{i})(any(isnan(site_list_monthly.(site_names{i}).RECO), 2), :) = [];
    site_list_monthly.(site_names{i})(any(isnan(site_list_monthly.(site_names{i}).TS_5cm), 2), :) = [];
    reco_monthly = site_list_monthly.(site_names{i}).RECO;
    temp_monthly = site_list_monthly.(site_names{i}).TS_5cm;
    [f_monthly,gof1_monthly]=fit(temp_monthly,reco_monthly,'exp1');

    ci_monthly = confint(f_monthly);
    a_monthly = f_monthly.a;
    b_monthly = f_monthly.b;
    r2_monthly = gof1_monthly.rsquare;
    temp_sort_monthly = sort(temp_monthly);
    y_fit_monthly = a_monthly.*exp(b_monthly.*temp_sort_monthly);
    y_cbl_monthly = ci_monthly(1,1).*exp(ci_monthly(1,2).*temp_sort_monthly);
    y_cbh_monthly = ci_monthly(2,1).*exp(ci_monthly(2,2).*temp_sort_monthly);
    q10_monthly = exp(10.*f_monthly.b);
    cil_monthly = exp(10.*ci_monthly(1,2));
    ciu_monthly = exp(10.*ci_monthly(2,2));

    fs_monthly{i} = f_monthly;
    temp_sorts_monthly{i} = temp_sort_monthly;
    temps_monthly{i} = temp_monthly;
    recos_monthly{i} = reco_monthly;
    cbls_monthly{i} = y_cbl_monthly;
    cbhs_monthly{i} = y_cbh_monthly;
    q10s_monthly{i} = q10_monthly;
    cils_monthly{i} = cil_monthly;
    cius_monthly{i} = ciu_monthly;
    r2s_monthly{i} = r2_monthly;

    %construct yearly relationships using 5cm soil temp
    site_list_yearly.(site_names{i})(any(isnan(site_list_yearly.(site_names{i}).RECO), 2), :) = [];
    site_list_yearly.(site_names{i})(any(isnan(site_list_yearly.(site_names{i}).TS_5cm), 2), :) = [];
    reco_yearly = site_list_yearly.(site_names{i}).RECO;
    temp_yearly = site_list_yearly.(site_names{i}).TS_5cm;
    [f_yearly,gof1_yearly]=fit(temp_yearly,reco_yearly,'exp1');
    
    ci_yearly = confint(f_yearly);
    a_yearly = f_yearly.a;
    b_yearly = f_yearly.b;
    r2_yearly = gof1_yearly.rsquare;
    temp_sort_yearly = sort(temp_yearly);
    y_fit_yearly = a_yearly.*exp(b_yearly.*temp_sort_yearly);
    y_cbl_yearly = ci_yearly(1,1).*exp(ci_yearly(1,2).*temp_sort_yearly);
    y_cbh_yearly = ci_yearly(2,1).*exp(ci_yearly(2,2).*temp_sort_yearly);
    q10_yearly = exp(10.*f_yearly.b);
    cil_yearly = exp(10.*ci_yearly(1,2));
    ciu_yearly = exp(10.*ci_yearly(2,2));

    fs_yearly{i} = f_yearly;
    temp_sorts_yearly{i} = temp_sort_yearly;
    temps_yearly{i} = temp_yearly;
    recos_yearly{i} = reco_yearly;
    cbls_yearly{i} = y_cbl_yearly;
    cbhs_yearly{i} = y_cbh_yearly;
    q10s_yearly{i} = q10_yearly;
    cils_yearly{i} = cil_yearly;
    cius_yearly{i} = ciu_yearly;
    r2s_yearly{i} = r2_yearly;
end

%calculate min and max range in q10 across sites
q10s_temporal_avg = horzcat(cell2mat(q10s_og)',cell2mat(q10s_daily)',cell2mat(q10s_weekly)',...
    cell2mat(q10s_monthly)',cell2mat(q10s_yearly)');
q10s_temporal_avg_min = min(q10s_temporal_avg(:,1:4),[],2);
q10s_temporal_avg_max = max(q10s_temporal_avg(:,1:4),[],2);
q10s_temporal_avg_range = q10s_temporal_avg_max-q10s_temporal_avg_min;
min(q10s_temporal_avg_range) %0.0968 at Oas
max(q10s_temporal_avg_range) %1.1584 at SCC

%create table of temporal averaged q10s and their 95% ci
cils_temporal_avg = horzcat(cell2mat(cils_og)',cell2mat(cils_daily)',cell2mat(cils_weekly)',...
    cell2mat(cils_monthly)',cell2mat(cils_yearly)');
cius_temporal_avg = horzcat(cell2mat(cius_og)',cell2mat(cius_daily)',cell2mat(cius_weekly)',...
    cell2mat(cius_monthly)',cell2mat(cius_yearly)');
r2s_temporal_avg = horzcat(cell2mat(r2s_og)',cell2mat(r2s_daily)',cell2mat(r2s_weekly)',...
    cell2mat(r2s_monthly)',cell2mat(r2s_yearly)');

%find the number of observations used to construct each functional relationship when temporally averging 
ns_temporal_avg = [];
for i = 1:16
n_temporal_avg = horzcat(size(temps_og{i},1),size(temps_daily{i},1),size(temps_weekly{i},1),...
    size(temps_monthly{i},1),size(temps_yearly{i},1));
ns_temporal_avg = vertcat(ns_temporal_avg,n_temporal_avg);
end

temporal_avg_q10s_95ci=string(round(q10s_temporal_avg,2))+'['+string(round(cils_temporal_avg,2))+', '+string(round(cius_temporal_avg,2))+']'+string(round(r2s_temporal_avg,2))+','+string(ns_temporal_avg);

%calculate the percentage of days where measurements are available at each site (for table 1)
for i = 1:numel(site_names)
    current_num_years = unique(site_list_daily.(site_names{i}).Year);
    available_days = size(current_num_years,1).*153;
    days_measured = size(site_list_daily.(site_names{i}),1);

    percent_days_avail{i} = (days_measured./available_days).*100;
end

%calculate the average temp range at each site
for i = 1:size(site_names,1)
    current_temp_sort_daily = temp_sorts_daily{i};
    current_temp_sort_monthly = temp_sorts_monthly{i};
    current_temp_sort_yearly = temp_sorts_yearly{i};
    temp_diff_daily = max(current_temp_sort_daily) - min(current_temp_sort_daily);
    temp_diff_monthly = max(current_temp_sort_monthly) - min(current_temp_sort_monthly);
    temp_diff_yearly = max(current_temp_sort_yearly) - min(current_temp_sort_yearly);
    temp_range_reduced = temp_diff_daily - temp_diff_monthly;
    temp_range_reduced_yearly = temp_diff_daily - temp_diff_yearly;
    temp_diffs_daily(i) = temp_diff_daily;
    temp_diffs_monthly(i) = temp_diff_monthly;
    temp_diffs_yearly(i) = temp_diff_yearly;
    temp_ranges_reduced(i) = temp_range_reduced;
    temp_ranges_reduced_yearly(i) = temp_range_reduced_yearly;
end
temp_diff_daily_avg = mean(temp_diffs_daily); %18.3915
temp_diff_monthly_avg = mean(temp_diffs_monthly); %12.2952
temp_range_reduced_avg = mean(temp_ranges_reduced); %6.0963
temp_range_reduced_yearly_avg = mean(temp_ranges_reduced_yearly); %15.4237

%colors for plotting
yellow_str = '#FDBD3C';
yellow = sscanf(yellow_str(2:end),'%2x%2x%2x',[1 3])/255;
green_str = '#80CB58';
green = sscanf(green_str(2:end),'%2x%2x%2x',[1 3])/255;
blue_str = '#2696EB';
blue = sscanf(blue_str(2:end),'%2x%2x%2x',[1 3])/255;
navy_str = '#3D26A7';
navy = sscanf(navy_str(2:end),'%2x%2x%2x',[1 3])/255;

plot_letters = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)'};

%plot impacts of temporal averaging at each site
fig = figure(); 
tlo = tiledlayout(4,4,"TileSpacing","compact");
h = gobjects(1,16);
for i = 1:16
    ax = nexttile(tlo);
    hold on 
    p_og = plot(ax,fs_og{i},temps_og{i},recos_og{i},'.');
    set(p_og,'color','k','LineWidth',3)
    delete(p_og(1))
    fill(ax,[temp_sorts_og{i}; flipud(temp_sorts_og{i})], [cbls_og{i}; flipud(cbhs_og{i})],'k', ...
        'EdgeColor', 'none', 'facealpha', 0.20)
    p_daily = plot(ax,fs_daily{i},temps_daily{i},recos_daily{i},'.');
    set(p_daily,'color',navy,'LineWidth',3,'LineStyle','--')
    delete(p_daily(1))
    % fill(ax,[temp_sorts_daily{i}; flipud(temp_sorts_daily{i})], [cbls_daily{i}; flipud(cbhs_daily{i})],'b', ...
    %     'EdgeColor', 'none', 'facealpha', 0.10)
    p_weekly = plot(ax,fs_weekly{i},temps_weekly{i},recos_weekly{i},'.');
    set(p_weekly,'color',blue,'LineWidth',3,'LineStyle','-.')
    delete(p_weekly(1))
    % fill(ax,[temp_sorts_weekly{i}; flipud(temp_sorts_weekly{i})], [cbls_weekly{i}; flipud(cbhs_weekly{i})],'b', ...
    %     'EdgeColor', 'none', 'facealpha', 0.10)
    p_monthly = plot(ax,fs_monthly{i},temps_monthly{i},recos_monthly{i},'.');
    set(p_monthly,'color',green,'LineWidth',3,'LineStyle',':')
    delete(p_monthly(1))
    % fill(ax,[temp_sorts_monthly{i}; flipud(temp_sorts_monthly{i})], [cbls_monthly{i}; flipud(cbhs_monthly{i})],'g', ...
    %     'EdgeColor', 'none', 'facealpha', 0.10)
    p_yearly = plot(ax,fs_yearly{i},temps_yearly{i},recos_yearly{i},'.');
    set(p_yearly,'color',yellow,'LineWidth',3)
    delete(p_yearly(1))
    box on;
    ax2 = gca; 
    ax2.FontSize = 14;
    Str = strrep(site_names_full{i},'Interpreter','none');
    title(Str,'Color','black','FontSize',14,'Interpreter','none');
    if i == 1
        legend("Half-hourly","","Daily",...
            "Weekly","Monthly",...
            "Seasonal",'location','northwest','FontSize',13)
        text(0.9,0.9,plot_letters{i}, 'FontWeight', 'bold', 'FontSize', 14,'Units','normalized')
    else
        b = gca; legend(b,'off');
        text(0.05,0.9,plot_letters{i}, 'FontWeight', 'bold', 'FontSize', 14,'Units','normalized')
    end
    xlabel(""); ylabel("");
    if i == 1
        ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    elseif i == 5
        ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    elseif i == 9
        ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    elseif i == 13
        xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
        ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    elseif i>13
        xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
    end
end

%AGU Poster
for i = 15%1:16
    %ax = nexttile(tlo);
    hold on 
    p_og = plot(fs_og{i},temps_og{i},recos_og{i},'.');
    set(p_og,'color','k','LineWidth',3)
    delete(p_og(1))
    fill([temp_sorts_og{i}; flipud(temp_sorts_og{i})], [cbls_og{i}; flipud(cbhs_og{i})],'k', ...
        'EdgeColor', 'none', 'facealpha', 0.20)
    p_daily = plot(fs_daily{i},temps_daily{i},recos_daily{i},'.');
    set(p_daily,'color',navy,'LineWidth',3,'LineStyle','--')
    delete(p_daily(1))
    % fill(ax,[temp_sorts_daily{i}; flipud(temp_sorts_daily{i})], [cbls_daily{i}; flipud(cbhs_daily{i})],'b', ...
    %     'EdgeColor', 'none', 'facealpha', 0.10)
    p_weekly = plot(fs_weekly{i},temps_weekly{i},recos_weekly{i},'.');
    set(p_weekly,'color',blue,'LineWidth',3,'LineStyle','-.')
    delete(p_weekly(1))
    % fill(ax,[temp_sorts_weekly{i}; flipud(temp_sorts_weekly{i})], [cbls_weekly{i}; flipud(cbhs_weekly{i})],'b', ...
    %     'EdgeColor', 'none', 'facealpha', 0.10)
    p_monthly = plot(fs_monthly{i},temps_monthly{i},recos_monthly{i},'.');
    set(p_monthly,'color',green,'LineWidth',3,'LineStyle',':')
    delete(p_monthly(1))
    % fill(ax,[temp_sorts_monthly{i}; flipud(temp_sorts_monthly{i})], [cbls_monthly{i}; flipud(cbhs_monthly{i})],'g', ...
    %     'EdgeColor', 'none', 'facealpha', 0.10)
    p_yearly = plot(fs_yearly{i},temps_yearly{i},recos_yearly{i},'.');
    set(p_yearly,'color',yellow,'LineWidth',3)
    delete(p_yearly(1))
    box on;
    ax2 = gca; 
    ax2.FontSize = 20;
    %Str = strrep(site_names_full{i},'Interpreter','none');
    %title(Str,'Color','black','FontSize',14,'Interpreter','none');
    legend("Half-hourly","","Daily",...
        "Weekly","Monthly",...
        "Seasonal",'location','northwest','FontSize',20)
    %text(0.9,0.9,plot_letters{i}, 'FontWeight', 'bold', 'FontSize', 14,'Units','normalized')
    xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',20)
    ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',20)
    xlim([-1 19])
    ylim([0 8])
end

%%%%%%%%%%
%How many daily measurements need to be taken each month? Construt
%functional relationships using 1 through 30 days per month for 500
%iterations in each set
%This part takes a long time to run, so outputs were saved in a .mat file
[man_fs,man_q10s,man_cils,man_cius,man_q10_counts,man_rdm_temps,man_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'Man',1);
[oas_fs,oas_q10s,oas_cils,oas_cius,oas_q10_counts,oas_rdm_temps,oas_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'Oas',2);
[obs_fs,obs_q10s,obs_cils,obs_cius,obs_q10_counts,obs_rdm_temps,obs_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'Obs',3);
[ojp_fs,ojp_q10s,ojp_cils,ojp_cius,ojp_q10_counts,ojp_rdm_temps,ojp_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'Ojp',4);
[scb_fs,scb_q10s,scb_cils,scb_cius,scb_q10_counts,scb_rdm_temps,scb_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'SCB',5);
[scc_fs,scc_q10s,scc_cils,scc_cius,scc_q10_counts,scc_rdm_temps,scc_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'SCC',6);
[sj2_fs,sj2_q10s,sj2_cils,sj2_cius,sj2_q10_counts,sj2_rdm_temps,sj2_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'SJ2',7);
[bzf_fs,bzf_q10s,bzf_cils,bzf_cius,bzf_q10_counts,bzf_rdm_temps,bzf_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'BZF',8);
[bzs_fs,bzs_q10s,bzs_cils,bzs_cius,bzs_q10_counts,bzs_rdm_temps,bzs_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'BZS',9);
[eml_fs,eml_q10s,eml_cils,eml_cius,eml_q10_counts,eml_rdm_temps,eml_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'EML',10);
[fcr_fs,fcr_q10s,fcr_cils,fcr_cius,fcr_q10_counts,fcr_rdm_temps,fcr_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'Fcr',11);
[ich_fs,ich_q10s,ich_cils,ich_cius,ich_q10_counts,ich_rdm_temps,ich_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'ICh',12);
[ics_fs,ics_q10s,ics_cils,ics_cius,ics_q10_counts,ics_rdm_temps,ics_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'ICs',13);
[prr_fs,prr_q10s,prr_cils,prr_cius,prr_q10_counts,prr_rdm_temps,prr_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'Prr',14);
[rpf_fs,rpf_q10s,rpf_cils,rpf_cius,rpf_q10_counts,rpf_rdm_temps,rpf_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'Rpf',15);
[uaf_fs,uaf_q10s,uaf_cils,uaf_cius,uaf_q10_counts,uaf_rdm_temps,uaf_rdm_reco] = data_randomization(site_list_daily,...
    site_list_monthly,cils_daily,cius_daily,'Uaf',16);

save("site_data_randomization_5cm_0703.mat","man_fs","man_q10s","man_cils","man_cius","man_q10_counts",...
    "oas_fs","oas_q10s","oas_cils","oas_cius","oas_q10_counts",...
    "obs_fs","obs_q10s","obs_cils","obs_cius","obs_q10_counts",...
    "ojp_fs","ojp_q10s","ojp_cils","ojp_cius","ojp_q10_counts",...
    "scb_fs","scb_q10s","scb_cils","scb_cius","scb_q10_counts",...
    "scc_fs","scc_q10s","scc_cils","scc_cius","scc_q10_counts",...
    "sj2_fs","sj2_q10s","sj2_cils","sj2_cius","sj2_q10_counts",...
    "bzf_fs","bzf_q10s","bzf_cils","bzf_cius","bzf_q10_counts",...
    "bzs_fs","bzs_q10s","bzs_cils","bzs_cius","bzs_q10_counts",...
    "eml_fs","eml_q10s","eml_cils","eml_cius","eml_q10_counts","eml_rdm_temps","eml_rdm_reco",...
    "fcr_fs","fcr_q10s","fcr_cils","fcr_cius","fcr_q10_counts",...
    "ich_fs","ich_q10s","ich_cils","ich_cius","ich_q10_counts",...
    "ics_fs","ics_q10s","ics_cils","ics_cius","ics_q10_counts",...
    "prr_fs","prr_q10s","prr_cils","prr_cius","prr_q10_counts",...
    "rpf_fs","rpf_q10s","rpf_cils","rpf_cius","rpf_q10_counts",...
    "uaf_fs","uaf_q10s","uaf_cils","uaf_cius","uaf_q10_counts")

%Plot degraded data with full data relationship
gray_str = '#3F3F3F';
gray = sscanf(gray_str(2:end),'%2x%2x%2x',[1 3])/255;

mylinestyles = repelem([{'-'},{'--'},{':'},{'-.'},{'--o'}], [7 7 7 7 7]);%{'-','--',':'};

%plot 1, 5, 10, 15, 20, 25 days per month for EML
tiledlayout(3,2);
nexttile
for n = 1:500
    current_iter_temp = eml_rdm_temps{n};
    current_iter_reco = eml_rdm_reco{n};
    hold on
    p1 = plot(eml_fs{n,1},current_iter_temp{1},current_iter_reco{1},'.');
    set(p1,'color',blue,'LineWidth',3)
    delete(p1(1))
end
hold on
p_og = plot(fs_daily{10},temps_daily{10},recos_daily{10},'.');
set(p_og,'color','k','LineWidth',3)
delete(p_og(1))
fill([temp_sorts_daily{10}; flipud(temp_sorts_daily{10})], [cbls_daily{10}; flipud(cbhs_daily{10})],'k', ...
     'EdgeColor', 'none', 'facealpha', 0.20)
xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
title('1 day per month','FontSize',14)
b = gca; legend(b,'off');
box on;
ax = gca;
ax.FontSize = 14; 
nexttile
for n = 1:500
    current_iter_temp = eml_rdm_temps{n};
    current_iter_reco = eml_rdm_reco{n};
    hold on
    p5 = plot(eml_fs{n,5},current_iter_temp{5},current_iter_reco{5},'.');
    set(p5,'color',blue,'LineWidth',3)
    delete(p5(1))
end
hold on
p_og = plot(fs_daily{10},temps_daily{10},recos_daily{10},'.');
set(p_og,'color','k','LineWidth',3)
delete(p_og(1))
fill([temp_sorts_daily{10}; flipud(temp_sorts_daily{10})], [cbls_daily{10}; flipud(cbhs_daily{10})],'k', ...
     'EdgeColor', 'none', 'facealpha', 0.20)
xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
title('5 days per month','FontSize',14)
b = gca; legend(b,'off');
box on;
ax = gca;
ax.FontSize = 14; 
nexttile
for n = 1:500
    current_iter_temp = eml_rdm_temps{n};
    current_iter_reco = eml_rdm_reco{n};
    hold on
    p10 = plot(eml_fs{n,10},current_iter_temp{10},current_iter_reco{10},'.');
    set(p10,'color',blue,'LineWidth',3)
    delete(p10(1))
end
hold on
p_og = plot(fs_daily{10},temps_daily{10},recos_daily{10},'.');
set(p_og,'color','k','LineWidth',3)
delete(p_og(1))
fill([temp_sorts_daily{10}; flipud(temp_sorts_daily{10})], [cbls_daily{10}; flipud(cbhs_daily{10})],'k', ...
     'EdgeColor', 'none', 'facealpha', 0.20)
xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
title('10 days per month','FontSize',14)
b = gca; legend(b,'off');
box on;
ax = gca;
ax.FontSize = 14; 
nexttile
for n = 1:500
    current_iter_temp = eml_rdm_temps{n};
    current_iter_reco = eml_rdm_reco{n};
    hold on
    p15 = plot(eml_fs{n,15},current_iter_temp{15},current_iter_reco{15},'.');
    set(p15,'color',blue,'LineWidth',3)
    delete(p15(1))
end
hold on
p_og = plot(fs_daily{10},temps_daily{10},recos_daily{10},'.');
set(p_og,'color','k','LineWidth',3)
delete(p_og(1))
fill([temp_sorts_daily{10}; flipud(temp_sorts_daily{10})], [cbls_daily{10}; flipud(cbhs_daily{10})],'k', ...
     'EdgeColor', 'none', 'facealpha', 0.20)
xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
title('15 days per month','FontSize',14)
b = gca; legend(b,'off');
box on;
ax = gca;
ax.FontSize = 14; 
nexttile
for n = 1:500
    current_iter_temp = eml_rdm_temps{n};
    current_iter_reco = eml_rdm_reco{n};
    hold on
    p20 = plot(eml_fs{n,20},current_iter_temp{20},current_iter_reco{20},'.');
    set(p20,'color',blue,'LineWidth',3)
    delete(p20(1))
end
hold on
p_og = plot(fs_daily{10},temps_daily{10},recos_daily{10},'.');
set(p_og,'color','k','LineWidth',3)
delete(p_og(1))
fill([temp_sorts_daily{10}; flipud(temp_sorts_daily{10})], [cbls_daily{10}; flipud(cbhs_daily{10})],'k', ...
     'EdgeColor', 'none', 'facealpha', 0.20)
xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
title('20 days per month','FontSize',14)
b = gca; legend(b,'off');
box on;
ax = gca;
ax.FontSize = 14; 
nexttile
for n = 1:500
    current_iter_temp = eml_rdm_temps{n};
    current_iter_reco = eml_rdm_reco{n};
    hold on
    p25 = plot(eml_fs{n,25},current_iter_temp{25},current_iter_reco{25},'.');
    set(p25,'color',blue,'LineWidth',3)
    delete(p25(1))
end
hold on
p_og = plot(fs_daily{10},temps_daily{10},recos_daily{10},'.');
set(p_og,'color','k','LineWidth',3)
delete(p_og(1))
fill([temp_sorts_daily{10}; flipud(temp_sorts_daily{10})], [cbls_daily{10}; flipud(cbhs_daily{10})],'k', ...
     'EdgeColor', 'none', 'facealpha', 0.20)
xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
title('25 days per month','FontSize',14)
b = gca; legend(b,'off');
box on;
ax = gca;
ax.FontSize = 14; 

%find the number of days where each site reaches 100 percent
sites_100_perc = [find(man_q10_counts == 100,1),find(oas_q10_counts == 100,1),find(obs_q10_counts == 100,1),...
    find(ojp_q10_counts == 100,1),find(scb_q10_counts == 100,1),find(scc_q10_counts == 100,1),...
    find(sj2_q10_counts == 100,1),find(bzf_q10_counts == 100,1),find(bzs_q10_counts == 100,1),...
    find(eml_q10_counts == 100,1),find(fcr_q10_counts == 100,1),find(ich_q10_counts == 100,1),...
    find(ics_q10_counts == 100,1),find(prr_q10_counts == 100,1),find(rpf_q10_counts == 100,1),find(uaf_q10_counts == 100,1)]';
min(sites_100_perc) %8
max(sites_100_perc) %21
mean(sites_100_perc) %17.56

avg_perc_15_days = (man_q10_counts(15)+oas_q10_counts(15)+obs_q10_counts(15)+ojp_q10_counts(15)+...
    scb_q10_counts(15)+scc_q10_counts(15)+sj2_q10_counts(15)+bzf_q10_counts(15)+bzs_q10_counts(15)+...
    eml_q10_counts(15)+fcr_q10_counts(15)+ich_q10_counts(15)+ics_q10_counts(15)+prr_q10_counts(15)+rpf_q10_counts(15)+uaf_q10_counts(15))/16;

avg_perc_10_days = (man_q10_counts(10)+oas_q10_counts(10)+obs_q10_counts(10)+ojp_q10_counts(10)+...
    scb_q10_counts(10)+scc_q10_counts(10)+sj2_q10_counts(10)+bzf_q10_counts(10)+bzs_q10_counts(10)+...
    eml_q10_counts(10)+fcr_q10_counts(10)+ich_q10_counts(10)+ics_q10_counts(10)+prr_q10_counts(10)+rpf_q10_counts(10)+uaf_q10_counts(10))/16;

avg_perc_5_day = (man_q10_counts(5)+oas_q10_counts(5)+obs_q10_counts(5)+ojp_q10_counts(5)+...
    scb_q10_counts(5)+scc_q10_counts(5)+sj2_q10_counts(5)+bzf_q10_counts(5)+bzs_q10_counts(5)+...
    eml_q10_counts(5)+fcr_q10_counts(5)+ich_q10_counts(5)+ics_q10_counts(5)+prr_q10_counts(5)+rpf_q10_counts(5)+uaf_q10_counts(5))/16;

avg_perc_1_day = (man_q10_counts(1)+oas_q10_counts(1)+obs_q10_counts(1)+ojp_q10_counts(1)+...
    scb_q10_counts(1)+scc_q10_counts(1)+sj2_q10_counts(1)+bzf_q10_counts(1)+bzs_q10_counts(1)+...
    eml_q10_counts(1)+fcr_q10_counts(1)+ich_q10_counts(1)+ics_q10_counts(1)+prr_q10_counts(1)+rpf_q10_counts(1)+uaf_q10_counts(1))/16;

%plot number of days required per month for each site
CM = slanCM('parula',19);
CM = CM(1:18,:);
CM = flipud(CM);
CT = {':','--','-'};
CT = horzcat(CT,CT,CT,CT,CT,CT,'-');
CM = vertcat(CM,[0 0 0]);
nan_leg = nan(19,1);

npoints_month = (1:30);
hold on
plot(npoints_month,ics_q10_counts,'color',CM(1,:),'LineWidth',2,'LineStyle',CT{:,1}) %58.10
plot(npoints_month,obs_q10_counts,'color',CM(2,:),'LineWidth',2,'LineStyle',CT{:,2}) %59.54
plot(npoints_month,ich_q10_counts,'color',CM(3,:),'LineWidth',2,'LineStyle',CT{:,3}) %59.69
plot(npoints_month,ojp_q10_counts,'color',CM(4,:),'LineWidth',2,'LineStyle',CT{:,4}) %67.84
plot(npoints_month,scc_q10_counts,'color',CM(5,:),'LineWidth',2,'LineStyle',CT{:,5}) %69.77
plot(npoints_month,sj2_q10_counts,'color',CM(6,:),'LineWidth',2,'LineStyle',CT{:,6}) %70.20
plot(npoints_month,scb_q10_counts,'color',CM(7,:),'LineWidth',2,'LineStyle',CT{:,7}) %80.39
plot(npoints_month,oas_q10_counts,'color',CM(8,:),'LineWidth',2,'LineStyle',CT{:,8}) %81.63
plot(npoints_month,eml_q10_counts,'color',CM(9,:),'LineWidth',2,'LineStyle',CT{:,9}) %82.48
plot(npoints_month,rpf_q10_counts,'color',CM(10,:),'LineWidth',2,'LineStyle',CT{:,10}) %89.07
plot(npoints_month,prr_q10_counts,'color',CM(11,:),'LineWidth',2,'LineStyle',CT{:,11}) %89.74
plot(npoints_month,fcr_q10_counts,'color',CM(13,:),'LineWidth',2,'LineStyle',CT{:,13}) %94.12
plot(npoints_month,man_q10_counts,'color',CM(12,:),'LineWidth',2,'LineStyle',CT{:,12}) %95.10
plot(npoints_month,uaf_q10_counts,'color',CM(14,:),'LineWidth',2,'LineStyle',CT{:,14}) %95.68
plot(npoints_month,bzs_q10_counts,'color',CM(15,:),'LineWidth',2,'LineStyle',CT{:,15}) %96.95
plot(npoints_month,bzf_q10_counts,'color',CM(16,:),'LineWidth',2,'LineStyle',CT{:,16}) %97.93
ylabel("Percentage iterations within original Q_{10} uncertainty")
xlabel("Days per month")
legend("US-Ics  (58%)","CA-Obs  (60%)","US-Ich  (60%)","CA-Ojp  (68%)","CA-SCC  (70%)","CA-SJ2  (70%)",...
    "CA-SCB  (80%)","CA-Oas  (82%)","US-EML  (82%)","US-Rpf  (89%)","US-Prr  (90%)",...
    "US-Fcr  (94%)","CA-Man  (95%)","US-Uaf  (96%)","US-BZS  (97%)","US-BZF  (98%)",'location','southeast','NumColumns',1)
legend('boxoff')
ylim([20 100])
xlim([1 30])
box on;
fontsize(14,"points")

%%%%%%%%%%
%To what temporal extent can functional relationships be applied across time?
%Construct functional relationships from a moving window of years
for i = 1:size(site_names,1)
    current_site = site_list_daily.(site_names{i});
    current_site_years = unique(current_site.Year);
    n=1;
    for j = 1:size(current_site_years,1)
        if j == 1
            f_yrs_all = cell(size(current_site_years,1),1);
            temp_yrs_all = cell(size(current_site_years,1),1);
            reco_yrs_all = cell(size(current_site_years,1),1);
            temp_sort_yrs_all = cell(size(current_site_years,1),1);
            y_cbl_yrs_all = cell(size(current_site_years,1),1);
            y_cbh_yrs_all = cell(size(current_site_years,1),1);
            q10_yrs_all = cell(size(current_site_years,1),1);
            yrs_subset_k_j = cell(size(current_site_years,1),1);
        end
    for k = 1:size(current_site_years,1)
        if k == 1
            f_yrs = cell(1,1);
            temp_yrs = cell(1,1);
            reco_yrs = cell(1,1);
            temp_sort_yrs = cell(1,1);
            y_cbl_yrs = cell(1,1);
            y_cbh_yrs = cell(1,1);
            q10_yrs = cell(1,1);
            yrs_subset_k = cell(1,1);
        end

        if k >= n
            years_subset = current_site_years(n:k);
            first_year = years_subset(1,:);
            current_year = years_subset(end);

            reco_yr_subset = current_site((current_site.Year >= first_year),:);
            temp_yr_subset = current_site((current_site.Year >= first_year),:);
            reco_yr_subset2 = reco_yr_subset.RECO((reco_yr_subset.Year <= current_year),:);
            temp_yr_subset2 = temp_yr_subset.TS_5cm((temp_yr_subset.Year <= current_year),:);
            [f_yr_subset,gof1_yr_subset]=fit(temp_yr_subset2,reco_yr_subset2,'exp1');
            ci_yr_subset = confint(f_yr_subset);
            a_yr_subset = f_yr_subset.a;
            b_yr_subset = f_yr_subset.b;
            temp_sort_subset = sort(temp_yr_subset2);
            y_fit_yr_subset = a_yr_subset.*exp(b_yr_subset.*temp_sort_subset);
            y_cbl_yr_subset = ci_yr_subset(1,1).*exp(ci_yr_subset(1,2).*temp_sort_subset);
            y_cbh_yr_subset = ci_yr_subset(2,1).*exp(ci_yr_subset(2,2).*temp_sort_subset);
            q10_yr_subset = exp(10.*f_yr_subset.b);
            cil_yr_subset = exp(10.*ci_yr_subset(1,2));
            ciu_yr_subset = exp(10.*ci_yr_subset(2,2));
        else
            years_subset = nan;
            first_year = nan;
            current_year = nan;
            reco_yr_subset = nan;
            temp_yr_subset = nan;
            reco_yr_subset2 = nan;
            temp_yr_subset2 = nan;
            f_yr_subset=nan;
            gof1_yr_subset=nan;
            ci_yr_subset = nan;
            a_yr_subset = nan;
            b_yr_subset = nan;
            temp_sort_subset = nan;
            y_fit_yr_subset = nan;
            y_cbl_yr_subset = nan;
            y_cbh_yr_subset = nan;
            q10_yr_subset = nan;
            cil_yr_subset = nan;
            ciu_yr_subset = nan;
        end

        f_yrs{k} = f_yr_subset;
        temp_yrs{k} = temp_yr_subset2;
        reco_yrs{k} = reco_yr_subset2;
        temp_sort_yrs{k} = temp_sort_subset;
        y_cbl_yrs{k} = y_cbl_yr_subset;
        y_cbh_yrs{k} = y_cbh_yr_subset;
        q10_yrs{k} = q10_yr_subset;
        yrs_subset_k{k} = years_subset;

        f_yrs_all{j} = f_yrs;
        temp_yrs_all{j} = temp_yrs;
        reco_yrs_all{j} = reco_yrs;
        temp_sort_yrs_all{j} = temp_sort_yrs;
        y_cbl_yrs_all{j} = y_cbl_yrs;
        y_cbh_yrs_all{j} = y_cbh_yrs;
        q10_yrs_all{j} = q10_yrs;
        yrs_subset_k_j{j}=yrs_subset_k;

        f_yrs_sites{i} = f_yrs_all;
        temp_yrs_sites{i} = temp_yrs_all;
        reco_yrs_sites{i} = reco_yrs_all;
        temp_sort_yrs_sites{i} = temp_sort_yrs_all;
        y_cbl_yrs_sites{i} = y_cbl_yrs_all;
        y_cbh_yrs_sites{i} = y_cbh_yrs_all;
        q10_yrs_sites{i} = q10_yrs_all;
        yrs_subset_all{i} = yrs_subset_k_j;
    end
    n=n+1;
    end
end

%remove nans from moving window
for i = 1:size(site_names,1)
    current_site_subset = yrs_subset_all{i};
    current_site_f = f_yrs_sites{i};
    current_site_temp = temp_yrs_sites{i};
    current_site_reco = reco_yrs_sites{i};
    current_site_temp_sort = temp_sort_yrs_sites{i};
    current_site_y_cbl = y_cbl_yrs_sites{i};
    current_site_y_cbh = y_cbh_yrs_sites{i};
    current_site_q10 = q10_yrs_sites{i};
    for j = 1:size(current_site_subset,1)
        current_yr_subset = current_site_subset{j};
        current2_site_f = current_site_f{j};
        current2_site_temp = current_site_temp{j};
        current2_site_reco = current_site_reco{j};
        current2_site_temp_sort = current_site_temp_sort{j};
        current2_site_y_cbl = current_site_y_cbl{j};
        current2_site_y_cbh = current_site_y_cbh{j};
        current2_site_q10 = current_site_q10{j};
        if j == 1
            yrs_subset_new = cell(size(current_site_subset,1),1);
            q10_subset_new = cell(size(current_site_subset,1),1);
            f_subset_new = cell(size(current_site_subset,1),1);
            reco_subset_new = cell(size(current_site_subset,1),1);
            temp_subset_new = cell(size(current_site_subset,1),1);
        end
        current_iter_subset = current_yr_subset;
        current_iter_subset(cellfun(@(current_iter_subset) any(isnan(current_iter_subset)),current_iter_subset)) = [];
        yrs_subset_new{j} = current_iter_subset;
        sites_yrs_subset_new{i} = yrs_subset_new;

        current_f_subset = current2_site_f;
        current_f_subset(cellfun(@(current2_site_q10) any(isnan(current2_site_q10)),current2_site_q10)) = [];
        f_subset_new{j} = current_f_subset;
        sites_yrs_f_new{i} = f_subset_new;

        current_reco_subset = current2_site_reco;
        current_reco_subset(cellfun(@(current_reco_subset) any(isnan(current_reco_subset)),current_reco_subset)) = [];
        reco_subset_new{j} = current_reco_subset;
        sites_yrs_reco_new{i} = reco_subset_new;

        current_temp_subset = current2_site_temp;
        current_temp_subset(cellfun(@(current_temp_subset) any(isnan(current_temp_subset)),current_temp_subset)) = [];
        temp_subset_new{j} = current_temp_subset;
        sites_yrs_temp_new{i} = temp_subset_new;

        current_temp_sort_subset = current2_site_temp_sort;
        current_temp_sort_subset(cellfun(@(current_temp_sort_subset) any(isnan(current_temp_sort_subset)),current_temp_sort_subset)) = [];
        temp_sort_subset_new{j} = current_temp_sort_subset;
        sites_yrs_temp_sort_new{i} = temp_sort_subset_new;

        current_y_cbl_subset = current2_site_y_cbl;
        current_y_cbl_subset(cellfun(@(current_y_cbl_subset) any(isnan(current_y_cbl_subset)),current_y_cbl_subset)) = [];
        y_cbl_subset_new{j} = current_y_cbl_subset;
        sites_yrs_y_cbl_new{i} = y_cbl_subset_new;

        current_y_cbh_subset = current2_site_y_cbh;
        current_y_cbh_subset(cellfun(@(current_y_cbh_subset) any(isnan(current_y_cbh_subset)),current_y_cbh_subset)) = [];
        y_cbh_subset_new{j} = current_y_cbh_subset;
        sites_yrs_y_cbh_new{i} = y_cbh_subset_new;

        current_q10_subset = current2_site_q10;
        current_q10_subset(cellfun(@(current_q10_subset) any(isnan(current_q10_subset)),current_q10_subset)) = [];
        q10_subset_new{j} = current_q10_subset;
        sites_yrs_q10_new{i} = q10_subset_new;
    end
end

%pull out all one year q10s at each site
for i = 1:size(site_names,1)
    current_site_q10 = sites_yrs_q10_new{i};
    current_site_subset = sites_yrs_subset_new{i};
    for j = 1:size(current_site_q10,1)
        if j == 1
            oneyear_site_yr_new = cell(1,1);
            oneyear_q10_new = cell(1,1);
        end
        current2_site_q10 = current_site_q10{j};
        current2_site_subset = current_site_subset{j};

        %grab first cell (year)
        oneyear_q10_subset = current2_site_q10{1};
        oneyear_site_subset = current2_site_subset{1};
        oneyear_q10_new{j} = oneyear_q10_subset;
        oneyear_site_yr_new{j} = oneyear_site_subset;
        oneyear_q10s_new{i} = oneyear_q10_new;
        oneyear_yrs_subset_new{i} = oneyear_site_yr_new;
    end
end

%for one year q10s at each site, how many fall within the 95% CI of the original q10?
for i = 1:size(oneyear_q10s_new,2)
    current_site_1yr_q10 = oneyear_q10s_new{i};
    current_site2_1yr_q10 = cell2mat(current_site_1yr_q10);
    current_site_yr_subset = oneyear_yrs_subset_new{i};
    current_site2_yr_subset = cell2mat(current_site_yr_subset);
    current_cil = cils_daily{i};
    current_ciu = cius_daily{i}; 

    count=current_site2_1yr_q10(current_site2_1yr_q10>=current_cil & current_site2_1yr_q10<=current_ciu);
    year_count=current_site2_yr_subset(current_site2_1yr_q10>=current_cil & current_site2_1yr_q10<=current_ciu);

    sites_1yr_count{i} = count;
    years_1yr_count{i} = year_count;
end
size_sites_1yr = cellfun(@length,years_1yr_count,'UniformOutput',false);
sum(cell2mat(size_sites_1yr)>0) %11 sites have at least 1 year
sum(cell2mat(size_sites_1yr)==0) %5 sites have no years

clear('current_site_f','current_site_temp','current_site_reco','current_site_temp_sort',...
    'current_site_y_cbl','current_site_y_cbh','current_site_q10','current2_site_f',...
    'current2_site_temp','current2_site_reco','current2_site_temp_sort',...
    'current2_site_y_cbl','current2_site_y_cbh','current2_site_q10')

CM = slanCM('parula',7); %parula, magma
CM = repelem(CM,3,1);
CM = CM(1:18,:);
CM = flipud(CM);
CT = {':','--','-'};
CT = horzcat(CT,CT,CT,CT,CT,CT,'-');
CM = vertcat(CM,[0 0 0]);
nan_leg = nan(19,1);

%plot functional relationships with varying temporal extents
fig2 = figure(); 
tlo2 = tiledlayout(4,4,"TileSpacing","compact");
h2 = gobjects(1,16);
for i = 1:size(site_names,1)
    ax3 = nexttile(tlo2);
    current_site_f = sites_yrs_f_new{i};
    current_site_temp = sites_yrs_temp_new{i};
    current_site_reco = sites_yrs_reco_new{i};
    current_site_temp_sort = sites_yrs_temp_sort_new{i};
    current_site_y_cbl = sites_yrs_y_cbl_new{i};
    current_site_y_cbh = sites_yrs_y_cbh_new{i};
    current_site_q10 = sites_yrs_q10_new{i};
    for j = 1:size(current_site_f,1)
        current2_site_f = current_site_f{j};
        current2_site_temp = current_site_temp{j};
        current2_site_reco = current_site_reco{j};
        current2_site_temp_sort = current_site_temp_sort{j};
        current2_site_y_cbl = current_site_y_cbl{j};
        current2_site_y_cbh = current_site_y_cbh{j};
        current2_site_q10 = current_site_q10{j};
        for k = 1:size(current2_site_f,2)
            hold on
            p_10yr = plot(ax3,current2_site_f{k},current2_site_temp{k},current2_site_reco{k},'.');
            set(p_10yr,'color',CM(k,:),'LineWidth',3,'LineStyle',CT{:,k})
            delete(p_10yr(1))
        end
    end
    p_daily = plot(ax3,fs_daily{i},temps_daily{i},recos_daily{i},'.');
    set(p_daily,'color','k','LineWidth',3)
    delete(p_daily(1))
    fill(ax3,[temp_sorts_daily{i}; flipud(temp_sorts_daily{i})], [cbls_daily{i}; flipud(cbhs_daily{i})],'k', ...
        'EdgeColor', 'none', 'facealpha', 0.20)
    xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
    ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    Str = strrep(site_names_full{i},'Interpreter','none');
    title(Str,'Color','black','FontSize',14,'Interpreter','none');
    box on;
    ax3 = gca; 
    ax3.FontSize = 14;
    text(0.65,0.07,size(site_list_yearly.(site_names{i}),1) + " years total", 'FontSize', 14,'Units','normalized')
    text(0.05,0.9,plot_letters{i}, 'FontWeight', 'bold', 'FontSize', 14,'Units','normalized')
    xlabel(""); ylabel("");

    if i == 1
        ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    elseif i == 5
        ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    elseif i == 9
        ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    elseif i == 13
        xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
        ylabel("Reco (gC m^{-2} d^{-1})",'FontSize',14)
    elseif i>13
        xlabel(['Soil Temperature (' char(176) 'C)'],'FontSize',14)
    end

    if i==16
        legend('','','','FontSize',14,'Location', 'eastoutside')
        leg_names = ["1 year","2 years","3 years","4 years","5 years","6 years","7 years",...
        "8 years","9 years","10 years","11 years","12 years","13 years","14 years",...
        "15 years","16 years","17 years","18 years","Full extent"];
        for m = 1:19
            plot(ax3,nan_leg, nan_leg,CT{:,m},'Color',CM(m,:),'LineWidth',3,'DisplayName',leg_names(m))
        end
    else
        b = gca; legend(b,'off');
    end
end

%plot swarmchart for temporal extents
fig2 = figure(); 
tlo2 = tiledlayout(4,4,"TileSpacing","compact");
h2 = gobjects(1,16);
for i = 1:size(sites_yrs_q10_new,2)
    ax3 = nexttile(tlo2);
    current_site_swarm = sites_yrs_q10_new{i};
    for j = 1:size(current_site_swarm,1)
        current_site_subset_swarm = current_site_swarm{j};
        for k = 1:size(current_site_subset_swarm,2)
            if current_site_subset_swarm{k} > 50
                current_site_subset_swarm{k} = NaN; %set q10 to nan if > 50
            end
            hold on
            s=swarmchart(ax3,k,current_site_subset_swarm{k},40,CM(k,:),'f');
            s.MarkerEdgeColor = 'k';
        end
    end
    errorbar(ax3,size(current_site_swarm,1),q10s_daily{i},q10s_daily{i}-cils_daily{i},cius_daily{i}-q10s_daily{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);  
    yline(cils_daily{i})
    yline(cius_daily{i})
    xlim([0.5 size(current_site_swarm,1)+0.5])
    xlabel("Years",'FontSize',14)
    ylabel("Q_{10}",'FontSize',14)
    Str = strrep(site_names_full{i},'Interpreter','none');
    title(Str,'Color','black','FontSize',14,'Interpreter','none');
    box on;
    ax3 = gca; 
    ax3.FontSize = 14;

    if i==16
        legend('','','','FontSize',14,'Location', 'eastoutside')
        leg_names = ["1 year","2 years","3 years","4 years","5 years","6 years","7 years",...
        "8 years","9 years","10 years","11 years","12 years","13 years","14 years",...
        "15 years","16 years","17 years","18 years","Full extent"];
        for m = 1:19
            plot(ax3,nan_leg, nan_leg,'-','Color',CM(m,:),'LineWidth',3,'DisplayName',leg_names(m))
        end
    else
        b = gca; legend(b,'off');
    end
end

%AGU poster
%fig2 = figure(); 
%tlo2 = tiledlayout(4,4,"TileSpacing","compact");
%h2 = gobjects(1,16);
for i = 2%1:size(sites_yrs_q10_new,2)
    %ax3 = nexttile(tlo2);
    current_site_swarm = sites_yrs_q10_new{i};
    for j = 1:size(current_site_swarm,1)
        current_site_subset_swarm = current_site_swarm{j};
        for k = 1:size(current_site_subset_swarm,2)
            if current_site_subset_swarm{k} > 50
                current_site_subset_swarm{k} = NaN; %set q10 to nan if > 50
            end
            hold on
            s=swarmchart(k,current_site_subset_swarm{k},100,CM(k,:),'f');
            s.MarkerEdgeColor = 'k';
        end
    end
    errorbar(size(current_site_swarm,1),q10s_daily{i},q10s_daily{i}-cils_daily{i},cius_daily{i}-q10s_daily{i},'o','color','k','linewidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10);  
    yline(cils_daily{i})
    yline(cius_daily{i})
    xlim([0.5 size(current_site_swarm,1)+0.5])
    xlabel("Years",'FontSize',20)
    ylabel("Q_{10}",'FontSize',20)
    %Str = strrep(site_names_full{i},'Interpreter','none');
    %title(Str,'Color','black','FontSize',14,'Interpreter','none');
    box on;
    ax3 = gca; 
    ax3.FontSize = 20;
    ylim([1.9 3.5])
        % legend('','','','FontSize',14,'Location', 'eastoutside')
        % leg_names = ["1 year","2 years","3 years","4 years","5 years","6 years","7 years",...
        % "8 years","9 years","10 years","11 years","12 years","13 years","14 years",...
        % "15 years","16 years","17 years","18 years","Full extent"];
        % for m = 1:19
        %     plot(nan_leg, nan_leg,'-','Color',CM(m,:),'LineWidth',3,'DisplayName',leg_names(m))
        % end
end

%%%%%%%%%%
%Look at US-Rpf over time
start_time_rpf = datetime(2009,05,01);
end_time_rpf = datetime(2020,09,30);
time_range_rpf = (start_time_rpf:days(1):end_time_rpf)';

rpf = site_list_daily.Rpf;
rpf = retime(rpf,time_range_rpf); %insert -9999 for missing data

dl = datetime('1-May-2009');
dr = datetime('30-Sep-2020');
hold on
plot(rpf.TIMESTAMP,rpf.RECO,'k')
ylabel("Reco (gC m^{-2} d^{-1})")
xlim([dl dr])
box on;
fontsize(14,"points")

%Construct functional relationships at Rpf until 2013 and after 2014
rpf_shortterm = site_list_daily.Rpf((site_list_daily.Rpf.Year < 2014),:);
rpf_longterm = site_list_daily.Rpf((site_list_daily.Rpf.Year >= 2014),:);

[f_short,gof1_short]=fit(rpf_shortterm.TS_5cm,rpf_shortterm.RECO,'exp1');
ci_short = confint(f_short);
a_short = f_short.a;
b_short = f_short.b;
temp_sort_short = sort(rpf_shortterm.TS_5cm);
y_fit_short = a_short.*exp(b_short.*temp_sort_short);
y_cbl_short = ci_short(1,1).*exp(ci_short(1,2).*temp_sort_short);
y_cbh_short = ci_short(2,1).*exp(ci_short(2,2).*temp_sort_short);
q10_short = exp(10.*f_short.b);
cil_short = exp(10.*ci_short(1,2));
ciu_short = exp(10.*ci_short(2,2));

[f_long,gof1_long]=fit(rpf_longterm.TS_5cm,rpf_longterm.RECO,'exp1');
ci_long = confint(f_long);
a_long = f_long.a;
b_long = f_long.b;
temp_sort_long = sort(rpf_longterm.TS_5cm);
y_fit_long = a_long.*exp(b_long.*temp_sort_long);
y_cbl_long = ci_long(1,1).*exp(ci_long(1,2).*temp_sort_long);
y_cbh_long = ci_long(2,1).*exp(ci_long(2,2).*temp_sort_long);
q10_long = exp(10.*f_long.b);
cil_long = exp(10.*ci_long(1,2));
ciu_long = exp(10.*ci_long(2,2));

blue2_str = '#2C90CB';
blue2 = sscanf(blue2_str(2:end),'%2x%2x%2x',[1 3])/255;
red2_str = '#DD443A';
red2 = sscanf(red2_str(2:end),'%2x%2x%2x',[1 3])/255;

hold on 
p_short = plot(f_short,rpf_shortterm.TS_5cm,rpf_shortterm.RECO,'.');
set(p_short,'color',red2,'LineWidth',4)
fill([temp_sort_short; flipud(temp_sort_short)], [y_cbl_short; flipud(y_cbh_short)],'r', ...
    'EdgeColor', 'none', 'facealpha', 0.10)
p_long = plot(f_long,rpf_longterm.TS_5cm,rpf_longterm.RECO,'.');
set(p_long,'color',blue2,'LineWidth',4)
fill([temp_sort_long; flipud(temp_sort_long)], [y_cbl_long; flipud(y_cbh_long)],'b', ...
    'EdgeColor', 'none', 'facealpha', 0.10)
xlabel(['Soil Temperature (' char(176) 'C)'])
ylabel("Reco (gC m^{-2} d^{-1})")
box on;
ylim([0 18])
xlim([-1 17])
legend("","Short-term response (2009-2012)","","","Long-term response (2014-2020)",'location','northwest')
text(12,17,'Short-term Q_{10}: 4.03')
text(12,16,'Long-term Q_{10}: 3.56')
fontsize(14,"points")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%How can functional relationships be applied across space?
%First plot Q10s for each site

%q10 at bzs and bzf
bzs_overall_q10 = q10s_daily{9}; %2.9073
bzf_overall_q10 = q10s_daily{8}; %1.5064

%sort based on landcover type
%deciduous broadleaf forest, evergreen needleleaf forest, open shrublands, permanent wetlands
ix2 = [2,15,1,3,4,6,7,9,14,16,10,11,12,5,8,13];
sites_sort2 = site_names_full(ix2);
m2 = q10s_daily(ix2);
cils_sort2=cell2mat(cils_daily(ix2)); cius_sort2=cell2mat(cius_daily(ix2));

%min and max q10 across sites
min(cell2mat(q10s_daily)) %1.5064
max(cell2mat(q10s_daily)) %3.9218

%plot q10s for 2cm, 5cm, and 10cm soil temps
hold on
errorbar(1:16,cell2mat(q10s_ta),cell2mat(q10s_ta)-cell2mat(cils_ta),cell2mat(cius_ta)-cell2mat(q10s_ta),'o','color',yellow,'linewidth',2,'MarkerEdgeColor',yellow,'MarkerFaceColor','white',...
        'MarkerSize',10);
errorbar(1:16,cell2mat(q10s_2cm),cell2mat(q10s_2cm)-cell2mat(cils_2cm),cell2mat(cius_2cm)-cell2mat(q10s_2cm),'o','color',green,'linewidth',2,'MarkerEdgeColor',green,'MarkerFaceColor','white',...
        'MarkerSize',10);
errorbar(1:16,cell2mat(q10s_daily),cell2mat(q10s_daily)-cell2mat(cils_daily),cell2mat(cius_daily)-cell2mat(q10s_daily),'o','color',blue,'linewidth',2,'MarkerEdgeColor',blue,'MarkerFaceColor','white',...
        'MarkerSize',10);
errorbar(1:16,cell2mat(q10s_10cm),cell2mat(q10s_10cm)-cell2mat(cils_10cm),cell2mat(cius_10cm)-cell2mat(q10s_10cm),'o','color',navy,'linewidth',2,'MarkerEdgeColor',navy,'MarkerFaceColor','white',...
        'MarkerSize',10);
xlim([0.25 16.75])
ylim([1.25 5.01])
set(gca,'xtick',1:1:16);
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[site_names_full],'fontsize',14)
xtickangle(45)
ylabel('Intercept','fontsize',14);
box on;
ylabel('Q_{10} (unitless)','fontsize',14);
legend("Air temperature","2 cm soil temperature","5 cm soil temperature","10 cm soil temperature",'location','northwest')

%Is q10 at various depth correlated with the mean temperature?
for i = 1:16
    current_temp_ta=temps_ta{i};
    mean_current_temp_ta = nanmean(current_temp_ta);
    temps_ta_mean{i} = mean_current_temp_ta;

    current_temp_2cm=temps_2cm{i};
    mean_current_temp_2cm = nanmean(current_temp_2cm);
    temps_2cm_mean{i} = mean_current_temp_2cm;

    current_temp_daily=temps_daily{i};
    mean_current_temp_daily = nanmean(current_temp_daily);
    temps_daily_mean{i} = mean_current_temp_daily;

    current_temp_10cm=temps_10cm{i};
    mean_current_temp_10cm = nanmean(current_temp_10cm);
    temps_10cm_mean{i} = mean_current_temp_10cm;
end

hold on
errorbar(1:16,cell2mat(temps_ta_mean),nan,nan,'o','color','k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','white',...
        'MarkerSize',10);
errorbar(1:16,cell2mat(temps_2cm_mean),nan,nan,'o','color','k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
errorbar(1:16,cell2mat(temps_daily_mean),nan,nan,'^','color','k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','white',...
        'MarkerSize',10);
errorbar(1:16,cell2mat(temps_10cm_mean),nan,nan,'^','color','k','linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerSize',10);
xlim([0.25 16.75])
set(gca,'xtick',1:1:16);
set(gca,'TickLabelInterpreter','none')
set(gca,'xticklabel',[site_names_full],'fontsize',14)
xtickangle(45)
ylabel(['Temperature (' char(176) 'C)'],'FontSize',14)
box on;
legend("Air temp","","","","","","","","","","","","","","","",...
    "2 cm soil temp","","","","","","","","","","","","","","","",...
    "5 cm soil temp","","","","","","","","","","","","","","","",...
    "10 cm soil temp",'location','northeast')

%Plot map of flux tower locations
tower_table = readtable("/Users/jmp838/Desktop/Research/Project 2 - Functional Benchmarks/Flux Towers/Ameriflux/flux_tower_lat_lon.csv","TreatAsMissing","NA");

hold on
m_proj('mercator','lon',[-170 -80],'lat',[50 72]);
m_grid('ytick',[50 55 60 65 70],'xtick',[-160 -140 -120 -100 -80],...
    'yticklabels',[50 55 60 65 70],'xticklabels',[-160 -140 -120 -100 -80],'linest','-')
m_coast('patch',[0.9 0.9 0.9]);
m_scatter(tower_table.Longitude,tower_table.Latitude,'filled','k')
for k=1:16
    m_text(round(tower_table.Longitude(k),2),round(tower_table.Latitude(k),2),tower_table.Site(k));
end
