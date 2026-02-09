function radar_processing(date,spc,slm,ROOTPATH)
tic
%spc = 100;, slm = 9;, ryear = 2016;

% loads each radar file in a directory and picks the surface, flattens the data,
% then averages over a user defined trace window

% requires "dist", "medfilt2", "lowess", and MERRA-2 P-E & 2-m T


root_dir = strcat(ROOTPATH,date,'/');
input_dir = strcat(root_dir,'mat_files/');
output_dir = strcat(root_dir,'flat100/');
ryear = str2double(date(1:4));

% if the output folder does not exist, create it!
if exist(output_dir,'dir') == 0 
    mkdir(output_dir);
end

if exist(strcat(root_dir,'fig_flat100/'),'dir') == 0 
    mkdir(strcat(root_dir,'fig_flat100/')); % also make folder to store the figures
end

fnames = dir(strcat(input_dir,'IRSNO1B*.mat'));
nm = size(fnames,1); 
nm

% empty matrices for carrying over traces to next file
catdat = [];
catLon = [];
catLat = [];
xdat = [];
xLon = [];
xLat = [];
xTM = [];
xTSD = [];
xGPS = [];
warning('off','all')

% let's add MERRA-2 P-E and temperature, then we can get a good
% approximation of the depth density profile

load(strcat('merra2_pme.mat'),'lat','lon','ann','yr');
ny = find(yr == ryear);
M2_yr = yr(1:ny);
pme = ann(:,:,1:ny);

pme_mn = mean(ann,3);
pme_sd = std(ann,[],3);

load(strcat('merra2_t2m.mat'),'ann');
t2m = ann(:,:,1:ny);

for j = 1:nm %:10 %:nm 
    
    outfilename = strcat(root_dir,'fig_flat100/',fnames(j,1).name(1:end-4),'.jpg');
    if exist(outfilename,'file') ~= 0
         continue
    end
    fname = strcat(input_dir,fnames(j,1).name);
    load(fname,'Data','Latitude','Longitude','Time','Truncate_Mean','Truncate_Std_Dev','GPS_Time');

    tr = size(Data,2);
    del_t = Time(2)-Time(1);    %two way travel time of each range bin
    srf = NaN(1,tr);    % empty array for surface pick
    d2 = NaN.*Data;
    oData = Data;   % keep original data
    oLon = Longitude;   % keep original Longitude values
    [range]=dist(Latitude,Longitude);   % get the distance between adjacent traces
    xval2 = cumsum([0 range]) ./ 1000;  %cumulative distance
    %Data = medfilt2(oData,[slm floor(slm/2)],'symmetric');  % apply a median filter to minimize noise
    Data = ordfilt2(oData,floor(slm*floor(slm/2)/8),ones(floor(slm/4),floor(slm/2)));
    % line up the traces by maximum slope for surface picking
    
    %-- start brooke version 
%         for k = 1:tr    % for each trace
%             dat = Data(:,k);    % grab return power
%             df = dat(slm:end) - dat(1:end-(slm-1)); % determine slope over the user defined window
%             df = movmean(df,slm);   % smooth to eliminate spurious data
%             q = find(isfinite(df) == 0);    % remove non-finite data
%             df(q) = NaN;
%             dmn = nanmean(df);  % mean slope
%             dsd = nanstd(df);   % standard deviation of slope
%             [~,v] = max(df);    % get the location of the maximum slope
%             if v <= numel(df)-1
%                 while v <= numel(df)-1 && df(v+1) > dmn + dsd*3
%                     v = v + 1;  % get the "deepest" slope above 3 sigma of slopes, and make the "surface"
%                 end
%             end
%             srf(k) = v + round(slm/2);    % account for the slope window by adding half the window size to surface pick
%             if v == numel(df)
%                 [~,v] = max(dat);         % if end, just make the surface the maximum return
%                 srf(k) = v;  
%             end
%         end
    %-- start brooke version 
    
    %-- start marissa version 
    
    df_array = Data(slm:end,:)-Data(1:end-(slm-1),:);
    df_array = movmean(df_array,slm,1);
    q = find(~isfinite(df_array));
    df_array(q) = NaN;
    dmn_array = nanmean(df_array,1); 
    dsd_array = nanstd(df_array,1); % not sure why this isn't the same as nanstd of a single column??? slightly off
    [~,v_array] = max(df_array);
    
    for k=1:tr
        if v_array(k) <= numel(df_array(:,k))
            while v_array(k) <= numel(df_array(:,k))-1 && df_array((v_array(k)+1),k) > dmn_array(k) + dsd_array(k)*3
                v_array(k) = v_array(k) + 1; 
            end
        end
    end
    
    srf = v_array + round(slm/2);    % account for the slope window by adding half the window size to surface pic        
    badind = v_array == length(Data);
    df_array_alt = df_array;
    df_array_alt(:,~badind) = NaN; 
    
    [~,v_arrayalt] = max(df_array_alt);  % if end, just make the surface the maximum return
    
    srf(badind) = NaN;
    v_good = nansum([srf;v_arrayalt]);
    srf = v_good;
   
    
    %--- End marissa version
    
    %-- Brooke's way
    dataout = lowess([(1:tr)' srf'],0.02,0); % generate a smooth surface
    %--
    
    %-- Marissa Alternate 
    %dataoutnew = smooth([1:tr],srf,0.02,'rloess');
    %dataout = [ (1:tr)' srf' dataoutnew ];
    %--
    
    sm = round(dataout(:,3))'; % round to bin index

    wdw = 15;   % window for finding surface outliers, can be modified to user preference
    bp = abs(sm - srf); % find difference between smooth surface and original surface
    q = find(bp > wdw); % if larger than window tolerance
    srf(q) = sm(q); % set equal to the smoothed surface at that trace

    % this loop shifts the orginal Data based on the surface picks, here I
    % artifically set the zero surface to range bin 200.  If there is no
    % surface pick, it does not shift it
    
%         -- Brooke Code
    for k = 1:tr
        if isnan(srf(k)) == 1
            d2(:,k) = oData(:,k);
        else
            if srf(k) > size(oData,1)
                srf(k) = size(oData,1);
            end
            d2(:,k) = circshift(oData(:,k),200-srf(k));
            if srf(k) > 200
                d2(end+(200-srf(k)+1):end,k) = NaN;
            end
        end
    end

%         -- end Brooke code 
    
%         %-- Marissa code -- FIX THIS -- 
%         good = isnan(srf);
%         d2(:,good) = oData(:,good);
%         
%         srflarge = srf > size(oData,1);
%         srf(srflarge) = size(oData,1);
%         
%         for k = 1:tr
%             if srflarge(tr)==1
%                 d2(:,k) = circshift(oData(:,k),200-srf(k));
%             end
%         end
%         
%         undersurface = srf > 200;
%   
%         for k = 1:tr
%             if undersurface == 1
%                 d2(end+(200-srf(k)+1):end,k) = NaN;
%             end
%         end
%         
%         %d2(end+(200-srf
    
%         for k = 1:tr
%             if isnan(srf(k)) == 1
%                 d2(:,k) = oData(:,k);
%             else
%                 if srf(k) > size(oData,1)
%                     srf(k) = size(oData,1);
%                 end
%                 d2(:,k) = circshift(oData(:,k),200-srf(k));
%                 if srf(k) > 200
%                     d2(end+(200-srf(k)+1):end,k) = NaN;
%                 end
%             end
%         end
    %-- end Marissa code 



    % combine current file with the data left over from the previous (if
    % any)
    temp = NaN(size(d2,1),size(xdat,2));
    v = min([size(d2,1) size(xdat,1)]);
    temp(1:v,:) = xdat(1:v,:);
    d2 = [temp d2];% combine actual retrun power
    clear temp
    Latitude = [xLat; Latitude];   % combine all ancillary data as well
    Longitude = [xLon; Longitude];
    Truncate_Mean = [xTM; Truncate_Mean];
    Truncate_Std_Dev = [xTSD; Truncate_Std_Dev];
    GPS_Time = [xGPS; GPS_Time];

    % data are flattened (and combined), now we need to average to get close to the
    % desired spacing
    tr = size(d2,2);
    range = dist(Latitude,Longitude);   % get trace spacing
    tdist = mean(range);    % approximate spacing distance
    clear range
    tavg = round(spc / tdist);  % find the approx. number of traces needed to equal the user defined spacing
    cols = floor(tr / tavg);    % number of columns in final product, based on the spacing and number of traces
    xt = rem(tr, tavg);     % number of traces left over to save for next file
    i = 1;
    
%-- Brooke code       
%         while i <= cols         % loop over the desired number of traces
%             ndat(:,i) = mean(d2(:,1:tavg),2);   % average data together
%             nLat(i) = mean(Latitude(1:tavg));
%             nLon(i) = mean(Longitude(1:tavg));
%             nTM(i) = mean(Truncate_Mean(1:tavg));
%             nTSD(i) = mean(Truncate_Std_Dev(1:tavg));
%             nGPS(i) = mean(GPS_Time(1:tavg));
%             d2(:,1:tavg) = [];  % remove the data used, loop until all traces complete
%             Latitude(1:tavg) = [];
%             Longitude(1:tavg) = [];
%             Truncate_Mean(1:tavg) = [];
%             Truncate_Std_Dev(1:tavg) = [];
%             GPS_Time(1:tavg) = [];
%             i = i+1
%         end

%         % extra data for next file is stored
%         xdat = d2;
%         xLon = Longitude;
%         xLat = Latitude;
%         xTM = Truncate_Mean;
%         xTSD = Truncate_Std_Dev;
%         xGPS = GPS_Time;
%-- Brooke code ending     

%-- Marissa Coding 

    [r,c] = size(d2);
    nlay = floor(c/tavg);
    d2_new = d2(:,1:int16(floor(nlay*tavg)));
    out = permute(reshape(d2_new,[r,int16(tavg),int16(nlay)]),[3,1,2]);
    ndat = mean(out,3)';
    
    [rl,~] = size(Longitude);
    twodsize = 1:int16(floor(rl/tavg)*tavg);
    long = Longitude(twodsize);
    lati = Latitude(twodsize);
    Truncate_Mean_new = Truncate_Mean(twodsize);
    Truncate_Std_Dev_new = Truncate_Std_Dev(twodsize);
    GPS_Time_new = GPS_Time(twodsize);
  
    newshape = [int16(floor(rl/nlay)),int16(nlay)];
    nLon = mean(reshape(long,newshape));
    nLat = mean(reshape(lati,newshape));
    nTM = mean(reshape(Truncate_Mean_new,newshape));
    nTSD = mean(reshape(Truncate_Std_Dev_new,newshape));
    nGPS = mean(reshape(GPS_Time_new,newshape));
   
    
    % Remove excess variables 
    clear d2_new Truncate_Mean_new Truncate_Std_Dev_new GPS_Time_new
    
    % Excess for later 
    start = int16(nlay*tavg);
    if nlay == 0 || tavg == 0
        continue
    end
    xdat = d2(:,start:end);
    xLon = Longitude(start:end);
    xLat = Latitude(start:end);
    xTM = Truncate_Mean(start:end);
    xTSD = Truncate_Std_Dev(start:end);
    xGPS = GPS_Time(start:end);  
   
% -- End marissa coding 



    % approximate depths
    ttt = ((1:size(Data,1))'-200).*del_t ./ 10^6;   % two-way travel time relative to surface at bin 200
    rho = 0.50; % guess at material density, used can change, but this is not important
    e = (1+0.845.*rho).^2;  % use Kovacs relationship between density and dielectric constant
    c = 299792458;  % speed of light
    d = (c .* ttt) / (2 *sqrt(e));  % approximate depths based on the simple density assumption
    [range] = dist(nLat,nLon);  % trace separation in meters
    xval = cumsum([0 range]) ./ 1000;   % cumulative distance in kilometers


    %% BEGIN PLOTTING SECTION

    mnL = min(nLon);
    mxL = max(nLon);

    %qr = sprintf('_nLon_%.2f_xLon_%.2f',mnL,mxL);  % for pole hole flight only

    zz = oLon(end) - oLon(1);
    % make plot of surface pick
    figure('Visible','off')
    if size(ndat,1) < 2000
        kl = size(ndat,1);
    else
        kl = 2000;
    end

    d1 = nanmedian(ndat(200:end,:),2);
    mp2(1) = nanmedian(d1);
    mp2(2) = nanstd(nanstd(~isnan(ndat(200:end,:))));

    ll = mp2(1) - mp2(2).*75;    % color scale, could be changed....
    ul = mp2(1) + mp2(2).*350;

    clim = [ll, ul];

    if isnan(ll)|| isnan(ul) || ~isfinite(ll) || ~isfinite(ul)
        do_clim = 'No';
    else
        do_clim = 'Yes';
    end

    do_clim = 'No';



    if strcmp(do_clim,'No')
        imagesc(xval,d(100:kl),ndat(100:kl,:))
        ax1 = gca;
        ax1.Box = 'off';
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');

        imagesc(ax2,nLon,d(100:kl),ndat(100:kl,:))
        ax2.XAxisLocation = 'top';
        ax2.YAxisLocation = 'right';        
    else

        imagesc(xval,d(100:kl),ndat(100:kl,:),clim)
        %imagesc(xval,d(100:kl),ndat(100:kl,:))
        ax1 = gca;
        ax1.Box = 'off';
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');

        imagesc(ax2,nLon,d(100:kl),ndat(100:kl,:),clim)
        %imagesc(ax2,nLon,d(100:kl),ndat(100:kl,:))
        ax2.XAxisLocation = 'top';
        ax2.YAxisLocation = 'right';   

    end
    ax2.Box = 'off';

    colormap(flipud(brewermap([],'Blues')))
    colorbar

    ax1.XLabel.String = 'Distance (km)';
    ax1.XLabel.FontSize = 14;
    ax1.YLabel.String = 'Approx. Depth (m)';
    ax1.YLabel.FontSize = 14;
    ax2.XLabel.String = 'Longitude';
    ax2.XLabel.FontSize = 14;
    outerpos1 = ax1.OuterPosition;
    ti1 = ax1.TightInset; 
    outerpos2 = ax2.OuterPosition;
    ti2 = ax2.TightInset; 
    left = outerpos1(1) + ti1(1);
    bottom = outerpos1(2) + ti1(2);
    ax_width = outerpos1(3) - ti1(1) - ti2(3);
    ax_height = outerpos1(4) - ti1(2) - ti2(4);
    ax1.Position = [left bottom ax_width ax_height];
    ax2.Position = [left bottom ax_width ax_height];
    ax3 = axes('Position',ax2.Position);
    plot(ax3,xval,repmat(0,size(xval)),'r');
    ax3.Box = 'off';
    ax3.Color = 'none';
    ax3.XLim = ax1.XLim;
    ax3.YLim = ax1.YLim;
    ax3.YColor = 'none';
    ax3.YDir = 'reverse';
    if zz < 0
        ax2.XDir = 'reverse';
    end
    %saveas(gcf,strcat(output_dir,'fig\',fnames(j,1).name(1:end-4),qr,'.jpg'))   % Pole Hole flight only
    saveas(gcf,strcat(root_dir,'fig_flat100/',fnames(j,1).name(1:end-4),'.jpg'))

    % plot the orginal data with surface pick
    figure('Visible','off')

    uy = nanmin(srf)-50;
    ly = nanmax(srf)+50;
    if isnan(uy) == 1 || uy < 1
        uy = 1;
    end
    if isnan(ly) == 1 || ly > size(Data,1)
        ly = size(Data,1);
    end

    if strcmp(do_clim,'No')
        imagesc(xval2,1:size(Data,1),Data)
        ax1 = gca;
        ax1.Box = 'off';
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');

        imagesc(ax2,oLon,1:size(Data,1),Data);
        ax2.XAxisLocation = 'top';
        ax2.YAxisLocation = 'right';  
    else
        %imagesc(xval2,1:size(Data,1),Data)
        imagesc(xval2,1:size(Data,1),Data,clim)
        ax1 = gca;
        ax1.Box = 'off';
        ax1_pos = ax1.Position; % position of first axes
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none');

        imagesc(ax2,oLon,1:size(Data,1),Data,clim);
        %imagesc(ax2,oLon,1:size(Data,1),Data);
        ax2.XAxisLocation = 'top';
        ax2.YAxisLocation = 'right';  
    end
    ax2.Box = 'off';
    colormap(flipud(brewermap([],'Blues')))
    colorbar
    hold on

    ax1.YLim = [uy ly];
    ax2.YLim = [uy ly];

    ax1.XLabel.String = 'Distance (km)';
    ax1.XLabel.FontSize = 14;
    ax1.YLabel.String = 'Range Bin';
    ax1.YLabel.FontSize = 14;
    ax2.XLabel.String = 'Longitude';
    ax2.XLabel.FontSize = 14;
    outerpos1 = ax1.OuterPosition;
    ti1 = ax1.TightInset; 
    outerpos2 = ax2.OuterPosition;
    ti2 = ax2.TightInset; 
    left = outerpos1(1) + ti1(1);
    bottom = outerpos1(2) + ti1(2);
    ax_width = outerpos1(3) - ti1(1) - ti2(3);
    ax_height = outerpos1(4) - ti1(2) - ti2(4);
    ax1.Position = [left bottom ax_width ax_height];
    ax2.Position = [left bottom ax_width ax_height];
    ax3 = axes('Position',ax2.Position);



    if strcmp(do_clim,'No')
        imagesc(Data)
    else
        imagesc(Data,clim)
    end

    hold on
    plot(ax3,srf,'r');
    ax3.Box = 'off';
    ax3.Color = 'none';
    ax3.YLim = ax1.YLim;
    ax3.YColor = 'none';
    ax3.XColor = 'none';
    ax3.YDir = 'reverse';
    if zz < 0
        ax2.XDir = 'reverse';
    end
    ax4 = axes('Position',ax2.Position);
    ax4.XLim = ax1.XLim;
    ax4.YLim = ax1.YLim;
    ax4.Color = 'none';
    ax4.XTickLabel = [];
    ax4.YTickLabel = [];
    %saveas(gcf,strcat(output_dir,'fig\',fnames(j,1).name(1:end-4),qr,'_orig.jpg'))  % Pole Hole only
    saveas(gcf,strcat(root_dir,'fig_flat100/',fnames(j,1).name(1:end-4),'_orig.jpg'))
    Data = ndat;
    del_t = Time(2) - Time(1);
    % save flat data
    Latitude = nLat;
    Longitude = nLon;
    Truncate_Mean = nTM;
    Truncate_Std_Dev = nTSD;
    GPS_Time = nGPS;
    % find the nearest MERRA-2 pixels
    M2_Latitude = NaN(1,size(Data,2));
    M2_Longitude = NaN(1,size(Data,2));
    M2_PmE = NaN(size(pme,3),size(Data,2));
    M2_T2m = NaN(size(pme,3),size(Data,2));
    for w = 1:size(Data,2)
        [~,delLa] = min(abs(lat - Latitude(w)));
        [~,delLo] = min(abs(lon - Longitude(w)));
        M2_Latitude(w) = lat(delLa);
        M2_Longitude(w) = lon(delLo);
        M2_PmE(:,w);
        pme(delLo,delLa,:);
        M2_PmE(:,w) = reshape(squeeze(pme(delLo,delLa,:)),[size(pme,3) 1]);
        M2_T2m(:,w) = reshape(squeeze(t2m(delLo,delLa,:)),[size(pme,3) 1]);
    end
    save(strcat(output_dir,'/',fnames(j,1).name),'Data','Latitude','Longitude','Truncate_Mean','Truncate_Std_Dev','GPS_Time','srf','tavg','spc','del_t','M2_Latitude','M2_Longitude','M2_PmE','M2_T2m','M2_yr')
    clear ndat nLon nLat nTM nTSD nGPS srf srf2 d2 Latitude Longitude Truncate_Mean Truncate_Std_Dev GPS_Time Data
    close all

toc
end
    
    