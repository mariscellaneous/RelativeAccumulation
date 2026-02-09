function comb_echo(date,ROOTPATH)


% quickly combine the surface-flattened and stacked ~ 5km sgements into
% chunks of user defined length


spc = 25;   % make 25 km segments


root_dir = strcat(ROOTPATH,date,'/');
input_dir = strcat(root_dir,'flat100/');
output_dir = strcat(root_dir,'comb/');
rws = 10000;    % number of range bins to store, can change if you want

datediff = strrep(date,'.','');

%yname = strcat(input_dir,'IRSNO1B_20161114_02_1*'); %specific 
yname = strcat(input_dir,'IRSNO1B_*'); % all points 

fnames = dir(yname);
nm = size(fnames,1); 
nm

% if the output folder does not exist, create it!
if exist(output_dir,'dir') == 0 
    mkdir(output_dir);
end

if exist(strcat(root_dir,'fig_comb/'),'dir') == 0 
    mkdir(strcat(root_dir,'fig_comb/'));
end

fr = spc ./ 5;  % each echogram is ~5 km
nf = floor(nm ./ fr);   % number of files to make

elm = 1;
while nm > 0
    % if there are not enough to fit desired spacing, use whats left
    if nm >= fr
        t = fr;
    else
        t = nm;
    end
    % determine the number of traces that we need
    nb = 0;
    nlat = [];
    nlon = [];
    delt = 0;
    for ind = 1:t            
        fname = strcat(input_dir,fnames(ind,1).name);
        load(fname,'Latitude','Longitude','M2_yr');
        if ind > 1
            rng = dist([nlat;Latitude(1)],[nlon;Longitude(1)]);
        end
        if delt > 1000      % if the following file is more than a kilometer from the end of the previous file 
            t = ind - 1;  %stop concatenation
        else
            nb = nb + numel(Latitude);  % otherwise, add to dat
            nlat = Latitude(end);
            nlon = Longitude(end);
        end
    end
   % create nan matrices to fill
   cdat = NaN(rws,nb);
   clat = NaN(1,nb); clon = clat; mlat = clat; mlon = clat; mid = clat; cgt = clat; ctm = clat; ctsd = clat;
   mpme = NaN(numel(M2_yr),nb); mt2m = mpme;
   s1 = 1;
   for ind = 1:t
       fname = strcat(input_dir,fnames(ind,1).name);
       load(fname);
       for g = 1:size(Data,2)
           mval(g) = nanmean(nanmean(Data(1:150,max([1 g-5]):min([g+5 size(Data,2)]))));
           rval(g) = nanmean(nanmean(Data(200:220,max([1 g-5]):min([g+5 size(Data,2)]))));
       end

       qq = find(isfinite(mval) == 0);
       mval(qq) = 0;
       rval = rval - mval;
       Data = (Data - repmat(mval,size(Data,1),1)) ./ repmat(rval,size(Data,1),1);
       q = min([rws size(Data,1)]);
       s2 = s1 + size(Data,2) - 1;
       cdat(1:q,s1:s2) = Data(1:q,:);
       clat(:,s1:s2) = Latitude;
       clon(:,s1:s2) = Longitude;
       mlat(:,s1:s2) = M2_Latitude;
       mlon(:,s1:s2) = M2_Longitude;
       %mid(:,s1:s2) = M2_id;
       mpme(:,s1:s2) = M2_PmE;
       mt2m(:,s1:s2) = M2_T2m;
       cgt(:,s1:s2) = GPS_Time;
       ctm(:,s1:s2) = Truncate_Mean;
       ctsd(:,s1:s2) = Truncate_Std_Dev;
       s1 = s1 + size(Data,2);
       mval = [];
       rval = [];
    end
   
    % get the approximate depth
    ttt = ((1:size(cdat,1))'-200).*del_t ./ 10^6;
    rho = 0.50;
    e = (1+0.845.*rho).^2;
    c = 299792458;
    d = (c .* ttt) / (2 *sqrt(e));
    [range]=dist(clat,clon);
    xval = cumsum([0 range]) ./ 1000;

    % make a plot

    zz = clon(end) - clon(1);
    % make plot of surface pick
    f = figure('Visible','off','Units','normalized');
    f.Position = [0 0 1 0.75];
    if size(cdat,1) < 5000
        kl = size(cdat,1);
    else
        kl = 5000;
    end

    %imagesc(xval,d(100:kl),cdat(100:kl,:),[0.1 0.9])
    imagesc(xval,d(100:kl),cdat(100:kl,:))
    ax1 = gca;
   
    ax1.Box = 'off';
    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none');
    %imagesc(ax2,clon,d(100:kl),cdat(100:kl,:),[0.1 0.9])
    imagesc(ax2,clon,d(100:kl),cdat(100:kl,:))
    ax2.XAxisLocation = 'top';
    ax2.YAxisLocation = 'right'; 
    colormap(flipud(brewermap([],'Blues')))
    ax2.Box = 'off';

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

    saveas(gcf,strcat(root_dir,'fig_comb/',fnames(t,1).name(1:end-10),num2str(elm,'%03.f'),'.jpg'))

    close all
   
   % now save.
   Data = cdat;
   Latitude = clat;
   Longitude = clon;
   M2_Latitude = mlat;
   M2_Longitude = mlon;
   M2_id = mid;
   M2_PmE = mpme;
   M2_T2m = mt2m;
   Truncate_Mean = ctm;
   GPS_Time = cgt;
   Truncate_Std_Dev = ctsd;
     
   save(strcat(output_dir,fnames(t,1).name(1:end-10),num2str(elm,'%03.f'),'.mat'),'Data','Longitude','Latitude','del_t','GPS_Time',...
            'M2_id','M2_Latitude','M2_Longitude','M2_PmE','M2_T2m','M2_yr','Truncate_Mean','Truncate_Std_Dev')
   elm = elm + 1;
   fnames(1:t) = [];
    nlat = [];
    nlon = [];
    nm = numel(fnames);
end

end
   