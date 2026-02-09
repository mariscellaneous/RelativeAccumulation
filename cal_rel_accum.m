function cal_rel_accum(date,ROOTPATH)


    % calculate accumulation rates from layer picks and MERRA-2 climate

    % fnames = dir('G:\SnowRadar\20161109\mat_files\flat100\comb\IRSNO1B*.mat');
    % input_dir = 'G:\SnowRadar\20161109\mat_files\flat100\comb\';

    %date = input('Date Range? ');
    % load the flattened and stacked echograms


    fnames = dir(strcat(ROOTPATH,date,'/comb/IRSNO1B*.mat'));
    input_dir = strcat(ROOTPATH,date,'/comb/');

    ryear = str2double(date(1:4));
    nm = size(fnames,1); 
    rho0 = 0.35;
    m2_acc = [];
    sr_acc = [];
    nlat = [];
    nlon = [];
    fid = [];
    rel_acc = [];
    per_acc = [];
    load('merra2_pme.mat','lat','lon','ann','yr');
    ny = find(yr == ryear);
    pme = ann(:,:,1:ny);
    pme_mn = mean(pme,3);
    clear pme

    load('merra2_t2m.mat','ann');
    t2m = ann(:,:,1:ny);
    t2m_mn = mean(t2m,3);
    clear t2m
    Y = repmat(lon,[1 numel(lat)]);
    X = repmat(lat',[numel(lon) 1]);
    for i = 1:nm
        load(strcat(input_dir,fnames(i).name))
        a = exist('lay');
        if a == 1
            % first let's get an interpolated view of M-2 climate
            mpme = interp2(X,Y,pme_mn,Latitude,Longitude);
            mt2m = interp2(X,Y,t2m_mn,Latitude,Longitude);
            
            % loop over each pick to estimate the age needed at each trace to
            % match the MERRA-2 accumulation
            ages = NaN .* lay;
            for q = 1:size(Data,2)
                if isfinite(lay(q)) == 1
                    acc = mpme(q) ./ 1000;  % mean accumulation rate across the entire file
                    temp = mt2m(q);   % mean temperature across the entire file
                    [depth, rho, age] = HL(acc,temp,rho0);  % get density, depth, and age from density model

                    % what is the mean "depth: of the layer in bins
                    twt = (lay(q) - 193) .* del_t ./  10^6;
                    % determine the permittivity
                    [e] = perm(rho);
                    % then convert depths to twtt based on the permittivity
                    tt = d2twtt(e);   % seconds
                    nD = interp1(tt,depth,twt);       % get mean depth of layer in meters
                    % now get mean age from density model
                    nA = interp1(depth,age,nD);       % load in g/cm^2
                    ages(q) = nA;
                end
            end
            nA = nanmean(ages);
            
    %         % loop over each pick, get density profile based on M2 climate,
    %         % claculate depths in m w.e., then estimate the age based on the M2
    %         % accumulation rate
    %         p = find(isfinite(lay) == 1);
    %         acc = mean(mean(M2_PmE(:,p),1)) ./ 1000;  % mean accumulation rate across the entire file
    %         temp = mean(mean(M2_T2m(:,p),1));   % mean temperature across the entire file
    %         [depth, rho, age] = HL(acc,temp,rho0);  % get density, depth, and age from density model
    % 
    %         % what is the mean "depth: of the layer in bins
    %         twt = (mean(lay(p)) - 193) .* del_t ./  10^6;
    %         % determine the permittivity
    %         [e] = perm(rho);
    %         % then convert depths to twtt based on the permittivity
    %         tt = d2twtt(e);   % seconds
    %         nD = interp1(tt,depth,twt);       % get mean depth of layer in meters
    %         % now get mean age from density model
    %         nA = interp1(depth,age,nD);       % load in g/cm^2

            % now run the interative procedure to get self-consistent density
            % profile and accumulation rate
            accR = NaN .* lay;
            accM = NaN .* lay;
            for q = 1:size(Data,2)
                iacc = mpme(q) ./ 1000;  % mean accumulation rate at site
                itemp = mt2m(q);   % mean temperature at site
                if isfinite(lay(q)) == 1
                    itwt = (lay(q) - 193) .* del_t ./  10^6;
                    nacc = convacc(rho0,itemp,iacc,itwt,nA);
                    accR(q) = nacc;
                    accM(q) = iacc;
                else
                    accR(q) = NaN;
                    accM(q) = iacc;
                end
            end
            fid = [fid repmat(i,[1 size(Data,2)])];
            m2_acc = [m2_acc accM];
            sr_acc = [sr_acc accR];
            rel_acc = [rel_acc accR - acc];
            per_acc = [rel_acc (accR - acc) ./ acc];
            nlat = [nlat Latitude];
            nlon = [nlon Longitude];
            save(strcat(input_dir,fnames(i).name),'-append','accR')
        end
        clear lay

            
    end          
end