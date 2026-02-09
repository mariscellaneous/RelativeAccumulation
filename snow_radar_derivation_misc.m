df_array = Data(slm:end,:)-Data(1:end-(slm-1),:);
df_array = movmean(df_array,slm,1);
q = find(~isfinite(df_array));
df_array(q) = NaN;
dmn_array = nanmean(df_array,1); 
dsd_array = nanstd(df_array,1);
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

dataoutnew = smooth([1:tr],srf,0.02,'rloess');
dataout = [ (1:tr)' srf' dataoutnew ];

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
   

good = isnan(srf);
d2(:,good) = oData(:,good);

srflarge = srf > size(oData,1);
srf(srflarge) = size(oData,1);

for k = 1:tr
    if srflarge(tr)==1
        d2(:,k) = circshift(oData(:,k),200-srf(k));
    end
end

undersurface = srf > 200;

for k = 1:tr
    if undersurface == 1
        d2(end+(200-srf(k)+1):end,k) = NaN;
    end
end

%d2(end+(200-srf

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
