% start at the beginning of a flight, flip through echograms, and pick a
% single horizon either automatically or manually.
% When moving to the next echogram, try to use the same layer if possible,
% but its not critical

% load the flattened and stacked echograms
date = '2017.11.16'; % FINISH DOING THIS ONE, goes to ~10, just ~20 more!

% CHANGE ROOTPATH BELOW!!!!
ROOTPATH = '/Users/mada4363/Desktop/Computing_Accum/SAMPLE_DATA/'


fnames = dir(strcat(ROOTPATH,date,'/comb/IRSNO1B*.mat'));
input_dir = strcat(ROOTPATH,date,'/comb/');

% fnames = dir('G:\Pole_Hole\Snow_Radar\2016\mat_files\flat100\comb\IRSNO1B*.mat');
% input_dir = 'G:\Pole_Hole\Snow_Radar\2016\mat_files\flat100\comb\';

nm = size(fnames,1); 
%filt = [3 10 3;1 2 1; 0 0 0; -1 -2 -1; -3 -10 -3];  % Sobel filter from
%Brooke
filt = [-1 0 1;-2 0 2; -1 0 -1; 1 2 1; -8 -15 0];

lp = [];
%positionval = [500 500 400 300];
positionval = [0 0 0.5 1];
wdw = 25;
% loop over each file, you will be asked if you want to pick a layer or
% skip to the next one

for i = 2:nm
    load(strcat(input_dir,fnames(i).name))
    fnames(i).name
    % take a median filter\
    Data = movmean(Data,[5 3]);
    % get the approximate depth
    ttt = ((1:size(Data,1))'-200).*del_t ./ 10^6;
    rho = 0.50;
    e = (1+0.845.*rho).^2;
    c = 299792458;
    d = (c .* ttt) / (2 *sqrt(e));
    [range]=dist(Latitude,Longitude);
    xval = cumsum([0 range]) ./ 1000;
    nval = sum(isnan(Data),2);
    i1 = find(nval == 0,1,'first');
    i2 = find(nval == 0,1,'last');
    % normalize every trace to mean and standard deviation
    for g = 1:size(Data,2)
        mval(g) = nanmean(nanmean(Data(1:150,max([1 g-5]):min([g+5 size(Data,2)]))));
        rval(g) = nanmean(nanmean(Data(200:220,max([1 g-5]):min([g+5 size(Data,2)]))));
    end

    qq = find(isfinite(mval) == 0);
    mval(qq) = 0;
    rval = rval - mval;
    Data = (Data - repmat(mval,size(Data,1),1)) ./ repmat(rval,size(Data,1),1);
    % use the filtered data for picking, it enhances edges
    dat_ = conv2(Data,filt,'same');
    %dat_ = Data;
    %dat_ = edge3(dat_,'Sobel',0.0001);
    %dat_ = imsharpen(dat_,'Radius',2,'Amount',1);
    %dat_ = imadjust(dat_,[0.001 0.999],[]);
    %bob = 0:size(dat_,1)-1;
    %cat = zeros(size(dat_,2),1);
    %[bobs,cats]=meshgrid(bob,cat);
    %dat_= dat_.*bobs'/1000;
    %dat_ = imreducehaze(dat1,0.9);
    %dat_ = conv2(dat_,filt2,'same');
    % display a figure with the echogram and decide if you want to pick a
    % layer
    h1 = figure('units','normalized','position',positionval);
    %imagesc(dat_,[-0.15 0.5])
    imagesc(dat_)
    %set(gca,'Logarithmic');
    colormap(flipud(brewermap([],'Blues')))
    %colormap(bone) %brooke method
    %caxis([-0.15 0.5]) %brooke method
    
    colorbar();
    
    
    %caxis([0 0.4])
    ylim([100 2000])
    xlabel('Trace Number','FontSize',14)
    ylabel('Range Bin','FontSize',14)
    % now do you want to pick?
    kk = input('Do you want to pick layers (0) automatically, (1) manually, or (2) not at all?: ');
    
    % start the picking procedure
    if kk == 2
        close all
    end
    
    while kk <= 1

        if kk == 0

            % are you ready to pick?
            v = input('How do you want to pick a point? (0)use last pick or (1)select point: ');

                if v == 0
                    x = 1;
                    y = lp(end);
                else
                    [x, y] = ginput(1);
                    % get the peak close to that pick
                    x = round(x);
                    y = round(y);
                end
                z = y-wdw:y+wdw;
                dat = dat_(z,x);

                [peakLoc,peakMag] = peakfinder(dat*-1,(max(dat*-1)-min(dat*-1))/4);
                [~,l] = max(peakMag*-1);
                lay(x) = z(peakLoc(l));
    %         end
            x1 = x + 1;
            nm = size(dat_,2);
            % pick to the right
            for j = x1:nm
                z = lay(j-1)-wdw:lay(j-1)+wdw;
                dat = dat_(z,j);
                [peakLoc,peakMag] = peakfinder(dat*-1,(max(dat*-1)-min(dat*-1))/4);
                [~,l] = max(peakMag*-1);
                lay(j) = z(peakLoc(l));
            end
            x1 = x - 1;

            while x1 >= 1
                z = lay(x1+1)-wdw:lay(x1+1)+wdw;
                dat = dat_(z,x1);
                [peakLoc,peakMag] = peakfinder(dat,(max(dat)-min(dat))/4);
                [~,l] = max(peakMag);
                lay(x1) = z(peakLoc(l));
                x1 = x1 - 1;
            end

            imagesc(dat_)
            hold on
            e1 = plot(lay,'r');
            colormap(flipud(brewermap([],'Blues')))
            colorbar
            ylim([100 2000])
            xlabel('Trace Number','FontSize',14)
            ylabel('Range Bin','FontSize',14)

            rr = input('Are you satisfied with pick? (1) yes, (0) no: ');
            if rr == 0
                delete(e1);
                gg = input('Do you want to change the range window?: (0)no, (1)yes: ');
                if gg == 1
                    wdw = input('Enter new range bin window: ');
                end
                kk = input('Do you want to pick layers (0) automatically, (1) manually, or (2) not at all?: ');
            else
                lp = lay;
                save(strcat(input_dir,fnames(i).name),'-append','lay')
                kk = 2;
            end
        % if you want to spline pick....    
        elseif kk == 1

            yy(1) = 100;
            yy(2) = 2000;
            
            hold on
            % now we need to pick the layers 
            ww = input('Do you want to pick a layer - Yes (1), No (0): ');
            cnt = 1;
            if ww == 0
                close all
                kk = 2;
                continue
            else

                while ww == 1
                    [x,y] = getline;
                    Vq = interp1(x,y,1:size(Data,2),'pchip');
                    q1 = plot(Vq,'r');
                    qq = input('Are you satisfied - Yes (1), No (0): ');
                    while qq == 0
                        delete(q1)
                        [x,y] = getline;
                        Vq = interp1(x,y,1:size(Data,2),'pchip');
                        q1 = plot(Vq,'r');
                        qq = input('Are you satisfied - Yes (1), No (0): ');
                    end
                    % using the intial pick, we now will pick peaks
                    wdw = 15;
                    for j = 1:size(Data,2)
                        pk = round(Vq(j));
                        ny  = (pk-wdw:pk+wdw);
                        dat = dat_(ny,j);
                        %dat = smo(ny,j);
                        [pLoc,pMag] = peakfinder(dat*-1,(max(dat*-1)-min(dat*-1))/4,[],1);
                        [~,oo] = min(abs(ny(pLoc)-pk));
                        %[~,oo] = max(pMag);
                        if numel(oo) == 0
                            lay(cnt,j) = NaN;
                        else
                            lay(cnt,j) = ny(pLoc(oo));
                        end
                    end
                    delete(q1)
                    e1 = plot(lay(cnt,:),'y');
                    %ww = input('Do you want to pick a layer - Yes (1), No (0): ');
                    ww = 0;     % only pick one layer
                    %cnt = cnt + 1;
                    ylim([yy(1) yy(2)]);
                end
                ww = input('Do you want to remove segments of the picks - Yes (1), No (0): ');
                while ww == 1
                    vv = input('Do you want to start from the left (1), start from the right (2), or choose a segment in the middle (3)?: ');
                    if vv == 1
                        [xx,~] = ginput(1);
                        lay(:,1:xx) = NaN;
                    else if vv == 2
                            [xx,~] = ginput(1);
                            lay(:,xx:end) = NaN;
                        else if vv == 3
                            [xx,~] = ginput(2);
                            lay(:,xx(1):xx(2)) = NaN;
                        end
                        end
                    end
                    delete(e1)
                    e1 = plot(lay(cnt,:),'y');
                    ww = input('?Do you want to remove segments of the picks - Yes (1), No (0): ');
                end
                
                rr = input('Are you satisfied with pick? (1) yes, (0) no: ');
                if rr == 0
                    delete(e1);
                    kk = input('Do you want to pick layers (0) automatically, (1) manually, or (2) not at all?: ');
                else
                    lp = lay;
                    save(strcat(input_dir,fnames(i).name),'-append','lay')
                    kk = 2;
                end
            end
        end
    end
       
    clear lay Data Latitude Longitude dat_ mval rval
    close all
end