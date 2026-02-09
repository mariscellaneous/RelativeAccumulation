function rad_nc2mat(date,ROOTPATH)

    % read a directory that contains one to several OIB snow radar netcdf files
    % and convert them to Matlab matrix files, similar to the ones available on
    % the CReSIS website

    % set directory where the .nc files are stored

    root_dir = strcat(ROOTPATH,date,'/');
    input_dir = strcat(root_dir,'Originals/');
    output_dir = strcat(root_dir,'mat_files');

    fnames = dir(strcat(input_dir,'*.nc'));
    % if the output folder does not exist, create it!
    if exist(output_dir,'dir') == 0 
        mkdir(output_dir);
    end

    % now loop over each file, load the variables of interest, and save as a
    % matrix file in the output directory
     nm = size(fnames,1);
     v_names = {'fasttime','lat','lon','altitude','roll','pitch','heading',...
         'time','amplitude','Elevation_Correction','Truncate_Bins',...
         'Truncate_Mean','Truncate_Median','Truncate_Std_Dev'};   %snow radar
     nv_names = {'Time','Latitude','Longitude','Elevation','Roll','Pitch',...
         'Heading','GPS_Time','Data','Elevation_Correction','Truncate_Bins',...
         'Truncate_Mean','Truncate_Median','Truncate_Std_Dev'};   % snow radar

    %  v_names = {'fasttime','lat','lon','altitude','roll','pitch','heading',...
    %      'time','amplitude'};   % accumulation radar
    %  nv_names = {'Time','Latitude','Longitude','Elevation','Roll','Pitch',...
    %      'Heading','GPS_Time','Data'};  % accumulation radar
     for i = 1:nm
         fname = strcat(input_dir,fnames(i).name);
         % get information about the netcdf file
         outfilename = strcat(output_dir,'/',fnames(i).name(1:end-2),'mat');
         if exist(outfilename,'file') ~= 0
             fname;
             ' Already exists';
          
             continue
         end
         
         try 
            finfo = ncinfo(fname);
         catch
            'Broken file';
            fname;
            continue
         end
         % get a list of variable names to find the ones we are interested in!
         varNames = {finfo.Variables.Name};
         % loop over the variables we want, match their names to the original
         for k = 1:numel(v_names)
             dimMatch = find(strcmpi(varNames,v_names(k)) == 1);
             if isempty(dimMatch) == 0
                 dat = ncread(fname,char(v_names(k)));
                 eval([char(nv_names(k)),' = dat;']);
                 clear dat
             end
         end
         % save only the variables we want in the output dir
         save(strcat(output_dir,'/',fnames(i).name(1:end-2),'mat'),nv_names{:});
     end

     clear all
     close all
 
end