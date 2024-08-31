function [repeat,obj]=epsi_class_MetaData_alreadyexist(data_path,obj)

                    %fprintf('Initializing epsi_class with previously created Meta_Data \n')

                    % saving the previous meta_data in case it will be
                    % modified. Mostly for backward compatibility
                    MDfile=dir(fullfile(data_path,'Meta_Data.mat'));

                    str_savefilename= ...
                        fullfile(MDfile.folder,...
                                ['Meta_Data_' datestr(MDfile.datenum,"yyyymmdd") '.mat']);

                    load(fullfile(data_path,'Meta_Data'),'Meta_Data')
                    save(str_savefilename,'Meta_Data');

                    obj.Meta_Data = Meta_Data;

                    % Meta_Data includes the path to the epsi processing
                    % library, but if the data was processed on another
                    % machine, you won't have access to it.
                    % Find the epsi library on your machine and add it as process path
                    spltpath=strsplit(path,':');
                    epsilib_path=spltpath{~cellfun(@isempty, ...
                        cellfun(@(x) ...
                        strfind(x,'epsilib'),spltpath, ...
                        'UniformOutput',false))};
                    obj.Meta_Data.paths.process_library=fileparts(epsilib_path);
                    obj.Meta_Data.paths.calibration = fullfile(obj.Meta_Data.paths.process_library,'CALIBRATION','ELECTRONICS');

                    % Always redefine the data path as the current
                    % directory or the directory you input
                    obj.Meta_Data.paths.data=data_path;
                    obj.Meta_Data.paths.raw_data = fullfile(data_path,'raw');
                    obj.Meta_Data.paths.mat_data = fullfile(data_path,'mat');
                    obj.Meta_Data.paths.profiles = fullfile(data_path,'profles');
                    obj.Meta_Data.paths.figures = fullfile(data_path,'figs');

                    %ALB I do not think this is needed.
                    % obj.Meta_Data = epsiSetup_get_raw_suffix(obj.Meta_Data);


                    % Check that Meta_Data has everything you need. If it's
                    % missing something, set repeat=0 and checkMD=[]
                    if  isdir(obj.Meta_Data.paths.process_library) && ...
                            isdir(obj.Meta_Data.paths.data) && ...
                            isdir(obj.Meta_Data.paths.calibration) && ...
                            isclassfield(obj.Meta_Data.paths,'raw_data') && ...
                            isclassfield(obj.Meta_Data.PROCESS,'nfft')
                        Meta_Data = obj.Meta_Data;
                        %                         save(fullfile(obj.Meta_Data.paths.data,'Meta_Data'),'Meta_Data');

                        repeat = 0; %Stop repeating. Keep this Meta_Data.
                    else
                        repeat  =  1;
                    end

end