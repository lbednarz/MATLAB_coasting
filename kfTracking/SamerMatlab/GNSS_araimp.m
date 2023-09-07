classdef GNSS_araimp
    properties(Constant)
        ftp = 'cddis.nasa.gov';  % ftp address for ephemeris
        
        thresh_CN0_L1 = 25;                                 % threshold for L1 C/N0, dBHz
        thresh_CN0_L2 = 17;                                 % threshold for L1 C/N0, dBHz
        lambda_L1 = (utilities.c/(utilities.fL1_GPS*10^6)); % L1 carrier wavelength, m
        lambda_L2 = (utilities.c/(utilities.fL2_GPS*10^6)); % L2 carrier wavelength, m
        cSlipThresh = 10;                                   % this is threshold for checking cycle ambiguity jumps when using iono free observables (that is why I picked this large), m
        mask = 7.5;                                         % satellite elevation mask, deg
    end
    properties(GetAccess=public)
        obs           % GPS/GLO observations
        ephGPS        % GPS ephemeris (rinex format)
        ephGLO        % Glonass ephemeris (igs precise SV pos and clock correction data?interpolated)
        startGPSTime  % time at first GPS observation
        endGPSTime    % time at last GPS observation
        x0_ecef       % receiver position estimate at first GPS observation
        PRN
        GLO_slots     % These are used to compute frequency of each slot
        
        % GPS/GLO observation and ephemeris file locations
        filePath     
    end
    methods
        function obj = GNSS_araimp()
            
            obj.filePath=utilities.filePath;
            
            % read RINEX 3.02 format:
            obj = obj.readRinex302();
            
            % read GPS ephemeris data (from Rinex)
            if exist(['Navfest/Day1_SpoofData/GPSeph3.19n']) ~= 0
                ephFile2 = [obj.filePath,'GPSeph2.19n'];
            else
                ephFile2 = [];
            end
            obj.ephGPS = GNSS_araimp.readEphGPS([obj.filePath,'GPSeph.19n'], ephFile2, obj.startGPSTime(1,2));
            
            % read GLO ephemeris data (from .mat obtained from IGS) (50+GLOPRN)
            obj.ephGLO = load([obj.filePath, 'GLOeph.mat']);
            % re-organizing GLO eph compatible with Ryans notation (37+GLOPRN)
            % [GPSweek TOW PRN x y z clock_corr]
            %if ~strcmp(utilities.scenario,'home')
                obj.ephGLO.SatOrbit(:,3) = obj.ephGLO.SatOrbit(:,3)-13; % PRN adjustment
            %end
        end
        function obj = readRinex302(obj)
            % This function opens RINEX 3.02 observation files.
            % Follows RINEX 3.02 standard. Reads Multiconstellation observables and
            % generates an output matrix
            
            idFile  = fopen ([obj.filePath,'obs.19O']);
            
            generalHeader = {'week', 'epoch', 'flag', 'prn'};
            
            %Initialize values
            measurementsInterval = -1;
            GLOslot = [];
            % Read header
            while (true)
                line = fgetl(idFile);                                                   %get line
                splitLine = strsplit(line);                                             %Line splited by spaces
                
                if strfind(line,'APPROX POSITION XYZ')                                  % Receiver aprox position
                    x = str2double(splitLine(2));
                    y = str2double(splitLine(3));
                    z = str2double(splitLine(4));
                    r_ecef_0 = [x;y;z];
                end
                XYZ_station=[197035.9650 -4754995.8064  4232146.4680];
                
                
                if strfind(line, 'TIME OF FIRST OBS')
                    %Read time
                    years = str2double(splitLine(2));
                    months = str2double(splitLine(3));
                    days = str2double(splitLine(4));
                    hours = str2double(splitLine(5));
                    minutes = str2double(splitLine(6));
                    seconds = str2double(splitLine(7));
                    times = [years, months, days, hours, minutes, seconds];
                    timest=real(times);
                    startTime=utc2gps(timest,0); %Transform date to seconds of week
                    
                elseif strfind(line, 'TIME OF LAST OBS')
                    %Read time
                    yeare = str2double(splitLine(2));
                    monthe = str2double(splitLine(3));
                    daye = str2double(splitLine(4));
                    houre = str2double(splitLine(5));
                    minutee = str2double(splitLine(6));
                    seconde = str2double(splitLine(7));
                    timee = [yeare, monthe, daye, houre, minutee, seconde];
                    timeed=real(timee);
                    endTime=utc2gps(timeed,0); %Transform date to seconds of week
                
                    
                elseif ~isempty(strfind(line,'SYS / # / OBS TYPES'))                    % Observation types for the different constellations (C1C, D1 and S1 only  )
                    constellation = line(1);
                    if constellation        == 'G'
                        hasGps = 1;
                    elseif constellation    == 'R'
                        hasGlonass = 1;
                    elseif constellation    == 'C'
                        hasBeidou = 1;
                    end
                    
                    nObservables = str2double(line(2:7));                                  % Number of observables
                    observables = splitLine(3:end - 7);                                     % Take the observables only (Not the text regions)
                    observables = [generalHeader, observables];                             % Append the header, as the data will be ordered like specified in obsrvables now.
                    
                    if nObservables >13 %Two line case
                        line2 = fgetl(idFile);
                        splitLine2 = strsplit(line2);
                        observables = [observables, splitLine2(2:end - 7) ];
                    end
                    
                    observablesHeader{uint8(constellation)} = observables;                  % Contains all the observables for the constellations.
                    %use constellation letter for indexation
                 
                  % CT: Gets GLO slots    
                elseif strfind(line, 'GLONASS SLOT / FRQ #')
                    for i = 1:length(splitLine)
                        try
                            if strcmp(splitLine{i}(1),'R')
                                  GLOslot =  [GLOslot;str2double(splitLine{i}(2:3))+37 str2double(splitLine{i+1})];
                            else
                               continue; 
                            end                              
                        end
                    end  
                    
                elseif strfind(line,'INTERVAL')
                    measurementsInterval=str2double(line(5:10));                       % Measurement intervals (Default 1)
                elseif strfind(line,'END OF HEADER')
                    break;                                                              % End of header loop
                end
            end
            
            if measurementsInterval == -1                                               %If itnerval not set interval = 1
                measurementsInterval = 1;
            end
            % Read body
            
            epoch = 0;                                                                  % Epoch counter
            
            while(~feof(idFile))                                                        % Until end of file
                line = fgetl(idFile);                                                   % Get line
                splitLine = strsplit(line);                                             % Split line by spaces
                
                if strfind(line, '>')                                                   % New epoch
                    epoch = epoch + 1;
                    %Read time
                    year = str2double(splitLine(2));
                    month = str2double(splitLine(3));
                    day = str2double(splitLine(4));
                    hour = str2double(splitLine(5));
                    minute = str2double(splitLine(6));
                    second = str2double(splitLine(7));
                    time = [year, month, day, hour, minute, second];
                    time=real(time);
                    [gpsWeek,tow]=utc2gps(time,0); %Transform date to seconds of week
                    
                    currentEpoch = tow;                                                 % Output
                    currentSatellites = str2double(splitLine(9));                      % Satellite information
                    currentFlag = str2double(splitLine(8));                            % flag (use/not use)
                    
                else
                    error 'Loading not correct, satellites skiped'                     % Simple check, it should never jump if using the right rinex version
                end
                
                if currentSatellites == 0
                    'No satellites in epoch'
                end
                nObs = 1;
                obj.obs{epoch} = GPSdata();
                for i = 1:real(currentSatellites)% Read the epoch satellites
                    %disp(currentSatellites);
                    line = fgetl(idFile);
                    %disp(line)
                    constellation = line(1);                                            % First character indicates de constellation
                    
                    % welcomes to GPS and GLO only
                    if ~strcmp(constellation,'G') && ~strcmp(constellation,'R')
                        continue;
                    end
                    if strcmp(constellation,'R') && ~utilities.useGLO
                        continue;
                    end
                    if strcmp(constellation,'G') && ~utilities.useGPS
                        continue;
                    end
                    prn = str2double([line(2) line(3)]);                              % Satellites PRN number
                    if strcmp(constellation,'R')
                        prn = prn+37;
                    end
                    %disp(prn);
                    nObservables =cellfun('length',observablesHeader(uint8(constellation))) - 4; %The header also includes things that are not measurements
                    %disp(nObservables)
                    measurementsPosition = (4:16:16*nObservables+4);                    %Vector containing the columns of the measurements. Each 16 columns theres a measurement
                    %disp(measurementsPosition);
                    if measurementsPosition(end) > length(line)
                        temp=ones(1,size(length(line)+1:measurementsPosition(end),2))*char(32);     % adding more whitespace
                        temp=char(temp);
                        line=[line,char(temp)];
                        measurementsPosition(end) = length(line);       %Correction of a wierd bug
                    end
                    
                    measurementsValue = zeros(1,nObservables); %Initialize vector to store data
                    for m = 1:nObservables % Number of observables in the line (Generally 3)
                        value = line(measurementsPosition(m):measurementsPosition(m+1)); % Column position of measurement. Measurements take 16 columns
                        measurementsValue(m) = str2double(value(1:14));                      % convert string to double
                    end
                    
                    measurementsValue = real(measurementsValue);                        % output of str2double is imaginary
                    if measurementsValue(1) == 0                                        % if PSR value equals 0
                        continue;                                                       % Skip line (Satellite has no information on L1)
                    end
                    
                    switch constellation
                        case 'G'
                            if (any(strfind(utilities.scenario,'SkyDel')))
                                % CT: measurementsValue is in the order of
                                % SkyDel (no **W)
                                % [C1C L1C D1C S1C C1W S1W C2W L2W D2W S2W C2L L2L D2L S2L ]
                                obj.obs{epoch}.val = [obj.obs{epoch}.val; real([gpsWeek,currentEpoch,prn,measurementsValue(1:4),nan(1,6),measurementsValue(5:end)])];
                            elseif ~(any(strfind(utilities.scenario,'2022')) || any(strfind(utilities.scenario,'2023')))
                                % CT: measurementsValue is in the order of
                                % [C1C L1C D1C S1C C1W S1W C2W L2W D2W S2W C2L L2L D2L S2L ]
                                obj.obs{epoch}.val = [obj.obs{epoch}.val; real([gpsWeek,currentEpoch,prn,measurementsValue])];
                            else
                                % SK change due to different RINEX out of
                                % new septentrio: (no C1W or S1W)
                                % [C1C L1C D1C S1C C2W L2W D2W S2W C2L L2L D2L S2L]
                                obj.obs{epoch}.val = [obj.obs{epoch}.val; real([gpsWeek,currentEpoch,prn,measurementsValue(1:4),nan,nan,measurementsValue(5:12)])];
                            end
                        case 'R'
                            % CT: measurementsValue is in the order of 
                            % [C1C L1C D1C S1C  0  0  C2C L2C D2C S2C  0  0  0 0]
                            obj.obs{epoch}.val = [obj.obs{epoch}.val; real([gpsWeek,currentEpoch,prn,measurementsValue(1:4),nan,nan,measurementsValue(5:8),nan,nan, nan, nan])];
                    end
                    
                    nObs= nObs+1;
                end
            end
            
            obj.x0_ecef = r_ecef_0;
            obj.startGPSTime = startTime;
            obj.endGPSTime = endTime;
            obj.GLO_slots = GLOslot;
            fclose(idFile);
        end
        
    end
    methods(Static)
        function [prns, pr_L1, pr_L2, cr_L1, cr_L2, isChannelWSelected, CN0_L1, CN0_L2] = getBestDualFreqObs(obs, isChannelWSelected)
            % returns an optimum dual frequency code range and carrier
            % cycles for a given instance of rinex observables.
            %
            % "isChannelWSelected" indicates prior L2 channel choice. 
            % -1 = That PRN did not exist in the previous
            % epoch, so we are free to set an optimum channel for L2
            % signal. 
            %
            % 0 = That PRN did exist, and L channel was selected, so the
            % algorithm keeps picking L2 signal from W channel 
            %
            % 1 = That PRN did exist, and W channel was selected, so the
            % algorithm keeps picking L2 signal from W channel
            
             prns = obs.PRN;
             % Use C channel for selecting L1 signals
             pr_L1 = obs.C1C;
             cr_L1 = obs.L1C;
             CN0_L1 = obs.S1C;
             
             % Select best C/N0 out of available L2 signals through W and L
             % channels 
             %
             % W channel
             pr_L2_W = obs.C2W;
             cr_L2_W = obs.L2W;
             CN0_L2_W = obs.S2W;
             % L channel
             pr_L2_L = obs.C2L;
             cr_L2_L = obs.L2L;
             CN0_L2_L = obs.S2L;
             %
             % replace NaN's with zeros so it is never selected
             CN0_L2_W(isnan(CN0_L2_W)) = 0;
             pr_L2_W(isnan(pr_L2_W)) = 0;
             cr_L2_W(isnan(cr_L2_W)) = 0;
             CN0_L2_L(isnan(CN0_L2_L)) = 0;
             pr_L2_L(isnan(pr_L2_L)) = 0;
             cr_L2_L(isnan(cr_L2_L)) = 0;
             cr_L1(isnan(cr_L1)) = 0;
             pr_L1(isnan(pr_L1)) = 0;
             
             %
             % select highest CN0 L2 signals between W and L channels
             %
             % Below array is to keep which channel is picked for each PRN:
             % e.g. 1=W 0=L. As long as maintaining lock to a specific PRN,
             % this L2 choice doesnot change.
             for i = 1:length(isChannelWSelected)
                 if isChannelWSelected(i) < 0 % newly joined PRN
                     isChannelWSelected(i) = CN0_L2_W(i) >= CN0_L2_L(i);
                 end
             end   
             pr_L2 = diag(isChannelWSelected)*pr_L2_W + diag(~isChannelWSelected)*pr_L2_L;
             cr_L2 = diag(isChannelWSelected)*cr_L2_W + diag(~isChannelWSelected)*cr_L2_L;
             CN0_L2 = diag(isChannelWSelected)*CN0_L2_W + diag(~isChannelWSelected)*CN0_L2_L;
                 
             if any(isnan(cr_L1))
                 disp('');
             end
             % eliminate L1 and L2 signals with no data (NaN)  or data with
             % zero values
             [indNaN,~] = find( (~isnan(pr_L1) .* ~(pr_L1==0)).* ...
                                (~isnan(cr_L1) .* ~(cr_L1==0)).* ...
                                (~isnan(pr_L2) .* ~(pr_L2==0)).* ...
                                (~isnan(cr_L2) .* ~(cr_L2==0)) == 0);             
             prns(indNaN) = []; 
             pr_L1(indNaN) = [];
             pr_L2(indNaN) = [];
             cr_L1(indNaN) = [];
             cr_L2(indNaN) = [];
             isChannelWSelected(indNaN) = [];
             CN0_L1(indNaN) = [];
             CN0_L2(indNaN) = [];
             
             % eliminate L1 and L2 signals with low CN0
             temp3 = (CN0_L1 > GNSS_araimp.thresh_CN0_L1) & ...
                     (CN0_L2 > GNSS_araimp.thresh_CN0_L2);
             [indToRemove,~] = find(temp3==0);
             %
             prns(indToRemove) = []; 
             pr_L1(indToRemove) = [];
             pr_L2(indToRemove) = [];
             cr_L1(indToRemove) = [];
             cr_L2(indToRemove) = [];
             isChannelWSelected(indToRemove) = [];
             CN0_L1(indToRemove) = [];
             CN0_L2(indToRemove) = [];
        end     
        function obs_IF = getIonoFreeObs(obs_f1, obs_f2, f1, f2)
            % Computes iono free combination for given two observables
            % (obs1, obs2) at different frequencies (f1, f2)
            
            % iono-free combination factors
            alpha = f1^2/(f1^2-f2^2);
            beta = f2^2/(f1^2-f2^2);
           
            obs_IF = alpha*obs_f1-beta*obs_f2;

        end
        
        function [pr_IF, cr_IF] = getIonoFreeObs_GLO(pr_L1, pr_L2, cr_L1, cr_L2, clockOffset, GLOslot)
            
            if ~isempty(pr_L1)
                clockOffset(:,1) = []; 
                for i = 1:length(pr_L1)  
                    GLOslot(i) = GLOslot(i);
                    % Frequency correction
                    L1freq_nom = utilities.fL1_0_GLO + GLOslot(i)*utilities.dfL1;        %nominal carrier frequency (Hz)
                    L2freq_nom = utilities.fL2_0_GLO + GLOslot(i)*utilities.dfL2; 

                    % Frequency correction is ignored because GLONASS satellite
                    % frequency correction (s/s) is not available in IGS, it is
                    % available only in NovAtel's message format
                    gamma = 0;
                    L1freq = L1freq_nom*(gamma+1);            %true carrier frequency (Hz)
                    L2freq = L2freq_nom*(gamma+1);
                    L1_lambda = utilities.c/(L1freq*1E+6);    %carrier wavelength (m)
                    L2_lambda = utilities.c/(L2freq*1E+6);            
                    k = (L1freq/L2freq)^2;                      %ionosphere-free factor

                    % Pseudoranges and carrier phases corrected for satellite clock bias
                    pr_L1(i) = pr_L1(i);%+utilities.c*(clockOffset(i));              %L1 code (m)
                    cr_L1(i) = cr_L1(i);%+utilities.c*(clockOffset(i))/L1_lambda;    %L1 carrier (cycles)
                    pr_L2(i) = pr_L2(i);%+utilities.c*(clockOffset(i));              %L2 code (m)
                    cr_L2(i) = cr_L2(i);%+utilities.c*(clockOffset(i))/L2_lambda;    %L2 carrier (cycles)


                    pr_IF(i,1) = (k/(k-1))*(pr_L1(i)-(1/k)*pr_L2(i));                       % code, m
                    cr_IF(i,1) = (k/(k-1))*(L1_lambda*cr_L1(i)-L2_lambda*(1/k)*cr_L2(i));   % carr, m  
                end
            else
                pr_IF = [];
                cr_IF = [];
            end
        end
        function EPH_t = readEphGPS(ephDay1, ephDay2, first_sec)
            % read ephemeris
            if isempty(ephDay2)
                EPH_t = GNSS_araimp.read_navrinex(ephDay1);
            else
                EPH_t=[GNSS_araimp.read_navrinex(ephDay1) ;GNSS_araimp.read_navrinex(ephDay2)];
            end
        end
        function ephemeris = read_navrinex(filename)
            fid = fopen(filename);
            if fid == -1
                errordlg(['The file ''' filename ''' does not exist.']);
                return;
            end
            % skip through header
            current_line = fgetl(fid);
            if num2str(current_line(5:12))>3
                RnxVer=str2num(current_line(5:12));
            end
            end_of_header = 0;
            while end_of_header == 0
                current_line = fgetl(fid);
                if strfind(current_line,'END OF HEADER')
                    end_of_header=1;
                end
            end
            j = 0;
            while feof(fid) ~= 1
                j = j+1;
                
                current_line = fgetl(fid);
                % parse epoch line (ignores SV clock bias, drift, and drift rate)
                if RnxVer>3
                    current_line=current_line(2:end);
                    [PRN, Y, M, D, H, min, sec,af0,af1,af2] = GNSS_araimp.parsef(current_line, {'I2' 'I5' 'I3' 'I3' 'I3' 'I3' ...
                    'F5.1','D19.12','D19.12','D19.12'});%GNSS_araimp.parsef304(current_line);
                else
                    [PRN, Y, M, D, H, min, sec,af0,af1,af2] = GNSS_araimp.parsef(current_line, {'I2' 'I3' 'I3' 'I3' 'I3' 'I3' ...
                    'F5.1','D19.12','D19.12','D19.12'});
                end
                % Broadcast orbit line 1
                current_line = fgetl(fid);
                [IODE Crs delta_n M0] = GNSS_araimp.parsef(current_line, {'D22.12' 'D19.12' 'D19.12' 'D19.12'});
                % Broadcast orbit line 2
                current_line = fgetl(fid);
                [Cuc e Cus sqrtA] = GNSS_araimp.parsef(current_line, {'D22.12' 'D19.12' 'D19.12' 'D19.12'});
                % Broadcast orbit line 3
                current_line = fgetl(fid);
                [toe Cic OMEGA Cis] = GNSS_araimp.parsef(current_line, {'D22.12' 'D19.12' 'D19.12' 'D19.12' 'D19.12'});
                % Broadcast orbit line 4
                current_line = fgetl(fid);
                [i0 Crc omega OMEGA_dot] = GNSS_araimp.parsef(current_line, {'D22.12' 'D19.12' 'D19.12' 'D19.12' 'D19.12'});
                % Broadcast orbit line 5
                current_line = fgetl(fid);
                [i_dot L2_codes GPS_wk L2_dataflag ] = GNSS_araimp.parsef(current_line, {'D22.12' 'D19.12' 'D19.12' 'D19.12' 'D19.12'});
                % Broadcast orbit line 6
                current_line = fgetl(fid);
                [URA SV_health TGD IODC] = GNSS_araimp.parsef(current_line, {'D22.12' 'D19.12' 'D19.12' 'D19.12' 'D19.12'});
                % Broadcast orbit line 7
                current_line = fgetl(fid);
                [msg_trans_t fit_int ] = GNSS_araimp.parsef(current_line, {'D22.12' 'D19.12' 'D19.12'});
                
                varargin=[Y, M, D, H, min, sec];
                [gps_week, toc] = GNSS_araimp.cal2gpstime(varargin);
                % ephemeris(j,:) = [PRN, M0, delta_n, e, sqrtA, OMEGA, i0, omega, OMEGA_dot, i_dot, Cuc, Cus, Crc, Crs, Cic, Cis, toe, IODE, GPS_wk,0, af0,af1,af2,0, 0, toc, TGD];
                
                % CT: This is in the same order in Ryan's code: Zeros are
                % those Ryan never uses in finding satellite position and
                % they are not even available in NASA's website.
                ephemeris(j,:) = [PRN, 0, 0, IODE, 0, GPS_wk, 0, toe, sqrtA^2, ...
                    delta_n, M0, e, omega, Cuc, Cus, Crc, Crs, ...
                    Cic, Cis, i0, i_dot, OMEGA, OMEGA_dot, 0, ...
                    toc, TGD, af0, af1, af2, 0, 0, URA]';
            end
        end
        function [gps_week, gps_seconds] = cal2gpstime(varargin)
            % Unpack
            if nargin == 1
                cal_time = varargin{1};
                year = cal_time(:,1);
                month = cal_time(:,2);
                day = cal_time(:,3);
                hour = cal_time(:,4);
                min = cal_time(:,5);
                sec = cal_time(:,6);
                clear cal_time
            else
                year = varargin{1};
                month = varargin{2};
                day = varargin{3};
                hour = varargin{4};
                min = varargin{5};
                sec = varargin{6};
            end
            % Seconds in one week
            secs_per_week = 604800;
            % Converts the two digit year to a four digit year.
            % Two digit year represents a year in the range 1980-2079.
            if (year >= 80 & year <= 99)
                year = 1900 + year;
            end
            if (year >= 0 & year <= 79)
                year = 2000 + year;
            end
            % Calculates the 'm' term used below from the given calendar month.
            if (month <= 2)
                y = year - 1;
                m = month + 12;
            end
            if (month > 2)
                y = year;
                m = month;
            end
            % Computes the Julian date corresponding to the given calendar date.
            JD = floor( (365.25 * y) ) + floor( (30.6001 * (m+1)) ) + ...
                day + ( (hour + min / 60 + sec / 3600) / 24 ) + 1720981.5;
            % Computes the GPS week corresponding to the given calendar date.
            gps_week = floor( (JD - 2444244.5) / 7 );
            % Computes the GPS seconds corresponding to the given calendar date.
            gps_seconds=round(((((JD-2444244.5)/7)-gps_week)*secs_per_week)/0.5)*0.5;
        end
        function varargout = parsef304(input)
            %GNSS_araimp.parsef   parse string value using FORTRAN formatting codes
            %   [val1,val2,...valn] = GNSS_araimp.parsef(input, format)
            %   input is string input value
            %   format is cell array of format codes
            if isletter(input(1))
                %input=input(2:end);
                input = input(2:end);
            end
            varargout = num2cell(str2num(input));
        end
        function varargout = parsef(input, format)
            %GNSS_araimp.parsef   parse string value using FORTRAN formatting codes
            %   [val1,val2,...valn] = GNSS_araimp.parsef(input, format)
            %   input is string input value
            %   format is cell array of format codes
            global input_
            input_ = input;
            %if isletter(input(1))
            %    input_ = input_(2:end);
            %end
            varargout = GNSS_araimp.getvals(1, format, 1);
            clear global input_
        end
        function [output, idx] = getvals(idx, format, reps)
            global input_
            count = 1;
            output = {};
            for k = 1:reps
                odx = 1;
                for i = 1:length(format)
                    fmt = format{i};
                    switch class(fmt)
                        case 'double'
                            count = fmt;
                        case 'char'
                            type = fmt(1);
                            if type == 'X'
                                idx = idx+count;
                            else
                                [len,cnt] = sscanf(fmt,'%*c%d',1);
                                if cnt ~= 1
                                    error(['Invalid format specifier: ''',fmt,'''']);
                                end
                                switch type
                                    case {'I','i'}
                                        for j = 1:count
                                            [val,cnt] = sscanf(input_(idx:min(idx+len-1,end)),'%d',1);
                                            if cnt == 1
                                                output{odx}(j,k) = val;
                                            else
                                                output{odx}(j,k) = NaN;
                                            end
                                            idx = idx+len;
                                        end
                                    case {'F','f'}
                                        for j = 1:count
                                            [val,cnt] = sscanf(input_(idx:min(idx+len-1,end)),'%f',1);
                                            if cnt == 1
                                                output{odx}(j,k) = val;
                                            else
                                                output{odx}(j,k) = NaN;
                                            end
                                            idx = idx+len;
                                        end
                                    case {'E','D','G'}
                                        for j = 1:count
                                            [val,cnt] = sscanf(input_(idx:min(idx+len-1,end)),'%f%*1[DdEe]%f',2);
                                            if cnt == 2
                                                output{odx}(j,k) = val(1) * 10^val(2); %#ok<AGROW>
                                            elseif cnt == 1
                                                output{odx}(j,k) = val;
                                            else
                                                output{odx}(j,k) = NaN;
                                            end
                                            idx = idx+len;
                                        end
                                    case 'A'
                                        for j = 1:count
                                            output{odx}{j,k} = input_(idx:min(idx+len-1,end));
                                            idx = idx+len;
                                        end
                                    otherwise
                                        error(['Invalid format specifier: ''',fmt,'''']);
                                end
                                odx = odx+1;
                            end
                            count = 1;
                        case 'cell'
                            [val, idx] = GNSS_araimp.getvals(idx, fmt, count);
                            if length(val) == 1
                                output(odx) = val;
                            else
                                output{odx} = val;
                            end
                            odx = odx+1;
                            count = 1;
                    end
                end
            end
        end
        function joinedData = innerJoin( data1, data2, cols1, cols2)
            %UNTITLED4 Summary of this function goes here
            %   Detailed explanation goes here
            
            [~,r1,r2] = intersect( data1(:,cols1) , data2(:,cols2) , 'rows');
            joinedData = [data1(r1,:) data2(r2,:)];
            
        end
        function eph_PRN = getInstantEph(prns, tow, ephTable)
            
            eph_PRN = [];
            for i = 1:length(prns)
                 eph_PRNi = ephTable(find(ephTable(:,1)==prns(i)),:);
                 dt = tow - eph_PRNi(:,8);
                 if ~isempty(min(dt(find(dt>0))))
                     tow2read = tow-min(dt(find(dt>0)));
                 else
                     % TODO: CT: this part is approximation in case an ephemeris data
                     % from previous days ephemeris file is needed and it
                     % is not available. It picks the ephemeris data of the
                     % nearest possible future time:
                     tow2read = tow-dt(1);  
                 end
                 
                 eph_PRNi = eph_PRNi(find(eph_PRNi(:,8)==tow2read),:);
                 eph_PRN = [eph_PRN;eph_PRNi];
            end
            
        end
        function [SVpos, clock] = getGLOSVPos(prns, tow, ephTable)
            % SK DEBUG: This here is faster that the commented portion
            % below.  NOTE: 18 leap sec. are hard coded.
            tempEph=ephTable.SatOrbit;
            tempEph=tempEph(tempEph(:,2)*20-18*20==round(tow(1)*20),:);
            [~,ia,~]=intersect(tempEph(:,3),prns);
            tempEph=tempEph(ia,:);
            tempEph(:,2)=tempEph(:,2)-18;
            SVpos = tempEph(:,3:6);
            clock = [tempEph(:,3) tempEph(:,7)];
            %{
            for i = 1:length(prns)
                 eph_PRNi = ephTable.SatOrbit(find(ephTable.SatOrbit(:,3)==prns(i)),:);
                 eph_PRNi(:,2) = eph_PRNi(:,2)-18; % TODO : leap sec is hard coded
                 dt = tow(i) - eph_PRNi(:,2);
                 if ~isempty(min(dt(find(dt>0))))
                     tow2read = tow(i)-min(dt(find(dt>0)));
                 else
                     % TODO: CT: this part is approximation in case an ephemeris data
                     % from previous days ephemeris file is needed and it
                     % is not available. It picks the ephemeris data of the
                     % nearest possible future time:
                     tow2read = tow(i)-dt(1);  
                 end
                 
                 SVpos = [SVpos; prns(i) eph_PRNi(find(eph_PRNi(:,2)==tow2read),4:6)];
                 clock = [clock; prns(i) eph_PRNi(find(eph_PRNi(:,2)==tow2read),7)];
            end
            %}
        end
    end
end