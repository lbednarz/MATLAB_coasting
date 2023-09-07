classdef utilities
    properties (Constant)
        % Constants
        c = 2.99792458e8; % Speed of light in m/s
        mu = 3.986004418e14; % Earth gravitational parameter (m^3/s^2)
        omegaWGS = 7.2921151467e-5; % WGS84 rotation rate of Earth (rad/s)
        
        % GNSS parameters
        fL1_GPS = 1575.42; % MHz
        fL2_GPS = 1227.60; % MHz
        useGLO = false; % Boolean value: true or false
        useGPS = true; % Boolean value: true or false
        scenario = 'test'; % Define the value here (it seems to be a string, e.g., 'SkyDel2023')
        fL1_0_GLO = 1602; % MHz center frequency for GLO L1
        dfL1 = 562.5e3; % FDMA frequency step for L1
        fL2_0_GLO = 1246; % MHz center frequency for GLO L2
        dfL2 = 562.5e3; % FDMA frequency step for L2
        
        % File path
        filePath = 'C:\Users\logan\Documents\Repos\GNSS_SDR-master\brdc3110.19n'; % Define the value here
    end

    methods (Static)
        function lambda = getLambda(frequency)
            lambda = utilities.c / (frequency * 10^6);
        end
        function constellation = getConstel(prns)
            % This function returns the constellation based on PRNs
            % 1 - GPS, 2 - GLONASS, etc...
        
            % Example PRN ranges (adjust as needed):
            GPS_PRNs = 1:32;
            GLONASS_PRNs = 65:96;
            
            if any(ismember(prns, GPS_PRNs))
                constellation = 1; % GPS
            elseif any(ismember(prns, GLONASS_PRNs))
                constellation = 2; % GLONASS
            else
                constellation = 0; % Unknown or another constellation
            end
        end
    end
end
