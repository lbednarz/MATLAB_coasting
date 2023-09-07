function [GPSweek, GPSsec] = convertToGPSTime(days_diff)
    % Calculate full weeks and the remainder in days
    GPSweek = floor(days_diff / 7);
    remainderDays = rem(days_diff, 7);
    
    % Convert remainder days into seconds
    GPSsec = remainderDays * 86400;
    
end  