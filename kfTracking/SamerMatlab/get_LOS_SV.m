function [SV_pos_PRN, LOS_ecef_PRN] = get_LOS_SV(pseudoranges, prn,  time, x_ecef)
% TODO: "utilities" should be passed in through the function runtime

    % Read ephemeris data
    ephemeris = read_navrinex(utilities.filePath);
    prns = ephemeris(:, 1);
     
    % Convert time
    [week, sec] = convertToGPSTime(14562);
    if isempty(time)
        time = [week, sec];
    else 
        time = [week, 0];
    end

    % Get Satellite positions and biases
    [SV_pos, SV_bias] = findSVpos(ephemeris, time, 0);
    SV_pos_PRN = SV_pos(SV_pos(:,1) == prn, :);
    LOS_ecef_PRN = [];
    
    % If no SV positions, return early
    if isempty(SV_pos)
        LOS_ecef_PRN = [];
        return;
    end
    
    % Apply corrections if user position is provided
    if ~isempty(x_ecef)
        rot = [];
        if any(utilities.getConstel(prns) == 1)
            rot = utilities.omegaWGS .* (pseudoranges - x_ecef) ./ utilities.c; % GPS
        end
        
        x = SV_pos_PRN(:,2) .* cos(rot) + SV_pos_PRN(:,3) .* sin(rot);
        y = -SV_pos_PRN(:,2) .* sin(rot) + SV_pos_PRN(:,3) .* cos(rot);
        z = SV_pos_PRN(:,4);
        SV_pos_PRN = [SV_pos_PRN(:,1) x y z]; % final satellite positions

        % Calculate LOS for each satellite
        LOS_ecef = zeros(size(SV_pos_PRN,1), 3);
        for i = 1:size(SV_pos_PRN, 1)
            LOS_ecef_PRN(i,:) = (SV_pos_PRN(i,2:4) - x_ecef') / norm(SV_pos_PRN(i,2:4) - x_ecef);
        end
    end 
    return;
end

%% older version below

%function [SV_pos, LOS_ecef] = get_LOS_SV(prs, x_ecef)
%     ephemeris = read_navrinex(utilities.filePath);
%     prns = ephemeris(:, 1);
%      
%     [week, sec] = convertToGPSTime(14562);
%     time = [week, sec];
%     disp(['GPSweek: ', num2str(week), ' GPSsec: ', num2str(sec)]);
%     
%     [SV_pos,SV_bias] = findSVpos(ephemeris,time,0);
%     LOS_ecef = [];
%     
%     if ~isempty(x_ecef)
%         if ~isempty(SV_pos)
%             counter = 1;
%             rot = [];
%             if any(utilities.getConstel(prns)==1)
%                 rot = utilities.omegaWGS.*(prs - x_ecef(3+counter))./utilities.c; % GPS
%                 counter = counter + 1;
%             end
%             % not programming GLONASS yet
% 
%             x = SV_pos(:,2).*cos(rot)+SV_pos(:,3).*sin(rot);
%             y = -SV_pos(:,2).*sin(rot)+SV_pos(:,3).*cos(rot);
%             z = SV_pos(:,4);
%             SV_pos =[SV_pos(:,1) x y z]; %final satellite positions
%         else
%             SV_pos = [];
%         end
% 
%         for i=1:size(SVpos,1)
%             % LOS ECEF will have form [PRN x y z]
%             LOS_ecef(i,:) = (SV_pos(i,2:4)-x_ecef(1:3)')/norm(SV_pos(i,2:4)-x_ecef(1:3)');
%         end
%     end
