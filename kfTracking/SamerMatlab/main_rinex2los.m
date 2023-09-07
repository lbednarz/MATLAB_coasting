function [SV_pos, LOS_ecef] = get_LOS_SV(prs, x_ecef)
    ephemeris = read_navrinex(utilities.filePath);
    prns = ephemeris(:, 1);
     
    [week, sec] = convertToGPSTime(14562);
    time = [week, sec];
    disp(['GPSweek: ', num2str(week), ' GPSsec: ', num2str(sec)]);
    
    [SV_pos,SV_bias] = findSVpos(ephemeris,time,0);
    LOS_ecef = [];
    
    if ~isempty(x_ecef)
        if ~isempty(SV_pos)
            counter = 1;
            rot = [];
            if any(utilities.getConstel(prns)==1)
                rot = utilities.omegaWGS.*(prs - x_ecef(3+counter))./utilities.c; % GPS
                counter = counter + 1;
            end
            % not programming GLONASS yet
    %         if any(utilities.getConstel(prns)==2)
    %             rot = [rot; utilities.omegaWGS.*(pr_IFF(iGLObegin:end) - x_ecef(3+counter))./utilities.c]; % GLO
    %         end
            x = SV_pos(:,2).*cos(rot)+SV_pos(:,3).*sin(rot);
            y = -SV_pos(:,2).*sin(rot)+SV_pos(:,3).*cos(rot);
            z = SV_pos(:,4);
            SV_pos =[SV_pos(:,1) x y z]; %final satellite positions
        else
            SV_pos = [];
        end

        for i=1:size(SVpos,1)
            LOS_ecef(i,:) = (SV_pos(i,2:4)-x_ecef(1:3)')/norm(SV_pos(i,2:4)-x_ecef(1:3)');
        end
    end