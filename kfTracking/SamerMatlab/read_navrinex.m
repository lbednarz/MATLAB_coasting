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

