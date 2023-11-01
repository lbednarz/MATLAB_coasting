function R_EN = ecefToEnu(lat, long)
    R_EN = [-sin(long), cos(long), 0;
            -sin(lat)*cos(long), -sin(lat)*sin(long), cos(lat);
            cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)];
end