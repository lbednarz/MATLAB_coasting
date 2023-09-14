function plot_CN0_vs_Pc(sysc, Qc, CN0dBHz_range)

Pc_plot = zeros(1, length(CN0dBHz_range));
c = 2.99792458e8;
L5 = 1176.45e6;
lam = c/L5;

for i = 1:length(CN0dBHz_range)
    CN0dBHz = CN0dBHz_range(i);
    CN0 = 10.^(CN0dBHz/10);
    [n,~] = size(sysc.C);
    R = 1/CN0;

    [~,~,Pc] = kalman(sysc,Qc,R,[]);
    Pc_pos = 2*pi*Pc + 1;
    Pc(1,1) = lam*(1*sqrt(Pc(1,1))); 
    Pc_plot(i) = Pc(1,1);
end

    figure;
    subplot(1,1,1);
    plot(CN0dBHz_range, Pc_plot);
    xlabel('Recieved CN0 (dBHz)', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('Position error (m)', 'Interpreter', 'latex', 'FontSize', 16);
    title('Position error vs. $\frac{C}{N_0}$', 'Interpreter', 'latex');
    grid on; 
    set(gca, 'FontSize', 14);
    
    % Save as a PNG with 300 dots per inch (DPI)
    print('ssP', '-dpng', '-r300');

end



