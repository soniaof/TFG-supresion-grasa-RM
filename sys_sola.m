function Sys=sys_sola

    Sys.dt=10e-6;                           % system dwell time [sec]
    Sys.gamma=4.2576*1e7;                 	% gyromagnetic ratio [Hz/T]
    Sys.Gmax=45;                            % max grad amp [mT/m]
    Sys.Smax=200;                           % max slew rate [T/m/sec]
    Sys.field=1.5;                          % field strength [T]
    Sys.fat=3.5;                            % fat shift [ppm]

return