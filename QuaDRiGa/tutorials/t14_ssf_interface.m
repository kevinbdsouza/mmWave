%% How to manually set SSPs in QuaDRiGa (Satellite Scenario)
%
% Sometimes it is necessary to adjust the small-scale fading directly. This tutorial creates a
% simple satellite scenario. Then, a list of large-scale parameters (LSPs) and their corresponding
% small-scale parameters (SSPs) is given. It is also demonstrated how those parameters are adjusted
% and how they influence the resulting channel coefficients. 

%% Setting general parameters
% We set up some basic parameters such as center frequency and sample density.

set(0,'defaultTextFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontName','Times')                    % Default Font Type
set(0,'defaultTextFontName','Times')                    % Default Font Type
set(0,'defaultFigurePaperPositionMode','auto')          % Default Plot position
set(0,'DefaultFigurePaperType','<custom>')              % Default Paper Type

s = qd_simulation_parameters;                           % Basic simulation parameters
s.center_frequency      = 2.185e9;                      % Center Frequency
s.sample_density        = 4;                            % 4 samples / half wave length
s.show_progress_bars    = 0;                            % Disable progress bars

%% Defining the layout
% We set up our array antennas for the transmitter at the satellite and the receiver. We use
% synthetic patch antennas for this case. Two elements are crossed by an angle of 90 degree. The
% signal is then split and fed with a 90 degree phase shift to both elements generating RHCP and
% LHCP signals. We further create a track with one segment. We then assign LSPs and automatically
% calculate SSPs. Those parameters can then be edited. 

sat_el      = 28.4;                                     % Satellite elevation angle
sat_az      = 161.6;                                    % Satellite azimuth angle (South = 180 deg)
rx_latitude = 51;                                       % Latitude of the Rx

l = qd_layout( s );                                     % Create a new layout
l.set_satellite_pos( rx_latitude, sat_el, sat_az );     % Set satellite position

a = qd_arrayant('custom',120,120,0);                    % Ppatch antenna with 120 degree opening
a.copy_element(1,2);                                    % Copy element 1 to element 2
a.rotate_pattern(90,'x',2);                             % Rotate second pattern by 90 degrees
a.coupling = 1/sqrt(2) * [1 1;1j -1j];                  % LHCP and RHCP
a.rotate_pattern( sat_el , 'y' );                       % Point satellite antenna to the receiver
a.rotate_pattern( 270-sat_az , 'z' );
l.tx_array = a;                                         % Assign sattelite antenna

b = a.copy;                                             % Create a copy for the receiver.
b.rotate_pattern(-90,'y');                              % Rotate to face sky-wards.
l.rx_array = b;                                         % Assign MT antenna

l.randomize_rx_positions( 500,0,0,0 );                  % Random Rx position
t = qd_track('linear',5);                               % Linear track, 5 m length
t.interpolate_positions( s.samples_per_meter);          % Interpolate to sample density
t.initial_position = l.rx_position;                     % Set rx position
t.scenario = 'MIMOSA_10-45_LOS';                        % Set LOS sattelite scenario
l.track = t;                                            % Assign track

l.gen_lsf_parameters;                                   % Generate LSPs

cb = l.init_builder;                                    % Generate builder objects
cb.gen_ssf_parameters;                                  % Generate SSF paraemters


%%
% The LSPs are stored in "t.par". The SSPs are stored in "cb". 
% The parameter are as follows:
% 
% * *t.par.ds* : 
%   The delay spread in [s]. This value is used to assign the powers "cb.pow" and the delays
%   "cb.taus" of each SC. You can only change the value of "t.par.ds" before calling "cb =
%   l.init_builder;".
%
% * *cb.taus* : 
%   The individual delays in [s] for each scattering cluster.
%
% * *t.par.kf* :
%   A vector of values for temporal evolution of the Ricean K-factor in [dB] for each snapshot
%   (sample point). The first value "t.par.kf(1)" is used to scale the power of the first SC. The
%   vector is applied right before applying the antenna patterns. Hence, it overwrites any changes
%   you make to the SC powers in "cb.pow". It should always be edited before scaling the individual
%   SC powers, if needed.  
%
% * *cb.pow* :
%   The power of each SC normalized to a sum-power of 1. The initial Ricean K-factor is already
%   applied. Any changes you make will be overwritten by "t.par.kf". Hence, if you edit "cb.pow",
%   you must make sure to also correct the values in"t.par.kf" and "cb.par.kf".

%% 
% You can edit power values in cb.pow as you wish. Here, we scale the first and second SC.

cb.pow(1:2) = [0.5, 0.2];                               % Set power ration of first two clusters
cb.pow = cb.pow / sum(cb.pow);                          % Make sure that the sum-power is 1
kf = 10*log10( cb.pow(1) / sum( cb.pow(2:end) ) );      % Calcualte new KF
t.par.kf = t.par.kf - ( t.par.kf(1) - kf );             % Correct KF in "t.par.kf" and "cb.par.kf"
cb.kf = 10.^( 0.1*t.par.kf(1) );

%%
% * *t.par.pg* :
%   A vector of values for total path gain in [dB] for each snapshot. This is equal to the sum-power
%   (e.g. in [dBm]) minus the transmit power (also in [dBm]). This vector is applied after the
%   KF-scaling and right before applying the antenna patterns.
%
% * *t.par.asA* :
%   Azimuth angle spread of arrival in [deg]. You can only change the value of "t.par.asA" before
%   calling "[~,~,cb] = l.get_channels;". 
%
% * *cb.AoA* :
%   Azimuth arrival angles for each SC in global spherical coordinates given in [rad]. The first
%   angle corresponds to the satellite  azimuth position. 
%
% * *t.par.esA* : 
%   Elevation angle spread of arrival in [deg]. You can only change the value of "t.par.esA" before
%   generating the builder objects.
%
% * *cb.EoA* : 
%    Elevation arrival angles for each SC in global spherical coordinates given in [rad]. The first
%    angle corresponds to the satellite elevation position. 
%       
% * *t.par.asD* :
%   Azimuth angle spread of departure in [deg]. This value is irrelevant for satellite setups where
%   the departure angle is almost the same for each SC. You can only change the value of "t.par.asD"
%   before generating the builder objects.
%
% * *cb.AoD* : 
%   Azimuth departure angles for each SC in global spherical coordinates given in [rad]. All angles
%   should point to the terminal position as seen from the satellite.
%
% * *t.par.esD* : 
%   Elevation angle spread of departure in [deg]. This value is irrelavant for satellite setups
%   where the departure angle is almost the same for each SC.You can only change the value of
%   "t.par.esA" before generating the builder objects. 
%
% * *cb.EoA* :
%   Elevation departure angles for each SC in global spherical coordinates given in [rad]. All
%   angles correspond roughly to the  negative satellite elevation position. 
%
% After making the required changes, you can combine the SC and the antenna patterns to calculate
% the channel coefficient by: 

c = cb.get_channels;

%%
% A narrow band time-series is calculated by:

h = squeeze( sum(c.coeff,3) );

%%
% Plot for the power on each MIMO sub-channel.

power = abs( reshape( h,4,[] ) ).^2;
power = 10*log10(power).';

figure
[~,dist] = t.get_length;
plot(dist,t.par.pg,'--k','Linewidth',2)
hold on
plot(dist,power)
hold off
title('Simulated Power')
xlabel('Distance from start point [m]')
ylabel('Received Power [dBm]')
axis([0,5,min( power(:) )-5,max( power(:) )+2])
grid on
legend('LSP: t.par.pg','LHCP-LHCP','LHCP-RHCP','RHCP-LHCP','RHCP-RHCP','Location','SouthEast');

