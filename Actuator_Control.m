%% ========================================================================
%  Automated Cyclic Stretch Control Script for VCA Stage
%  Author: Gia Kang
%  Institution: CTE Lab, Carleton University, Canada
%  Date: Oct 2025
%  Version: 1.0

%  Description:
%  This MATLAB script controls a voice-coil actuator (VCA) stage 
%  via serial communication to perform cyclic uniaxial stretching 
%  of PDMS cell culture wells. The program applies a predefined number 
%  of loading/unloading cycles with user-specified displacement, 
%  hold duration, and interval time between cycles. Written and tested for
%  VCA model Zaber X-DMQ12L-AE55D12.
%
%  Key Parameters:
%    numCycles    - Total number of loading cycles (1 cycle = loading + unloading)
%    intervalTime - Pause duration (in seconds) between consecutive cycles
%    holdTime     - Duration (in seconds) to hold at peak strain before unloading
%    relMove      - Total relative displacement (in millimeters) applied per cycle
%
%  Usage Notes:
%  1. Manually “home” the actuator in Zaber Launcher before running the script to ensure 
%     accurate positioning.
%  2. Verify that the COM port number (e.g., 'COM3') matches the connected device.
%  3. Adjust relMove and timing parameters according to desired strain 
%     magnitude and rate calibration for your PDMS sample.
%  4. Instantaneous/single-cycle tests may be performed outside the incubator, 
%     but minimize the duration cells spend outside controlled conditions.
%
%  Dependencies:
%    - Zaber Motion Library for MATLAB
%
%% ========================================================================

%% ----------------------- USER-DEFINED PARAMETERS -------------------------
numCycles = 10;       % Total number of loading cycles (1 cycle = loading + unloading)
intervalTime = 30;    % Duration (in seconds) between consecutive cycles
holdTime = 0;         % Hold time (in seconds) at maximum stretch before unloading
relMove = -3;         % Total relative displacement (in millimeters) per cycle
comPort = 'COM3';     % COM port of connected device

%% ----------------------- INITIALIZE CONNECTION --------------------------
% Open serial communication with the VCA device
connection = Connection.openSerialPort(comPort); 

try
    % Enable device alerts (optional, for communication safety)
    connection.enableAlerts();

    % Detect connected devices
    deviceList = connection.detectDevices();
    device = deviceList(1);   % Select the first detected device

    % Recommended: Home the stage manually in Zaber Launcher before running the script.
    % Uncomment if automatic homing is desired.
    % device.getAxis(1).home();  
    % Uncomment if required to move the actuator to an absolute position.
    % device.getAxis(1).moveAbsolute(9.5, Units.LENGTH_MILLIMETRES); 

    %% ----------------------- EXECUTE LOADING ---------------------
    for i = 1:numCycles
        fprintf("Starting experiment\n");
        fprintf("Cycle %d\n", i);

        % Move stage to apply tensile displacement (relative movement)
        device.getAxis(1).moveRelative(relMove, Units.LENGTH_MILLIMETRES);

        % Hold at peak strain (if holdTime > 0)
        if holdTime > 0
            pause(holdTime);
        else
        end

        % Return stage to the starting position (reverse displacement)
        device.getAxis(1).moveRelative(-relMove, Units.LENGTH_MILLIMETRES);

        % Pause between cycles to allow recovery or cell relaxation
        pause(intervalTime);
    end

    % Close the serial connection after completion
    connection.close();

catch exception
    % Ensure connection closes even if an error occurs
    connection.close();
    rethrow(exception);
end

%% ----------------------- FINAL MESSAGE ----------------------------------
fprintf("Experiment completed!\n");
