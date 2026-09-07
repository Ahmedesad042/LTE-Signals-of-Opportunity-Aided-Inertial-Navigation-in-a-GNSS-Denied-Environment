% LTE Signals-of-Opportunity Aided Inertial Navigation
% in a GNSS-Denied Environment
%
% This proof-of-concept simulation maps LTE base stations while GNSS is
% available. During GNSS denial, inter-receiver TDOAs are converted into
% base-station bearings. Multiple mapped-station bearings produce an LTE
% position fix, which is fused with the INS using an extended Kalman filter.

clear;
clc;
close all;

%% PART 1 — Configuration

% Timing and sampling
dt = 0.1;
T = 1000;
t = (0:dt:T).';
N = numel(t);

imuRate = 1/dt;
gnssRate = 1;
cellularRate = 1;

gnssOutageStart = 100;
gnssOutageEnd = T;

% Independent random seeds keep subsystems reproducible.
randomSeed.heading = 47;
randomSeed.imu = 48;
randomSeed.gnss = 49;
randomSeed.tdoa = 50;
randomSeed.lteIQ = 51;

% Sensor errors
accelerometerBias = [0.015, -0.010];       % [m/s^2]
accelerometerNoiseStd = 0.03;              % [m/s^2]

gyroscopeBias = deg2rad(0.02);             % [rad/s]
gyroscopeNoiseStd = deg2rad(0.05);         % [rad/s]

gnssPositionNoiseStd = 2.0;                % [m]
gnssVelocityNoiseStd = 0.10;               % [m/s]
headingSensorNoiseStd = deg2rad(0.1);      % [rad]

% Receiver positions: [forward, left] relative to aircraft centre [m]
receiverBody = [
     5,  0;
    -5, -4;
    -5,  4
    ];

numberOfReceivers = size(receiverBody,1);
referenceReceiverIndex = 1;

nonreferenceReceiverIndices = 1:numberOfReceivers;
nonreferenceReceiverIndices(referenceReceiverIndex) = [];

% Base-station positions are simulation truth: [East, North] [m]
trueBaseStationPosition = 1000 * [
     8.0, -0.5;
    11.0,  3.0;
     4.0,  4.0;
     6.0,  7.0;
    12.0,  8.0;
     9.5, 10.0;
    14.1,  5.0;
    13.5,  7.0
    ];

% PCI's must be between 0 and 503
baseStationCellID = [
     10;
     71;
    184;
    307;
    402;
    499;
     25;
    138
    ];

numberOfBaseStations = size(trueBaseStationPosition,1);
speedOfLight = 299792458;

% TDOA measurement settings
usePerfectTDOA = true;
tdoaSimulationNoiseStd = 0.5e-9;
tdoaAssumedNoiseStd = 0.5e-9;

% LTE PSS/SSS settings
lteFFTSize = 2048;
lteSampleRate = 30.72e6;

pssCyclicPrefixLength = 144;
sssCyclicPrefixLength = 144;

pssRoot = [25, 29, 34];
numberOfSSSCandidates = 168;

lteTimingOffsetStart = 250;
lteTimingOffsetStep = 100;
lteTrailingGuardSamples = 100;

lteTestSNRdB = 5;
pssAcceptanceThreshold = 0.50;
sssAcceptanceThreshold = 0.50;

% Mapping settings
numberOfMappingObservations = 10;

mappingOptimizationOptions = optimset( ...
    'Display','off', ...
    'TolX',1e-6, ...
    'TolFun',1e-6, ...
    'MaxIter',2000, ...
    'MaxFunEvals',5000);

% Navigation-filter settings
numberOfNavigationStates = 8;

initialHeadingUncertainty = deg2rad(1);
initialAccelerometerBiasUncertainty = 0.05;
initialGyroscopeBiasUncertainty = deg2rad(0.1);

accelerometerBiasRandomWalkStd = 1e-4;
gyroscopeBiasRandomWalkStd = deg2rad(0.001);

minimumBaseStationsForPositionFix = 3;

minimumBearingInformationRcond = 1e-8;
lteInnovationGateThreshold = 9.21;

%% PART 2 — Generate the true aircraft trajectory

position = zeros(N,2);            % [East, North] [m]
speed = zeros(N,1);               % [m/s]
heading = zeros(N,1);             % Counterclockwise from East [rad]

accelerationCommand = zeros(N,1); % Forward acceleration [m/s^2]
turnRateCommand = zeros(N,1);     % Positive means left turn [rad/s]

position(1,:) = [0,0];
speed(1) = 50;
heading(1) = 0;

% Motion commands
accelerationCommand(t >= 30  & t < 50)  = 1.0;
turnRateCommand(t >= 70      & t < 100) = deg2rad(2);
accelerationCommand(t >= 130 & t < 150) = -0.75;
turnRateCommand(t >= 170     & t < 220) = deg2rad(-1.5);
accelerationCommand(t >= 240 & t < 260) = 0.5;
turnRateCommand(t >= 240     & t < 280) = deg2rad(1);

% Propagate position, speed and heading
for k = 1:N-1

    speed(k+1) = speed(k) + accelerationCommand(k)*dt;
    heading(k+1) = heading(k) + turnRateCommand(k)*dt;

    intervalSpeed = 0.5*(speed(k) + speed(k+1));
    intervalHeading = 0.5*(heading(k) + heading(k+1));

    position(k+1,:) = position(k,:) ...
        + intervalSpeed ...
        * [cos(intervalHeading),sin(intervalHeading)]*dt;
end

% Global velocity
velocity = [
    speed.*cos(heading), ...
    speed.*sin(heading)
    ];

% Analytic global acceleration
globalAcceleration = [
    accelerationCommand.*cos(heading) ...
        - speed.*turnRateCommand.*sin(heading), ...
    accelerationCommand.*sin(heading) ...
        + speed.*turnRateCommand.*cos(heading)
    ];

truth.time = t;
truth.position = position;
truth.velocity = velocity;
truth.speed = speed;
truth.heading = heading;
truth.acceleration = globalAcceleration;
truth.turnRate = turnRateCommand;
truth.accelerationCommand = accelerationCommand;

%% PART 3 — Simulate heading, IMU and GNSS measurements

% Heading sensor
rng(randomSeed.heading,'twister');

headingSensorNoise = ...
    headingSensorNoiseStd*randn(N,1);

measuredHeading = truth.heading + headingSensorNoise;
measuredHeading = atan2(sin(measuredHeading),cos(measuredHeading));

headingSensor.time = t;
headingSensor.measurement = measuredHeading;
headingSensor.noise = headingSensorNoise;
headingSensor.noiseStd = headingSensorNoiseStd;

% Ideal IMU quantities in body coordinates
trueAccelerationBody = [
    truth.accelerationCommand, ...
    truth.speed.*truth.turnRate
    ];

trueGyroscopeZ = truth.turnRate;

% Noisy IMU measurements
rng(randomSeed.imu,'twister');

accelerometerNoise = ...
    accelerometerNoiseStd*randn(N,2);

gyroscopeNoise = ...
    gyroscopeNoiseStd*randn(N,1);

measuredAccelerationBody = ...
    trueAccelerationBody ...
    + accelerometerBias ...
    + accelerometerNoise;

measuredGyroscopeZ = ...
    trueGyroscopeZ ...
    + gyroscopeBias ...
    + gyroscopeNoise;

imu.time = t;
imu.rate = imuRate;

imu.accelerationBody = measuredAccelerationBody;
imu.gyroscopeZ = measuredGyroscopeZ;

imu.trueAccelerationBody = trueAccelerationBody;
imu.trueGyroscopeZ = trueGyroscopeZ;

imu.accelerometerBias = accelerometerBias;
imu.gyroscopeBias = gyroscopeBias;
imu.accelerometerNoise = accelerometerNoise;
imu.gyroscopeNoise = gyroscopeNoise;

% GNSS measurement times
gnssSampleStep = round(1/(gnssRate*dt));
gnssSampleIndices = (1:gnssSampleStep:N).';
gnssTime = t(gnssSampleIndices);
numberOfGnssSamples = numel(gnssTime);

trueGnssPosition = truth.position(gnssSampleIndices,:);
trueGnssVelocity = truth.velocity(gnssSampleIndices,:);

% Noisy GNSS measurements
rng(randomSeed.gnss,'twister');

gnssPositionNoise = ...
    gnssPositionNoiseStd*randn(numberOfGnssSamples,2);

gnssVelocityNoise = ...
    gnssVelocityNoiseStd*randn(numberOfGnssSamples,2);

measuredGnssPosition = ...
    trueGnssPosition + gnssPositionNoise;

measuredGnssVelocity = ...
    trueGnssVelocity + gnssVelocityNoise;

% Remove measurements during the outage
gnssAvailable = ...
    gnssTime < gnssOutageStart ...
    | gnssTime > gnssOutageEnd;

measuredGnssPosition(~gnssAvailable,:) = NaN;
measuredGnssVelocity(~gnssAvailable,:) = NaN;

gnss.time = gnssTime;
gnss.rate = gnssRate;
gnss.sampleIndices = gnssSampleIndices;
gnss.available = gnssAvailable;

gnss.position = measuredGnssPosition;
gnss.velocity = measuredGnssVelocity;

gnss.truePosition = trueGnssPosition;
gnss.trueVelocity = trueGnssVelocity;

gnss.positionNoise = gnssPositionNoise;
gnss.velocityNoise = gnssVelocityNoise;

gnss.outageStart = gnssOutageStart;
gnss.outageEnd = gnssOutageEnd;

%% PART 4 — Generate receiver geometry and TDOA measurements

% Transform receiver positions from aircraft coordinates to global coordinates
receiverGlobal = zeros(N,2,numberOfReceivers);

for k = 1:N

    psi = truth.heading(k);

    rotationBodyToGlobal = [
        cos(psi), -sin(psi);
        sin(psi),  cos(psi)
        ];

    for receiverIndex = 1:numberOfReceivers

        receiverOffsetGlobal = rotationBodyToGlobal ...
            * receiverBody(receiverIndex,:).';

        receiverGlobal(k,:,receiverIndex) = ...
            truth.position(k,:) + receiverOffsetGlobal.';
    end
end

truth.receiverBody = receiverBody;
truth.receiverGlobal = receiverGlobal;

% Cellular observation times
cellularSampleStep = round(1/(cellularRate*dt));

if abs(cellularSampleStep*dt-1/cellularRate) > 1e-12
    error('The cellular rate is not compatible with the simulation time step.');
end

cellularSampleIndices = (1:cellularSampleStep:N).';
cellularTime = t(cellularSampleIndices);
numberOfCellularSamples = numel(cellularTime);

% Validate receiver geometry
receiverBaselineBody = ...
    receiverBody(nonreferenceReceiverIndices,:) ...
    - receiverBody(referenceReceiverIndex,:);

if rank(receiverBaselineBody) < 2
    error('The receiver baselines cannot estimate a 2D arrival direction.');
end

numberOfTDOAPairs = numel(nonreferenceReceiverIndices);

% Allocate true cellular geometry
trueCellularRange = zeros( ...
    numberOfCellularSamples, ...
    numberOfReceivers, ...
    numberOfBaseStations);

truePropagationDelay = zeros( ...
    numberOfCellularSamples, ...
    numberOfReceivers, ...
    numberOfBaseStations);

trueInterReceiverTDOA = zeros( ...
    numberOfCellularSamples, ...
    numberOfTDOAPairs, ...
    numberOfBaseStations);

% Calculate exact ranges, delays and TDOAs
for cellularIndex = 1:numberOfCellularSamples

    trajectoryIndex = cellularSampleIndices(cellularIndex);

    for baseStationIndex = 1:numberOfBaseStations

        for receiverIndex = 1:numberOfReceivers

            receiverPosition = reshape( ...
                receiverGlobal(trajectoryIndex,:,receiverIndex),1,2);

            currentRange = norm( ...
                trueBaseStationPosition(baseStationIndex,:) ...
                - receiverPosition);

            trueCellularRange( ...
                cellularIndex,receiverIndex,baseStationIndex) = ...
                currentRange;

            truePropagationDelay( ...
                cellularIndex,receiverIndex,baseStationIndex) = ...
                currentRange/speedOfLight;
        end

        referenceDelay = truePropagationDelay( ...
            cellularIndex,referenceReceiverIndex,baseStationIndex);

        for pairIndex = 1:numberOfTDOAPairs

            receiverIndex = ...
                nonreferenceReceiverIndices(pairIndex);

            trueInterReceiverTDOA( ...
                cellularIndex,pairIndex,baseStationIndex) = ...
                truePropagationDelay( ...
                    cellularIndex,receiverIndex,baseStationIndex) ...
                - referenceDelay;
        end
    end
end

% Add measurement noise unless the perfect-TDOA test is enabled
if usePerfectTDOA

    interReceiverTDOANoise = ...
        zeros(size(trueInterReceiverTDOA));

    appliedTDOANoiseStd = 0;

else

    rng(randomSeed.tdoa,'twister');

    interReceiverTDOANoise = ...
        tdoaSimulationNoiseStd ...
        * randn(size(trueInterReceiverTDOA));

    appliedTDOANoiseStd = tdoaSimulationNoiseStd;
end

measuredInterReceiverTDOA = ...
    trueInterReceiverTDOA + interReceiverTDOANoise;

% Check the physical TDOA limits
receiverBaselineLength = vecnorm(receiverBaselineBody,2,2);
maximumPossibleTDOA = receiverBaselineLength/speedOfLight;

for pairIndex = 1:numberOfTDOAPairs

    pairTDOA = trueInterReceiverTDOA(:,pairIndex,:);

    if any(abs(pairTDOA(:)) ...
            > maximumPossibleTDOA(pairIndex) + 1e-15)

        error('A calculated TDOA exceeds its receiver-baseline limit.');
    end
end

% Store cellular observations
cellular.time = cellularTime;
cellular.rate = cellularRate;
cellular.sampleIndices = cellularSampleIndices;

cellular.referenceReceiverIndex = referenceReceiverIndex;
cellular.nonreferenceReceiverIndices = nonreferenceReceiverIndices;

cellular.trueRange = trueCellularRange;
cellular.truePropagationDelay = truePropagationDelay;
cellular.trueTDOA = trueInterReceiverTDOA;

cellular.measuredTDOA = measuredInterReceiverTDOA;
cellular.tdoaNoise = interReceiverTDOANoise;
cellular.appliedTDOANoiseStd = appliedTDOANoiseStd;
cellular.assumedTDOANoiseStd = tdoaAssumedNoiseStd;

cellular.maximumPossibleTDOA = maximumPossibleTDOA;

%% PART 5 — Generate LTE PSS and SSS references

% The 62 synchronization subcarriers surrounding, but excluding, DC
dcIndex = lteFFTSize/2 + 1;

synchronizationSubcarrierIndices = [
    (dcIndex + (-31:-1)).';
    (dcIndex + (1:31)).'
    ];

% Generate the three possible PSS waveforms
numberOfPSSCandidates = numel(pssRoot);

pssReferenceWaveform = zeros( ...
    lteFFTSize + pssCyclicPrefixLength, ...
    numberOfPSSCandidates);

for candidateNID2 = 0:2

    root = pssRoot(candidateNID2+1);
    pssSequence = zeros(62,1);

    n = (0:30).';
    pssSequence(1:31) = ...
        exp(-1i*pi*root.*n.*(n+1)/63);

    n = (31:61).';
    pssSequence(32:62) = ...
        exp(-1i*pi*root.*(n+1).*(n+2)/63);

    frequencyGrid = zeros(lteFFTSize,1);
    frequencyGrid(synchronizationSubcarrierIndices) = pssSequence;

    timeDomainSymbol = sqrt(lteFFTSize) ...
        * ifft(ifftshift(frequencyGrid),lteFFTSize);

    pssReferenceWaveform(:,candidateNID2+1) = [
        timeDomainSymbol(end-pssCyclicPrefixLength+1:end);
        timeDomainSymbol
        ];
end

% Generate the three fundamental SSS m-sequences
initialMSequence = [0;0;0;0;1];

xS = zeros(31,1);
xC = zeros(31,1);
xZ = zeros(31,1);

xS(1:5) = initialMSequence;
xC(1:5) = initialMSequence;
xZ(1:5) = initialMSequence;

for sequenceIndex = 0:25

    xS(sequenceIndex+6) = mod( ...
        xS(sequenceIndex+3) + xS(sequenceIndex+1),2);

    xC(sequenceIndex+6) = mod( ...
        xC(sequenceIndex+4) + xC(sequenceIndex+1),2);

    xZ(sequenceIndex+6) = mod( ...
        xZ(sequenceIndex+5) ...
        + xZ(sequenceIndex+3) ...
        + xZ(sequenceIndex+2) ...
        + xZ(sequenceIndex+1),2);
end

sTilde = 1-2*xS;
cTilde = 1-2*xC;
zTilde = 1-2*xZ;

% Precompute all 168 × 3 possible SSS sequences
nSSS = (0:30).';

sssReferenceSequence = zeros( ...
    62,numberOfSSSCandidates,numberOfPSSCandidates);

for candidateNID2 = 0:2

    c0 = cTilde(mod(nSSS+candidateNID2,31)+1);
    c1 = cTilde(mod(nSSS+candidateNID2+3,31)+1);

    for candidateNID1 = 0:numberOfSSSCandidates-1

        qPrime = floor(candidateNID1/30);

        q = floor( ...
            (candidateNID1 + qPrime*(qPrime+1)/2)/30);

        mPrime = candidateNID1 + q*(q+1)/2;

        m0 = mod(mPrime,31);
        m1 = mod(m0 + floor(mPrime/31) + 1,31);

        s0 = sTilde(mod(nSSS+m0,31)+1);
        s1 = sTilde(mod(nSSS+m1,31)+1);
        z1 = zTilde(mod(nSSS+mod(m0,8),31)+1);

        sssSequence = zeros(62,1);
        sssSequence(1:2:end) = s0.*c0;
        sssSequence(2:2:end) = s1.*c1.*z1;

        sssReferenceSequence( ...
            :,candidateNID1+1,candidateNID2+1) = ...
            sssSequence;
    end
end

%% PART 6 — Simulate LTE waveforms and detect each PCI

trueNID1FromIQ = floor(baseStationCellID/3);
trueNID2FromIQ = mod(baseStationCellID,3);

detectedCellIDFromIQ = nan(numberOfBaseStations,1);
detectedNID1FromIQ = nan(numberOfBaseStations,1);
detectedNID2FromIQ = nan(numberOfBaseStations,1);

truePssTimingOffsetFromIQ = nan(numberOfBaseStations,1);
detectedPssTimingOffsetFromIQ = nan(numberOfBaseStations,1);

pssCorrelationMetricFromIQ = nan(numberOfBaseStations,1);
sssCorrelationMetricFromIQ = nan(numberOfBaseStations,1);

pssDetectionAvailableFromIQ = false(numberOfBaseStations,1);
cellIDDetectionCorrectFromIQ = false(numberOfBaseStations,1);

sssWaveformLength = lteFFTSize + sssCyclicPrefixLength;
pssWaveformLength = lteFFTSize + pssCyclicPrefixLength;

rng(randomSeed.lteIQ,'twister');

for baseStationIndex = 1:numberOfBaseStations

    transmittedCellID = baseStationCellID(baseStationIndex);
    transmittedNID1 = trueNID1FromIQ(baseStationIndex);
    transmittedNID2 = trueNID2FromIQ(baseStationIndex);

    % Generate the SSS OFDM waveform
    transmittedSSSSequence = sssReferenceSequence( ...
        :,transmittedNID1+1,transmittedNID2+1);

    sssFrequencyGrid = zeros(lteFFTSize,1);
    sssFrequencyGrid(synchronizationSubcarrierIndices) = ...
        transmittedSSSSequence;

    sssTimeDomainSymbol = sqrt(lteFFTSize) ...
        * ifft(ifftshift(sssFrequencyGrid),lteFFTSize);

    transmittedSSSWaveform = [
        sssTimeDomainSymbol(end-sssCyclicPrefixLength+1:end);
        sssTimeDomainSymbol
        ];

    transmittedPSSWaveform = ...
        pssReferenceWaveform(:,transmittedNID2+1);

    transmittedSynchronizationWaveform = [
        transmittedSSSWaveform;
        transmittedPSSWaveform
        ];

    % Insert the synchronization waveform at an artificial sample offset
    synchronizationTimingOffset = ...
        lteTimingOffsetStart ...
        + lteTimingOffsetStep*(baseStationIndex-1);

    receivedBufferLength = ...
        synchronizationTimingOffset ...
        + length(transmittedSynchronizationWaveform) ...
        + lteTrailingGuardSamples;

    receivedSynchronizationWaveform = ...
        zeros(receivedBufferLength,1);

    synchronizationStartIndex = ...
        synchronizationTimingOffset + 1;

    synchronizationEndIndex = ...
        synchronizationStartIndex ...
        + length(transmittedSynchronizationWaveform) - 1;

    receivedSynchronizationWaveform( ...
        synchronizationStartIndex:synchronizationEndIndex) = ...
        transmittedSynchronizationWaveform;

    truePssTimingOffset = ...
        synchronizationTimingOffset + sssWaveformLength;

    % Add complex Gaussian noise
    signalPower = mean( ...
        abs(transmittedSynchronizationWaveform).^2);

    noisePower = signalPower/10^(lteTestSNRdB/10);

    complexNoise = sqrt(noisePower/2) * ( ...
        randn(receivedBufferLength,1) ...
        + 1i*randn(receivedBufferLength,1));

    receivedSynchronizationWaveform = ...
        receivedSynchronizationWaveform + complexNoise;

    % Search for the three possible PSS waveforms
    numberOfCorrelationPositions = ...
        receivedBufferLength-pssWaveformLength+1;

    candidatePeakMagnitude = zeros(numberOfPSSCandidates,1);
    candidateTimingOffset = zeros(numberOfPSSCandidates,1);

    for candidateIndex = 1:numberOfPSSCandidates

        currentReference = ...
            pssReferenceWaveform(:,candidateIndex);

        currentCorrelation = abs(conv( ...
            receivedSynchronizationWaveform, ...
            conj(flipud(currentReference)), ...
            'valid'));

        if numel(currentCorrelation) ~= numberOfCorrelationPositions
            error('Unexpected PSS correlation length.');
        end

        [candidatePeakMagnitude(candidateIndex),peakIndex] = ...
            max(currentCorrelation);

        candidateTimingOffset(candidateIndex) = peakIndex-1;
    end

    [detectedPSSPeak,detectedCandidateIndex] = ...
        max(candidatePeakMagnitude);

    detectedPssTimingOffset = ...
        candidateTimingOffset(detectedCandidateIndex);

    detectedPssStartIndex = detectedPssTimingOffset+1;
    detectedPssEndIndex = ...
        detectedPssStartIndex+pssWaveformLength-1;

    detectedPSSReference = ...
        pssReferenceWaveform(:,detectedCandidateIndex);

    receivedPSSWindow = receivedSynchronizationWaveform( ...
        detectedPssStartIndex:detectedPssEndIndex);

    pssCorrelationMetric = detectedPSSPeak ...
        / max(norm(receivedPSSWindow) ...
        * norm(detectedPSSReference),eps);

    truePssTimingOffsetFromIQ(baseStationIndex) = ...
        truePssTimingOffset;

    detectedPssTimingOffsetFromIQ(baseStationIndex) = ...
        detectedPssTimingOffset;

    pssCorrelationMetricFromIQ(baseStationIndex) = ...
        pssCorrelationMetric;

    if pssCorrelationMetric < pssAcceptanceThreshold
        continue;
    end

    pssDetectionAvailableFromIQ(baseStationIndex) = true;

    detectedNID2 = detectedCandidateIndex-1;
    detectedNID2FromIQ(baseStationIndex) = detectedNID2;

    % Locate and OFDM-demodulate the preceding SSS symbol
    detectedSssStartIndex = ...
        detectedPssStartIndex-sssWaveformLength;

    receivedSssUsefulStartIndex = ...
        detectedSssStartIndex+sssCyclicPrefixLength;

    receivedSssUsefulEndIndex = ...
        receivedSssUsefulStartIndex+lteFFTSize-1;

    if detectedSssStartIndex < 1 ...
            || receivedSssUsefulEndIndex > receivedBufferLength

        error('The detected SSS symbol lies outside the received buffer.');
    end

    receivedSssTimeDomainSymbol = ...
        receivedSynchronizationWaveform( ...
            receivedSssUsefulStartIndex:receivedSssUsefulEndIndex);

    receivedSssFrequencyGrid = fftshift( ...
        fft(receivedSssTimeDomainSymbol,lteFFTSize) ...
        / sqrt(lteFFTSize));

    receivedSSSSequence = ...
        receivedSssFrequencyGrid( ...
            synchronizationSubcarrierIndices);

    % Compare the received SSS against all 168 NID1 candidates
    candidateSSSSequences = ...
        sssReferenceSequence(:,:,detectedNID2+1);

    candidateSSSNorm = sqrt( ...
        sum(abs(candidateSSSSequences).^2,1)).';

    receivedSSSNorm = norm(receivedSSSSequence);

    sssCorrelationMetric = abs( ...
        candidateSSSSequences' * receivedSSSSequence) ...
        ./ max(candidateSSSNorm*receivedSSSNorm,eps);

    [detectedSSSPeak,detectedNID1Index] = ...
        max(sssCorrelationMetric);

    detectedNID1 = detectedNID1Index-1;
    detectedCellID = 3*detectedNID1+detectedNID2;

    detectedNID1FromIQ(baseStationIndex) = detectedNID1;
    detectedCellIDFromIQ(baseStationIndex) = detectedCellID;
    sssCorrelationMetricFromIQ(baseStationIndex) = detectedSSSPeak;

    cellIDDetectionCorrectFromIQ(baseStationIndex) = ...
        detectedCellID == transmittedCellID;
end

%% PART 7 — Associate detected cells with measurement slots

numberOfMeasurementSlots = ...
    size(cellular.measuredTDOA,3);

associationStatusCode.noDetection = uint8(0);
associationStatusCode.associated = uint8(1);
associationStatusCode.unknownCellID = uint8(2);
associationStatusCode.ambiguousCellID = uint8(3);
associationStatusCode.rejected = uint8(4);

% One IQ detection was performed per station and is repeated at every epoch.
observedCellID = repmat( ...
    detectedCellIDFromIQ.', ...
    numberOfCellularSamples,1);

cellDetectionAvailable = repmat( ...
    pssDetectionAvailableFromIQ.', ...
    numberOfCellularSamples,1);

singleSnapshotIdentificationAccepted = ...
    pssDetectionAvailableFromIQ ...
    & isfinite(detectedCellIDFromIQ) ...
    & sssCorrelationMetricFromIQ >= sssAcceptanceThreshold;

cellIdentificationAccepted = repmat( ...
    singleSnapshotIdentificationAccepted.', ...
    numberOfCellularSamples,1);

associatedBaseStationIndex = nan( ...
    numberOfCellularSamples,numberOfMeasurementSlots);

associationStatus = repmat( ...
    associationStatusCode.noDetection, ...
    numberOfCellularSamples,numberOfMeasurementSlots);

associationCorrect = false( ...
    numberOfCellularSamples,numberOfMeasurementSlots);

% TDOAs reorganized by associated catalogue entry
associatedTDOA = nan( ...
    numberOfCellularSamples, ...
    numberOfTDOAPairs, ...
    numberOfBaseStations);

% Simulation truth used only to evaluate association
trueSourceBaseStationIndex = repmat( ...
    1:numberOfBaseStations, ...
    numberOfCellularSamples,1);

for cellularIndex = 1:numberOfCellularSamples

    for measurementIndex = 1:numberOfMeasurementSlots

        if ~cellDetectionAvailable(cellularIndex,measurementIndex)
            continue;
        end

        if ~cellIdentificationAccepted(cellularIndex,measurementIndex)

            associationStatus(cellularIndex,measurementIndex) = ...
                associationStatusCode.rejected;

            continue;
        end

        currentCellID = ...
            observedCellID(cellularIndex,measurementIndex);

        matchingCatalogueIndices = find( ...
            baseStationCellID == currentCellID);

        if isempty(matchingCatalogueIndices)

            associationStatus(cellularIndex,measurementIndex) = ...
                associationStatusCode.unknownCellID;

        elseif numel(matchingCatalogueIndices) > 1

            associationStatus(cellularIndex,measurementIndex) = ...
                associationStatusCode.ambiguousCellID;

        else

            matchedBaseStationIndex = ...
                matchingCatalogueIndices(1);

            associatedBaseStationIndex( ...
                cellularIndex,measurementIndex) = ...
                matchedBaseStationIndex;

            associationStatus(cellularIndex,measurementIndex) = ...
                associationStatusCode.associated;

            associatedTDOA( ...
                cellularIndex,:,matchedBaseStationIndex) = ...
                cellular.measuredTDOA( ...
                    cellularIndex,:,measurementIndex);

            associationCorrect(cellularIndex,measurementIndex) = ...
                matchedBaseStationIndex ...
                == trueSourceBaseStationIndex( ...
                    cellularIndex,measurementIndex);
        end
    end
end

associationSuccessful = ...
    associationStatus == associationStatusCode.associated;

cellular.numberOfMeasurementSlots = numberOfMeasurementSlots;
cellular.observedCellID = observedCellID;
cellular.cellDetectionAvailable = cellDetectionAvailable;
cellular.cellIdentificationAccepted = cellIdentificationAccepted;

cellular.associatedBaseStationIndex = associatedBaseStationIndex;
cellular.associatedTDOA = associatedTDOA;
cellular.associationStatus = associationStatus;
cellular.associationStatusCode = associationStatusCode;
cellular.associationSuccessful = associationSuccessful;

cellular.trueSourceBaseStationIndex = trueSourceBaseStationIndex;
cellular.associationCorrect = associationCorrect;

%% PART 8 — Select observations for base-station mapping

% Match cellular observations to GNSS samples using trajectory indices
[hasMatchingGnssSample,gnssIndexAtCellularTime] = ismember( ...
    cellular.sampleIndices,gnss.sampleIndices);

gnssAvailableAtCellularTime = ...
    false(numberOfCellularSamples,1);

gnssAvailableAtCellularTime(hasMatchingGnssSample) = ...
    gnss.available( ...
        gnssIndexAtCellularTime(hasMatchingGnssSample));

headingAvailableAtCellularTime = isfinite( ...
    headingSensor.measurement(cellular.sampleIndices));

% Require associated TDOAs for every receiver pair and base station
allAssociatedTDOAAvailable = all(isfinite(reshape( ...
    cellular.associatedTDOA, ...
    numberOfCellularSamples,[])),2);

mappingCandidateMask = ...
    cellular.time < gnssOutageStart ...
    & gnssAvailableAtCellularTime ...
    & headingAvailableAtCellularTime ...
    & allAssociatedTDOAAvailable;

mappingCandidateIndices = find(mappingCandidateMask);

if numel(mappingCandidateIndices) < 2
    error('Fewer than two valid mapping observations are available.');
end

% Select observations distributed across the pre-outage interval
numberOfSelectedMappingObservations = min( ...
    numberOfMappingObservations, ...
    numel(mappingCandidateIndices));

mappingSelectionPositions = round(linspace( ...
    1, ...
    numel(mappingCandidateIndices), ...
    numberOfSelectedMappingObservations)).';

mappingCellularIndices = ...
    mappingCandidateIndices(mappingSelectionPositions);

mappingGnssIndices = ...
    gnssIndexAtCellularTime(mappingCellularIndices);

mappingTrajectoryIndices = ...
    cellular.sampleIndices(mappingCellularIndices);

% Operational mapping inputs
mappingTime = cellular.time(mappingCellularIndices);

mappingAircraftPosition = ...
    gnss.position(mappingGnssIndices,:);

mappingAircraftHeading = ...
    headingSensor.measurement(mappingTrajectoryIndices);

mappingTDOA = ...
    cellular.associatedTDOA(mappingCellularIndices,:,:);

% Calculate receiver positions using measured aircraft position and heading
mappingReceiverGlobal = zeros( ...
    numberOfSelectedMappingObservations, ...
    2, ...
    numberOfReceivers);

for observationIndex = 1:numberOfSelectedMappingObservations

    psi = mappingAircraftHeading(observationIndex);

    rotationBodyToGlobal = [
        cos(psi), -sin(psi);
        sin(psi),  cos(psi)
        ];

    for receiverIndex = 1:numberOfReceivers

        receiverOffsetGlobal = rotationBodyToGlobal ...
            * receiverBody(receiverIndex,:).';

        mappingReceiverGlobal( ...
            observationIndex,:,receiverIndex) = ...
            mappingAircraftPosition(observationIndex,:) ...
            + receiverOffsetGlobal.';
    end
end

%% PART 9 — Estimate and freeze the base-station map

numberOfMappingSamples = ...
    numberOfSelectedMappingObservations;

mappingDirection = nan( ...
    numberOfMappingSamples,2,numberOfBaseStations);

mappingBearing = nan( ...
    numberOfMappingSamples,numberOfBaseStations);

mappingDirectionMagnitude = nan( ...
    numberOfMappingSamples,numberOfBaseStations);

% Estimate a bearing from each selected aircraft position
for observationIndex = 1:numberOfMappingSamples

    referencePosition = reshape( ...
        mappingReceiverGlobal( ...
            observationIndex,:,referenceReceiverIndex),1,2);

    receiverGeometryMatrix = zeros(numberOfTDOAPairs,2);

    for pairIndex = 1:numberOfTDOAPairs

        receiverIndex = ...
            nonreferenceReceiverIndices(pairIndex);

        receiverPosition = reshape( ...
            mappingReceiverGlobal( ...
                observationIndex,:,receiverIndex),1,2);

        receiverGeometryMatrix(pairIndex,:) = ...
            -(receiverPosition-referencePosition);
    end

    for baseStationIndex = 1:numberOfBaseStations

        rangeDifference = speedOfLight*reshape( ...
            mappingTDOA( ...
                observationIndex,:,baseStationIndex),[],1);

        directionEstimate = ...
            receiverGeometryMatrix\rangeDifference;

        directionMagnitude = norm(directionEstimate);

        if ~isfinite(directionMagnitude) ...
                || directionMagnitude <= sqrt(eps)
            continue;
        end

        directionUnit = directionEstimate/directionMagnitude;

        mappingDirection( ...
            observationIndex,:,baseStationIndex) = ...
            directionUnit.';

        mappingDirectionMagnitude( ...
            observationIndex,baseStationIndex) = ...
            directionMagnitude;

        mappingBearing( ...
            observationIndex,baseStationIndex) = ...
            atan2(directionUnit(2),directionUnit(1));
    end
end

% Obtain an initial position by intersecting the bearing lines
mappingInitialBaseStationPosition = ...
    nan(numberOfBaseStations,2);

mappingBearingGeometryRank = ...
    zeros(numberOfBaseStations,1);

referenceReceiverPositions = ...
    mappingReceiverGlobal(:,:,referenceReceiverIndex);

for baseStationIndex = 1:numberOfBaseStations

    bearingLineMatrix = zeros(numberOfMappingSamples,2);
    bearingLineRightHandSide = zeros(numberOfMappingSamples,1);

    for observationIndex = 1:numberOfMappingSamples

        direction = reshape( ...
            mappingDirection( ...
                observationIndex,:,baseStationIndex),1,2);

        lineNormal = [-direction(2),direction(1)];

        bearingLineMatrix(observationIndex,:) = lineNormal;

        bearingLineRightHandSide(observationIndex) = ...
            lineNormal ...
            * referenceReceiverPositions(observationIndex,:).';
    end

    mappingBearingGeometryRank(baseStationIndex) = ...
        rank(bearingLineMatrix);

    if mappingBearingGeometryRank(baseStationIndex) < 2
        continue;
    end

    mappingInitialBaseStationPosition(baseStationIndex,:) = ...
        (bearingLineMatrix\bearingLineRightHandSide).';
end

% Prepare receiver coordinates for exact TDOA refinement
referenceReceiverEast = ...
    referenceReceiverPositions(:,1);

referenceReceiverNorth = ...
    referenceReceiverPositions(:,2);

nonreferenceReceiverEast = reshape( ...
    mappingReceiverGlobal( ...
        :,1,nonreferenceReceiverIndices), ...
    numberOfMappingSamples,numberOfTDOAPairs);

nonreferenceReceiverNorth = reshape( ...
    mappingReceiverGlobal( ...
        :,2,nonreferenceReceiverIndices), ...
    numberOfMappingSamples,numberOfTDOAPairs);

mappingEstimatedBaseStationPosition = ...
    nan(numberOfBaseStations,2);

mappingBaseStationPositionCovariance = ...
    nan(2,2,numberOfBaseStations);

mappingFinalCost = nan(numberOfBaseStations,1);
mappingResidualRMS = nan(numberOfBaseStations,1);
mappingOptimizationExitFlag = zeros(numberOfBaseStations,1);
mappingInformationRank = zeros(numberOfBaseStations,1);

rangeDifferenceVariance = ...
    (speedOfLight*tdoaAssumedNoiseStd)^2;

if rangeDifferenceVariance <= 0
    error('The assumed TDOA uncertainty must be greater than zero.');
end

% Refine each position using the exact range-difference equations
for baseStationIndex = 1:numberOfBaseStations

    initialPosition = ...
        mappingInitialBaseStationPosition(baseStationIndex,:);

    if any(~isfinite(initialPosition))
        continue;
    end

    measuredRangeDifference = ...
        speedOfLight*mappingTDOA(:,:,baseStationIndex);

    rangeToReference = @(stationPosition) sqrt( ...
        (referenceReceiverEast-stationPosition(1)).^2 ...
        + (referenceReceiverNorth-stationPosition(2)).^2);

    rangeToNonreference = @(stationPosition) sqrt( ...
        (nonreferenceReceiverEast-stationPosition(1)).^2 ...
        + (nonreferenceReceiverNorth-stationPosition(2)).^2);

    tdoaCostFunction = @(stationPosition) sum(sum( ...
        (rangeToNonreference(stationPosition) ...
        - rangeToReference(stationPosition) ...
        - measuredRangeDifference).^2));

    [estimatedPosition,finalCost,exitFlag] = fminsearch( ...
        tdoaCostFunction, ...
        initialPosition, ...
        mappingOptimizationOptions);

    mappingEstimatedBaseStationPosition( ...
        baseStationIndex,:) = estimatedPosition;

    mappingFinalCost(baseStationIndex) = finalCost;
    mappingResidualRMS(baseStationIndex) = sqrt( ...
        finalCost/(numberOfMappingSamples*numberOfTDOAPairs));

    mappingOptimizationExitFlag(baseStationIndex) = exitFlag;

    % Jacobian of all exact range-difference residuals
    mappingJacobian = zeros( ...
        numberOfMappingSamples*numberOfTDOAPairs,2);

    jacobianRow = 0;

    for observationIndex = 1:numberOfMappingSamples

        referenceVector = estimatedPosition ...
            - referenceReceiverPositions(observationIndex,:);

        referenceDerivative = ...
            referenceVector/norm(referenceVector);

        for pairIndex = 1:numberOfTDOAPairs

            receiverIndex = ...
                nonreferenceReceiverIndices(pairIndex);

            receiverPosition = reshape( ...
                mappingReceiverGlobal( ...
                    observationIndex,:,receiverIndex),1,2);

            receiverVector = ...
                estimatedPosition-receiverPosition;

            receiverDerivative = ...
                receiverVector/norm(receiverVector);

            jacobianRow = jacobianRow+1;

            mappingJacobian(jacobianRow,:) = ...
                receiverDerivative-referenceDerivative;
        end
    end

    mappingInformationMatrix = ...
        (mappingJacobian.'*mappingJacobian) ...
        / rangeDifferenceVariance;

    mappingInformationRank(baseStationIndex) = ...
        rank(mappingInformationMatrix);

    if mappingInformationRank(baseStationIndex) < 2
        continue;
    end

    mappingBaseStationPositionCovariance( ...
        :,:,baseStationIndex) = ...
        pinv(mappingInformationMatrix) ...
        + gnssPositionNoiseStd^2*eye(2);
end

mappingUsableBaseStation = ...
    mappingBearingGeometryRank == 2 ...
    & mappingOptimizationExitFlag > 0 ...
    & mappingInformationRank == 2 ...
    & all(isfinite(mappingEstimatedBaseStationPosition),2);

landmarkMap.cellID = ...
    baseStationCellID(mappingUsableBaseStation);

landmarkMap.catalogueIndex = ...
    find(mappingUsableBaseStation);

landmarkMap.position = ...
    mappingEstimatedBaseStationPosition( ...
        mappingUsableBaseStation,:);

landmarkMap.positionCovariance = ...
    mappingBaseStationPositionCovariance( ...
        :,:,mappingUsableBaseStation);

landmarkMap.residualRMS = ...
    mappingResidualRMS(mappingUsableBaseStation);

landmarkMap.count = nnz(mappingUsableBaseStation);
landmarkMap.frozen = true;

%% PART 10 — Configure the navigation EKF

% State: [East, North, vEast, vNorth, heading, bAx, bAy, bGyro]'
stateIndex.position = 1:2;
stateIndex.velocity = 3:4;
stateIndex.heading = 5;
stateIndex.accelerometerBias = 6:7;
stateIndex.gyroscopeBias = 8;

if ~gnss.available(1) ...
        || any(~isfinite(gnss.position(1,:))) ...
        || any(~isfinite(gnss.velocity(1,:)))
    error('An initial GNSS measurement is required.');
end

initialNavigationState = [
    gnss.position(1,1);
    gnss.position(1,2);
    gnss.velocity(1,1);
    gnss.velocity(1,2);
    headingSensor.measurement(1);
    0;
    0;
    0
    ];

initialPositionUncertainty = gnssPositionNoiseStd;
initialVelocityUncertainty = gnssVelocityNoiseStd;

initialNavigationCovariance = diag([
    initialPositionUncertainty^2;
    initialPositionUncertainty^2;
    initialVelocityUncertainty^2;
    initialVelocityUncertainty^2;
    initialHeadingUncertainty^2;
    initialAccelerometerBiasUncertainty^2;
    initialAccelerometerBiasUncertainty^2;
    initialGyroscopeBiasUncertainty^2
    ]);

% GNSS measures position and velocity
gnssMeasurementMatrix = zeros(4,numberOfNavigationStates);
gnssMeasurementMatrix(:,1:4) = eye(4);

gnssMeasurementNoiseCovariance = diag([
    gnssPositionNoiseStd^2;
    gnssPositionNoiseStd^2;
    gnssVelocityNoiseStd^2;
    gnssVelocityNoiseStd^2
    ]);

% The independent heading sensor measures state element 5
headingMeasurementMatrix = zeros(1,numberOfNavigationStates);
headingMeasurementMatrix(stateIndex.heading) = 1;

headingMeasurementNoiseVariance = ...
    headingSensorNoiseStd^2;

% IMU noise and bias-random-walk intensities
processNoiseSourceCovariance = diag([
    accelerometerNoiseStd^2;
    accelerometerNoiseStd^2;
    gyroscopeNoiseStd^2;
    accelerometerBiasRandomWalkStd^2;
    accelerometerBiasRandomWalkStd^2;
    gyroscopeBiasRandomWalkStd^2
    ]);

navigationIdentityMatrix = ...
    eye(numberOfNavigationStates);

%% PART 11 — Generate LTE bearings during the GNSS outage

bodyReceiverGeometryMatrix = -receiverBaselineBody;

if rank(bodyReceiverGeometryMatrix) < 2
    error('The receiver geometry cannot estimate a 2D direction.');
end

bodyReceiverGeometryPseudoInverse = ...
    pinv(bodyReceiverGeometryMatrix);

rangeDifferenceCovariance = ...
    (speedOfLight*tdoaAssumedNoiseStd)^2 ...
    * eye(numberOfTDOAPairs);

rawBodyDirectionCovariance = ...
    bodyReceiverGeometryPseudoInverse ...
    * rangeDifferenceCovariance ...
    * bodyReceiverGeometryPseudoInverse.';

numberOfLandmarks = landmarkMap.count;

lteBearingAvailable = false( ...
    numberOfCellularSamples,numberOfLandmarks);

lteBearingHistory = nan( ...
    numberOfCellularSamples,numberOfLandmarks);

lteBearingVarianceHistory = nan( ...
    numberOfCellularSamples,numberOfLandmarks);

lteBodyBearingVarianceHistory = nan( ...
    numberOfCellularSamples,numberOfLandmarks);

lteDirectionMagnitudeHistory = nan( ...
    numberOfCellularSamples,numberOfLandmarks);

lteGlobalDirectionHistory = nan( ...
    numberOfCellularSamples,2,numberOfLandmarks);

lteHeadingUsedHistory = nan(numberOfCellularSamples,1);

outageCellularIndices = find( ...
    cellular.time >= gnssOutageStart ...
    & cellular.time <= gnssOutageEnd);

for cellularIndex = outageCellularIndices.'

    trajectoryIndex = ...
        cellular.sampleIndices(cellularIndex);

    availableHeading = ...
        headingSensor.measurement(trajectoryIndex);

    if ~isfinite(availableHeading)
        continue;
    end

    lteHeadingUsedHistory(cellularIndex) = availableHeading;

    rotationBodyToGlobal = [
        cos(availableHeading), -sin(availableHeading);
        sin(availableHeading),  cos(availableHeading)
        ];

    for landmarkIndex = 1:numberOfLandmarks

        catalogueIndex = ...
            landmarkMap.catalogueIndex(landmarkIndex);

        currentTDOA = reshape( ...
            cellular.associatedTDOA( ...
                cellularIndex,:,catalogueIndex),[],1);

        if any(~isfinite(currentTDOA))
            continue;
        end

        rangeDifference = speedOfLight*currentTDOA;

        rawDirectionBody = ...
            bodyReceiverGeometryPseudoInverse ...
            * rangeDifference;

        directionMagnitude = norm(rawDirectionBody);

        if ~isfinite(directionMagnitude) ...
                || directionMagnitude <= sqrt(eps)
            continue;
        end

        directionBody = ...
            rawDirectionBody/directionMagnitude;

        directionNormalizationJacobian = ...
            (eye(2)-directionBody*directionBody.') ...
            / directionMagnitude;

        normalizedDirectionCovariance = ...
            directionNormalizationJacobian ...
            * rawBodyDirectionCovariance ...
            * directionNormalizationJacobian.';

        bodyBearingGradient = [
            -directionBody(2), ...
             directionBody(1)
            ];

        bodyBearingVariance = ...
            bodyBearingGradient ...
            * normalizedDirectionCovariance ...
            * bodyBearingGradient.';

        bodyBearingVariance = max(bodyBearingVariance,eps);

        directionGlobal = ...
            rotationBodyToGlobal*directionBody;

        globalBearing = atan2( ...
            directionGlobal(2),directionGlobal(1));

        globalBearingVariance = ...
            bodyBearingVariance ...
            + headingSensorNoiseStd^2;

        lteBearingAvailable( ...
            cellularIndex,landmarkIndex) = true;

        lteBearingHistory( ...
            cellularIndex,landmarkIndex) = ...
            globalBearing;

        lteBearingVarianceHistory( ...
            cellularIndex,landmarkIndex) = ...
            globalBearingVariance;

        lteBodyBearingVarianceHistory( ...
            cellularIndex,landmarkIndex) = ...
            bodyBearingVariance;

        lteDirectionMagnitudeHistory( ...
            cellularIndex,landmarkIndex) = ...
            directionMagnitude;

        lteGlobalDirectionHistory( ...
            cellularIndex,:,landmarkIndex) = ...
            directionGlobal.';
    end
end

%% PART 12 — Calculate standalone LTE position measurements

ltePositionHistory = nan(numberOfCellularSamples,2);
ltePositionAvailable = false(numberOfCellularSamples,1);

ltePositionCovarianceHistory = ...
    nan(2,2,numberOfCellularSamples);

lteStationsUsedHistory = ...
    zeros(numberOfCellularSamples,1);

lteBearingResidualRMSHistory = ...
    nan(numberOfCellularSamples,1);

lteGeometryRcondHistory = ...
    nan(numberOfCellularSamples,1);

headingVariance = headingSensorNoiseStd^2;
quarterTurnMatrix = [0,-1;1,0];

for cellularIndex = outageCellularIndices.'

    availableLandmarkIndices = find( ...
        lteBearingAvailable(cellularIndex,:));

    numberOfAvailableLandmarks = ...
        numel(availableLandmarkIndices);

    if numberOfAvailableLandmarks ...
            < minimumBaseStationsForPositionFix
        continue;
    end

    stationPositions = ...
        landmarkMap.position(availableLandmarkIndices,:);

    stationPositionCovariance = ...
        landmarkMap.positionCovariance( ...
            :,:,availableLandmarkIndices);

    globalDirections = squeeze( ...
        lteGlobalDirectionHistory( ...
            cellularIndex,:,availableLandmarkIndices)).';

    bodyBearingVariances = ...
        lteBodyBearingVarianceHistory( ...
            cellularIndex,availableLandmarkIndices).';

    bearingLineMatrix = [
        -globalDirections(:,2), ...
         globalDirections(:,1)
        ];

    bearingLineRightHandSide = sum( ...
        bearingLineMatrix.*stationPositions,2);

    if rank(bearingLineMatrix) < 2
        continue;
    end

    % Initial unweighted estimate of reference-receiver position
    estimatedReferenceReceiverPosition = ...
        bearingLineMatrix\bearingLineRightHandSide;

    geometryValid = true;
    bearingPositionInformationMatrix = nan(2,2);

    % Update distance-dependent bearing weights
    for weightingIteration = 1:2

        stationDistance = vecnorm( ...
            stationPositions ...
            - estimatedReferenceReceiverPosition.',2,2);

        independentLineVariance = zeros( ...
            numberOfAvailableLandmarks,1);

        for landmarkCounter = 1:numberOfAvailableLandmarks

            lineNormal = ...
                bearingLineMatrix(landmarkCounter,:);

            mapNormalVariance = ...
                lineNormal ...
                * stationPositionCovariance( ...
                    :,:,landmarkCounter) ...
                * lineNormal.';

            independentLineVariance(landmarkCounter) = ...
                stationDistance(landmarkCounter)^2 ...
                * bodyBearingVariances(landmarkCounter) ...
                + mapNormalVariance ...
                + 1e-6;
        end

        % The same heading error rotates every bearing.
        headingLineSensitivity = stationDistance;

        bearingLineCovariance = ...
            diag(independentLineVariance) ...
            + headingVariance ...
            * (headingLineSensitivity ...
            * headingLineSensitivity.');

        bearingPositionInformationMatrix = ...
            bearingLineMatrix.' ...
            * (bearingLineCovariance\bearingLineMatrix);

        informationRcond = ...
            rcond(bearingPositionInformationMatrix);

        if rank(bearingPositionInformationMatrix) < 2 ...
                || informationRcond ...
                < minimumBearingInformationRcond

            geometryValid = false;
            break;
        end

        estimatedReferenceReceiverPosition = ...
            bearingPositionInformationMatrix ...
            \ (bearingLineMatrix.' ...
            * (bearingLineCovariance ...
            \ bearingLineRightHandSide));
    end

    if ~geometryValid ...
            || any(~isfinite( ...
            estimatedReferenceReceiverPosition))
        continue;
    end

    referenceReceiverCovariance = ...
        pinv(bearingPositionInformationMatrix);

    % Convert reference-receiver position to aircraft-centre position
    availableHeading = ...
        lteHeadingUsedHistory(cellularIndex);

    rotationBodyToGlobal = [
        cos(availableHeading), -sin(availableHeading);
        sin(availableHeading),  cos(availableHeading)
        ];

    referenceReceiverOffsetGlobal = ...
        rotationBodyToGlobal ...
        * receiverBody(referenceReceiverIndex,:).';

    estimatedAircraftPosition = ...
        estimatedReferenceReceiverPosition ...
        - referenceReceiverOffsetGlobal;

    offsetHeadingJacobian = ...
        -rotationBodyToGlobal ...
        * quarterTurnMatrix ...
        * receiverBody(referenceReceiverIndex,:).';

    aircraftPositionCovariance = ...
        referenceReceiverCovariance ...
        + offsetHeadingJacobian ...
        * headingVariance ...
        * offsetHeadingJacobian.';

    aircraftPositionCovariance = ...
        0.5*(aircraftPositionCovariance ...
        + aircraftPositionCovariance.') ...
        + 1e-6*eye(2);

    lineResidual = ...
        bearingLineMatrix ...
        * estimatedReferenceReceiverPosition ...
        - bearingLineRightHandSide;

    ltePositionHistory(cellularIndex,:) = ...
        estimatedAircraftPosition.';

    ltePositionCovarianceHistory( ...
        :,:,cellularIndex) = ...
        aircraftPositionCovariance;

    ltePositionAvailable(cellularIndex) = true;

    lteStationsUsedHistory(cellularIndex) = ...
        numberOfAvailableLandmarks;

    lteBearingResidualRMSHistory(cellularIndex) = ...
        sqrt(mean(lineResidual.^2));

    lteGeometryRcondHistory(cellularIndex) = ...
        rcond(bearingPositionInformationMatrix);
end

numberOfLTEPositionFixes = ...
    nnz(ltePositionAvailable);

%% PART 13 — Run reference and LTE-aided navigation filters

referenceSolutionIndex = 1;
lteAidedSolutionIndex = 2;
numberOfNavigationSolutions = 2;

navigationState = repmat( ...
    initialNavigationState,1,numberOfNavigationSolutions);

navigationCovariance = repmat( ...
    initialNavigationCovariance,1,1,numberOfNavigationSolutions);

navigationStateHistory = nan( ...
    N,numberOfNavigationStates,numberOfNavigationSolutions);

for solutionIndex = 1:numberOfNavigationSolutions
    navigationStateHistory(1,:,solutionIndex) = ...
        initialNavigationState.';
end

gnssIndexAtTrajectory = zeros(N,1);
gnssIndexAtTrajectory(gnss.sampleIndices) = ...
    (1:numberOfGnssSamples).';

cellularIndexAtTrajectory = zeros(N,1);
cellularIndexAtTrajectory(cellular.sampleIndices) = ...
    (1:numberOfCellularSamples).';

ltePositionMeasurementMatrix = zeros( ...
    2,numberOfNavigationStates);

ltePositionMeasurementMatrix(:,stateIndex.position) = eye(2);

lteInnovationNISHistory = nan(numberOfCellularSamples,1);
lteUpdateAccepted = false(numberOfCellularSamples,1);

for k = 1:N-1

    for solutionIndex = 1:numberOfNavigationSolutions

        currentState = ...
            navigationState(:,solutionIndex);

        currentCovariance = ...
            navigationCovariance(:,:,solutionIndex);

        currentPosition = ...
            currentState(stateIndex.position);

        currentVelocity = ...
            currentState(stateIndex.velocity);

        currentHeading = ...
            currentState(stateIndex.heading);

        correctedAccelerationBody = ...
            imu.accelerationBody(k,:).' ...
            - currentState(stateIndex.accelerometerBias);

        correctedTurnRate = ...
            imu.gyroscopeZ(k) ...
            - currentState(stateIndex.gyroscopeBias);

        rotationBodyToGlobal = [
            cos(currentHeading), -sin(currentHeading);
            sin(currentHeading),  cos(currentHeading)
            ];

        correctedAccelerationGlobal = ...
            rotationBodyToGlobal*correctedAccelerationBody;

        predictedState = currentState;

        predictedState(stateIndex.position) = ...
            currentPosition ...
            + currentVelocity*dt ...
            + 0.5*correctedAccelerationGlobal*dt^2;

        predictedState(stateIndex.velocity) = ...
            currentVelocity ...
            + correctedAccelerationGlobal*dt;

        predictedState(stateIndex.heading) = ...
            currentHeading+correctedTurnRate*dt;

        predictedState(stateIndex.heading) = atan2( ...
            sin(predictedState(stateIndex.heading)), ...
            cos(predictedState(stateIndex.heading)));

        accelerationHeadingDerivative = [
            -sin(currentHeading), -cos(currentHeading);
             cos(currentHeading), -sin(currentHeading)
            ] * correctedAccelerationBody;

        stateTransitionMatrix = navigationIdentityMatrix;

        stateTransitionMatrix( ...
            stateIndex.position,stateIndex.velocity) = ...
            dt*eye(2);

        stateTransitionMatrix( ...
            stateIndex.position,stateIndex.heading) = ...
            0.5*accelerationHeadingDerivative*dt^2;

        stateTransitionMatrix( ...
            stateIndex.velocity,stateIndex.heading) = ...
            accelerationHeadingDerivative*dt;

        stateTransitionMatrix( ...
            stateIndex.position,stateIndex.accelerometerBias) = ...
            -0.5*rotationBodyToGlobal*dt^2;

        stateTransitionMatrix( ...
            stateIndex.velocity,stateIndex.accelerometerBias) = ...
            -rotationBodyToGlobal*dt;

        stateTransitionMatrix( ...
            stateIndex.heading,stateIndex.gyroscopeBias) = ...
            -dt;

        processNoiseMapping = zeros( ...
            numberOfNavigationStates,6);

        processNoiseMapping( ...
            stateIndex.position,1:2) = ...
            0.5*rotationBodyToGlobal*dt^2;

        processNoiseMapping( ...
            stateIndex.velocity,1:2) = ...
            rotationBodyToGlobal*dt;

        processNoiseMapping(stateIndex.heading,3) = dt;

        processNoiseMapping( ...
            stateIndex.accelerometerBias,4:5) = ...
            sqrt(dt)*eye(2);

        processNoiseMapping( ...
            stateIndex.gyroscopeBias,6) = sqrt(dt);

        processNoiseCovariance = ...
            processNoiseMapping ...
            * processNoiseSourceCovariance ...
            * processNoiseMapping.';

        predictedCovariance = ...
            stateTransitionMatrix ...
            * currentCovariance ...
            * stateTransitionMatrix.' ...
            + processNoiseCovariance;

        predictedCovariance = ...
            0.5*(predictedCovariance+predictedCovariance.');

        currentState = predictedState;
        currentCovariance = predictedCovariance;

        % Heading update
        headingMeasurement = ...
            headingSensor.measurement(k+1);

        headingInnovation = atan2( ...
            sin(headingMeasurement ...
            - headingMeasurementMatrix*currentState), ...
            cos(headingMeasurement ...
            - headingMeasurementMatrix*currentState));

        headingInnovationVariance = ...
            headingMeasurementMatrix ...
            * currentCovariance ...
            * headingMeasurementMatrix.' ...
            + headingMeasurementNoiseVariance;

        headingKalmanGain = ...
            currentCovariance ...
            * headingMeasurementMatrix.' ...
            / headingInnovationVariance;

        currentState = ...
            currentState ...
            + headingKalmanGain*headingInnovation;

        currentState(stateIndex.heading) = atan2( ...
            sin(currentState(stateIndex.heading)), ...
            cos(currentState(stateIndex.heading)));

        covarianceCorrectionMatrix = ...
            navigationIdentityMatrix ...
            - headingKalmanGain*headingMeasurementMatrix;

        currentCovariance = ...
            covarianceCorrectionMatrix ...
            * currentCovariance ...
            * covarianceCorrectionMatrix.' ...
            + headingKalmanGain ...
            * headingMeasurementNoiseVariance ...
            * headingKalmanGain.';

        % GNSS update when available
        gnssIndex = gnssIndexAtTrajectory(k+1);

        if gnssIndex > 0 && gnss.available(gnssIndex)

            gnssMeasurement = [
                gnss.position(gnssIndex,:).';
                gnss.velocity(gnssIndex,:).'
                ];

            gnssInnovation = ...
                gnssMeasurement ...
                - gnssMeasurementMatrix*currentState;

            gnssInnovationCovariance = ...
                gnssMeasurementMatrix ...
                * currentCovariance ...
                * gnssMeasurementMatrix.' ...
                + gnssMeasurementNoiseCovariance;

            gnssKalmanGain = ...
                currentCovariance ...
                * gnssMeasurementMatrix.' ...
                / gnssInnovationCovariance;

            currentState = ...
                currentState+gnssKalmanGain*gnssInnovation;

            covarianceCorrectionMatrix = ...
                navigationIdentityMatrix ...
                - gnssKalmanGain*gnssMeasurementMatrix;

            currentCovariance = ...
                covarianceCorrectionMatrix ...
                * currentCovariance ...
                * covarianceCorrectionMatrix.' ...
                + gnssKalmanGain ...
                * gnssMeasurementNoiseCovariance ...
                * gnssKalmanGain.';
        end

        % LTE update for the aided solution only
        cellularIndex = ...
            cellularIndexAtTrajectory(k+1);

        if solutionIndex == lteAidedSolutionIndex ...
                && cellularIndex > 0 ...
                && ltePositionAvailable(cellularIndex)

            lteMeasurement = ...
                ltePositionHistory(cellularIndex,:).';

            lteMeasurementCovariance = ...
                ltePositionCovarianceHistory( ...
                    :,:,cellularIndex);

            lteInnovation = ...
                lteMeasurement ...
                - ltePositionMeasurementMatrix*currentState;

            lteInnovationCovariance = ...
                ltePositionMeasurementMatrix ...
                * currentCovariance ...
                * ltePositionMeasurementMatrix.' ...
                + lteMeasurementCovariance;

            if all(isfinite(lteInnovationCovariance(:))) ...
                    && rcond(lteInnovationCovariance) > eps

                lteInnovationNIS = ...
                    lteInnovation.' ...
                    * (lteInnovationCovariance\lteInnovation);
            else
                lteInnovationNIS = Inf;
            end

            lteInnovationNISHistory(cellularIndex) = ...
                lteInnovationNIS;

            if isfinite(lteInnovationNIS) ...
                    && lteInnovationNIS ...
                    <= lteInnovationGateThreshold

                lteKalmanGain = ...
                    currentCovariance ...
                    * ltePositionMeasurementMatrix.' ...
                    / lteInnovationCovariance;

                currentState = ...
                    currentState ...
                    + lteKalmanGain*lteInnovation;

                covarianceCorrectionMatrix = ...
                    navigationIdentityMatrix ...
                    - lteKalmanGain ...
                    * ltePositionMeasurementMatrix;

                currentCovariance = ...
                    covarianceCorrectionMatrix ...
                    * currentCovariance ...
                    * covarianceCorrectionMatrix.' ...
                    + lteKalmanGain ...
                    * lteMeasurementCovariance ...
                    * lteKalmanGain.';

                lteUpdateAccepted(cellularIndex) = true;
            end
        end

        currentState(stateIndex.heading) = atan2( ...
            sin(currentState(stateIndex.heading)), ...
            cos(currentState(stateIndex.heading)));

        currentCovariance = ...
            0.5*(currentCovariance+currentCovariance.');

        navigationState(:,solutionIndex) = currentState;
        navigationCovariance(:,:,solutionIndex) = ...
            currentCovariance;

        navigationStateHistory(k+1,:,solutionIndex) = ...
            currentState.';
    end
end

referenceNavigationStateHistory = ...
    navigationStateHistory(:,:,referenceSolutionIndex);

lteAidedNavigationStateHistory = ...
    navigationStateHistory(:,:,lteAidedSolutionIndex);

numberOfLTEUpdatesAccepted = ...
    nnz(lteUpdateAccepted);

numberOfLTEUpdatesRejected = ...
    numberOfLTEPositionFixes-numberOfLTEUpdatesAccepted;

%% PART 14 — Print results and plot navigation performance

referencePosition = ...
    referenceNavigationStateHistory(:,stateIndex.position);

lteAidedPosition = ...
    lteAidedNavigationStateHistory(:,stateIndex.position);

referencePositionError = ...
    referencePosition-truth.position;

lteAidedPositionError = ...
    lteAidedPosition-truth.position;

referenceHorizontalError = vecnorm( ...
    referencePositionError,2,2);

lteAidedHorizontalError = vecnorm( ...
    lteAidedPositionError,2,2);

preOutageIndices = t < gnssOutageStart;

outageIndices = ...
    t >= gnssOutageStart ...
    & t <= gnssOutageEnd;

preOutagePositionRMSE = sqrt(mean( ...
    referenceHorizontalError(preOutageIndices).^2));

outagePositionRMSE = sqrt(mean( ...
    referenceHorizontalError(outageIndices).^2));

maximumOutagePositionError = max( ...
    referenceHorizontalError(outageIndices));

finalPositionError = ...
    referenceHorizontalError(end);

lastPreOutageIndex = find( ...
    t < gnssOutageStart,1,'last');

%% Calculate LTE-aided outage metrics

lteAidedOutageRMSE = sqrt(mean( ...
    lteAidedHorizontalError(outageIndices).^2));

lteAidedMaximumOutageError = max( ...
    lteAidedHorizontalError(outageIndices));

lteAidedFinalError = ...
    lteAidedHorizontalError(end);

% Positive values mean that LTE aiding improved the result.
if outagePositionRMSE > eps

    rmseImprovementPercent = 100 * ...
        (outagePositionRMSE-lteAidedOutageRMSE) ...
        / outagePositionRMSE;
else
    rmseImprovementPercent = NaN;
end

if finalPositionError > eps

    finalErrorImprovementPercent = 100 * ...
        (finalPositionError-lteAidedFinalError) ...
        / finalPositionError;
else
    finalErrorImprovementPercent = NaN;
end

%% Evaluate standalone LTE position fixes

lteFixCellularIndices = find( ...
    ltePositionAvailable ...
    & cellular.time >= gnssOutageStart ...
    & cellular.time <= gnssOutageEnd);

if isempty(lteFixCellularIndices)

    lteOnlyPositionRMSE = NaN;
    lteOnlyMaximumError = NaN;

else

    lteFixTrajectoryIndices = ...
        cellular.sampleIndices(lteFixCellularIndices);

    lteOnlyPositionError = ...
        ltePositionHistory(lteFixCellularIndices,:) ...
        - truth.position(lteFixTrajectoryIndices,:);

    lteOnlyHorizontalError = ...
        vecnorm(lteOnlyPositionError,2,2);

    lteOnlyPositionRMSE = sqrt(mean( ...
        lteOnlyHorizontalError.^2));

    lteOnlyMaximumError = ...
        max(lteOnlyHorizontalError);
end

%% Find the first outage cellular observation

firstOutageCellularIndex = find( ...
    cellular.time >= gnssOutageStart ...
    & cellular.time <= gnssOutageEnd, ...
    1,'first');



%% Print final base-station landmark map

fprintf('\nFinal base-station landmark map\n\n');

fprintf(['BS    PCI    Estimated East [m]    ', ...
    'Estimated North [m]    ', ...
    'True East [m]    True North [m]\n']);

fprintf(['--    ---    ------------------    ', ...
    '-------------------    ', ...
    '-------------    --------------\n']);

for landmarkIndex = 1:landmarkMap.count

    baseStationIndex = ...
        landmarkMap.catalogueIndex(landmarkIndex);

    fprintf(['%2d    %3d    %18.2f    ', ...
        '%19.2f    %13.2f    %14.2f\n'], ...
        baseStationIndex, ...
        landmarkMap.cellID(landmarkIndex), ...
        landmarkMap.position(landmarkIndex,1), ...
        landmarkMap.position(landmarkIndex,2), ...
        trueBaseStationPosition(baseStationIndex,1), ...
        trueBaseStationPosition(baseStationIndex,2));
end

fprintf('\nUsable landmarks: %d of %d\n', ...
    landmarkMap.count,numberOfBaseStations);

%% Print IQ-based PCI detection results

numberOfCorrectCellDetections = ...
    nnz(cellIDDetectionCorrectFromIQ);

cellDetectionAccuracyFromIQ = ...
    100*numberOfCorrectCellDetections ...
    / numberOfBaseStations;

fprintf('\nLTE PCI detection from IQ waveforms\n\n');

fprintf(['BS    True PCI    Detected PCI    ', ...
    'True timing    Det. timing    ', ...
    'PSS metric    SSS metric    Correct\n']);

fprintf(['--    --------    ------------    ', ...
    '-----------    -----------    ', ...
    '----------    ----------    -------\n']);

for baseStationIndex = 1:numberOfBaseStations

    fprintf(['%2d    %8d    %12.0f    ', ...
        '%11.0f    %11.0f    ', ...
        '%10.3f    %10.3f    %7d\n'], ...
        baseStationIndex, ...
        baseStationCellID(baseStationIndex), ...
        detectedCellIDFromIQ(baseStationIndex), ...
        truePssTimingOffsetFromIQ(baseStationIndex), ...
        detectedPssTimingOffsetFromIQ(baseStationIndex), ...
        pssCorrelationMetricFromIQ(baseStationIndex), ...
        sssCorrelationMetricFromIQ(baseStationIndex), ...
        cellIDDetectionCorrectFromIQ(baseStationIndex));
end

fprintf('\nIQ-based PCI detection accuracy: %.2f %%\n', ...
    cellDetectionAccuracyFromIQ);

%% Print cell-association summary

numberOfAssociationAttempts = ...
    numel(cellular.associationStatus);

numberAssociated = nnz( ...
    cellular.associationStatus ...
    == associationStatusCode.associated);

numberNoDetection = nnz( ...
    cellular.associationStatus ...
    == associationStatusCode.noDetection);

numberUnknown = nnz( ...
    cellular.associationStatus ...
    == associationStatusCode.unknownCellID);

numberAmbiguous = nnz( ...
    cellular.associationStatus ...
    == associationStatusCode.ambiguousCellID);

numberRejected = nnz( ...
    cellular.associationStatus ...
    == associationStatusCode.rejected);

correctAssociatedMeasurements = ...
    cellular.associationCorrect( ...
        cellular.associationSuccessful);

if isempty(correctAssociatedMeasurements)
    associationAccuracy = NaN;
else
    associationAccuracy = ...
        100*mean(double(correctAssociatedMeasurements));
end

fprintf('\nCell-association summary\n');

fprintf('Association attempts: %d\n', ...
    numberOfAssociationAttempts);

fprintf('Successfully associated: %d\n', ...
    numberAssociated);

fprintf('No detection: %d\n', ...
    numberNoDetection);

fprintf('Unknown cell ID: %d\n', ...
    numberUnknown);

fprintf('Ambiguous cell ID: %d\n', ...
    numberAmbiguous);

fprintf('Rejected detections: %d\n', ...
    numberRejected);

fprintf('Association accuracy: %.2f %%\n', ...
    associationAccuracy);

%% Print INS/GNSS verification metrics

fprintf('\nINS/GNSS navigation verification\n\n');

fprintf('Pre-outage horizontal RMSE: %.2f m\n', ...
    preOutagePositionRMSE);

fprintf('Outage horizontal RMSE: %.2f m\n', ...
    outagePositionRMSE);

fprintf('Maximum outage error: %.2f m\n', ...
    maximumOutagePositionError);

fprintf('Final horizontal error: %.2f m\n', ...
    finalPositionError);

fprintf('\nEstimated biases immediately before outage\n');

fprintf('Forward accelerometer bias: %.5f m/s^2\n', ...
    referenceNavigationStateHistory( ...
        lastPreOutageIndex,stateIndex.accelerometerBias(1)));

fprintf('True forward accelerometer bias: %.5f m/s^2\n', ...
    accelerometerBias(1));

fprintf('Left accelerometer bias: %.5f m/s^2\n', ...
    referenceNavigationStateHistory( ...
        lastPreOutageIndex,stateIndex.accelerometerBias(2)));

fprintf('True left accelerometer bias: %.5f m/s^2\n', ...
    accelerometerBias(2));

fprintf('Gyroscope bias: %.5f deg/s\n', ...
    rad2deg(referenceNavigationStateHistory( ...
        lastPreOutageIndex,stateIndex.gyroscopeBias)));

fprintf('True gyroscope bias: %.5f deg/s\n', ...
    rad2deg(gyroscopeBias));

%% Print LTE availability and update counts

if isempty(firstOutageCellularIndex)

    fprintf('\nFirst outage cellular observation: none\n');

else

    fprintf('\nFirst outage cellular observation: %.1f s\n', ...
        cellular.time(firstOutageCellularIndex));
end

fprintf('LTE bearing position fixes generated: %d\n', ...
    numberOfLTEPositionFixes);

fprintf('LTE bearing position fixes applied: %d\n', ...
    numberOfLTEUpdatesAccepted);

fprintf('LTE fixes rejected by innovation gate: %d\n', ...
    numberOfLTEUpdatesRejected);

%% Print GNSS-outage navigation comparison

fprintf('\nGNSS-outage navigation comparison\n\n');

fprintf(['Method              Outage RMSE [m]    ', ...
    'Maximum error [m]    Final error [m]\n']);

fprintf(['------------------  -----------------    ', ...
    '-----------------    ---------------\n']);

fprintf('INS only            %17.2f    %17.2f    %15.2f\n', ...
    outagePositionRMSE, ...
    maximumOutagePositionError, ...
    finalPositionError);

fprintf('LTE-bearing aided   %17.2f    %17.2f    %15.2f\n', ...
    lteAidedOutageRMSE, ...
    lteAidedMaximumOutageError, ...
    lteAidedFinalError);

fprintf('\nOutage RMSE improvement: %.2f %%\n', ...
    rmseImprovementPercent);

fprintf('Final-error improvement: %.2f %%\n', ...
    finalErrorImprovementPercent);

fprintf('LTE bearing position fixes: %d\n', ...
    numberOfLTEUpdatesAccepted);

fprintf('Standalone LTE position RMSE: %.2f m\n', ...
    lteOnlyPositionRMSE);

fprintf('Standalone LTE maximum error: %.2f m\n', ...
    lteOnlyMaximumError);

%% Figure 1 — Trajectories, LTE fixes and base stations

outageStartIndex = find( ...
    t >= gnssOutageStart,1,'first');

figure;

plot(truth.position(:,1)/1000, ...
    truth.position(:,2)/1000, ...
    'k','LineWidth',1.5);

hold on;

plot(referencePosition(:,1)/1000, ...
    referencePosition(:,2)/1000, ...
    'r--','LineWidth',1.5);

plot(lteAidedPosition(:,1)/1000, ...
    lteAidedPosition(:,2)/1000, ...
    'b-.','LineWidth',1.5);

plot(ltePositionHistory(ltePositionAvailable,1)/1000, ...
    ltePositionHistory(ltePositionAvailable,2)/1000, ...
    'co','MarkerSize',3);

plot(trueBaseStationPosition(:,1)/1000, ...
    trueBaseStationPosition(:,2)/1000, ...
    'r^','MarkerSize',8, ...
    'MarkerFaceColor','r');

plot(landmarkMap.position(:,1)/1000, ...
    landmarkMap.position(:,2)/1000, ...
    'ms','MarkerSize',8, ...
    'LineWidth',1.5);

plot(truth.position(1,1)/1000, ...
    truth.position(1,2)/1000, ...
    'go','MarkerSize',8, ...
    'MarkerFaceColor','g');

plot(truth.position(outageStartIndex,1)/1000, ...
    truth.position(outageStartIndex,2)/1000, ...
    'yo','MarkerSize',8, ...
    'MarkerFaceColor','y', ...
    'MarkerEdgeColor','k');

plot(truth.position(end,1)/1000, ...
    truth.position(end,2)/1000, ...
    'ko','MarkerSize',8, ...
    'MarkerFaceColor','k');

axis equal;
grid on;

xlabel('East [km]');
ylabel('North [km]');

title('INS and LTE-Bearing-Aided Navigation');

legend( ...
    'True trajectory', ...
    'INS only', ...
    'LTE-bearing aided', ...
    'Standalone LTE fixes', ...
    'True base stations', ...
    'Mapped base stations', ...
    'Start', ...
    'GNSS outage begins', ...
    'End', ...
    'Location','best');

%% Figure 2 — Horizontal position error

figure;

plot(t,referenceHorizontalError, ...
    'r','LineWidth',1.5);

hold on;

plot(t,lteAidedHorizontalError, ...
    'b','LineWidth',1.5);

currentVerticalLimits = ylim;

plot([gnssOutageStart,gnssOutageStart], ...
    currentVerticalLimits, ...
    'k--','LineWidth',1.5);

grid on;

xlabel('Time [s]');
ylabel('Horizontal position error [m]');

title('Horizontal Position Error Comparison');

legend( ...
    'INS only', ...
    'LTE-bearing aided', ...
    'GNSS outage begins', ...
    'Location','best');

%% Figure 3 — East and North position errors

figure;

subplot(2,1,1);

plot(t,referencePositionError(:,1), ...
    'r','LineWidth',1.2);

hold on;

plot(t,lteAidedPositionError(:,1), ...
    'b','LineWidth',1.2);

currentVerticalLimits = ylim;

plot([gnssOutageStart,gnssOutageStart], ...
    currentVerticalLimits, ...
    'k--','LineWidth',1);

grid on;
ylabel('East error [m]');
title('Position-Error Components');

legend( ...
    'INS only', ...
    'LTE-bearing aided', ...
    'GNSS outage begins', ...
    'Location','best');

subplot(2,1,2);

plot(t,referencePositionError(:,2), ...
    'r','LineWidth',1.2);

hold on;

plot(t,lteAidedPositionError(:,2), ...
    'b','LineWidth',1.2);

currentVerticalLimits = ylim;

plot([gnssOutageStart,gnssOutageStart], ...
    currentVerticalLimits, ...
    'k--','LineWidth',1);

grid on;
xlabel('Time [s]');
ylabel('North error [m]');

legend( ...
    'INS only', ...
    'LTE-bearing aided', ...
    'GNSS outage begins', ...
    'Location','best');

%% Figure 4 — Estimated IMU biases

figure;

%% Forward accelerometer bias

subplot(3,1,1);

plot(t,referenceNavigationStateHistory( ...
    :,stateIndex.accelerometerBias(1)), ...
    'r--','LineWidth',1.4);

hold on;

plot(t,lteAidedNavigationStateHistory( ...
    :,stateIndex.accelerometerBias(1)), ...
    'b','LineWidth',1.2);

plot(t,accelerometerBias(1)*ones(N,1), ...
    'k--','LineWidth',1.2);

currentVerticalLimits = ylim;

plot([gnssOutageStart,gnssOutageStart], ...
    currentVerticalLimits, ...
    'm--','LineWidth',1.2);

grid on;
xlim([t(1),t(end)]);

ylabel('Bias [m/s^2]');
title('Forward Accelerometer Bias Estimate');

legend( ...
    'Heading-aided INS', ...
    'LTE-aided INS', ...
    'True bias', ...
    'GNSS outage begins', ...
    'Location','best');

%% Left accelerometer bias

subplot(3,1,2);

plot(t,referenceNavigationStateHistory( ...
    :,stateIndex.accelerometerBias(2)), ...
    'r--','LineWidth',1.4);

hold on;

plot(t,lteAidedNavigationStateHistory( ...
    :,stateIndex.accelerometerBias(2)), ...
    'b','LineWidth',1.2);

plot(t,accelerometerBias(2)*ones(N,1), ...
    'k--','LineWidth',1.2);

currentVerticalLimits = ylim;

plot([gnssOutageStart,gnssOutageStart], ...
    currentVerticalLimits, ...
    'm--','LineWidth',1.2);

grid on;
xlim([t(1),t(end)]);

ylabel('Bias [m/s^2]');
title('Left Accelerometer Bias Estimate');

legend( ...
    'Heading-aided INS', ...
    'LTE-aided INS', ...
    'True bias', ...
    'GNSS outage begins', ...
    'Location','best');

%% Gyroscope bias

subplot(3,1,3);

plot(t,rad2deg(referenceNavigationStateHistory( ...
    :,stateIndex.gyroscopeBias)), ...
    'r--','LineWidth',1.4);

hold on;

plot(t,rad2deg(lteAidedNavigationStateHistory( ...
    :,stateIndex.gyroscopeBias)), ...
    'b','LineWidth',1.2);

plot(t,rad2deg(gyroscopeBias)*ones(N,1), ...
    'k--','LineWidth',1.2);

currentVerticalLimits = ylim;

plot([gnssOutageStart,gnssOutageStart], ...
    currentVerticalLimits, ...
    'm--','LineWidth',1.2);

grid on;
xlim([t(1),t(end)]);

xlabel('Time [s]');
ylabel('Bias [deg/s]');
title('Gyroscope Bias Estimate');

legend( ...
    'Heading-aided INS', ...
    'LTE-aided INS', ...
    'True bias', ...
    'GNSS outage begins', ...
    'Location','best');

sgtitle('EKF Bias Estimates');