%CDLChannel TR 38.900 CDL channel model

%  Copyright 2016-2017 The MathWorks, Inc.

function [pathGains,sampleTimes,AoDSpec,AoASpec,pdp] = CDLChannelNlos(model,info,pdp,insize,coupling,Phi)
    
    %Bandwidth 
    fc = model.CarrierFrequency;
    if fc == 6e9
    bw = 20e6; 
    interCarSpace = 0.5e6;
    else 
    bw = 100e6;
    interCarSpace = 2e6;
    end 
    
    nSub = (bw/interCarSpace) + 1;
    
    % named column subscripts
    [~,power] = nr5g.internal.namedDelayProfileColumns();
    
    %----------------------------------------------------------------------
    % channel matrix creation for LOS clusters, sub-clustered (split) NLOS
    % clusters and other NLOS clusters (step 11 and ray offset angles of
    % step 7)

    [HNLOS,AoDSpec,AoASpec,pdp] = generateClusterGains(model,insize,pdp,coupling,Phi,info.ClusterTypes,'NLOS',fc,bw,nSub,interCarSpace);  % size T_H-by-(Lin-LOS-NSC)-by-P-by-R
    
            
    %pathGains = cat(2,HLOS,HSplitNLOS,HNLOS); % size T_H-by-L-by-P-by-R
    pathGains = HNLOS;         %%change for LOS and NLOS
    sampleTimes = getSampleTimes(model,insize);
    
    %----------------------------------------------------------------------
    %  normalisations
    
    if (model.NormalizeChannelOutputs)
        pathGains = pathGains / sqrt(size(pathGains,4));
    end
    if (model.NormalizePathGains)
        p = sum(10.^(pdp(:,power)/10));
        pathGains = pathGains / sqrt(p);
    end
    
end

% Generation of cluster gains, covering:
%    * Section 7.5 Step 7: ray offset angles
%    * Section 7.5 Step 11: Generate channel coefficients
% Note that the phases produced by "Step 10: Draw initial random phases"
% are provided as an input 'Phi' here and used to created the polarization
% term (also, 'Phi' can be externally specified rather than randomly
% generated depending on MODEL.InitPhase). Ray coupling indices produced by
% "Step 8: coupling of rays" are provided as an input 'coupling' here. Note
% that steps 5 (Generate delays), 6 (Generate cluster powers), 7 (Generate
% arrival and departure angles) and 9 (Generate XPRs) are not performed and
% the values of these small-scale parameters for a single link are provided
% either from the tabulated CDL profiles or from the custom CDL option
% (MODEL.DelayProfile='Custom').
function [H,AoDSpec,AoASpec,pdp] = generateClusterGains(model,insize,pdp,coupling,Phi,allClusterTypes,thisClusterType,fc,bw,nSub,interCarSpace)
    
    % named column subscripts
    [~,power,AoD,AoA,ZoD,ZoA] = nr5g.internal.namedDelayProfileColumns();
    
    % get dimensionality information for a single cluster
    [D,t] = getClusterDimensionality(model,insize);
    T_H = D(1); % number of time-domain channel samples
    S = D(3);   % number of transmit antennas
    U = D(4);   % number of receive antennas
    
    % get part of delay profile corresponding to this cluster type
    pdp = pdp(strcmpi(allClusterTypes,thisClusterType),:);
    
    % get number of clusters and set up the corresponding dimension of the
    % channel matrix output 'H' ('H' will have dimension sizes 'D')
    N = size(pdp,1);
    D(2) = N;
    
    % create empty channel matrix 
    H = zeros(D);
    
    % if no clusters are active for this cluster type, return empty channel
    % matrix
    if (N==0)
        return;
    end
    
    % convert cluster powers to a linear scale
    P = 10.^(pdp(:,power).'/10);
    
    % get part of Phi matrix corresponding to this cluster type
    Phi = Phi(strcmpi(allClusterTypes,thisClusterType),:,:);
    
    % get part of ray coupling matrix corresponding to this cluster type
    coupling = coupling(strcmpi(allClusterTypes,thisClusterType),:,:);
    
    % get number of rays (note that some clusters are generated with a
    % subset of this number of rays, see 'R_i' below)
    M = nr5g.internal.getNumberOfRays(model);

    % get per-cluster parameters
    C_ASA = model.DesiredASA;
    C_ASD = model.DesiredASD;
    C_ZSA = model.DesiredZSA;
    C_ZSD = model.DesiredZSD;
    XPR = model.XPR;    
     
    % apply ray offset angles to AOA and AOD angles according to Section
    % 7.5 Equation 7.5-13
    ray_offset_alpha = unifrnd(-2,2,[1 20]); %kron([0.0447 0.1413 0.2492 0.3715 0.5129 0.6797 0.8844 1.1481 1.5195 2.1551],[2 -2]); % size 1-by-M
    if (strcmpi(thisClusterType,'LOS'))
        ray_offset_alpha(:) = 0; % no ray offseting for LOS path (only first ray will be used anyway)
    end
    phi_AOA = pdp(:,AoA) + C_ASA*ray_offset_alpha;   % size N-by-M
    phi_AOD = pdp(:,AoD) + C_ASD*ray_offset_alpha;   % size N-by-M
    
    % apply ray offset angles to ZOA and ZOD angles according to Section
    % 7.5 Equation 7.5-18
    theta_ZOA = pdp(:,ZoA) + C_ZSA*ray_offset_alpha; % size N-by-M
    theta_ZOD = pdp(:,ZoD) + C_ZSD*ray_offset_alpha; % size N-by-M
    
    % scale angles according to Section 7.7.5.1
    if(model.AngleScaling)
        phi_AOA = scaleCDLAngles(model,pdp,phi_AOA,'DesiredASA','DesiredMeanAoA',@(x)mod(x+180,360)-180); % wrap to -180...180
        phi_AOD = scaleCDLAngles(model,pdp,phi_AOD,'DesiredASD','DesiredMeanAoD',@(x)mod(x+180,360)-180); % wrap to -180...180
        theta_ZOA = scaleCDLAngles(model,pdp,theta_ZOA,'DesiredZSA','DesiredMeanZoA',@(x)max(min(x,180),0)); % clip to 0...180
        theta_ZOD = scaleCDLAngles(model,pdp,theta_ZOD,'DesiredZSD','DesiredMeanZoD',@(x)max(min(x,180),0)); % clip to 0...180
    end
    
    % wrap zenith angles to [0,360] and map [180,360] to [180,0]
    theta_ZOA = limitZenithAngles(theta_ZOA);
    theta_ZOD = limitZenithAngles(theta_ZOD);
    
    % perform coupling of rays (step 8)
    phi_AOA = phi_AOA(coupling(:,:,1));     % AoD to AoA coupling
    theta_ZOA = theta_ZOA(coupling(:,:,2)); % ZoD to ZoA coupling
    theta_ZOD = theta_ZOD(coupling(:,:,3)); % AoD to ZoD coupling
    
    % spherical unit vectors of departure for each cluster and each ray
    rhat_tx = phat(phi_AOD,theta_ZOD); % size 3-by-N-by-M
    AoDSpec = rhat_tx;
    %disp(size(AoDSpec));
    
    % transmit antenna location vector
    c = 299792458;                          % light speed, m/s
    lambda_0 = c / model.CarrierFrequency;  % carrier wavelength, m
    dbar_tx = model.TransmitAntennaArray.ElementPositions*lambda_0 + repmat(model.TransmitAntennaArray.Position,1,S); % size 3-by-S
    
    transmitFieldTerm = zeros(2,N*M,S);
    transmitLocationTerm = zeros(N*M,S);
    % for each transmit antenna:
    for s = 1:S

        % get transmit field term
        transmitFieldTerm(:,:,s) = getFieldTerm(model.TransmitAntennaArray,theta_ZOD,phi_AOD,s);

        fStart = fc - bw/2;
        for Nsb = 1:nSub 
          f = fStart + (Nsb-1)*interCarSpace;  
        % get transmit location term
        transmitLocationTerm(:,s,Nsb) = getLocationTerm(rhat_tx,dbar_tx,s,f);
        end 

    end
    
    % spherical unit vectors of arrival for each cluster and each ray
    rhat_rx = phat(phi_AOA,theta_ZOA); % size 3-by-N-by-M
    AoASpec = rhat_rx;
    disp(size(AoASpec));
    
    % receive antenna location vectors
    dbar_rx = model.ReceiveAntennaArray.ElementPositions*lambda_0 + repmat(model.ReceiveAntennaArray.Position,1,U); % size 3-by-U
    
    receiveFieldTerm = zeros(2,N*M,U);
    receiveLocationTerm = zeros(N*M,U);
    % for each receive antenna:
    for u = 1:U

        % get receive field term
        receiveFieldTerm(:,:,u) = getFieldTerm(model.ReceiveAntennaArray,theta_ZOA,phi_AOA,u);
        
        fStart = fc - bw/2;
        for Nsb = 1:nSub 
          f = fStart + (Nsb-1)*interCarSpace;  
        
        % get receive location term
        receiveLocationTerm(:,u,Nsb) = getLocationTerm(rhat_rx,dbar_rx,u,f);
        
        end 
        
    end
    
    % get polarization term, size 4-by-(N*M)
    if (strcmpi(thisClusterType,'LOS'))
        polarizationTerm = getLOSPolarizationMatrix(Phi);
    else
        polarizationTerm = getNLOSPolarizationMatrix(XPR,Phi);
    end
    
    % get Doppler term
    dopplerTerm = getDopplerTerm(model,rhat_rx,t); % size T_H-by-(N*M)

    % establish set of active rays for each cluster
    R = getActiveRays(Phi,T_H);
    
    % for each transmit antenna:
    %fStart = fc - bw/2;
    for Nsb = 1:nSub 
       %f = fStart + (nSub-1)*interCarSpace;
        
    for s = 1:S

        % for each receive antenna:
        for u = 1:U

            % calculate the term of the coefficients combining the transmit
            % and receive antenna fields and the polarization matrix;
            % output is of size (N*M)-by-1
            fieldTerms = ((receiveFieldTerm(1,:,u).*polarizationTerm(1,:) + receiveFieldTerm(2,:,u).*polarizationTerm(2,:)) .* transmitFieldTerm(1,:,s) + ...
                          (receiveFieldTerm(1,:,u).*polarizationTerm(3,:) + receiveFieldTerm(2,:,u).*polarizationTerm(4,:)) .* transmitFieldTerm(2,:,s)).';

            % combine the field terms with the location terms for the
            % current transmit and receive antenna; output is of size
            % (N*M)-by-1
            allTerms = fieldTerms .* receiveLocationTerm(:,u,Nsb) .* transmitLocationTerm(:,s,Nsb);

            % zero out any entries corresponding to unused rays for the
            % current set of clusters; output is of size (N*M)-by-1
            allTerms(R(:)==0) = 0;

            % apply the Doppler term; output is of size T_H-by-(N*M)
            allTerms = dopplerTerm .* allTerms.';

            % combine the rays and then assign into the appropriate part of
            % the overall cluster gain array; assignment is of size
            % T_H-by-N
            H(:,:,s,u,Nsb) = sum(reshape(allTerms,[T_H N M]),3);

        end

    end
    
    delay = pdp(:,1);
    expDelay = exp(-1i*2*pi*f*delay);
    delayTerm = expDelay;    %ignored because ray wise delays arent mentioned and cluster wise delay spread is ignored
    % -(rhat_rx.*dbar_rx)/c - (rhat_tx.*dbar_tx)/c;
    for nClus = 1:N
        H(:,nClus,:,:,:) = H(:,nClus,:,:,:) * delayTerm(nClus);
    end 
    
    
    end 

    % calculate number of active rays for each cluster, |R_i|
    modR_i = sum(R,2).';
    
    % apply final scaling term to H
    % the scaling here includes |R_i| because the desired cluster power
    % |R_i|/M for sub-clusters has already been applied to 'pdp', during
    % its initial construction
    H = H .* sqrt(P./modR_i);
    
end

% establish set of active rays for each cluster in matrix 'R'
% each row of R corresponds to one cluster, with the columns corresponding
% to the ray number - the 'i'th element of a row is set to 1 if ray R_i is
% active as follows:
% * R_i = 1 for LOS
% * R_i according to Table 7.5-5 for sub-clustered NLOS
% * R_i = 1...20 for other NLOS
% (R_i is only defined in TR 38.900 for sub-clustered NLOS clusters, but
% the variable is given the values above to allow uniform processing for
% the three cases)
% note that if T_H, the number of time samples, is zero then an
% appropriately sized empty matrix is created to facilitate creation of a
% correctly-shaped empty channel output
function R = getActiveRays(Phi,T_H)

    [N,M,~] = size(Phi);
    if (T_H==0)
        R = zeros(N,0);
    else
        R = ones(N,M);
        R(Phi(:,:,1)==-Inf) = 0;
    end
    
end

% TR 38.900 Equation 7.5-22/7.5-28/7.5-29, field term which represents the
% field pattern of a transmit or receive antenna element in the GCS,
% including any transformation from the LCS due to the antenna orientation
function fieldTerm_ant = getFieldTerm(array,theta,phi,ant)

    % reshape theta and phi to combine cluster and ray dimensions into a
    % single row
    theta = theta(:).';
    phi = phi(:).';
    
    % alpha: bearing angle of LCS w.r.t. GCS
    % beta: downtilt angle of LCS w.r.t. GCS
    % gamma: slant angle of LCS w.r.t GCS
    % note that the bearing (alpha) and downtilt (beta) angles in the LCS
    % are formed from the combination of the overall array orientation and
    % the orientation of the current element
    alpha = array.Orientation(1) + array.ElementOrientations(1,ant);
    beta = array.Orientation(2) + array.ElementOrientations(2,ant);
    gamma = array.Orientation(3);
    
    % TR 38.900 Section 7.1 Equations 7.1-7 and 7.1-8
    % LCS zenith (theta) and azimuth (phi) angles for a given LCS
    % orientation alpha/beta/gamma and Cartesian representation (rhat) of
    % GCS zenith and azimuth angles theta and phi
    theta_prime = acosd((cosd(beta).*cosd(gamma).*cosd(theta)) + (((sind(beta).*cosd(gamma).*cosd(phi-alpha)) - (sind(gamma).*sind(phi-alpha))).*sind(theta)));
    phi_prime = angle(    ((cosd(beta).*sind(theta).*cosd(phi-alpha)) - (sind(beta).*cosd(theta))) + ...
                       1i*((cosd(beta).*sind(gamma).*cosd(theta)) + (((sind(beta).*sind(gamma).*cosd(phi-alpha)) + (cosd(gamma).*sind(phi-alpha))).*sind(theta)))) * 180/pi;
    
    % TR 38.900 Section 7.1 Equation 7.1-15
    % Rotation of the spherical basis vectors theta_hat and phi_hat by an
    % angle Psi due to the orientation of the LCS w.r.t. the GCS
    Psi = angle(    ((sind(gamma).*cosd(theta).*sind(phi-alpha)) + cosd(gamma).*((cosd(beta).*sind(theta)) - (sind(beta).*cosd(theta).*cosd(phi-alpha)))) + ...
                 1i*((sind(gamma).*cosd(phi-alpha)) + (sind(beta).*cosd(gamma).*sind(phi-alpha)))) * 180/pi;
    
    % TR 36.873 Section 7.1.1
    % Polarization slant angle
    zeta = array.ElementOrientations(3,ant);
    
    % F = [F_theta_prime; F_phi_prime], where F_theta_prime and F_phi_prime
    % are the field components in the direction of the LCS spherical unit
    % vectors theta_hat_prime and phi_hat_prime respectively.
    F = feval(array.ElementPattern,theta_prime(:).',phi_prime(:).',zeta);
    
    % TR 38.900 Section 7.1 Equation 7.1-11
    % Rotation of element field pattern from LCS to GCS
    fieldTerm_ant = [F(1,:).*cosd(Psi) - F(2,:).*sind(Psi); F(1,:).*sind(Psi) + F(2,:).*cosd(Psi)];
    
end

% TR 38.900 Equation 7.5-22/7.5-28/7.5-29, exponential term which is a
% function of the spherical unit vector 'rhat' for arrival (receive) or
% departure (transmit) azimuth and zenith angles and the location of the
% receive or transmit antenna element 'dbar'
function locationTerm_ant = getLocationTerm(rhat,dbar,ant,f)

    % reshape rhat to combine cluster and ray dimensions into a single row
    rhat = reshape(rhat,[3 numel(rhat)/3]);

    c = 299792458;                          % light speed, m/s
    lambda_f = c / f;  % carrier wavelength, m
    locationTerm_ant = exp(1i*2*pi*rhat.'*dbar(:,ant)/lambda_f);

end

% TR 38.900 Equation 7.5-22/7.5-28/7.5-29, exponentional term which is a
% function of the Doppler due to UT velocity
function dopplerTerm = getDopplerTerm(model,rhat_rx,t)
    
    % reshape rhat_rx to combine cluster and ray dimensions into a single
    % row
    rhat_rx = reshape(rhat_rx,[3 numel(rhat_rx)/3]);

    c = 299792458;                            % light speed, m/s
    lambda_0 = c / model.CarrierFrequency;    % carrier wavelength, m
    v = model.MaximumDopplerShift * lambda_0; % UT speed, m/s
    theta_v = model.UTDirectionOfTravel(2);   % UT zenith angle of travel, degrees
    phi_v = model.UTDirectionOfTravel(1);     % UT azimuth angle of travel, degrees
    
    vbar = v * phat(phi_v,theta_v);    % UT velocity
    
    dopplerTerm = exp(kron(t,(1i*2*pi*rhat_rx.'*vbar/lambda_0).'));
    
end

% TR 38.900 Section 7.1 Equation 7.1-6
% Cartesian representation of points on the unit sphere with zenith angles
% 'theta' and azimuth angles 'phi'
function out = phat(phi,theta)

    phi = permute(phi,[3 1 2]);
    theta = permute(theta,[3 1 2]);
    
    sintheta = sind(theta);
    x = sintheta.*cosd(phi);
    y = sintheta.*sind(phi);
    z = cosd(theta);
    out = [x; y; z];
    
end

% TR 38.900 Section 7.5 step 11
% polarization matrix in Equations 7.5-22 and 7.5-28
function polarizationMatrix = getNLOSPolarizationMatrix(XPR,Phi)

    % reshape Phi to combine cluster and ray dimensions into a single row
    Phi = reshape(permute(Phi,[3 1 2]),[4 numel(Phi)/4]);
    
    % get cross polarization power ratio in linear scale, kappa
    kappa = 10^(XPR/10);
    
    % get polarization matrix in Equations 7.5-22 and 7.5-28
    % note that the 1st dimension of Phi, when interpreted as a 2-by-2
    % matrix, is transposed during this operation (1st dimension element
    % order [1 3 2 4]), as the order of the polarization terms in Phi
    % described by the standard and the MODEL.InitPhase parameter
    % (theta/theta, theta/phi, phi/theta, phi/phi) is effectively row major
    polarizationMatrix = [1; sqrt(1./kappa); sqrt(1./kappa); 1] .* exp(1i*Phi([1 3 2 4],:,:));

end

% TR 38.900 Section 7.5 step 11
% polarization matrix in Equation 7.5-29
function polarizationMatrix = getLOSPolarizationMatrix(Phi)

    polarizationMatrix = [1; 0; 0; -1;] .* exp(1i*Phi(:,:,1));
    
end

function theta = limitZenithAngles(theta)

    theta = mod(theta,360);
    theta(theta>180) = 360 - theta(theta>180);

end

function [D,t] = getClusterDimensionality(model,insize)

    t = getSampleTimes(model,insize); % sample times
    
    T = numel(t);                                            % number of time-domain channel samples
    P = insize(2);                                           % number of transmit antennas
    R = size(model.ReceiveAntennaArray.ElementPositions,2);  % number of receive antennas
    
    D = [T 1 P R];
    
end

function [t,F_cg] = getSampleTimes(model,insize)

    if (isfield(model,'SampleTimes'))
        
        % if filtering is off, or filtering is on and SampleTimes are
        % present, use SampleTimes to define time points
        t = model.SampleTimes(:);
        F_cg = [];
        
    else
        
        if (model.SampleDensity == Inf)
            
            % set coefficient generation sampling rate F_cg equal to input
            % sampling rate (SampleRate is mandatory if filtering is on)
            % and set the number of samples to the input waveform length
            F_cg = model.SampleRate;
            T = insize(1);
            
        else
            
            % set coefficient generation sampling rate F_cg according to
            % Doppler frequency and SampleDensity parameter
            F_cg = model.MaximumDopplerShift * 2 * model.SampleDensity;
            
            % calculate the number of time points at which to take channel
            % samples, including at least an extra half sample period
            % beyond the end of the waveform to allow for zero order hold
            % "balancing" (see zohChannelMatrixTimes in
            % applyCDLChannelMatrix)
            T = ceil(insize(1) * F_cg / model.SampleRate + 0.5);
            
        end
        
        % calculate vector of time points at sampling rate F_cg
        % (with special case for F_cg=0 for zero Doppler (T=1,t=0))
        if (T==1)
            t = 0;
        else
            t = (0:(T-1)).'/F_cg;
        end
        
        % add coefficients generation time offset
        t = t + model.InitialTime;

    end
    
end

% get field s.(f) or default value 'v' is 'f' is not a field of 's'
function v = getfieldordefault(s,f,v)
    
    if (isfield(s,f))
        v = s.(f);    
    end
    
end

% TR 38.900 Section 7.7.5.1 Scaling of angles
function phi = scaleCDLAngles(model,pdp,phi,AS_field,mu_field,postfn)

    AS_model = calculateRMSAngleSpread(pdp,phi);
    mu_model = mean(phi(:));

    AS_desired = getfieldordefault(model,AS_field,AS_model);
    mu_desired = getfieldordefault(model,mu_field,mu_model);

    if (AS_model==0)
        phi = postfn(phi - mu_model + mu_desired);
    else
        phi = postfn((AS_desired / AS_model * (phi - mu_model)) + mu_desired);
    end

end

% Calculate RMS angle spread
function AS = calculateRMSAngleSpread(pdp,phi)

    % named column subscripts
    [~,power] = nr5g.internal.namedDelayProfileColumns();

    % linear cluster powers
    P = 10.^(pdp(:,power)/10);
    
    % establish ray powers from cluster powers
    P = repmat(P,1,size(phi,2))/size(phi,2);
    
    % calculate RMS angle spread via circular standard deviation technique
    AS = circularStandardDeviation(phi,P);
    
end

function AS = circularStandardDeviation(phi,P)

    AS = sqrt(-2*log(abs(sum(exp(1i*phi(:)).*P(:))/sum(P(:)))));

end
