%% Introduction to MIMO Systems
% This example shows Multiple-Input-Multiple-Output (MIMO) systems, which
% use multiple antennas at the transmitter and receiver ends of a wireless
% communication system. MIMO systems are increasingly being adopted in
% communication systems for the potential gains in capacity they realize
% when using multiple antennas. Multiple antennas use the spatial dimension
% in addition to the time and frequency ones, without changing the
% bandwidth requirements of the system.
%
% For a generic communications link, this example focuses on transmit
% diversity in lieu of traditional receive diversity. Using the flat-fading
% Rayleigh channel, it illustrates the concept of Orthogonal Space-Time
% Block Coding, which is employable when multiple transmitter antennas are
% used. It is assumed here that the channel undergoes independent fading
% between the multiple transmit-receive antenna pairs.
%
% For a chosen system, it also provides a measure of the performance
% degradation when the channel is imperfectly estimated at the receiver, 
% compared to the case of perfect channel knowledge at the receiver.

% Copyright 2006-2014 The MathWorks, Inc.


%% PART 1: Transmit Diversity vs. Receive Diversity
%
% Using diversity reception is a well-known technique to mitigate the
% effects of fading over a communications link. However, it has mostly been
% relegated to the receiver end. In [ <#13 1> ], Alamouti proposes a
% transmit diversity scheme that offers similar diversity gains, using
% multiple antennas at the transmitter. This was conceived to be more
% practical as, for example, it would only require multiple antennas at the
% base station in comparison to multiple antennas for every mobile in a
% cellular communications system.
%
% This section highlights this comparison of transmit vs. receive diversity
% by simulating coherent binary phase-shift keying (BPSK) modulation over
% flat-fading Rayleigh channels. For transmit diversity, we use two
% transmit antennas and one receive antenna (2x1 notationally), while for
% receive diversity we employ one transmit antenna and two receive antennas
% (1x2 notationally).
%
% The simulation covers an end-to-end system showing the encoded and/or
% transmitted signal, channel model, and reception and demodulation of the
% received signal. It also provides the no-diversity link (single transmit-
% receive antenna case) and theoretical performance of second-order
% diversity link for comparison. It is assumed here that the channel is
% known perfectly at the receiver for all systems. We run the simulation
% over a range of Eb/No points to generate BER results that allow us to
% compare the different systems.

%%
% We start by defining some common simulation parameters
frmLen = 1000;       % frame length
numPackets = 1000;  % number of packets
EbNo = -10:2:10;      % Eb/No varying to 20 dB
N = 2;              % maximum number of Tx antennas
M = 2;              % maximum number of Rx antennas

%%
% and set up the simulation.

% Create comm.BPSKModulator and comm.BPSKDemodulator System objects
P = 2;				% modulation order
hMod   = comm.BPSKModulator;
hDemod = comm.BPSKDemodulator('OutputDataType','double');

% Create comm.OSTBCEncoder and comm.OSTBCCombiner System objects
hAlamoutiEnc = comm.OSTBCEncoder;
hAlamoutiDec = comm.OSTBCCombiner;
hA2amoutiDec = comm.OSTBCCombiner(2,2);

% Create two comm.AWGNChannel System objects for one and two receive
% antennas respectively. Set the NoiseMethod property of the channel to
% 'Signal to noise ratio (Eb/No)' to specify the noise level using the
% energy per bit to noise power spectral density ratio (Eb/No). The output
% of the BPSK modulator generates unit power signals; set the SignalPower
% property to 1 Watt.
hAWGN1Rx = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)',...
                            'SignalPower', 1);
hAWGN2Rx = clone(hAWGN1Rx);

% Create comm.ErrorRate calculator System objects to evaluate BER.
hErrorCalc1 = comm.ErrorRate;
hErrorCalc2 = comm.ErrorRate;
hErrorCalc3 = comm.ErrorRate;

% Since the comm.AWGNChannel System objects as well as the RANDI function
% use the default random stream, the following commands are executed so
% that the results will be repeatable, i.e., same results will be obtained
% for every run of the example. The default stream will be restored at the
% end of the example.
s = RandStream.create('mt19937ar', 'seed',55408);
prevStream = RandStream.setGlobalStream(s);

% Pre-allocate variables for speed
H = zeros(frmLen, N, M);
ber_noDiver  = zeros(3,length(EbNo));
ber_Alamouti = zeros(3,length(EbNo));
ber_MaxRatio = zeros(3,length(EbNo));
ber_thy2     = zeros(1,length(EbNo));
%%

% Set up a figure for visualizing BER results
h = gcf; 
grid on; 
hold on;
ax = gca;
ax.YScale = 'log';
xlim([EbNo(1), EbNo(end)]);
ylim([1e-4 1]);
xlabel('Eb/No (dB)');
ylabel('BER'); 
h.NumberTitle = 'off';
h.Renderer = 'zbuffer';
h.Name = 'Transmit vs. Receive Diversity';
title('Transmit vs. Receive Diversity');

% Loop over several EbNo points
for idx = 1:length(EbNo)
    reset(hErrorCalc1);
    reset(hErrorCalc2);
    reset(hErrorCalc3);
    % Set the EbNo property of the AWGNChannel System objects
    hAWGN1Rx.EbNo = EbNo(idx); 
    hAWGN2Rx.EbNo = EbNo(idx); 
    % Loop over the number of packets
    for packetIdx = 1:numPackets
        % Generate data vector per frame 
        data = randi([0 P-1], frmLen, 1); 
        
        % Modulate data
        modData = step(hMod, data);     

        % Alamouti Space-Time Block Encoder
        encData = step(hAlamoutiEnc, modData);
        
        % Create the Rayleigh distributed channel response matrix
        %   for two transmit and two receive antennas
        H(1:N:end, :, :) = (randn(frmLen/2, N, M) + ...
                         1i*randn(frmLen/2, N, M))/sqrt(2);
        %   assume held constant for 2 symbol periods
        H(2:N:end, :, :) = H(1:N:end, :, :);
        
        % Extract part of H to represent the 1x1, 2x1 and 1x2 channels
        H11 = H(:,1,1);
        H21 = H(:,:,1)/sqrt(2);
        H12 = squeeze(H(:,1,:));
        
        % Pass through the channels
        chanOut11 = H11 .* modData;
        
        chanOut22_11 = H(:,1,1) .* encData(:,1);
        chanOut22_21 = H(:,2,1) .* encData(:,2);
        chanOut22_12 = H(:,1,2) .* encData(:,1);
        chanOut22_22 = H(:,2,2) .* encData(:,2);
        chanOut22(:,1) = chanOut22_11 + chanOut22_21;
        chanOut22(:,2) = chanOut22_12 + chanOut22_22; 
        
        chanOut21 = sum(H21.* encData, 2);
%         chanOut12 = H12 .* repmat(modData, 1, 2);
        chanOut12(:,1) = chanOut22(:,1) + chanOut22(:,2);
        chanOut12(:,2) = chanOut12(:,1);
        
        % Add AWGN
        rxSig11 = step(hAWGN1Rx, chanOut11);
        rxSig21 = step(hAWGN1Rx, chanOut21);
        rxSig12 = step(hAWGN2Rx, chanOut12);
        rxSig22 = step(hAWGN2Rx, chanOut22);
              
        % Alamouti Space-Time Block Combiner
        decData = step(hAlamoutiDec, rxSig21, H21);
        decData1 = step(hA2amoutiDec, rxSig22, H);
        decData2 = step(hA2amoutiDec, rxSig12, H);

        % ML Detector (minimum Euclidean distance)
        demod11 = step(hDemod, rxSig11.*conj(H11));
        demod21 = step(hDemod, decData);
%         demod12 = step(hDemod, sum(rxSig12.*conj(H12), 2));
        demod12 = step(hDemod, decData2);
        demod22 = step(hDemod, decData1);
        
        % Calculate and update BER for current EbNo value
        %   for uncoded 1x1 system
        ber_noDiver(:,idx)  = step(hErrorCalc1, data, demod11);
        %   for Alamouti coded 2x2 system
        ber_Alamouti(:,idx) = step(hErrorCalc2, data, demod22);
        
        %   for Maximal-ratio combined 1x2 system
        ber_MaxRatio(:,idx) = step(hErrorCalc3, data, demod12);
        
    end % end of FOR loop for numPackets

    % Calculate theoretical second-order diversity BER for current EbNo
    ber_thy2(idx) = berfading(EbNo(idx), 'psk', 2, 2);
    
    % Plot results
    semilogy(EbNo(1:idx), ber_noDiver(1,1:idx), 'r*', ...
             EbNo(1:idx), ber_Alamouti(1,1:idx), 'go', ...
             EbNo(1:idx), ber_MaxRatio(1,1:idx), 'bs', ...
             EbNo(1:idx), ber_thy2(1:idx), 'm');
    legend('No Diversity (1Tx, 1Rx)', 'Alamouti (2Tx, 1Rx)',...
           'Maximal-Ratio Combining (1Tx, 2Rx)', ...
           'Theoretical 2nd-Order Diversity');
    
    drawnow;
end  % end of for loop for EbNo
 
% Perform curve fitting and replot the results
fitBER11 = berfit(EbNo, ber_noDiver(1,:));
fitBER21 = berfit(EbNo, ber_Alamouti(1,:));
fitBER12 = berfit(EbNo, ber_MaxRatio(1,:));
semilogy(EbNo, fitBER11, 'r', EbNo, fitBER21, 'g', EbNo, fitBER12, 'b');
hold off;

% Restore default stream
RandStream.setGlobalStream(prevStream);

%%
% The transmit diversity system has a computation complexity very similar
% to that of the receive diversity system.
%
% The resulting simulation results show that using two transmit antennas
% and one receive antenna provides the same diversity order as the
% maximal-ratio combined (MRC) system of one transmit antenna and two
% receive antennas.
%
% Also observe that transmit diversity has a 3 dB disadvantage when
% compared to MRC receive diversity. This is because we modelled the total
% transmitted power to be the same in both cases. If we calibrate the
% transmitted power such that the received power for these two cases is the
% same, then the performance would be identical. The theoretical
% performance of second-order diversity link matches the transmit diversity
% system as it normalizes the total power across all the diversity
% branches.
%
% The accompanying functional scripts, <matlab:edit('mrc1m.m') MRC1M.m> and
% <matlab:edit('ostbc2m.m') OSTBC2M.m> aid further exploration for the
% interested users.


%% PART 2: Space-Time Block Coding with Channel Estimation
%
% Building on the theory of orthogonal designs, Tarokh et al. [ <#13 2> ]
% generalized Alamouti's transmit diversity scheme to an arbitrary number
% of transmitter antennas, leading to the concept of Space-Time Block
% Codes. For complex signal constellations, they showed that Alamouti's
% scheme is the only full-rate scheme for two transmit antennas.
%
% In this section, we study the performance of such a scheme with two
% receive antennas (i.e., a 2x2 system) with and without channel
% estimation. In the realistic scenario where the channel state information
% is not known at the receiver, this has to be extracted from the received
% signal. We assume that the channel estimator performs this using
% orthogonal pilot signals that are prepended to every packet [ <#13 3> ].
% It is assumed that the channel remains unchanged for the length of the
% packet (i.e., it undergoes slow fading).
%
% A simulation similar to the one described in the previous section is
% employed here, which leads us to estimate the BER performance for a
% space-time block coded system using two transmit and two receive
% antennas.

%%
% Again we start by defining the common simulation parameters
frmLen = 100;           % frame length
maxNumErrs = 300;       % maximum number of errors
maxNumPackets = 3000;   % maximum number of packets
EbNo = 0:2:12;          % Eb/No varying to 12 dB
N = 2;                  % number of Tx antennas
M = 2;                  % number of Rx antennas
pLen = 8;               % number of pilot symbols per frame
W = hadamard(pLen);
pilots = W(:, 1:N);     % orthogonal set per transmit antenna

%%
% and set up the simulation.

% Create a comm.MIMOChannel System object to simulate the 2x2 spatially
% independent flat-fading Rayleigh channel
hChan = comm.MIMOChannel('MaximumDopplerShift', 0, ...
                         'SpatialCorrelation',  false, ...
                         'NumTransmitAntennas', N, ...
                         'NumReceiveAntennas',  M, ...
                         'PathGainsOutputPort', true);

% Change the NumReceiveAntennas property value of the hAlamoutiDec System
% object to M that is 2
release(hAlamoutiDec);
hAlamoutiDec.NumReceiveAntennas = M;

% Release the hAWGN2Rx System object
release(hAWGN2Rx);

% Set the global random stream for repeatability
s = RandStream.create('mt19937ar', 'seed',55408);
prevStream = RandStream.setGlobalStream(s);

% Pre-allocate variables for speed
HEst = zeros(frmLen, N, M);
ber_Estimate = zeros(3,length(EbNo));
ber_Known    = zeros(3,length(EbNo));

%%

% Set up a figure for visualizing BER results
clf(h); 
grid on; 
hold on;
ax = gca;
ax.YScale = 'log';
xlim([EbNo(1), EbNo(end)]);
ylim([1e-4 1]);
xlabel('Eb/No (dB)');
ylabel('BER'); 
h.NumberTitle = 'off';
h.Name = 'Orthogonal Space-Time Block Coding';
h.Renderer = 'zbuffer';
title('Alamouti-coded 2x2 System');

% Loop over several EbNo points
for idx = 1:length(EbNo)
    reset(hErrorCalc1);
    reset(hErrorCalc2);
    hAWGN2Rx.EbNo = EbNo(idx); 

    % Loop till the number of errors exceed 'maxNumErrs'
    % or the maximum number of packets have been simulated
    while (ber_Estimate(2,idx) < maxNumErrs) && ...
          (ber_Known(2,idx) < maxNumErrs) && ...
          (ber_Estimate(3,idx)/frmLen < maxNumPackets)
        % Generate data vector per frame 
        data = randi([0 P-1], frmLen, 1);
        
        % Modulate data
        modData = step(hMod, data);           
        
        % Alamouti Space-Time Block Encoder
        encData = step(hAlamoutiEnc, modData);

        % Prepend pilot symbols for each frame
        txSig = [pilots; encData];

        % Pass through the 2x2 channel        
        reset(hChan);
        [chanOut, H] = step(hChan, txSig);

        % Add AWGN
        rxSig = step(hAWGN2Rx, chanOut);

        % Channel Estimation
        %   For each link => N*M estimates
        HEst(1,:,:) = pilots(:,:).' * rxSig(1:pLen, :) / pLen;
        %   assume held constant for the whole frame
        HEst = HEst(ones(frmLen, 1), :, :);

        % Combiner using estimated channel
        decDataEst = step(hAlamoutiDec, rxSig(pLen+1:end,:), HEst);

        % Combiner using known channel
        decDataKnown = step(hAlamoutiDec, rxSig(pLen+1:end,:), ...
                            squeeze(H(pLen+1:end,:,:,:)));
        
        % ML Detector (minimum Euclidean distance)
        demodEst   = step(hDemod, decDataEst);      % estimated
        demodKnown = step(hDemod, decDataKnown);    % known      
        
        % Calculate and update BER for current EbNo value
        %   for estimated channel
        ber_Estimate(:,idx) = step(hErrorCalc1, data, demodEst);
        %   for known channel
        ber_Known(:,idx)    = step(hErrorCalc2, data, demodKnown);
    
    end % end of FOR loop for numPackets

    % Plot results
    semilogy(EbNo(1:idx), ber_Estimate(1,1:idx), 'ro');
    semilogy(EbNo(1:idx), ber_Known(1,1:idx), 'g*');
    legend(['Channel estimated with ' num2str(pLen) ' pilot symbols/frame'],...
           'Known channel');
    drawnow;
end  % end of for loop for EbNo

% Perform curve fitting and replot the results
fitBEREst   = berfit(EbNo, ber_Estimate(1,:));
fitBERKnown = berfit(EbNo, ber_Known(1,:));
semilogy(EbNo, fitBEREst, 'r', EbNo, fitBERKnown, 'g'); 
hold off;

% Restore default stream
RandStream.setGlobalStream(prevStream);

%%
% For the 2x2 simulated system, the diversity order is different than that
% seen for either 1x2 or 2x1 systems in the previous section.
%
% Note that with 8 pilot symbols for each 100 symbols of data, channel
% estimation causes about a 1 dB degradation in performance for the
% selected Eb/No range. This improves with an increase in the number of
% pilot symbols per frame but adds to the overhead of the link. In this
% comparison, we keep the transmitted SNR per symbol to be the same in both
% cases.
%
% The accompanying functional script, <matlab:edit('ostbc2m_e.m')
% OSTBC2M_E.m> aids further experimentation for the interested users.


%% PART 3: Orthogonal Space-Time Block Coding and Further Explorations
%
% In this final section, we present some performance results for orthogonal
% space-time block coding using four transmit antennas (4x1 system) using a
% half-rate code, G4, as per [ <#13 4> ].
%
% We expect the system to offer a diversity order of 4 and will compare it
% with 1x4 and 2x2 systems, which have the same diversity order also. To
% allow for a fair comparison, we use quaternary PSK with the half-rate G4
% code to achieve the same transmission rate of 1 bit/sec/Hz.
%
% These results take some time to generate on a single core. If you do not
% have Parallel Computing Toolbox(TM) (PCT) installed, we load the results
% from a prior simulation. The functional script <matlab:edit('ostbc4m.m')
% ostbc4m.m> is included, which, along with  <matlab:edit('mrc1m.m')
% mrc1m.m> and <matlab:edit('ostbc2m.m') ostbc2m.m>, was used to generate
% these results. If PCT is installed, these simulations are performed in
% parallel. In this case the functional scripts
% <matlab:edit('ostbc4m_pct.m') ostbc4m_pct.m>, <matlab:edit('mrc1m_pct.m')
% mrc1m_pct.m> and <matlab:edit('ostbc2m_pct.m') ostbc2m_pct.m> are used.
% The user is urged to use these scripts as a starting point to study other
% codes and systems.

[licensePCT,~] = license( 'checkout' , 'Distrib_Computing_Toolbox');
if ( licensePCT && ~isempty(ver('distcomp')))
    EbNo = 0:2:20;
    [ber11, ber14, ber22, ber41] = mimoOSTBCWithPCT(100,4e3,EbNo);
else
    load ostbcRes.mat;
end


% Set up a figure for visualizing BER results
clf(h);
grid on; 
hold on; 
h.Renderer  = 'zbuffer';
ax = gca;
ax.YScale = 'log';
xlim([EbNo(1), EbNo(end)]);
ylim([1e-5 1]);
xlabel('Eb/No (dB)');
ylabel('BER');
h.NumberTitle = 'off';
h.Name = 'Orthogonal Space-Time Block Coding(2)';
title('G4-coded 4x1 System and Other Comparisons');

% Theoretical performance of fourth-order diversity for QPSK
BERthy4 = berfading(EbNo, 'psk', 4, 4);

% Plot results
semilogy(EbNo, ber11, 'r*', EbNo, ber41, 'ms', EbNo, ber22, 'c^', ...
         EbNo, ber14, 'ko', EbNo, BERthy4, 'g');
legend('No Diversity (1Tx, 1Rx), BPSK', 'OSTBC (4Tx, 1Rx), QPSK', ...
       'Alamouti (2Tx, 2Rx), BPSK', 'Maximal-Ratio Combining (1Tx, 4Rx), BPSK', ...
       'Theoretical 4th-Order Diversity, QPSK');

% Perform curve fitting
fitBER11 = berfit(EbNo, ber11);
fitBER41 = berfit(EbNo(1:9), ber41(1:9));
fitBER22 = berfit(EbNo(1:8), ber22(1:8));
fitBER14 = berfit(EbNo(1:7), ber14(1:7));
semilogy(EbNo, fitBER11, 'r', EbNo(1:9), fitBER41, 'm', ...
         EbNo(1:8), fitBER22, 'c', EbNo(1:7), fitBER14, 'k');
hold off;

%% 
% As expected, the similar slopes of the BER curves for the 4x1, 2x2 and
% 1x4 systems indicate an identical diversity order for each system.
%
% Also observe the 3 dB penalty for the 4x1 system that can be attributed
% to the same total transmitted power assumption made for each of the three
% systems. If we calibrate the transmitted power such that the received
% power for each of these systems is the same, then the three systems would
% perform identically. Again, the theoretical performance matches the
% simulation performance of the 4x1 system as the total power is normalized
% across the diversity branches.
   
%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('mrc1m.m') mrc1m.m>
% * <matlab:edit('ostbc2m.m') ostbc2m.m>
% * <matlab:edit('ostbc4m.m') ostbc4m.m>
% * <matlab:edit('mimoOSTBCWithPCT.m') mimoOSTBCWithPCT.m>
% * <matlab:edit('mrc1m_pct.m') mrc1m_pct.m>
% * <matlab:edit('ostbc2m_pct.m') ostbc2m_pct.m>
% * <matlab:edit('ostbc4m_pct.m') ostbc4m_pct.m>


%% References
% # S. M. Alamouti, "A simple transmit diversity technique for wireless
% communications", IEEE(R) Journal on Selected Areas in
% Communications, Vol. 16, No. 8, Oct. 1998, pp. 1451-1458.
% # V. Tarokh, H. Jafarkhami, and A.R. Calderbank, "Space-time block codes
% from orthogonal designs", IEEE Transactions on Information Theory,
% Vol. 45, No. 5, Jul. 1999, pp. 1456-1467.
% # A.F. Naguib, V. Tarokh, N. Seshadri, and A.R. Calderbank, "Space-time
% codes for high data rate wireless communication: Mismatch
% analysis", Proceedings of IEEE International Conf. on
% Communications, pp. 309-313, June 1997.
% # V. Tarokh, H. Jafarkhami, and A.R. Calderbank, "Space-time block codes
% for wireless communications: Performance results", IEEE Journal on
% Selected Areas in Communications, Vol. 17,  No. 3, Mar. 1999, pp.
% 451-460.

displayEndOfDemoMessage(mfilename)
