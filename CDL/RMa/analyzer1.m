  % Example 3:
  %   Create a resource grid for 64 antennas, fill it with QPSK symbols
  %   and perform LTE OFDM modulation. Pass the resulting waveform
  %   through a 64-by-4 channel and plot the received waveform spectrum:

  enb.NDLRB = 25;
  enb.CyclicPrefix = 'Normal';
  grid = lteDLResourceGrid(enb,64);
  grid(:) = lteSymbolModulate(randi([0 1],numel(grid)*2,1),'QPSK');
  [txWaveform,txInfo] = lteOFDMModulate(enb,grid);

  cdl = nr5gCDLChannelTry;
  cdl.SampleRate = txInfo.SamplingRate;
  cdl.TransmitAntennaArray.Size = [2 4 2 2 2];
  cdl.TransmitAntennaArray.ElementSpacing = [0.5 0.5 2.0 1.0];
  cdl.ReceiveAntennaArray.Size = [2 1 2 1 1];
  
  % The antenna array elements are mapped to the waveform channels
  % (columns) in the order that a 5-D array of size
  % TransmitAntennaArray.Size or ReceiveAntennaArray.Size is linearly
  % indexed i.e. across the dimensions first to last. See the
  % TransmitAntennaArray or ReceiveAntennaArray property help for more
  % details.
  
  rxWaveform = cdl(txWaveform);

  analyzer = dsp.SpectrumAnalyzer('SampleRate',cdl.SampleRate);
  analyzer.Title = ['Received signal spectrum for ' cdl.DelayProfile];
  analyzer(rxWaveform);