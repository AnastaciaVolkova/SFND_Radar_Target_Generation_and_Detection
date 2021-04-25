clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz

% Max Range = 200m
Rmax = 200;

% Range Resolution = 1 m
Rres = 1;

% Max velocity = 100 m/s
Vmax  = 100; % m/s;

% Velocity resoltion
Vres = 3; % m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3e8;

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
r0 = 86; % m
v = 16.4; % m/s => 59km/h


%% FMCW Waveform Generation

% *%TODO* :
% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.


%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

Tchirp = 5.5*2*Rmax/c;

%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

B = c / (2*Rres);

slope = B / Tchirp;

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)                 
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = r0 + v*t(i);
    td(i) = 2*r_t(i)/c;
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi()*(t(i)*(fc + slope*t(i)/2)));
    Rx(i) = cos(2*pi()*((t(i)-td(i))*(fc + slope*(t(i)-td(i))/2)));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).* Rx(i);
end

%% RANGE MEASUREMENT

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr*Nd, 1]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Mix_fft = fft(Mix, Nr)./Nr;

 % *%TODO* :
% Take the absolute value of FFT output
Mix_fft = abs(Mix_fft);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
P1 = Mix_fft(1:Nr/2+1);

%plotting the range
 % *%TODO* :
 % plot FFT output 
 
figure ('Name','Range from First FFT')
fs = Nr*Nd/(Nd*Tchirp); % Hz sampling frequency
f = (0:fs/Nr:fs/2);
subplot(2,1,1)
plot(f/10^9, P1); 
xlabel("Beet Frequency, GHz");
range = c * f * Tchirp / (2*B);
subplot(2,1,2)
axis ([0 Rmax 0 1]);
plot(range, P1); 
xlabel("Target position, m");

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr, Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-Vmax,Vmax,Nd);
range_axis = linspace(-Rmax,Rmax,Nr/2)*((Nr/2)/(2*Rmax));
figure,surf(doppler_axis,range_axis,RDM);
xlabel("Doppler");
ylabel("range");

Nr = Nr / 2;
%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions
Tr = 8; % Size of training band in range dimension
Td = 4; % Size of training band in doppler dimension 
% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4; % Size of guard band in range dimension 
Gd = 2; % Size of guard band in the doppler dimension

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 10*log10(5);

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(Nr, Nd);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR
total_train = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1)-(1+2*Gd)*(1+2*Gr);

threshold_block = zeros(Nr,Nd);

for i=Tr+Gr+1:(Nr-Tr-Gr)    
    for j=Td+Gd+1:(Nd-Td-Gd)
        % Get copy training cells and guard cells.
        training_cells = db2pow(RDM([i-Gr-Tr:i+Gr+Tr], [j-Gd-Td:j+Gd+Td]));
        
        % Set to zero copy of guard cells.
        training_cells( [i-Gr: i+Gr], [j-Gd: j+Gd] ) = 0;
        
        % Get noise level.
        noise_level(i, j) = sum(training_cells, 'all')/total_train;
        
        % Compute threshold.
        threshold = offset + pow2db(noise_level(i,j));
        
        if RDM(i, j) > threshold
            threshold_block(i, j) = 1;
        else
            threshold_block(i, j) = 0;
        end
    end
end

% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
% Already done in as threshold_block = zeros(Nr,Nd);

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis, threshold_block);
xlabel("doppler, Hz");
ylabel("range, m");
colorbar;


 
 