clear;
clc;
close all;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3e8;

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
R = 100;
v = -50;

%% FMCW Waveform Generation

%Operating carrier frequency of Radar 
fc = 77e9;             %carrier freq
maxRange = 200;
rangeResolution = 1;
maxVelocity = 100;

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
B = c / (2 * rangeResolution);  % aka Bsweep; formula from lesson 3.1

% Lesson 3.2 suggested that Tchirp should be "at least 5-6 times larger
% than 2 * maxRange / c.  The "at least" part is so that you don't have to
% worry about chirp boundary effects.  However, the actual value of Tchirp
% should presumably be related to the velocity range or resolution or
% something?  In short, I have NO IDEA why this is the unique right answer.
Tchirp = 5.5 * (2 * maxRange / c);

slope = B / Tchirp;  % Lesson 2.4
                                                      
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd = 128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr = 1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t = linspace(0, Nd*Tchirp, Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx = zeros(1,length(t)); % transmitted signal
Rx = zeros(1,length(t)); % received signal
Mix = zeros(1,length(t)); % beat signal

%Similar vectors for range_covered and time delay.
r_t = zeros(1,length(t));
td = zeros(1,length(t));

% amplitude; the units shouldn't matter for this exercise, since we're not 
% looking at noise
A = 1;

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)-1        
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R + v * t(i);
    
    %For each time sample we need update the transmitted and
    %received signal.
    
    % Note, we should actually be restarting at each chirp -- the frequency
    % shouldn't be linearly increasing through all time! However, that will
    % introduce artifacts at each chirp boundary, which this gets to
    % conveniently avoid.
    Tx(i) = A * cos(2*pi*(fc * t(i) + 1/2 * slope * t(i) * t(i)));
    
    % We're assuming that light travel is instantaneous for our purposes,
    % so that the delay is the round-trip time to the car at time t(i)
    tau = 2 * r_t(i) / c;
    Rx(i) = A * cos(2*pi*(fc * (t(i)-tau) + 1/2 * slope * (t(i)-tau) * (t(i)-tau)));
    
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    % Note, this would work much more cleanly if we were working with
    % complex exponentials; as it is, we're introducing high-frequency 
    % noise by doing this.  (Maybe we're doing it so that we have a noise
    % level to work with?)
    Mix(i) = Tx(i) * Rx(i);
end

%% RANGE MEASUREMENT
 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
% This seems to be a carry-over from the next task, we do NOT need to
% reshape here.
Mix_reshaped = reshape(Mix, [Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
% Not sure why we care about normalizing here, since we will want the argmax
% which doesn't care about the scale
% Note that fft(Mix,Nr) is the Fourier transform of the first chirp on its
% own
beat_fft = fft(Mix, Nr) / Nr;

 % *%TODO* :
% Take the absolute value of FFT output
beat_fft = abs(beat_fft);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
beat_fft = beat_fft(1:Nr/2 + 1);

% find the range estimate as the argmax of this signal, note we need to
% shift by 1 because beat_fft starts with frequency (range) of 0!
[~, range_estimate] = max(beat_fft); range_estimate = range_estimate - 1;

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1);
 % *%TODO* :
 % plot FFT output
plot(0:Nr/2, beat_fft);
xline(range_estimate);
title(sprintf('Range estimate %.02f', range_estimate));
axis ([0 200 0 1]);

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map. You will implement CFAR on the generated RDM

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix = reshape(Mix, [Nr, Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix, Nr, Nd);

% Taking just one side of signal from Range dimension.
% (This wasn't explained in the course, but the reason we need to do it 
% even in the noiseless case is that our mixing introduces high-frequency
% components! If we were working with complex signals, we wouldn't need to
% filter out high frequencies just yet.)
sig_fft2 = sig_fft2(1:Nr/2, 1:Nd);
% This is REALLY confusing. What seems to be happening is that we want to
% fftshift the velocity dimension (because we want 0 at the center), and we
% can afford to fftshift the range dimension because there's enough samples
% there that even after discarding the high-frequency half, the last half 
% is unphysical (we're not detecting ranges above 200, and Nr > 200*4).
% We're declaring these to correspond to the equally-unphysical negative
% ranges, instead of zooming in on the range from 0 to 200 as in the
% previous exercise (and as in the RDM in lesson 3.5!)
sig_fft2 = fftshift(sig_fft2);

RDM = abs(sig_fft2);
RDM = 10*log10(RDM);

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);

% I find the definition below easier to understand, in particular it
% highlights that the steps in the range axis are of size 1 (or, well,
% almost).
% range_axis = linspace(-200,200,Nr/2) * ((Nr/2)/400);
range_axis = linspace(-Nr/4,Nr/4,Nr/2);

figure,surf(doppler_axis,range_axis,RDM);
% figure, pcolor(doppler_axis, range_axis, RDM);
% shading interp; colormap jet;
xlabel('Doppler velocity estimate');
ylabel('Range estimate');
title('Range Doppler Map; please ignore any ranges <0 or >200 as unphysical');

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
training_r = 5;
training_d = 5;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
guard_r = 3;
guard_d = 5;

% *%TODO* :
% offset the threshold by SNR value in dB
% um, we don't have noise as such?? That said, the peak seems to be about
% 15dB above the noise, so 10dB should be a reasonable offset.
threshold_dB = 10;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
% why would I want it 1x1??
%noise_level = zeros(1,1);
noise_level = zeros(size(RDM));
cfar = zeros(size(RDM));

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

for r = 1 + training_r + guard_r : Nr/2 - training_r - guard_r
    for d = 1 + training_d + guard_d : Nd - training_d - guard_d
        patch = RDM(r-training_r-guard_r : r+training_r+guard_r, d-training_d-guard_d : d+training_d+guard_d);
        patch = db2pow(patch);
        % make sure not to include the central (2*guard_r+1)*(2*guard_d+1)
        % rectangle in noise measurement
        mask = true(size(patch));
        mask(1 + training_r : end-training_r, 1 + training_d : end-training_d) = false;
        noise_level(r,d) = pow2db(mean(patch(mask)));
        
        cfar(r,d) = RDM(r,d) > noise_level(r,d) + threshold_dB;
    end
end

% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 

% actually, I took care of that already in the set-up.

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,cfar);
colorbar;


 
 