# SFND_Radar_Target_Generation_and_Detection
Udacity project Radar Target Generation and Detection
## FMCW Waveform Design
|parameter|value|line|
|-|-|-|
|Chirp frquency bandwidth B|150 MHz|66|
|Chirp time Tchirp, us|7.3|43|
|slope Hz²|2e13|68|

```
Tchirp = 5.5*2*Rmax/c;
B = c / (2*Rres);
slope = B / Tchirp;
```
## Simulation Loop
- Simulate Target movement (line 76)
```
r_t(i) = r0 + v*t(i);
```
- Calculate the beat or mixed signal for every timestamp (line 82:89)
```
Tx(i) = cos(2*pi()*(t(i)*(fc + slope*t(i)/2)));
Rx(i) = cos(2*pi()*((t(i)-td(i))*(fc + slope*(t(i)-td(i))/2)));
Mix(i) = Tx(i).* Rx(i);
```
## Range FFT (1st FFT)
- Implement the Range FFT on the Beat or Mixed Signal and plot the result (line 97:127).
```
Mix = reshape(Mix, [Nr*Nd, 1]);
Mix_fft = fft(Mix, Nr)./Nr;
Mix_fft = abs(Mix_fft);
P1 = Mix_fft(1:Nr/2+1);
```
![Image of range](https://github.com/AnastaciaVolkova/SFND_Radar_Target_Generation_and_Detection/blob/master/images/range.jpg)

## 2D CFAR
- Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map (line 202:222)
Output of 2D FFT
<br>
![Image of range/doppler](https://github.com/AnastaciaVolkova/SFND_Radar_Target_Generation_and_Detection/blob/master/images/range_doppler.jpg)
### CFAR implementation steps:
- Choose parameters for CFAR:
  - Choose the size of training band in range dimension (Tr=8)
  - Choose the size of training band in doppler dimension  (Td=4) 
  - Choose the size of guard band in range dimension (Gr=4)
  - Choose the size of guard band in the doppler dimension (Gd=2)
  - Define **offset** in db as 10*log10(5). Power of signal should be 5 times higher than average of training cell power. Values are chosen as power of 2. Values for Doppler dimension are less twice than values from range dimension as Nr=512, Nd=128.
  - Initialize threshold_block as matrix NrxNd zeros. It helps to avoid special zeroing of corners. As CUT has a shift of T+G.
- Move a window inside range-doppler db matrix. The center of window - cell under test. It is surrounded by guard cells, which are surrounded by training cells. 
  - Window is defined (line 205)
    ```
    training_cells = db2pow(RDM([i-Gr-Tr:i+Gr+Tr], [j-Gd-Td:j+Gd+Td]));
    ```
  - In order to avoid contribution of  cell under test and guard cells to average value, make them zero (line 208)
    ```
    training_cells( [i-Gr: i+Gr], [j-Gd: j+Gd] ) = 0;
    ```
  - Compute average of power of training cells inside window (line 211)
    ```
    noise_level(i, j) = sum(training_cells, 'all')/total_train;
    ```
- Apply threshold. Remember only those cells, which have power greater than threshold (line 214)
    ```
    threshold = offset + pow2db(noise_level(i,j));

    if RDM(i, j) > threshold
        threshold_block(i, j) = 1;
    else
        threshold_block(i, j) = 0;
    end
    ```
- Threshold peak coincides with target relative speed 16.4 m/s and distance to target 86 m.
<br>
![Image of threshold](https://github.com/AnastaciaVolkova/SFND_Radar_Target_Generation_and_Detection/blob/master/images/range_doppler_threshold.jpg)





