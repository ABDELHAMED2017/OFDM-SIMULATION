## FILES

- **OFDM_SIM.m** --  The main function.
- **ofdm_parameters.m** -- Configure some simulation parameters which includes some key ones as below
  - `ifft_size`, IFFT/FFT size
  - `carrier_count`, sub carrier numbers
  - `symb_size`,  modulation constellation mapping bit size, it relats to modulation scheme such as **BPSK, QPSK, 16QPSK and 256PSK**
  - `clipping`, OFDM amplitude clipping factor
  - `SNR_dB`, SNR with dB as unit
  - `guard_time`, guard time for between symbols
  - `spacing`, data subcarrier spacing

## ABOUT

The project is created for OFDM simulation. 

We transfer an image via the wireless system and test the performance of image transportation quality with different parameters. Furthur to explore the project.

## HISTORY

|      DATE      |          MODIFICATION           |
| :------------: | :-----------------------------: |
| *Dec 25, 2016* |        Creat the project        |
| *Dec 27, 2016* | described every single function |