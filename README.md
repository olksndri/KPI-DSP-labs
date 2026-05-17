# KPI-DSP-labs

Switching between labs done by changing variable in main.cpp: 

```C++
int lab = 1;
```

## KPI-DSP-lab 1 
This project branch implements Fast Fourier Transform (FFT), windowing functions, and spectral visualization for audio signals.

## KPI-DSP-lab 2  
This project branch implements software for calculation and visualization the mel-frequency cepstral coefficients of an audio signal.

## KPI-DSP-lab 3
This project branch implements simple voice activity detector (VAD) based on calculated MFCC.

## KPI-DSP-lab 4
This project branch implements JPEG based image compression algorithm. 

## KPI-DSP-lab 5 (Calculation and graphical work) 
This project branch implements Haar DWT input samples decomposition and visualisation.  


## System Requirements

The project is designed for **Linux** environments.

### Prerequisites
* **CMake** (3.10+)
* **Python 3** (with `venv` support)
* **C/C++ Compiler** (GCC or Clang)



## Installation and Setup

Run the setup script once to initialize the Python virtual environment and install visualization dependencies:

```bash
./scripts/setup.sh
```

## Usage
### Clean
To clean build files, run: 

```Bash
./scripts/clean.sh
```

### Build
To build the project, run:

```Bash
./scripts/build.sh
```

### Lab 1-3 Execution
To analyze an audio file, provide the path to a .wav file as an argument:

```Bash
./scripts/run.sh path/to/audio.wav
```

### Lab 4 Execution
To compress an image, provide the path to a .tif file as an argument:

```Bash
./scripts/run.sh path/to/image.tif
```


### Lab 5 (CGW) Execution
As simple, as that 

```Bash
./scripts/run.sh
```



