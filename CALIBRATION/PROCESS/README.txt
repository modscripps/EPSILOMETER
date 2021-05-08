Python programme to accurately measure the average Peak-to-Peak voltage of shear probes.

Install Xcode command tools, home-brew, python, matplotlib, scipy, numpy onto computer
xcode-select --install
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
export PATH=/usr/local/bin:/usr/local/share/python:$PATH
brew install python
easy_install pip
pip install numpy
brew install gfortran
pip install scipy
brew install pkg-config
pip install matplotlib
pip install git+git://github.com/matplotlib/matplotlib.git

Within folder of calibrator.py are the following files:
config.txt - to input all parameters
test_data/record_{serial number}.txt - in which all data from a single test is saved. This file will be overwritten if shear probe serial number is not changed in config.txt
README.txt
dwfcontstants.py - contains essential variables for Diligent python framework

Connect shear probe output signal to Channel 1 of Digilent Analog Discovery device.
Make sure Waveforms software is closed.
Set all parameters in config.txt, with the following format, note the spaces on either side of the equals sign: {Parameter Name} = {Parameter Value}
Save config.txt
Open Terminal, change to programme directory (cd /dir/Calibrator), and type "python calibrator.py"
Depending on what value nSamples has, this could take up to a few minutes.
Graphs will be printed to screen, and data will be written to text file in test_data folder

Error in P2P measurement is at 0.1% with order of 8 and cutoff frequency of 10Hz.
