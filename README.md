SKIF_beamlines

--- Install Spyder on your Linux, Window or Mac OSX machine https://docs.spyder-ide.org/installation.html

--- Install SRW code on your computer. Follow the instruction on the GitHub page https://github.com/ochubar/SRW . There are separate instructions for Window, Linux and MAc OSX systems.

--- As soon as you install the code and able to run Chubar's SRW examples (try it!) you are ready to adjust the environment for SKIF lib.

--- The last step to be done is to add PYTHONPATH of the SKIF_lib.py.
	The easiest way: Use Spyder Python Path manager to add SKIF_lib.py as your environment. Find the sign of the python on the toolbar is Spyder. Click it. Add the path to the SKIF_beamline project. It must be "directory you downloaded your SKIF_beamline folder/SKIF_XAS_beamline/
	
	The way full of pain and suffering:

	On Linux machine: 
		* In the terminal type: gedit ~/.bashrc
		* Add to this file one of these two lines (try it both!) 
		export PYTHONPATH=$PYTHONPATH:"directory you downloaded your SKIF_beamline folder/SKIF_XAS_beamline/:
		export PATH="directory you downloaded your SKIF_beamline folder/SKIF_XAS_beamline/:$PATH"
	
	On Window machine: 
		Please google it or use the "easiest way".

SKIF_lib.py is a script where you can find some functions that extend standard SRW library and will help you a lot. Use it wisely for SKIF development.

