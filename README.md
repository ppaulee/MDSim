MolSim
===

The Molecular Dynamics teaching code.

https://github.com/ppaulee/MDSim

**Build Instructions:**
* *(without Doxygen)*  
mkdir buildDir  
cd ./buildDir  
cmake -S ../ -B .  
make  
./MolSim ../eingabe-sonne.txt <end_time> <delta_t>

* *(with Doxygen)*  
mkdir buildDir  
cd ./buildDir  
cmake -S ../ -B . -DBUILD_DOXY=ON   
make  
make doxyDoc

Group ID : B
Member: Paul Ungermann, Michael Ethan, Jonas Hägele

Simulation  :
https://youtu.be/Q5dOevq1hPk