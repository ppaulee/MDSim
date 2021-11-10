Group B {#mainpage}
====================
###Members:
* Paul Ungermann
* Michael Ethan
* Jonas Hägele

###Code
https://github.com/ppaulee/MDSim

###Build Instructions:
####(without Doxygen):
* mkdir buildDir  
  cd buildDir  
  cmake -S ../ -B .  
  make  
  ./MolSim ../eingabe-sonne.txt <end_time> <delta_t>
####(with Doxygen):
* mkdir buildDir  
  cd buildDir  
  cmake -S ../ -B . -DBUILD_DOXY=ON   
  make  
  make doxyDoc

###Simulation:
https://youtu.be/Q5dOevq1hPk