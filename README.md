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
  ./MolSim ../input.txt -e 5 -s 0.0002 -a lj
####(with Doxygen):
* mkdir buildDir  
  cd buildDir  
  cmake -S ../ -B . -DBUILD_DOXY=ON   
  make  
  make doxyDoc

###Simulation:
https://youtu.be/Cdi7o-XbQ14
