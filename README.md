Group B {#mainpage}
====================
###Members:
* Paul Ungermann
* Michael Ethan
* Jonas Hägele

###Code
https://github.com/ppaulee/MDSim

###Compiler:
GNU Compiler 9.3.0


###Build Instructions:
####(without Doxygen):
* mkdir buildDir  
  cd buildDir  
  cmake -S ../ -B .  
  make  
  ctest \
  ./MolSim -f ../input.txt -e 5 -s 0.0002 -a lj -w 100
####(with Doxygen):
* mkdir buildDir  
  cd buildDir  
  cmake -S ../ -B . -DBUILD_DOXY=ON   
  make  
  make doxyDoc

###Simulation:
https://youtu.be/Cdi7o-XbQ14

###FAQ:
Datei wird nicht richtig eingelesen. \
Lösung: Line Endings mit \r\n \
Alternative Lösung: In MolSim.cpp Zeile 128 auskommentieren und Zeilen 129-130 einkommentieren.
