# java_Rmath
java_Rmath is a Java package providing a Java Native Interface for Rmath library, which includes random number generators and probability functions.

Build Instructions
------------------

* java_Rmath requires the installation of Rmath library (The command in Ubuntu is "sudo apt-get install r-mathlib"). 

* To compile the package, enter src directory and execute "make". Notice that you may have to change the extension of generated libraries in the Makefile based on your operating system. On OS X you have to change all the extensions of dynamic library to .dylib while on Linux the corresponding extensions are .so

* To clean generated file, type “make clean” on the command line. 

Running the tests
-----------------
* For testing, enter test directory and execute “make” . If you want to clean testing results and all class files, type "make clean". 

Source Repository
-----------------
java_Rmath's source-code repository is hosted here on GitHub: https://github.com/pkuxchen/java_Rmath.git


Authors
---------

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| LiZhen Nie | lznie111@gmail.com   | Department of Statistics Chicago University |
| Xiang Chen (maintainer)| pkuchenxiang@pku.edu.cn   | Visiting student, Department of Biostatistics  UCLA |
| Lu Zhang | lu.zhang@ucla.edu    | PhD student, Department of Biostatistics UCLA  |                            
| Sudipto Banerjee | sudipto@ucla.edu   | Professor, Department of Biostatistics  UCLA |
<!--- --->
                             


Licensing
---------
java_Rmath is licensed under the Creative Commons Attribution License. 

