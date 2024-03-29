# ReCoDE Diffusion Code

# Table of Contents
- [Abstract](#abstract)
- [Structure of the Code](#structure-of-the-code)
- [User Guide](#user-guide)
- [Features of the Code](#features-of-the-code)
- [Exercises](#exercises)
- [PETSc Installation](#petsc-installation)
- [Theory](#theory)

# Abstract
This code is part of the Research Computing and Data Science Examples (ReCoDE) project. The code itself is a 1-dimensional neutron diffusion solver written in Fortran in an object oriented format. The example will focus on features of the code that can be used as a teaching aid to give readers some experience with the concepts such that they can implement them in the exercises or directly in their own codes. An understanding of neutron diffusion and reactor physics is not required for this example, but a discussion of the theory can be found in the bottom section of this readme.

This project aims to provide examples of:
- [Compiled Codes and Makefiles](#compiled-codes-and-makefiles)
- [Compiler Directives](#compiler-directives)
- [Object Oriented Programming](#object-oriented-programming)
- [Reading input data from file](#reading-input-data-from-file)
- [Generating output files](#generating-output-files)
- [Paraview](#paraview)
- [Solving mathematical problems](#solving-mathematical-problems)
- [Discretisation of a spatial dimension](#discretisation-of-a-spatial-dimension)
- [Optimised data storage](#optimised-data-storage)
- [Incorporating external libraries (PETSc)](#incorporating-external-libraries)


# Structure of the Code
Object Oriented Programming is a method of coding by which a piece of software is designed around the different required data types and classes which will be used in the programme. In it's simplest sense, this translates to a code which is made of multiple different modules which make use of one another, rather than one large block of code. This style of programming is particularly advantageous when working with large coding projects, such as those used in research or data science. Object oriented codes are also much easier to read once their structure has been understood. As such, this section of the readme will give an explanation of how the Diffusion Code project is structured.

At it's simplest level, the code reads a specified problem from an input file, converts that to a system of equations which can then be solved, and outputs the resulting data to a set of files which can be read by an external program such a GNUPlot or Paraview.
$$ \text{Input File } \rightarrow \text{Generate Equations } \rightarrow \text{Solve } \rightarrow \text{Output File } $$

This explanation can now be further expanded in terms of complexity, where the structure will be given in terms files and modules. For the sake of readability, this is given as a full flow chart below. The **Problem** module reads through the input file, storing relevant data or passing it to the **Materials** module. Data from these modules is the used by the **MatGen** module to generate the system of equations. If PETSc is used, this data is then passed to the **PETScMat** and **PETScVec** modules, which are wrappers for the data library. These are then passed into the **PETScKSP** module which solves the problem. If PETSc isn't used, this data is passed into the **CRS** module, which stores the data effieicently, such that it can be fed into the **Solver** module. The solved data is then passed into the **Output** module, which generates an output both in .txt and .vtu format.

![DiffCode drawio](https://user-images.githubusercontent.com/83182489/173537965-ac15206a-dc13-4659-89f0-5b7141eb3091.png)

This can also be represented through some blocks of pseudocode

First read through the input file
```fortran
Open( Input File)
Read Problem Data
Read Material Data
Close( Input File)

Set Problem Data
Set Material Data
```

Then generate the equations and solve for the flux
```fortran
Get Problem Data
Get Material Data

Do 1, Problem Size
  Calculate Matrix Value
  Calculate Vector Value
EndDo

If (PETSc Used) Then
  Do 1, Problem Size
    Fill PETSc Matrix
    Fill PETSc Vector
  EndDo
  Flux = PETScSolver( PETSc Matrix, PETSc Vector)
Else
  Do 1, Problem Size
    Fill CRS Matrix
    Fill Vector
  EndDo
  Flux = Solver( CRS Matrix, Vector)
EndIf
```

Finally generate the output files
```fortran
Open( Output File)
Do 1, Problem Size
  Write( Output File) Position, Flux
  Write( Output File) Region Number
  Write( Output File) Node Number
EndDo
Close( Output File)
```

# User Guide
  ## Compiling the Code
  As fortran is a compiled language, the code must first be compiled before it can be executed. To do so, navigate to the **src** directory within the code and enter the commands:
  ```bash
  make clean
  make
  ```
  The **make clean** command first leans the directory of any module files or executables created by any past makes. This should generally be done before you compile the code with a new make option. The **make** command will execute the makefile contained within the same directory, converting the fortran code to an optimised form that can be read by your computer. This then generates an executable named **diffusion**. Once this has been done, the code can be executed from the parent directory using the command:
  ```bash
  ./diffusion
  ```
  This command tells the executable to run and will generate relevant output files containing the solution to the problem. 


  ## Changing Input Options
  The code is designed such that the user can easily change the problem which the code is attempting to solve. The code uses the input file **Input.txt** to read details about the problem, such as the positions of boundaries or materials in the problem. An example of such an input can be seen below:
  ```
  ------------------------------------------
  Regions: - Integer number of regions in the problem
  2
  ------------------------------------------
  Boundaries: - Real number positions of the boundaries between the regions (one per line)
  0.0
  0.5
  1.0
  ------------------------------------------
  Nodes: - Integer number of nodes in each region (one per line)
  5
  5
  ------------------------------------------
  Materials: - Fuel, Water or Steel (one per line)
  Fuel
  Fuel
  ------------------------------------------
  Boundary_Conditions: - Zero or Reflective (two parameters - one per line)
  Zero
  Zero
  ------------------------------------------
  ```
   For this example problem, we are stating that we have a geometry ranging from *x = 0.0* to *1.0*, half fuel and half steel with a central boundary at *x = 0.5*. As seen from the above input, the code needs four different parameters to be described to it.
  - **Regions** - An integer number of regions that exists within the problem. We have 1 region from *0.0* to *0.5* and another from *0.5* to *1.0*, hence we give the code the integer number **2**.
  - **Boundaries** - The positions of the boundaries within the problem. Our first boundary is at the start of our geometry so we enter the number **0.0**. We then have an internal boundary halfway through the problem seperating the regions so we enter the number **0.5**. Finally we have the exterior boundary of our geometry, so we enter the number **1.0**. The code will always read one more value here than the number of regions in the problem.
  - **Nodes** - This desribes how refined we want the geometry in each region. For the example we want a quick solve with just enough nodes to see the flux profile. As we need to describe this for each region we enter the value **10** twice. The code will always read the same number of values here as the number of regions in the problem.
  - **Materials** - This described the materials that are present within the system. The first half of our geometry is fuel, with the latter half being Steel, so we enter **Fuel** and **Steel**. The code will always read the same number of values here as the number of regions in the problem.
  - **Boundary Conditions** - This tells the code what boundaries exist at the edges of our problem. Two boudnary conditions have been implemented in out code, that of 'Zero' and 'Reflective'. The former simply ensures that the flux will tend to zero at the boundary, while the latter ensures that the derivative of the flux will tend to zero at the boundary.
  ## Reading Output Files
  The code generates two output files, **Output.txt** and **Output.vtu**. The former is a simple text file containing the position and flux of the solution to the problem. These are simply given in two columns such that they can be read easily by someting like a GNUPlot or Python script. An example of such a flux profile can be seen below:
  ```
  0.000000E+00  0.135813E+01
  0.166667E+00  0.137306E+01
  0.333333E+00  0.141907E+01
  0.500000E+00  0.150000E+01
  0.666667E+00  0.158093E+01
  0.833333E+00  0.162694E+01
  0.100000E+01  0.164187E+01
  ```
  This can then be plotted using tools such as GNUPlot to produce the profile seen below:

  ![FluxProfile](https://github.com/ImperialCollegeLondon/ReCoDE_Diffusion_Code/blob/main/images/FluxProfile.png)

  GNUPlot was chosen here as the application is very simple and the program itself is very easy to install. On Linux, you only need the commands:
  ```bash
  sudo apt update
  sudo apt install gnuplot
  ```
  For installation on Windows, you will need to download the software from http://www.gnuplot.info/download.html. The GNUPlot script used to generate the plot is named **fluxplot** and can be found in the home directory of the project. This can be run with the command:
  ```bash
  gnuplot -p fluxplot
  ```
  The  **Output.vtu** file stores additional data such as the region and cell numbers in an XML format. This can be directly read by software such as Paraview, allowing for more detailed visualisations than that of the previous flux profile. For visualisation purposes, the data has been smeared in a second dimension, which should give users an idea of how multi-dimensional cell data can be viewed. An example output from paraview can be seen in the image below:

  ![FluxParaview](https://github.com/ImperialCollegeLondon/ReCoDE_Diffusion_Code/blob/main/images/FluxParaview.png)
  
  Instructions on how to install paraview can be found at: https://www.paraview.org/Wiki/ParaView:Build_And_Install

  A user guid for paraview can be found at: https://docs.paraview.org/en/latest/

# Features of the Code

## Object Oriented Programming

As discussed previously in the [Structure of the Code](#structure-of-the-code) section, this code utilises Object Oriented Programming. An example of such an abject orented structure can be seen in the code snippet below. This shows some of the **Materials** module, which handles the storage of data pertaining the material properties of the problem.

```Fortran
Module Materials_Mod

    use Constants_Mod
    !!Stores standard material data and material data explicitly set via an input file

    Implicit None

  type, public :: t_material
        private
        Real(kind=dp) :: Sig_a, S
        Character(len=20) :: Name
    contains
    !!Procedures which handle the storing, calculation and retrieval of material data
        procedure, public :: SetName => SetMaterialName
        procedure, public :: GetName => GetMaterialName
        procedure, public :: SetProps => SetMaterialProperties
        procedure, public :: GetSig_a
        procedure, public :: GetS
  end type
```

We first define the name of the module such that it can be used by other modules in our code.

```Fortran
Module Materials_Mod
```
We then tell our module to make use of the **Constants** module.

```Fortran
use Constants_Mod
```
This method of declaring the name of this specific module and the modules which this code can be seen throughout our various pieces of code. We would like to store some data within this module such that it can be called forth at a later date. As such, we set up a public type which contains our chosen private data.

```Fortran
type, public :: t_material
        private
        Real(kind=dp) :: Sig_a, S
        Character(len=20) :: Name
```
We can then also choose what public procedures this module should utilise, which other pieces of code can then utilise. The most common routines you will often utilise in a code are 'gets' and 'sets' which pull stored data and set stored data respectively. In fortran you can also utilise pointers to point to certain pieces of memory. A nice example of this can be seen in the code block below where **SetName** points to **SetMaterialName** allowing us to use quick subroutine names when repeatedly calling a routine, but verbose names within the actual module.

```Fortran
      procedure, public :: SetName => SetMaterialName
      procedure, public :: GetName => GetMaterialName
      procedure, public :: SetProps => SetMaterialProperties
      procedure, public :: GetSig_a
      procedure, public :: GetS
end type
```







## Compiled Codes and Makefiles

Programming languages can be differentiated into 'Interpreted' and 'Compiled' languages. Lagnuages such as Python are interpreted languages, where the code is read directly by the computer and the actions described are performed. In contrast, languages such as fortran must be compiled by a program before the code can be used, where a program known as a compiler reads through the code, convering it to machine language. While time must be spent compiling the code, this generally results in very noticeable performance increases and is the reason why most high perfromance codes are compiled. Readers familiar with Python may be familiar with the Numpy library that is often used for numerical calculations, a library which makes use of compiled codes to give Python codes some dramatic reductions in numerical computation time.

Makefiles are particularly useful in a compiled code that utilises OOP as it allows for easy compilation of the required files and handles the order of dependencies. The code snippet below comes from the makefile used in the project and shows how an ordered list of objects is used when compiling the code. This is written such that the compiler will work through the files in a way that each dependency is compiled before the modules that utilise it, ending with the main.

```Makefile
OBJS= stdlib/Constants.o \
 inoutproc/Materials.o \
 inoutproc/Problem.o \
 inoutproc/Output.o \
 centralproc/Matrix_Base.o \
 centralproc/CRS.o \
 centralproc/Solver.o \
 centralproc/MatGen.o \
 main.o \
```



## Compiler Directives

The makefile can also facilitate compilation in a number of different ways depending on the options provided. The snippet below compiles the code when the command ```make debug``` is used on the command line, such that it uses all common flags set by **$(FC_FLAGS_COMMON)**. In addition, the flag ```-DDBEUG``` creates the compiler flag 'DEBUG'.
```Makefile
debug:   FC_FLAGS = $(FC_FLAGS_COMMON) -DDEBUG
```
This flag can then be checked with C pre-processor language to allow for differnet parts of the code to be compiled dependent on the flags incorporated. In the below example, a more detailed output is given by the **Solve** module when the code is compiled with ```DEBUG```.
```fortran
# ifdef DEBUG 
  Write(*,*) "---CG Convergence Succeeded---"
  Write(*,'(g0)',advance='no') "Succeeded after iterations:  "
  Write(*,'(g0)',advance='no') BCG_Iterations
  Write(*,'(g0)',advance='no') "  with residual:"
  Write(*,'(E14.6)') rho1
  Write(*,*) "-------------------------------"
# endif
```

To test this we can enter the following commands to view the modified output with additional information:
```bash
make clean
make debug
```

Then in our main directory we can run the newly compiled 'debug' version of the code
```bash
./diffusion
```


## Reading input data from file

Scientific codes will often make use of an input file which contains the problem specification. These codes must therefore be able to open an input file and read through it, extracting the relevant data such that it can be used to solve the problem. This can be seen in the code snippet below taken from the **Problem** module, where an input file has been opened and the boundary conditions read in.

```Fortran
!!Open input file containing problem specification
Open(InputFile, File='input.in', Status='Old')

!!Read in the boundary conditions of the problem
String_Read = ''
Do While (String_Read .NE. 'Boundary_Conditions:')
  Read(InputFile,*) String_Read
EndDo
Do ii = 1, 2
  Read(InputFile,*) String_Read
  If (String_Read == 'Zero') Then 
    this%Boundary_Conditions(ii) = 0
  ElseIf (String_Read == 'Reflective') Then 
    this%Boundary_Conditions(ii) = 1
  Else 
    Write(*,*) "ERROR: Unrecognised Boundary Condition"
  EndIf
EndDo
```


When reading an input file, fortran needs an Integer ID, Filename and Status. The ID allows the file to be referred to easily in future, and is good practice to set through a parameter for readability. The Filename simply matches the name of the file, and the status describes the state of the file, in this case it is an 'Old' file which already exists in the directory.
```Fortran
Integer, parameter :: InputFile = 101
Open(InputFile, File='input.in', Status='Old')
```

We then wish to read through our file until we reach the name of the boundary conditions. To do so, we loop through our file until we reach the 'Boundary_Conditions:' string. We now know that the next line will contain the name of our boundary condition, and hence this can be read in.
```Fortran
String_Read = ''
String_Read = ''
Do While (String_Read .NE. 'Boundary_Conditions:')
  Read(InputFile,*) String_Read
EndDo
```

We finally need to check what type of boundary has been specified and store this information. Here we have an If statement which will loop over the known boundary condition names.

```Fortran
Read(InputFile,*) String_Read
If (String_Read == 'Zero') Then 
  this%Boundary_Conditions(ii) = 0
ElseIf (String_Read == 'Reflective') Then 
  this%Boundary_Conditions(ii) = 1
Else 
  Write(*,*) "ERROR: Unrecognised Boundary Condition"
EndIf
```

Finally we should close this file, achieved through the close command and the associated ID.

```Fortran
Close(InputFile)
```


## Generating output files

Most scientific codes will want to generate outputs that can be read by external software such as GNUPlot. To do so, data generated by the code must be written to an output file in some given format. This is done in a way similar to that of reading the input, with the read statements replaced with write statements.

First we tell the code to open the output file, specifying 'Replace' to tell the code that we wish to overwrite the file if it is already present.

```Fortran
!!Generate textfile
Open(textfile,File='OutputFile.txt',Status='Replace') 
```


We then loop over our solutions, writing the values of the position and fluxes to the textfile, with each row corresponding to a node. In code snippet below we write down the data for the first node, then loop over the rest of the nodes in each of the regions.

```Fortran
Position = 0._dp
!!First node
NodeID = 1
Write(textfile,'(2E14.6)') Position, Flux(NodeID)
Do ii = 1, N_Regions
    Do jj = 1, RegionNodes(ii)-1
        NodeID = NodeID + 1
        Position = Position + (Boundary_Pos(ii+1)-Boundary_Pos(ii))/Real(RegionNodes(ii)-1,dp)
        Write(textfile,'(2E14.6)') Position, Flux(NodeID)
    EndDo 
EndDo
```

Finally we close the file. We now have a complete set of output data stored in text format that can be plotted or used for later analysis.
```Fortran
Close(textfile)
```



## Paraview

One particular output viewing software that may be of interest to some readers is Paraview. Paraview is a widely used software for data analysis and visualisation within the scientific community. While paraview could read in our textfile outputs to generate a graph, this does not demonstrate much of the power of the software. Instead the code has been written such that it will generate a 2-Dimensional plot of the flux that can be manipulated in Paraview.

To do so, we must generate an output file which uses the VTK file format and translate our resulting data to this format. We have also made use of the more modern .VTU file type, a form of VTK file which makes use of the html file structure. In addition to the readability we gain from such a structure, it is also recommended that users avoid the .VTK legacy file type if possible. A snippet from **OutputFile.vtu** can be seen below, where we have specified that this is a 'VTKFile' using the 'Unstructured Grid' format. Contained within the Unstructured Grid, we define the number of points and cells within the problem. We then go on to describe the relevant point and cell data. The vtk file format is described here: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf.

```html
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="18" NumberOfCells="8">
```
## Solving Mathematical Problems

- Discuss how a mathematical problem can be converted into something solveable

- How you can write your own modules to solve systems of equations

- Brief aside on preconditioners

In order to generate a solution to our problem we must first convert the general equation that we wish to solve into a system of equations. Scientific problems will often consist of a series of known properties of a system which, when multiplied by an unknown parameter, generate known values. Such systems of equations can be converted into a matrix and vector problem of the form $Ax=b$, where $A$ is our matrix describing the properties of the system, $x$ is a vector containing our unknown parameter of interest and $b$ is a known parameter of the system. 

Simple matrix mathematics tells us that this could be trivially solved if we can invert the matrix $A$, such that:

$$ A x = b$$

$$ A^{-1} A x = A^{-1} b $$

$$ x = A^{-1} b $$

Where $A^{-1}$ is the inverse of the matrix $A$. Our system produces a realtively simple Tri-Diagonal matrix (where only the 3 main diagonals contain data), and as such can be inverted through a number of relatively simple algorithms (such as the Thomas Algorithm). While this would be an oprimal way of solving such a problem, this is only for a specific case and not generally applicable to a wide range of problems. As such, we have instead incorporated an iterative solver which progressively approximates the solution to the problem. As the algorithm iterates, the error between the values of $Ax$ and $b$ are compared, where the problem is 'solved' once this error passes below a certain threshold.

Such solvers are widely used throughout scientific computing as they are widely applicable to a range of different problems involving a range of different matrices. While the commonly used Conjugate Gradient (CG) method is servicable for a number of different problems, it does require the matrix $A$ to be symmetric, a limitation which we would ideally like to avoid for this code. As such we have implemented the Bi-Conjugate Gradient (BCG) method for our solver, a modified form of the CG algorithm which can handle asymmetric matrices.


## Discretisation of a spatial dimension

- Discretisation in space using data stored within other modules

While our code could solve problems that involve homogenous systems, real systems will have properties that vary in some dimensions. A 1-Dimensional code such as this example code can facilitate variations in properties in the x dimension. To do so, we must have some method of disretising our spatial x dimension which allows us to consistently generate relevant systems of equations to describe the system. The specific deisretisation scheme that we utilise in our code is known as 'finite difference', a scheme that generates very simple systems of equations, especially in a one-dimensional case.

This discretisation scheme can be seen in the image below, whre we split out domain into nodes from $0$ to $n$. Each non-boundary node $i$ will therefore have corresponding neighbouring nodes $i-1$ and $i+1$ which will need to be incorporated into a specific equation for the resulting solution to the node, described in further detail in the [Theory](#Theory) section. 

![FD](https://github.com/ImperialCollegeLondon/ReCoDE_Diffusion_Code/blob/main/images/FD.png)

This therefore generates a system of equations of the form:

$$ P_{i-1}\phi_{i-1} + P_{i}\phi_{i} + P_{i+1}\phi_{i+1} = S_{i} $$

Where $P$ is some property of the system given by the scientific equation of interest, $\phi$ are the properties of interest that are being solved for (in our case neutron flux), and $S$ is a known result of the system (in our case the source of neutrons). As can be seen from the above equation, the solution at one node $i$ only depends on its immediate neighbours, so any matrix created for the system will be tridiagonal, meaning that only the main and adjacecnt diagonals will be filled.


## Optimised data storage

In addition to using faster codes such as fortran, we can also achieve some significant speed increases by being careful with our use of memory. One important optimisation in our code will be how we store the data in our matrix $A$. The matrix $A$ for a 3-node and 5-node system can be seen below, where $a$, $b$ and $c$ are our lower, main and upper diagonals respectively. We can already begin to see that we are storing a large portion of $0$ 's in our system, a probelm which will only increase as we increase our node numbers by orders of magnitude. 

$$
\left(\begin{array}{ccc} 
b & c & 0\\
a & b & c\\
0 & a & b\\
\end{array}\right)
\rightarrow
\left(\begin{array}{ccccc} 
b & c & 0 & 0 & 0\\
a & b & c & 0 & 0\\
0 & a & b & c & 0\\
0 & 0 & a & b & c\\
0 & 0 & 0 & a & b\\
\end{array}\right)
$$ 

These matrices are known as 'sparse', as much of the data that they contain is of value $0$. We can therefore drastically reduce the amount of memory needed to store this matrix by only saving the values and the positions that they occupy. A number of efficient storage systems exist for this problem, but we have chosen to use Compressed Row Storage (CRS) in this case as it makes no assumptions about the structure of the matrix and can therefore be applied to a large range of problems. A good explanation of this methodology can be found at http://netlib.org/linalg/html_templates/node91.html.


## Incorporating External Libraries

While we can make many efforts to write our own optimised data storage and solver routines, a number of libraries exist which will do this far more optimally. In many cases, people will have spent huge amounts of time to write libraries that can perform an action as efficiently as possible, so it is sensible for us to expect their routines to provided a performance increase when compared to ours.

For some smaller routines, you will often be able to find optimal examples shared online, but for larger problems you will often need to incorporate a whole library. In our code, we have facilitated the use of the 'Portable, Extensible Toolkit for Scientific Computation' (PETSc). Installing and using this library is explained in the respective [Installation](##Installing-PETSc) and [Compilation](##Compiling-with-PETSc) sections.

When using an external library, we will generally have to write a 'wrapper' in our code, which provides an interface between our own routines and the routines of the library. An example of this is shown in the code snippets below for the **PETSc_Init** module, which initialises the PETSc library. We first create our module and use the **#include** statements to tell the code to use our main path specified in our makefile and then look for this specific file **petscsys.h**. We then tell it to **use petscsys**, allowing us to perform actions which require this file, such as initialising PETSc.

```Fortran
Module PETSc_Init_Mod
!!Initialize the PETSc Database
#include <petsc/finclude/petscsys.h>
use petscsys
```

We then write our subroutine which does the actual initialisation, **PETSc_Init**. We check if the routine has been called already and if not we call **PETScInitialize**, a routine of the PETSc library. This sets up PETSc such that we can now define our relevant vectors or matrices and then solve the system of equations.

```Fortran
Subroutine PETSc_Init
  PetscErrorCode ierr
  Logical :: Called = .FALSE.
  If (.NOT. Called) Then
    Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    Called = .TRUE.
    If (ierr .NE. 0) Then
      Error Stop "Failed to Initialize PETSc"
    EndIf
  EndIf
End Subroutine PETSc_Init
```


# Exercises
This section contains a number of suggested exercises that will give the user a deeper understanding of the topics covered in the above descriptions and utilised throughout this code.  The exercises will increase progressively in complexity, with jumps in challenge roughly seperated by the numeric ID of the exercise. 

Solutions to the exercises can be found in the 'solutions' directory, where a version of the main code exists, modified such that it solves the problem. Where possible, changes made to the main code have been appropriately commented with a double dollar symbol such that users can more readily find these modifications. For example, changes made for exercise 1a can be found by searching for '$$ Exercise 1a' within the appropriate solution directory.

## Exercise 1a

**Task:** Add a print statement to the **SetName** Subroutine in the Materials Module

**Aims:**

-  Gain an initial understanding of how a code can utilise an OOP structure

For this exercise, add a print (or write) statement to the **SetName** Subroutine contained within the Materials Module, which prints the name of the material being set to the terminal.

## Exercise 1b

**Task:** Write a new subroutine which prints all material data stored in the module

**Aims:** 

- Develop an understanding of Object Orientation
- Learn how to write a new subroutine and implement it in a class structure
- Learn how to utilise a type to call the generated routine
- Opportunity to gain experience with Fortran output formatting

For this exercise, you must make your own subroutine within the Materials Module called **PrintMaterial**. This subroutine, when called from the Main, should print all material data contained within the materials class to the terminal. This is also a good opportunity for users to gain some experience with Fortran formatting if they wish, as this can make the terminal output much easier to read.

## Exercise 1c

**Task:** Make a new material which can be used in the input deck

**Aims:** 

- Experience with fortran logical arguments, namely **If Statements**
- Initial practice adjusting some exisitng logic within a subroutine

For this exercise, a new material will be added to the list of those that the code can handle. The material 'Iron' will be added with an absoprtion cross section of *4.0* and source of *0.1*.




## Exercise 2a

**Task:** Get input deck to instead explicitly read in material data rather than associating with a name

**Aims:** 

- Experience making larger adjustments to the logic of the code
- Understanding of how to adjust various sections of a code to handle a modification
- Experience handling input files

Currently, the code reads in a material name and then assigns material properties to the region based on a set of data stored in an If statement. For this exercise, the code will instead read in the absorption and source terms directly, such that 'Fuel' would instead become '1.0 6.0'. To achieve this, the user will need to adjust the **Problem** module, as well as adjust how they have defined their problem in the **input.in** file.



## Exercise 2b

**Task:** Add another two sets of cell data to the .vtu file where the absorption and source material properties of each cell can be viewed 

**Aims:** 

- Gain experience writing data to an output file
- Learn how to utilise vtu files to add data sets to paraview

The current vtu file output containts the Flux, Cell Number and Region Number data sets. In this exercise, the user will add the Source term data set to the vtu output, such that the source value of each cell can be viewed in Paraview. An additional Cell Data set will need to be added which uses the known Region Numbers to pull data from the **Materials** Module.


## Exercise 3a

**Task:** Smear the problem into another dimension to produce a 3D paraview output

**Aims:** 

- Gain further experience writing data to an output file
- Learn how to utilise vtu files to generate more complex outputs in paraview

The main version of the code already smears a 1-dimensional flux profile into a second dimension for the sake of a visualisation example. This could be furthered by also smearing the results in a third dimension. To do so, the user will need to make an additional set of nodes that have been translated in the z axis and ensure that these are then associated with the correct cells. Users should utilise the paraview vtk format guide to do so, noting that they will now need to define the cells as a three dimensional shape. The vtk file format is described here: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf



## Exercise 3b

**Task:** Create a Compressed Diagonal Storage Module and switch it out for the existing CRS Module

**Aims:** 

- Gain an understanding of how to write a completely new module and incorporate it into the OOP structure
- Gain experience with memory efficient programming
- Gain experinece with abstract classes

Currently the non-PETSc version of the code utilises a CRS system for handling the matrices involved in the problem. For this task, the user should create their own module which performs the same task utilising Compressed Diagonal Storage (CDS). The CDS module should have all the features of the CRS, allowing it to be used with the polymorphic matrix constructor and solver that already exists within the code. The user should ensure that they have added their CDS module to the makefile, and allocated it as 't_cds' where this has been done for 't_crs'. As a stretch goal, the user could also compare the solve times with their CDS module against that of the original CRS. An explanation of the CDS methodology can be found at: http://netlib.org/linalg/html_templates/node94.html



## Possible Extension Exercise

**Task:** Make a wrapper for the blas/lapack library and use instead of PETSc/CG

**Aims:** 

- Gain experience installing external libraries
- Understand how one can incorporate an external library into the code

This extension exercise involves incorporating an additional external library into the code. BLAS and LAPACK are often used in professional codes as the utilise a number of highly optimised mathematical routines. To solve the problem, the user should utilise DGETRF and DGETRI to perform LU decomposition and inversion of the matrix, then multiply it by the source vector using DGEMV. The code can be installed easily on a linux OS with the commands:

```bash
sudo apt-get install libblas-dev liblapack-dev
```

An explanation of how to use the code and each routine can be found at:
http://www.netlib.org/lapack/explore-html/

# Installing PETSc 
  The Portable, Extensible Toolkit for Scientific Computation (PETSc) is an optional external library which can be utilised in this problem to solve the system of equations generated by the code. This section will give an explanation of how to download and compile PETSc. Note that the configuring and testing of PETSc may take some time. Installation instructions can also be found at https://petsc.org/release/install/

  - Download a copy of PETSc from:

  https://petsc.org/release/download/

  - Place the downloaded PETSc folder into the directory of your choice and navigate inside it through a terminal. Configure the code by running:
  ```bash
  ./configure --download-mpich --download-fblaslapack=1
  ```
  - Check that the installation was successful by running:
  ```bash
  make all check
  ```

# Compiling with PETSc

  With PETSc installed on your system, the last stage is to set up the environment variables that the code will use to find your PETSc directory. The text below shows how this can be done on a Linux OS, by entering the following commands into your bashrc file:

  ```bash
export PETSC_DIR="/home/jack/petsc-3.16.0"

export PETSC_BUILDS_DIR="/home/jack/petsc-3.16.0/builds"
```

  Once you have installed PETSc (see [installation instructions](##Installing-PETSc)) and set up the required environment variables, the code can also be compiled to use it instead of the custom storage and solvers. To do so, navigate to the **src** directory within the code and enter the commands:
  ```bash
  make clean
  make petsc
  ```
  This will execute the makefile contained within the same directory with the petsc options, converting the fortran code to an optimised form that utilises the PETSc library. This then generates an executable named **diffusion_petsc**. Once this has been done, the code can be executed from the parent directory using the command:
  ```bash
  ./diffusion_petsc
  ```


# Theory

In the design of any nuclear reactor, it is critical to understand the distribution of neutrons throughout the system. As neutrons move through various parts of the reactor such as the fuel or the moderator, their behaviour will greatly change based on the medium in which they find themselves. This neutron behaviour can be most simply approximated as a form of diffusion occurring within the reactor, producing a solvable equation which can be used to accurately describe this behaviour. The analysis and solving of this equation is known as neutron diffusion theory, and is a critical area of nuclear physics and engineering, used widely in the industry to calculate neutron flux profiles  and multiplication factors within a reactor.

The first step in generating this equation is to describe the spatial neutron balance within a volume $dV$, where $dV = dx dy dz$, centred at $r$, where $r = (x,y,z)$. Assuming steady state:

$$ \text{neutrons lost in } dV/s = \text{neutrons produced in } dV/s$$

Splitting these into the separate sources and sinks:

$$ \text{neutrons leaking from } dV/s + \text{neutrons absorbed in } dV/s = \text{neutrons emitted in } dV/s + \text{fission neutrons produced in } dV/s$$

Where: 

$$ \text{neutrons leaking from } dV/s =  \left(  \frac{\partial}{\partial x} J_x (x,y,z) + \frac{\partial}{\partial y} J_y (x,y,z)  + \frac{\partial}{\partial z} J_z (x,y,z)   \right) dV$$

$$\text{neutrons absorbed in } dV/s = \Sigma_a (x,y,z) \phi(x,y,z) dV$$

$$ \text{neutrons emitted in } dV/s = S (x,y,z) dV$$

$$\text{fission neutrons produced in } dV/s = \nu \Sigma_F (x,y,z) \phi(x,y,z) dV$$

Combining these and eliminating $dV$ gives:

$$ \frac{\partial}{\partial x} J_x (r) + \frac{\partial}{\partial y} J_y (r) + \frac{\partial}{\partial z} J_x (z) + \Sigma_{a} (r) \phi(r) = S(r) + \nu \Sigma_f (r) \phi (r) $$

Given the vector definitions: 

$$ J(r) = \hat{i} J_x (r) + \hat{j} J_y (r) + \hat{k} J_z (r) \text{  and   } \nabla = \hat{i} \frac{\partial}{\partial x} + \hat{j} \frac{\partial}{\partial y} + \hat{k} \frac{\partial}{\partial z} $$

We can write the balance equation as:

$$ \nabla \cdot J(r) + \Sigma_{a} (r) \phi(r) = S(r) + \nu \Sigma_f (r) \phi (r) $$

Given that we are using diffusion to describe the behaviour of neutrons within a reactor, we can therefore make use of Fick's Law, which states that any solute (the neutrons) will diffuse from an area of high concentration (high neutron flux) to an area of low concentration. This law gives an equation that will relate the current in a system to the concentration gradient multiplied by a diffusion coefficient. In the simplest 1-D reactor case the current of neutrons will therefore be given by the negative of a diffusion coefficient multiplied by the gradient of the flux. A general form of this can be seen below, where $J$ is the neutron current (or diffusion flux), $D$ is the diffusion coefficient and $\phi$ is the neutron flux. 

$$ J = - D(r)\nabla\phi(r) $$

This can be substituted into the full form of the balance equation to form the neutron diffusion equation: 

$$ - \nabla \cdot D(r)\nabla\phi(r) + \Sigma_{a} (r) \phi(r) = S(r) + \nu \Sigma_f (r) \phi (r) $$

It is therefore also important to have an accurate method of calculation for the diffusion coefficient. The diffusion coefficient is equal to 1/3 multiplied by the inverse of the transport cross section. For simple neutron diffusion cases involving isotropic scatter, the transport cross section can be set to be equal to the total cross section, which is equal to the sum of the absorption and scattering (including self-scatter) cross sections. This can be seen below where $\Sigma_{tr}$ is the transport cross section, $\Sigma_{t}$ is the total cross section, $\Sigma_{a}$ is the absorption cross section and $\Sigma_{s}$ is the scattering cross section.

$$ D =  \frac{1}{3\Sigma_{tr}} \simeq \frac{1}{3\Sigma_{t}} = \frac{1}{3(\Sigma_{a}+\Sigma_{s})} $$

This relation exists as the transport cross section is defined as $\Sigma_{tr} = \Sigma_{t} - \bar{\mu}\Sigma_{s}$ , where $\bar{\mu}$ is the average cosine of the scattering angle. This has a value of $0$ in the laboratory system for isotropic scatter, so the transport and total cross sections can be approximated to one another in this case. 

This equation effectively relates the rate of change of neutrons within a system to a number of material properties and the flux, and hence can be solved with knowledge of the materials involved and the use of mathematical solvers. An example of the neutron diffusion equation can be seen below, for a 1-D slab reactor at steady state, where $\lambda$ is the eigenvalue of the system, $\nu$ is the average neutrons produced per fission and $\Sigma_f$ is the fission cross section. In this example a fission source of neutrons is being used instead of a simple volumetric source.

$$ -\frac{d}{d x} D(x) \frac{d \phi(x)}{d x}+\Sigma_{a}(x) \phi(x)=\frac{1}{\lambda} \nu \Sigma_{f}(x) \phi(x) $$

Neutron diffusion codes will principally solve this equation to calculate the neutron flux over a specified geometry. Problems involving a fission neutron source can be solved for the eigenvalue $\lambda$ of the system, utilising a fission source normalisation.
