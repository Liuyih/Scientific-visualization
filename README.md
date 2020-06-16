# Scientific-visualization
=================================================
---------------------12252019--------------------
=================================================
Instructions of "Learnply" Program:

1. Program contains the following folders, please make sure they are all under same directory
    ./learnply/
    ./Opengl/
    ./quadmesh_2D/
   
====================Run the Program=============================

2. To open the program, you need to install Microsoft Visual Studio.

3. ./Opengl/ is for environment setting, you do not need to change anything in this folder.

4. To run the program, after you install the VS, please go the ./learnply/ folder and open the "learnply.sln"

5. After you open the "sln" file, on the  "Solution Explorer" bar, click the source file, and check any of the cpp file,
   
if you find the syntax error, that is probably because of the version setting. To solve this problem, you need to right-click project, 
   
go to the "properties setting". In the "General" list, you will see "Windows SDK Verison", click the list and replace with the "10.017134.0". 
You do not need to type this after you click, it will pop out automatically.  

6. quadmesh_2D is for "ply" files, you do not need to edit anything in this folder.   

====================Program Control=============================

The program support 3 display modes:
    display mode 0: keyboard button '0' to active, show the whole mesh. 
    display mode 1: keyboard button '1' to active, show the wireframe of the mesh.  
    display mode 2: keyboard button '2' to active, show the mesh with coloring based on checkerboard algorithm. 
    display mode 8: keyboard button '8' to active, show the contour tree graph for the minecraft.png
    left mouse button: translate the model
    right mouse button: rotate the model
    shift + left button: zoom in/out
    mouse scroll button down: check the ID of the element
 
    keyboard button 'r': reset the viewpoint
    keyboard button '|': write the viewpoint
    keyboard button '^': load the viewpoint
    
=================================================
