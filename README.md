Project 6: Sea Level Rise (slr)
--------------------

Code to view flooded area after sea level rise. To run, first
make. The code processes elevation grids saved in ASCII (.asc)
files. The program takes the arguments:

$ ./slr <file>.asc <max slr> <increment> <grid size reduction factor>

file:
.asc file name to parse

max slr:
Maximum sea level flooding to calculate. For example, max slr
= 50 would allow you to toggle through slr values from 0 to 50.

increment:
Step size for rendering intermediate flooding levels. For
example, increment = 3 would allow the user to render flooded areas
in increments of 3 meters per keypress.

grid size reduction factor:
How many cells to average together, to reduce the resolution of the
input grid. The density is the square root of the number of grid cells
to reduce into one cell. For example, a density parameter of 3 will
average out the values in 9 neighboring cells into a single cell.

Controls
----------
'+' increase sea level rise by increment
'-' decrease sea level rise by increment

'2' bird-eye's view (default)
'x' or 'X' rotate over x axis
'y' or 'Y' rotate over y axis
'z' or 'Z' rotate over z axis
