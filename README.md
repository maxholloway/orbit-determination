# orbit-determination
Python 3 program to determine orbital elements from four observation coordinates.
With the orbital elements, one has the tools to know the location of Near-Earth 
Asteroid for up to 25 million years.

In order to run this program to find the location of an asteroid other than
the one described by the current input file, one must edite the "Inputs.txt"
file with inputs in a very particular manner. Namely, it must be of the following
format:

{year} {month} {day} {Civil time in {hours} {minutes} {seconds}} {Right Ascension (RA) in {hours} {minutes} {seconds}} {Declination in {hours} {minutes} {seconds}} {x distance to sun} {y distance to sun} {z distance to the sun}}

All fields are separated by 1 space, and do not include the curly braces. An
example of correct formatting is in the "Inputs.txt" file. Additionally, the 
program requires that there be 4 observations; each row of the inputs file represents
an observation of an asteroid.