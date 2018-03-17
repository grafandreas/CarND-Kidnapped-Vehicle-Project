# Overview
This repository contains my submission for the Kidnapped Car (Particle Filter) project.

# Overall result
The simulator marks the exercise / solution as passed.

# Notes on implementation

## Code structure

The code structure follows the code template closely. No additional methods or other refactoring was necessary.

## STL

For dealing with the lists of particles, observations, etc. I decided not to use the classic for(;;) code
structure, but opted for the use of the STL iterators + functions (transform, find etc.). The main reasons for
that was to learn about these STL features.

## Possible Code improvements

In real life, the vectors might be replaced by arrays, to have fixed length and faster access. Also, memory handling
was no priority in this project. In a real embedded system we would check to not create objects on the heap. But the
focus of this exercise is on the algorithm.

## Major challenges
Coming from Java/Xtend, the major challenge of this exercise was the unexpected behavior of some STL methods.

