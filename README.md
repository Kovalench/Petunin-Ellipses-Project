# Petunin's Ellipses for Two Point Sets (C++ + Python)

This project implements the full pipeline for constructing **Petunin’s ellipses** for two labeled point sets enclosed by a convex polygon.  
The C++ module generates or reads input data, computes geometric transformations, identifies main axes and perpendiculars, scales the data into a square, and writes all necessary parameters into an output file.  
The Python module loads this output and visualizes three stages:  
1) original space,  
2) scaled space,  
3) Petunin’s ellipses.

The project demonstrates algorithmic geometry, convex regions, coordinate transformations, and visual analytics.

---

## Features

### C++ Module
- Reads point sets from file or generates them randomly inside a convex polygon  
- Parses linear inequalities to build the polygon  
- Splits points into two labeled sets (+1 and -1)  
- Finds the **main line** via the two furthest points  
- Computes perpendiculars and projections onto the main line  
- Rotates point sets to align the main axis horizontally  
- Scales the rotated space into a perfect square  
- Computes radii used to construct Petunin’s ellipses  
- Exports all processed data to `output_data.txt`

### Python Module
- Loads the full output format from the C++ program  
- Visualizes:
  - Original points and polygon  
  - Main line and perpendiculars  
  - Scaled points and circles  
  - Final Petunin ellipses in original coordinates  
- Automatically detects polygon type (quadrilateral or pentagon)

---

## Project Structure

/src
petunin.cpp # Full C++ implementation (geometry & preprocessing)

/python
visualize.py # Matplotlib visualization of original, scaled, and ellipse spaces

/data
input_data.txt # Example dataset
output_data.txt # Generated results

README.md

---

## How the Algorithm Works

### 1. Input Data
The program accepts:
- A set of 2D points with labels: `1` and `-1`
- A set of linear inequalities defining a convex polygon  
  Format: `A B C` meaning:  
  `A*x + B*y <= C`

### 2. Main Line
The algorithm finds the two **furthest points** from Set1.  
These define the main axis of Petunin’s ellipses.

### 3. Perpendiculars
For each side of the line, the algorithm locates the point with the largest perpendicular distance and computes its orthogonal projection.

### 4. Rotation
The entire coordinate system is rotated such that the main axis becomes horizontal.

### 5. Scaling to Square
All points and polygon vertices are scaled uniformly into a square to ensure proper ellipse construction.

### 6. Output Format
`output_data.txt` contains:
- scaling factors  
- circle and ellipse radii  
- original and scaled sets  
- polygon vertices  
- perpendiculars  
- main line  
- centers

The Python script reconstructs everything from this file.

---
