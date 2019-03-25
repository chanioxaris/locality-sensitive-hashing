## Overview
Experimental implementation of the paper 'Locality-sensitive hashing of curves' published by A. Driemel and F. Silvestri. More information you can find at [this PDF file](http://drops.dagstuhl.de/opus/volltexte/2017/7203/pdf/LIPIcs-SoCG-2017-37.pdf) or at [docs folder](https://github.com/chanioxaris/Hashing-Search-PolygonalCurves/tree/master/doc).


## Locality-sensitive hashing
One of the main applications of LSH is to provide a method for efficient approximate nearest neighbor search algorithms. Consider an LSH family ![F](https://latex.codecogs.com/gif.latex?F). The algorithm has two main parameters: the width parameter k and the number of hash tables L.

In the first step, we define a new family ![G](https://latex.codecogs.com/gif.latex?G) of hash functions g, where each function g is obtained by concatenating k functions ![](https://latex.codecogs.com/gif.latex?h_1,....,&space;h_k) from ![F](https://latex.codecogs.com/gif.latex?F), i.e., ![](https://latex.codecogs.com/gif.latex?g(p)&space;=&space;[h_1(p),....,&space;h_k(p)]). In other words, a random hash function g is obtained by concatenating k randomly chosen hash functions from F. The algorithm then constructs L hash tables, each corresponding to a different randomly chosen hash function g.

In the preprocessing step we hash all n points from the data set S into each of the L hash tables. Given that the resulting hash tables have only n non-zero entries, one can reduce the amount of memory used per each hash table to O(n) using standard hash functions.

Given a query point q, the algorithm iterates over the L hash functions g. For each g considered, it retrieves the data points that are hashed into the same bucket as q. The process is stopped as soon as a point within distance c*R from q is found.

## Metric Distance

### Discrete Frechet Distance (DFD)
The [Fréchet distance](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance) is a measure of similarity between curves that takes into account the location and ordering of the points along the curves. It is named after [Maurice Fréchet](https://en.wikipedia.org/wiki/Maurice_Fr%C3%A9chet).

Imagine a man traversing a finite curved path while walking his dog on a leash, with the dog traversing a separate one. Assume that the dog varies her speed to keep as much slack in her leash as possible: the Fréchet distance between the two curves is the length of the shortest leash sufficient for both to traverse their separate paths. Note that the definition is symmetric with respect to the two curves—the Frechet distance would be the same if the dog was walking her owner.

![Frechet Distance](https://github.com/chanioxaris/locality-sensitive-hashing-curves/blob/master/img/frechet_distance.jpg)

### Dynamic Time Warping (DFT)
In [time series analysis](https://en.wikipedia.org/wiki/Time_series), dynamic time warping (DTW) is one of the algorithms for measuring similarity between two temporal sequences, which may vary in speed. For instance, similarities in walking could be detected using DTW, even if one person was walking faster than the other, or if there were accelerations and decelerations during the course of an observation. 

In general, DTW is a method that calculates an optimal match between two given sequences (e.g. time series) with certain restriction and rules:
* Every index from the first sequence must be matched with one or more indices from the other sequence, and vice versa
* The first index from the first sequence must be matched with the first index from the other sequence (but it does not have to be its only match)
* The last index from the first sequence must be matched with the last index from the other sequence (but it does not have to be its only match)
* The mapping of the indices from the first sequence to indices from the other sequence must be monotonically increasing, and vice versa

The [optimal match](https://en.wikipedia.org/wiki/Optimal_matching) is denoted by the match that satisfies all the restrictions and the rules and that has the minimal cost, where the cost is computed as the sum of absolute differences, for each matched pair of indices, between their values.


## Input data files

### Input file 
The format of input text file, described by the following structure:
```
@dimension {2,3,4} 
curve_id1 m1 (x11,y11) (x12,y12) (x13,y13) ... (x1m1,y1m1)
curve_id2 m2 (x21,y21) (x22,y22) (x23,y23) ... (x2m2,y2m2)
curve_id3 m3 (x31,y31) (x32,y32) (x33,y33) ... (x3m3,y3m3)
...
curve_idN mN (xN1,yN1) (xN2,yN2) (xN3,yN3) ... (xNmN,yNmN)
```
where ```mi``` the total points included in i curve, ```(xij,yij)``` the coordinates of point j in i curve, when dimesion is equal to 2

### Query file 
The format of query text file, described by the following structure:
```
@R: 2.0
curve_idS1 mS1 (xS11,yS11) (xS12,yS12) (xS13,yS13) ... (xS1mS1,yS1mS1)
curve_idS2 mS2 (xS21,yS21) (xS22,yS22) (xS23,yS23) ... (xS2mS2,yS2mS2)
curve_idS3 mS3 (xS31,yS31) (xS32,yS32) (xS33,yS33) ... (xS3mS3,yS3mS3)
...
curve_idSN mSN (xSN1,ySN1) (xSN2,ySN2) (xSN3,ySN3) ... (xSNmSN,ySNmSN)
```
where ```mSi``` the total points included in i curve, ```(xSij,ySij)``` the coordinates of point j in Si curve and ```R``` the radius of neighbours search


## Compile

`./makefile`

## Usage

`./lsh –d [input file] –q [query file] –k [lsh functions] -L [hashtables] -ο [output file] –function {DFT, DTW} –hash {classic, probabilistic} -stats (optional)`
