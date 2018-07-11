## Overview
Experimental implementation of the paper with title 'Locality-sensitive hashing of curves' published by A. Driemel and F. Silvestri. More information you can find at [this PDF file](http://drops.dagstuhl.de/opus/volltexte/2017/7203/pdf/LIPIcs-SoCG-2017-37.pdf) or at [docs folder](https://github.com/chanioxaris/Hashing-Search-PolygonalCurves/tree/master/doc) containing the paper.


### Metric Distance

#### Discrete Frechet Distance (DFD)
The [Fréchet distance](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance) is a measure of similarity between curves that takes into account the location and ordering of the points along the curves. It is named after [Maurice Fréchet](https://en.wikipedia.org/wiki/Maurice_Fr%C3%A9chet).

Imagine a man traversing a finite curved path while walking his dog on a leash, with the dog traversing a separate one. Assume that the dog varies her speed to keep as much slack in her leash as possible: the Fréchet distance between the two curves is the length of the shortest leash sufficient for both to traverse their separate paths. Note that the definition is symmetric with respect to the two curves—the Frechet distance would be the same if the dog was walking her owner.

![Frechet Distance](https://github.com/chanioxaris/Hashing-Search-PolygonalCurves/blob/master/img/frechet_distance.jpg)

#### Dynamic Time Warping (DFT)

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
