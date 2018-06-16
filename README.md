## Overview
Experimental implementation of the paper with title 'Locality-sensitive hashing of curves' published by A. Driemel and F. Silvestri. More information you can find at [this PDF file](http://drops.dagstuhl.de/opus/volltexte/2017/7203/pdf/LIPIcs-SoCG-2017-37.pdf) or at [docs folder](https://github.com/chanioxaris/Hashing-Search-PolygonalCurves/tree/master/doc) containing the paper.

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
The format of input text file, described by the following structure:
```
@R: 2.0
curve_idS1 mS1 (xS11,y11) (xS12,y12) (xS13,y13) ... (xS1mS1,yS1mS1)
curve_idS2 mS2 (xS21,y21) (xS22,y22) (xS23,y23) ... (xS2mS2,yS2mS2)
curve_idS3 mS3 (xS31,y31) (xS32,y32) (xS33,y33) ... (xS3mS3,yS3mS3)
...
curve_idSN mSN (xN1,yN1) (xN2,yN2) (xN3,yN3) ... (xNmN,yNmN)
```
where ```mSi``` the total points included in i curve, ```(xSij,ySij)``` the coordinates of point j in Si curve and ```R``` the radius of neighbours search


## Compile

`./makefile`

## Usage

`./lsh –d [input file] –q [query file] –k [lsh functions] -L [hashtables] -ο [output file] –function {DFT, DTW} –hash {classic, probabilistic} -stats (optional)`
