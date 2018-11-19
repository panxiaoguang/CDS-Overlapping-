# CDS-Overlapping  
The purpose of the project is to find a consensus CDS sequence based on the ID of the gene, thereby facilitating the silencing, knockout, etc. of the gene.
## Program composition   
  - EOCE.py
  - download.py  
## Usage  
The process will be completed in two steps. First, you need to download the genbank sequence to the local, and then extract the common area locally through a series of steps to finally complete the target.  
### download.py
``` sh
python3 download.py -i idlist.txt -o download/
```   
- note:  
   - [idlist.txt]() Need to write the ID into the file in order, split with a newline
   - download/ this is a folder that your downloadfile in
### EOCE.py   
``` sh
python3 EOCE.py -o download/ -u position/ -s sequence/
```
- note:
   - download/ the file you first download
   - position/ this is overlapping position like [start,stop].....[19334,19782]
   - sequence/ this overlapping sequence
## Dependent environment
Mainly used [biopython](https://biopython.org/wiki/Download) module download and sequence operations, so you need to install python and biopython modules first.
```
pip[pip3] install biopython
```
