# HOR Detection for rice
This script detects HORs for dimer-compressed monomer sequences of rice centromeres.

## Requirements
- Python 3.6+
- networkx
- pygraphviz

## Usage
The script detects HOR for one dimer-compressed sequence at a time. 
```sh
#examples
python HORdetect_rice.py -fi data/all.anno2.cps.cps.cps_seq.fa -pos data/all.anno2.cps.cps.cps.bed -p Chr02_13-65 -O . -nodeThr 2 -edgeThr 2 -minTrav 2
```

## Used Code
This script includes some code from HORmon.
- `BuildMonomerGraph` & `drawGraph`: build monomer graph and visualization
- `genAllCycles` & `genCycleInner` ：get initial cycles

## Original Copyright Notice 
HORmon  
Copyright (c) 2021 Saint Petersburg State University  

HORmon is free software; you can redistribute it and/or modify  
it under the terms of the GNU General Public License, Version 2,  
dated June 1991, as published by the Free Software Foundation.  

HORmon is distributed in the hope that it will be useful, but  
WITHOUT ANY WARRANTY; without even the implied warranty of  
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
General Public License for more details.  

You should have received a copy of the GNU General Public License along  
with this program; if not, write to the Free Software Foundation, Inc.,  
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.  









