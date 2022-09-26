# Expansion_Spatial_Transcriptomics
Expansion Spatial Transcriptomics data analysis from:

Expansion Spatial Transcriptomics\
Yuhang Fan, Zaneta Andrusivova, Yunming Wu, Chew Chai, Ludvig Larsson, Menxiao He, Liqun Luo, Joakim Lundeberg, Bo Wang

Copyright (C) 2022  Yuhang Fan & Zaneta Andrusivova

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

This code will perform data analysis for RNA-seq results from Ex-ST.

The code requires Python (version 3.7.12) and R (version 4.0.5).
The code runs on Ubuntu 18.04.5 LTS (GNU/Linux 5.4.0-124-generic x86_64)
Processor Intel(R) Core(TM) i7-6700K CPU @ 4.00GHz
 
Data files are too large to include within this Github repo and must be
downloaded separately from:

https://data.mendeley.com/ (Ex-ST data and MOB standard ST data)

https://www.10xgenomics.com/resources/datasets/mouse-brain-section-coronal-1-standard-1-0-0 (mouse hippocampus standard ST data)

http://mousebrain.org/adolescent/downloads.html (Single cell mouse brain RNA-seq data)


Python Toolboxes required are:
- samalg (version 0.8.9)
- scanpy (version 1.8.2)
- pandas (version 1.3.4)
- numpy (version 1.20.3)
- seaborn (version 0.11.2)
- sciPy (version 1.9.1)

R Toolboxes required are:
- STUtility (version 0.1.0)
