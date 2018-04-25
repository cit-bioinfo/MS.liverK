# MS.liverK

  Cancer is a highly heterogeneous disease, with marked differences between patients. Even within a specific type of cancer, such as hepatocellular carcinoma (HCC), heterogeneity persists on several points: diverse genetic mutations, resulting in diverse cancer phenotypes. Therefore, characterizing these differences within a cancer type is highly important to better understand the mechanisms that undergo cancer maintenance and progression. A more thorough classification of cancer samples into different subtypes also has tremendous benefits in terms of precision medicine, as is allows clinicians to target therapies to the specific disease of each patient. A striking example of targeted therapy is the treatment of resistant patients by anti-PD1 antibody immunotherapy. To better understand cancer heterogeneity and identify patients more likely to respond to a specific treatment, the notion of molecular subgroups has emerged.
  
  In HCC, several teams have designed unsupervised or supervised classification strategies. MS.liverK provides easy-to-use functions to characterize HCC samples from transcriptomic data. It implements 6 different molecular subtypes classification and scores 45 genes signatures about molecular subtypes, prognosis or biologial pathways.

Installation
========
Install from the GitHub repository using devtools:

    install.packages("devtools")
    library(devtools)
    devtools::install_github("FPetitprez/MS.liverK")

Dependencies
========
The R packages "pamr" and "survival" are required for this package to run properly.

License
========

MS.liverK is a free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
