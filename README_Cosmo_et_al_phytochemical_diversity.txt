This README_Cosmo_et_al_phytochemical_diversity.txt file was generated on 2021-03-23 by Leandro Giacobelli Cosmo

GENERAL INFORMATION

1. Title of Dataset: 

Cosmo_et_al_phytochemical_diversity

2. Author information:

Leandro G. Cosmo: Programa de Pós-Graduação em Ecologia, Institute of Biology, P.O.Box: 6109, University of Campinas – UNICAMP, 13083-970, Campinas, SP, Brazil/Ecology Graduate School, Biosciences Institute, University of São Paulo, São Paulo, 05508-900, Brazil

Lydia F. Yamaguchi: Department of Fundamental Chemistry, Institute of Chemistry, University of São Paulo, 05508-000, São Paulo-SP, Brazil
 
Gabriel M.F. Felix: Programa de Pós-Graduação em Ecologia, Institute of Biology, P.O.Box: 6109, University of Campinas – UNICAMP, 13083-970, Campinas, SP, Brazil

Massuo J. Kato: Department of Fundamental Chemistry, Institute of Chemistry, University of São Paulo, 05508-000, São Paulo-SP, Brazil

Rodrigo Cogni: Department of Ecology, University of São Paulo, São Paulo, 05508-900, Brazil

Martín Pareja: Department of Animal Biology, Institute of Biology, P.O. Box: 6109, University of Campinas – UNICAMP, 13083-970, Campinas, SP, Brazil

Corresponding author: Martín Pareja, E-Mail: mpareja@unicamp.br, Phone: +55 19-3521-6324.

2. Date of data collection: 

From July 2017 to April 2018

3. Geographic location of data collection: 

Reserva Biologica Municipal Serra do Japi, Jundiaí, São Paulo, Brazil 

4. Information about funding sources that supported the collection of the data:

This study was financed by a FAPESP-NSF Dimensions of Biodiversity-Biota collaborative grant “Chemically mediated multi-trophic interaction diversity across tropical gradients” (2014/50316-7) and by CAPES-Finance Code 001. MP received support from FAEPEX-UNICAMP (grant n° 60412 and PAPDIC 2014/4715). RC is funded by FAPESP (2013/25991-0), CNPq, (307015/2015-7 and 307447/2018-9), and a Newton Advanced Fellowship from the Royal Society. LGC currently receives a FAPESP scholarship (2019/22146-3).

DATA & FILE OVERVIEW

1. File List: 

Data files:

Cosmo_et_al_fulldataset.txt
Cosmo_et_al_fulldataset.xlsx
Cosmo_et_al_network.txt
Cosmo_et_al_network.xlsx
dataset_hplc_ms.txt
dataset_hplc_ms.xlsx
dataset_nmr.xlsx
dataset_nmr.txt

Scripts:

script_ICC_SEM.R
script_network_analysis.R
script_PD_measures.R
katz_function.R

DATA-SPECIFIC INFORMATION:

Cosmo_et_al_fulldataset.txt/Cosmo_et_al_fulldataset.xlsx: full dataset with all the variables used for the analyses in the main text.

Cosmo_et_al_network.txt/Cosmo_et_al_network.xlsx: incidence interaction matrix of the plant-caterpillar network analyzed in the main text. Rows correspond to plant samples, while columns correspond to caterpillar species.

dataset_hplc_ms.txt/dataset_hplc_ms.xlsx: dataset with the pre-processed HPLC-HRESIMS data used to compute the compositional dimension of phytochemical diversity. Rows correspond to plant samples, while columns correspond to the mass to charge ratios (m/z) and retention times of tentative compounds and their intensities.

dataset_nmr.txt/dataset_nmr.xlsx: dataset with the pre-processed NMR data used to compute the structural dimension of phytochemical diversity. Rows correspond to plant samples and columns to chemical shifts and their intensities.

script_ICC_SEM.R: R script to calculate the Intra-class correlation coefficient (ICC) and perform the structural equation model analysis presented in the main text.

script_network_analysis.R: R script to calculate the networks metrics described in the main text and fit the linear mixed model used in the main text.

script_PD_measures.R: R script to compute the phytochemical diversity metrics used in the main text from pre-processed HPLC-HRESIMS and NMR data, and decompose these metrics into their within a between-subject components.

katz_function.R: R script with a function that calculates Katz centrality for a given adjacency matrix. Used within the file script_network_analysis.R.
