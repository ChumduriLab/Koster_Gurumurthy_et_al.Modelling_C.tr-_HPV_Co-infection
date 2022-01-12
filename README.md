# Modelling Chlamydia and HPV co-infection in patient-derived ectocervix organoids reveals distinct cellular reprogramming

Cervical mucosa is confronted by coinfections with pathogenic microbes, yet their implications in pathogenesis remain unclear. Lack of in-vitro models 
recapitulating cervical epithelium has been a bottleneck to study coinfections. Using patient-derived ectocervical organoids, we systematically modeled 
individual and coinfection dynamics of Human papillomavirus (HPV)16 E6E7 and Chlamydia associated with carcinogenesis. Ectocervical stem cells from 
HPV-negative donors were genetically manipulated to introduce HPV16 E6E7 oncogenes to mimic HPV integration. Organoids from these stem cells developed 
the characteristics of precancerous lesions while retaining the capacity to self-renew and organize into mature stratified epithelium similar to healthy 
organoids. HPV16 E6E7 interfered with Chlamydia development and induced persistence. Unique transcriptional and post-translational responses induced by 
Chlamydia and HPV led to distinct reprogramming of host cell processes. Strikingly, Chlamydia impedes HPV-induced mechanisms that maintain cellular and 
genome integrity, including mismatch repair in the stem cells. Together, our study employing organoids demonstrates the hazard of multiple infections and 
the unique cellular microenvironment they create, potentially contributing to neoplastic progression.





















#### For more information, please read our manuscript published in Nature Communication:
	
	






	
	
*Stefanie Koster, Rajendra Kumar Gurumurthy, Naveen Kumar, Pon Ganish Prakash, Jayabhuvaneshwari Dhanraj, Sofia Bayer, Hilmar Berger, 
Shilpa Mary Kurian, Marina Drabkina, Hans-Joachim Mollenkopf, Christian Goosmann, Volker Brinkmann, Zachary Nagel, Mandy Mangler, 
Thomas F Meyer, Cindrilla Chumduri*. (doi: https://doi.org/10.1101/2021.04.15.439996)



























## Overview
```
Folders:
└── Data: Pre-processed, Differential Exp., GSEA RData files are stashed here.
└── Code: Contains scripts for microarray analysis used in this manuscript.
        └── Human Ectocervical Organoids
            │       ├── 00_Preprocessing_and_QC.Rmd
            │       └── 01_Diff_Expression_analysis.Rmd
            │       ├── 02_Ctr_HPV_transcription_modules_v2.Rmd
            │       └── 03_GSEA_Analysis.Rmd
            │       ├── 04_GSEA_plots.Rmd
            │       └── 05_AdditionalVisualizations_2D.Rmd
            │       ├── 06_AdditionalVisualizations_3D.Rmd
            │       └── 07_Helper_functions.R	    
            └── Cervical Intra Epithelia Tissue
            │       └── 08_External_CIN_Analysis.R	
```































### File Description:

1. R scripts: 
     - ***00_Preprocessing_and_QC.Rmd***: Microarray data QC and Preprocessing for hCEcto E6E7 organoids with or without Ctr infection.
      	 
     - ***01_Diff_Expression_analysis.Rmd***: To perform Differential expression analysis on hCEcto E6E7 organoids with or without Ctr infection. 
      
     - ***02_Ctr_HPV_transcription_modules_v2.Rmd***: To analyse gene sets to infer TF activity, Gene Ontology analyses for modules of Ctr-HPV-coinfection.	 
     
     - ***03_GSEA_Analysis.Rmd***: GSEA analysis on DGE results obtained for hCEcto E6E7 organoids with or without Ctr infection.
      
     - ***04_GSEA_plots.Rmd***: Snippets to visualize the obtained GSEA results.
      
     - ***05_AdditionalVisualizations_2D.Rmd***: Generating heatmaps, venn diagrams and dotplots shown in the manuscript (Fig.3 & 4, Suppl.Fig.3 & 4).
     	 
     - ***06_AdditionalVisualizations_3D.Rmd***: Generating heatmaps, venn diagrams and dotplots shown in the manuscript (Fig.3 & 4, Suppl.Fig.3 & 4).
     	 
     - ***07_Helper_functions.R***:Custom functions that can be utilized for mapping agilent probe ID to gene symbols, symbols to entrez ID and also 
	 					for additional datarepresentation.	    					
						
     - ***08_External_CIN_Analysis.R***: Profiling cervical intraepithelial neoplasia (CIN) characteristics in coinfections by analyzing the microarray 
						dataset of laser-captured epithelium from 128 CIN samples downloaded from the GEO 
						(accession no. GSE63514; Johan A. den Boon et.al 2015, PMID: 26056290) 	



























## Contact
Please email the corresponding author in the main manuscript for other requests!

*For specific queries related to this repository, please contact:*
https://github.com/ChumduriLab


www.chumdurilab.org





