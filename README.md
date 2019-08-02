# HPV_vacc_costeffect_Switzerland
Scripts used for the study on cost-effectiveness of HPV vaccination in Switzerland (quadri- and nonavalent vaccine)

Manuscript title: Impact and cost-effectiveness of nonavalent HPV vaccination in Switzerland: 
insights from a dynamic transmission model
Authors: Maurane Riesen, Johannes A. Bogaards, Nicola Low, Christian L. Althaus

Core script of the HPV cost-effectiveness project in Switzerland: 00_corescript.R

Maurane Riesen, version May 2019
---------------------------------------------------------------------------------------------------

The overall script is composed by a core script which is linked to 
- R files to retrieve required parameters and a folder containing parameter tables (01_params_script_*, tab_params)
- R files meant to be run on the high performance computing cluster UBELIX, used to do the 
  calibration and sensitivity analysis, in two distinct steps (02_UBE_part*)
- R files containing the different functions (03_function_*)

 ------------------------------------------------------------------------------------------------------
