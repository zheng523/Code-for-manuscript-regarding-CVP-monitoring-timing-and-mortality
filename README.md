## Data Availability

The code files for data extraction, cleaning and analysis in our manuscript entitled "Association between early central venous pressure monitoring and 30-day mortality in critically ill patients with comorbid heart failure and sepsis: a retrospective cohort study" are provided below. The annotations in these code files were originally written in Chinese and later translated into English by AI. The translations may not be perfectly accurate, but they do not affect the data analyses. 

### Prerequisites
- Install MIMIC-IV (v3.1) and MIMIC-IV-ECHO (v1.0 or v1.0.1, only the structured-measurement table) databases in your local PostgreSQL system, including the `mimiciv_hosp`, `mimiciv_icu`, `mimiciv_derived`, and `mimiciv_echo` (with only one table named 'echo') schemas. 

### Execution Steps
1. Run `lvef.sql`, then `icustays_data.sql`, and then `HFsepsis.sql` in pgAdmin or Navicat to extract data needed for later procedures.
   > **Note:** If the names and paths of your database and schemas are different, the code may require slight modifications.

2. Run `CVP.R` line by line in RStudio or Positron to clean and analyse the data.
   > **Note:** Certain operational steps were not included in the code and had to be performed manually, such as variable inspection, figure formatting and saving, result exporting, and table creation. Not all analysis results are presented in the manuscript.
