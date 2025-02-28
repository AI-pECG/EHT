# EHT Analysis Overview

This page provides instructions on how to process EHT data. For a given EHT experiment, there will be a batch of EHTs with post movements over time. In our pipeline, all this data is stored within a single .txt file (as seen in Sample_data/Results.txt). This pipeline will identify all the individual EHTs, assess post deflections over time, calculate common EHT metrics, and compile all results within a file "EHT_results".

# Instructions
Execute "sbatch EHT_matlab.sh folder_name pixel_calibration pacing_hz". For the sample data, folder_name = Sample_data, pixel_calibration = 67 (units of pixels/mm), and pacing_hz = 1 (unit of Hz). For this Sample Data, you'd therefore execute "sbatch EHT_matlab.sh Sample_data 67 1"

# Results
When you go to the Sample_data folder, you'll see:
- AcquireEHT_ (input data for each tissue analyzed)
- AcquireEHT_fig (output MATLAB force development over time plot for each tissue)
- AcquireEHT_result.txt (output text results for each tissue)
- EHT_results (which compiles all the tissue results that are within the folder)
