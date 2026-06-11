# Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

# LOGIN TO SERVER:
ssh -X -p 25260 m033372@mtc-b.phys.hawaii.edu

# KILL SCREEN:
screen -X -S 1665 quit

# TO COPY FILES TO SERVER:
scp -r -P 25260 /Users/glennjocher/Google\ Drive/MATLAB/neutrinos/nViewGUI/TSresults/fcnrunTS.m  m033372@mtc-b.phys.hawaii.edu:/scratch/MATLAB/nViewGUI/TSresults

# TO COPY FILES FROM SERVER:  (DO NOT LOG INTO SERVER)
scp -r -P 25260 m033372@mtc-b.phys.hawaii.edu:/data3/skimmedData/specialRuns/exp_0001_run_2184.glenn /Users/glennjocher/Desktop/




# BEFORE CREATING A GITHUB REPO:
find functions-matlab/ -size +100M -ls
find . -name '.DS_Store' -type f -delete
find . -name 'Icon*' -type f -delete