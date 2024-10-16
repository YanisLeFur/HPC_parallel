import os
import time

#weak scaling
ntasks = [1, 2, 4, 8, 16, 32,64]

#N = [625,1250,2500,5000,10000,20000,40000]
N = [40000,40000,40000,40000,40000,40000,40000]

# this is for strong scaling
for i in range(len(ntasks)):
    if(ntasks[i]==64):
        os.system(f"sbatch --ntasks={N[i]} --qos=parallel script.run")
        time.sleep(1)
    
    else:
        os.system(f"sbatch --ntasks={N[i]} script.run")
        time.sleep(1)
      