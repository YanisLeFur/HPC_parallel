import os
import time

#weak scaling
ntasks = [1, 2, 4, 8, 16, 32,64]

#N = [625,1250,2500,5000,10000,20000,40000]
N = [40000,40000,40000,40000,40000,40000,40000]

# this is for strong scaling
for i in range(len(ntasks)):
    if(ntasks[i]==64):
        f = open(f"batches/N_{N[i]}_t_{ntasks[i]}.job", 'w')
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --reservation=Course-math-454-final\n")
        f.write("#SBATCH --account=math-454\n")
        f.write(f"#SBATCH --ntasks={ntasks[i]} --qos=parallel\n")
        f.write(f"#SBATCH --output=output/N_{N[i]}_t_{ntasks[i]}.out\n\n")
        f.write(f"ulimit -l 127590")
        f.write(f"srun ./nbody-code examples/galaxy_{N[i]}.txt")
        f.close()
        os.system(f"sbatch batches/N_{N[i]}_t_{ntasks[i]}.job")
        time.sleep(1)
    
    else:
        f = open(f"batches/N_{N[i]}_t_{ntasks[i]}.job", 'w')
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --reservation=Course-math-454-final\n")
        f.write("#SBATCH --account=math-454\n")
        f.write(f"#SBATCH --ntasks={ntasks[i]}\n")
        f.write(f"#SBATCH --output=output/N_{N[i]}_t_{ntasks[i]}.out\n\n")
        f.write(f"ulimit -l 127590")
        f.write(f"srun ./nbody-code examples/galaxy_{N[i]}.txt")
        f.close()
        os.system(f"sbatch batches/N_{N[i]}_t_{ntasks[i]}.job")
        time.sleep(1)