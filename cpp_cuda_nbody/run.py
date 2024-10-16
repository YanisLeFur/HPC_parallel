import os
import time

nthreads = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
for nt in nthreads:
    
    f = open(f"batches/row_nt={nt}.job", 'w')
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --nodes=1\n")
    f.write("#SBATCH --time=0:5:0\n")
    f.write(f"#SBATCH --partition=gpu\n")
    f.write(f"#SBATCH --gres=gpu:1\n")
    f.write(f"#SBATCH --qos=gpu_free\n")
    f.write(f"#SBATCH --account=math-454\n")
    f.write(f"#SBATCH --reservation=Course-math-454-final\n")
    f.write(f"#SBATCH --output=output/row_nt={nt}.out\n\n")
    f.write(f"srun nvprof ./nbody-code-row examples/galaxy.txt {nt}")
    f.close()
    os.system(f"sbatch batches/row_nt={nt}.job")
    time.sleep(1)






# nthreads = [1, 2, 4, 8, 16, 32]

# # One thread per entrt
# for nt in nthreads:
    
#     f = open(f"batches/entry_nt={nt}.job", 'w')
#     f.write("#!/bin/bash\n")
#     f.write("#SBATCH --nodes=1\n")
#     f.write("#SBATCH --time=0:5:0\n")
#     f.write(f"#SBATCH --partition=gpu\n")
#     f.write(f"#SBATCH --gres=gpu:1\n")
#     f.write(f"#SBATCH --qos=gpu_free\n")
#     f.write(f"#SBATCH --account=math-454\n")
#     f.write(f"#SBATCH --reservation=Course-math-454-final\n")
#     f.write(f"#SBATCH --output=output/entry_nt={nt}.out\n\n")
#     f.write(f"srun nvprof ./nbody-code-entry examples/galaxy.txt {nt} {nt}")
#     f.close()
#     os.system(f"sbatch batches/entry_nt={nt}.job")
#     time.sleep(1)