#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G 
#SBATCH --time=10:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.pixy_n90.%J.out"
#SBATCH --job-name="pixy_caculate"

awk 'FNR>1 || NR==1' MKGDPRYuan_n121.*h*.subset*_dxy.txt > output_temp/MKGDPRYuan_n121.dxy.txt
awk 'FNR>1 || NR==1' MKGDPRYuan_n121.*h*.subset*_pi.txt > output_temp/MKGDPRYuan_n121.pi.txt
awk 'FNR>1 || NR==1' MKGDPRYuan_n121.*h*.subset*_fst.txt > output_temp/MKGDPRYuan_n121.fst.txt

cd output_temp

cut -f 1,2,5 MKGDPRYuan_n121.pi.txt > MKGDPRYuan_n121.pi2.txt
cut -f 1,2,3,6 MKGDPRYuan_n121.dxy.txt > MKGDPRYuan_n121.dxy2.txt
cut -f 1,2,3,6 MKGDPRYuan_n121.fst.txt > MKGDPRYuan_n121.fst2.txt

awk 'FNR > 1 {
    group = $1"\t"$2        # Create a key for the group using the first and second columns
    count[group]++          # Increment the count for this group
    sum[group] += $3        # Accumulate the sum of the third column for this group
    sumsq[group] += $3*$3   # Accumulate the square of the third column (for variance calculation)
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]                     # Calculate mean
        variance = (sumsq[group] / count[group]) - (mean * mean) # Calculate variance
        sd = (variance > 0) ? sqrt(variance) : 0             # Calculate standard deviation
        print group, mean, sd                                # Print group, mean, and SD
    }
}' MKGDPRYuan_n121.pi2.txt | sort  > MKGDPRYuan_n121.pi3.txt


awk 'FNR > 1 {
    group = $1"\t"$2"\t"$3   # Create a key for the group using the first, second, and third columns
    count[group]++           # Increment the count for this group
    sum[group] += $4         # Accumulate the sum of the fourth column for this group
    sumsq[group] += $4*$4    # Accumulate the square of the fourth column (for variance calculation)
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]                      # Calculate mean
        variance = (sumsq[group] / count[group]) - (mean * mean) # Calculate variance
        sd = (variance > 0) ? sqrt(variance) : 0              # Calculate standard deviation
        print group, mean, sd                                 # Print group, mean, and SD
    }
}' MKGDPRYuan_n121.dxy2.txt | sort  > MKGDPRYuan_n121.dxy3.txt

awk 'FNR > 1 {
    group = $1"\t"$2"\t"$3   # Create a key for the group using the first, second, and third columns
    count[group]++           # Increment the count for this group
    sum[group] += $4         # Accumulate the sum of the fourth column for this group
    sumsq[group] += $4*$4    # Accumulate the square of the fourth column (for variance calculation)
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]                      # Calculate mean
        variance = (sumsq[group] / count[group]) - (mean * mean) # Calculate variance
        sd = (variance > 0) ? sqrt(variance) : 0              # Calculate standard deviation
        print group, mean, sd                                 # Print group, mean, and SD
    }
}' MKGDPRYuan_n121.fst2.txt | sort  > MKGDPRYuan_n121.fst3.txt

