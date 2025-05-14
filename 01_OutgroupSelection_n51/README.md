### Selecting the Gb wild samples from (Yuan et al)[https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/advs.202003634] samples 

#### Reads trimming
```
module load  trimmomatic/0.39-zwxnnrx
trimmomatic PE -threads $thr $file1 $file2 $tDir/$name.R1.fq.gz $tDir/$name.U1.fq.gz $tDir/$name.R2.fq.gz $tDir/$name.U2.fq.gz ILLUMINACLIP:Adapters.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:75
```

