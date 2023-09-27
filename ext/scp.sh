input='/home/maa160/SnapHiC-D/ext/Lee_Astro_samples.txt'
while IFS= read -r line
do
	slurm_path="${line}"
	filename=$(basename ${line})
	local_path="/Users/megan/CompBio-Data/Lee2019/"
	sshpass -p "S!143netmegants" scp -P24 maa160@solar.cs.sfu.ca:${slurm_path} ${local_path}
done < "$input"
