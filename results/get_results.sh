SAMPLES=`cat ../data/samples_list.txt`

for SAMPLE in ${SAMPLES}
do
  aws s3 cp s3://cchauve-orchestration-ch/altera-pilot-study-canexia/${SAMPLE} . --exclude "*" --include "*indels*.vcf.tar.gz" --recursive
  tar xzvf ${SAMPLE}_indels_filtered_snpeff.vcf.tar.gz
  tar xzvf ${SAMPLE}_indels_unfiltered.vcf.tar.gz
  mv ${SAMPLE}_indels_filtered_snpeff.vcf.tar.gz ${SAMPLE}/
  mv ${SAMPLE}_indels_unfiltered.vcf.tar.gz ${SAMPLE}/
done
