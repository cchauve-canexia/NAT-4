awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 1.vcf > Natera_Tumor-7_C035_0122_010703.vcf; python filter_indels.py Natera_Tumor-7_C035_0122_010703.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 2.vcf > Natera_Tumor-6_C035_0121_010702.vcf; python filter_indels.py Natera_Tumor-6_C035_0121_010702.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 3.vcf > Natera_Tumor-2_C035_0117_010698.vcf; python filter_indels.py Natera_Tumor-2_C035_0117_010698.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 4.vcf > C035_0207_012541.vcf; python filter_indels.py C035_0207_012541.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 5.vcf > C035_0214_012548.vcf; python filter_indels.py C035_0214_012548.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 6.vcf > C035_0209_012543.vcf; python filter_indels.py C035_0209_012543.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 7.vcf > C035_0204_012538.vcf; python filter_indels.py C035_0204_012538.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 8.vcf > C035_0206_012540.vcf; python filter_indels.py C035_0206_012540.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 9.vcf > C035_0212_012546.vcf; python filter_indels.py C035_0212_012546.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 10.vcf > C035_0213_012547.vcf; python filter_indels.py C035_0213_012547.vcf
# awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 11.vcf > C035_0205_012539.vcf; python filter_indels.py C035_0205_012539.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 12.vcf > C035_0203_012537.vcf; python filter_indels.py C035_0203_012537.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 13.vcf > C035_0208_012542.vcf; python filter_indels.py C035_0208_012542.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 14.vcf > C035_0210_012544.vcf; python filter_indels.py C035_0210_012544.vcf
awk '{if (substr($1,1,1)=="#" || (length($4) + length($5)) >= 3) {print $0}}' 15.vcf > C035_0211_012545.vcf; python filter_indels.py C035_0211_012545.vcf
