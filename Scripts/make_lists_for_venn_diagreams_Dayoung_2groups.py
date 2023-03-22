# First all presence/absence venn diagrams
# 1. All samples 2
## Comment out the wanted path 
# LCM
dirpath = "C:/Users/pc/Desktop/Dayohari/220615_DEA_Dayoung_ToSend/LCM/WithoutPig4/files_2_groups/"
# MF
# dirpath = "/home/sieglinde/5_Projects/220615_DEA_Dayoung_ToSend/MF/files_2_groups/"

present_in_nose_LCM = []
present_in_lung_LCM = []

line_count = 0 
with open(dirpath + "DeSeq2_Full_Table_lungvsnose_nozeroexpression.tsv", "r") as full:
    for line in full:
        line = line.rstrip()
        splittedline = line.split('\t')
        ens_code = splittedline[2]
        nasal = splittedline[9:17]
        lung = splittedline[17:]
        if line_count == 0:
            print(nasal)
            print(lung)
        else:
            present_lung = False
            present_nose = False
            for count in nasal:
                if not count == '0':
                    present_nose = True
            for count in lung:
                if not count == '0':
                    present_lung = True
            
            if present_nose:
                present_in_nose_LCM.append(ens_code)
            if present_lung:
                present_in_lung_LCM.append(ens_code)
        line_count += 1

#print(present_in_nose)
print(len(present_in_nose_LCM))
#print(present_in_lung)
print(len(present_in_lung_LCM))

with open(dirpath + "present_nose.csv", "w") as LCM_nose:
    for ens_code in present_in_nose_LCM:
        LCM_nose.write(ens_code + "\n")

with open(dirpath + "present_lung.csv", "w") as LCM_lung:
    for ens_code in present_in_lung_LCM:
        LCM_lung.write(ens_code + "\n")

########################################################################################################################
########################################################################################################################
# DEG Lung vs Nose 
## LCM
DEG_NoseVsLung_LCM = []

line_count = 0 
with open(dirpath +  "DeSeq2_B_Lung_vs_A_Nasal_filtered.tsv", "r") as full:
    for line in full:
        line = line.rstrip()
        splittedline = line.split('\t')
        ens_code = splittedline[0]
        if line_count == 0:
            print(ens_code)
        else:
            DEG_NoseVsLung_LCM.append(ens_code)
        line_count += 1

print(len(DEG_NoseVsLung_LCM))

with open(dirpath + "DEG_NoseVsLung.csv", "w") as file:
    for ens_code in DEG_NoseVsLung_LCM:
        file.write(ens_code + "\n")