#!/bin/bash
echo "Enter the test file path: "
read variant_list
path_control=$variant_list"_variants.txt"
if [ ! -e "$path_control" ] 
then
  echo "file doesnt exists"
  exit
fi
#for all classes


#list of interested mutation from the first file
mutation_list="ARID1A FAT1 KMT2D PIK3CA COL4A4 KDM6A RYR2 TP53 FGFR3 XIST SCN11A ANK3 BRAF CSMD1 DNAH5 LRP1B SYNE1 MGAM PBRM1 CSMD2 CSMD3 ALB FLG MUC16 TBC1D12 APC KIAA2022 VHL KRAS RB1 CTNNB1 KEAP1 RELN DRD5 RP1 MXRA5 PTEN PPP2R1A USH2A DSCAM NRAS PIK3R1" 

#list of interested copy new list from the second file
copynew_list="RP11-170N11.1 FAM72C S100A16 SH2D2A PVRL4 LGR6 C4BPA VGLL4 LINC00620 LRRC3B LOC101927374 RNU6-69P PLXNA1 BCHE MECOM TP63 FAM43A SCARNA22 GRK4 NFKB1 LOC100506858 LOC102467074 PDE4D MARVELD2 DTWD2 PAIP2 DUSP22 CD83 ATXN1 LOC101928519 CASC15 LRRC16A OPN5 TBC1D32 SYNE1 MICALL2 RADIL FBXL18 DAGLB AC017060.1 JHDM1D-AS1 CTSB FGF20 SLC18A1 UBXN8 ADAM32 SNX31 MYC C9orf53 CDKN2A ADARB2 TAF3 ITGB1 CCNYL2 LINC00845 PGAP2 MMP10 NCAPD3 MIR1244-1|chr12 RP11-486A14.2 HPD MIR2276 GSX1 LINC00457 DGKH LINC00371 DHRS12 NEK3 TPTE2P3 SFTA3 MAPK1IP1L ULK4P1 ZSCAN10 RBFOX1 LINC00919 SMPD3 CHST4 WWOX MAF ABR WSCD1 LYRM9 RAD51L3-RFFL PTRF SLC4A1 HIGD1B PRPSAP1 GATA6-AS1 MAPK4 PALM GADD45B CHD6 FAM65C DYRK1A GGT5 EP300 ARFGAP3 PHF21B CELSR1 "

#generate associative array for mutations
declare -A FID1
count=0;
for m in $mutation_list; do
FID1[$m]=$count;
count=$((count+1));
done

#generate associative array for copy 
declare -A FID2
count=0;
for m in $copynew_list; do
FID2[$m]=$count;
count=$((count+1));
done
#I created two arrays since I do not know if these lists are intersecting or not

for file in $variant_list; do
    printf 'start %s\n' $file

    START_TIME=$SECONDS
    #get the line containing the patients
    patients_line=$(head -n 1 $file\_CNs.txt)

    #init the associative array for patients
    unset PAT
    declare -A PAT

    #get patient list and create an associative array from it
    ARR=($(echo $patients_line | tr " " "\n"))    
    count=-2; #2 is for the text "gene symbol"
    for pat in "${ARR[@]}"; do	
	PAT[$pat]=$count;
	count=$((count+1));
    done

    #get the patient list from the mutations file now
    cut -f 1 $file\_variants.txt | uniq > ../data/patients_temp.txt
    while read pat; do
	if [[ -z ${PAT[$pat]} ]]; then
	    PAT[$pat]=$count
	    count=$((count+1));
	fi
    done < ../data/patients_temp.txt
    
    no_patients=$count;
    no_f_values=$((no_patients*42));
    no_f_valuesm1=$((no_f_values-1));
    #create the whole feature array - did not measure time which is faster
    unset F
    for i in $(seq 0 $no_f_valuesm1); do
	F[$i]=0;
    done  
    ELAPSED_TIME=$((SECONDS - START_TIME))
    echo "prep time: " $ELAPSED_TIME

    #PART 1: handline the variants file
    #filter the file only with the useful features
    START_TIME=$SECONDS
    #get only the two columns
    cut -f 1,2 $file\_variants.txt | uniq > $file\_variants12.txt;

    grep 'ARID1A\|FAT1\|KMT2D\|PIK3CA\|COL4A4\|KDM6A\|RYR2\|TP53\|FGFR3\|XIST\|SCN11A\|ANK3\|BRAF\|CSMD1\|DNAH5\|LRP1B\|SYNE1\|MGAM\|PBRM1\|CSMD2\|CSMD3\|ALB\|FLG\|MUC16\|TBC1D12\|APC\|KIAA2022\|VHL\|KRAS\|RB1\|CTNNB1\|KEAP1\|RELN\|DRD5\|RP1\|MXRA5\|PTEN\|PPP2R1A\|USH2A\|DSCAM\|NRAS\|PIK3R1' $file\_variants12.txt > $file\_variants12.filtered.txt
    ELAPSED_TIME=$((SECONDS - START_TIME))
    echo "cut and grep time: " $ELAPSED_TIME

    prev_patient="null"
 
    #go over all the lines for the variants files
    START_TIME=$SECONDS
    while read l; do
	set -- $l
	#feature id for mutation
	if [[ ! -z ${FID1[$2]} ]]; then	    
	    patient_id=${PAT[$1]};
	    gene_feat_id=${FID1[$2]};

	    #compute the index
	    feat_index=$((patient_id*42+gene_feat_id));
	    F[$feat_index]=1;
	fi  	    
    done < $file\_variants12.filtered.txt
    ELAPSED_TIME=$((SECONDS - START_TIME))
    echo "variants time: " $ELAPSED_TIME

    START_TIME=$SECONDS
    grep 'RP11-170N11.1\|FAM72C\|S100A16\|SH2D2A\|PVRL4\|LGR6\|C4BPA\|VGLL4\|LINC00620\|LRRC3B\|LOC101927374\|RNU6-69P\|PLXNA1\|BCHE\|MECOM\|TP63\|FAM43A\|SCARNA22\|GRK4\|NFKB1\|LOC100506858\|LOC102467074\|PDE4D\|MARVELD2\|DTWD2\|PAIP2\|DUSP22\|CD83\|ATXN1\|LOC101928519\|CASC15\|LRRC16A\|OPN5\|TBC1D32\|SYNE1\|MICALL2\|RADIL\|FBXL18\|DAGLB\|AC017060.1\|JHDM1D-AS1\|CTSB\|FGF20\|SLC18A1\|UBXN8\|ADAM32\|SNX31\|MYC\|C9orf53\|CDKN2A\|ADARB2\|TAF3\|ITGB1\|CCNYL2\|LINC00845\|PGAP2\|MMP10\|NCAPD3\|MIR1244-1|chr12\|RP11-486A14.2\|HPD\|MIR2276\|GSX1\|LINC00457\|DGKH\|LINC00371\|DHRS12\|NEK3\|TPTE2P3\|SFTA3\|MAPK1IP1L\|ULK4P1\|ZSCAN10\|RBFOX1\|LINC00919\|SMPD3\|CHST4\|WWOX\|MAF\|ABR\|WSCD1\|LYRM9\|RAD51L3-RFFL\|PTRF\|SLC4A1\|HIGD1B\|PRPSAP1\|GATA6-AS1\|MAPK4\|PALM\|GADD45B\|CHD6\|FAM65C\|DYRK1A\|GGT5\|EP300\|ARFGAP3\|PHF21B\|CELSR1' $file\_CNs.txt > $file\_CNs.filtered.txt;
    ELAPSED_TIME=$((SECONDS - START_TIME))
    echo "grep time: " $ELAPSED_TIME

    #PART2: handle the CN file
    START_TIME=$SECONDS
    no_g_values=$((no_patients*99));
    no_g_valuesm1=$((no_g_values-1));
    unset G
    for i in $(seq 0 $no_g_valuesm1); do
	G[$i]=0;
    done  

   
    while read l; do
	#I do not know how to get the first word of a line faster other than this
	ARR=($(echo $l | tr " " "\n"))
	first=0;	
	for val in "${ARR[@]}"; do	
	    #first word is copy - check if it is really one
	    if [[ $first == 0 ]]; then
		first=1;
		if [[ ! -z ${FID2[$val]} ]]; then
		    feat_index=$((${FID2[$val]}*no_patients));
		else
		    #we can skip the rest of the line since we are not using this
		    break;
		fi
	    else
		if [[ "$val" -lt "0" ]]; then
		    G[$feat_index]=-1;
		elif [[ "$val" -gt "0" ]]; then
		    G[$feat_index]=1;
		fi
		feat_index=$((feat_index+1))
	    fi
	done
    done < $file\_CNs.filtered.txt
    ELAPSED_TIME=$((SECONDS - START_TIME))
    echo "copynew time: " $ELAPSED_TIME

    START_TIME=$SECONDS
    #PART3: Print the matrix
    for i in "${!PAT[@]}"; do
	patient_id=${PAT[$i]}
	if [[ $patient_id -gt -1 ]]; then
	    printf '%s ' $i 
	    index=$((patient_id*42));
	    for j in {0..41}; do
		printf '%d ' ${F[$index]}
		index=$((index+1));
	    done
	    index=$patient_id;
	    for j in {0..98}; do
		printf '%d ' ${G[$index]}
		index=$((index+no_patients));
	    done
	    printf '\n'
	fi 
    done > ../data/data.test.txt
#    sort $file.train.txt > $file.train.sorted.txt
    ELAPSED_TIME=$((SECONDS - START_TIME))
    echo "print time: " $ELAPSED_TIME
done
