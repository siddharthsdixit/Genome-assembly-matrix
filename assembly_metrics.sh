#!/bin/bash

USAGE="\nUSAGE: $0 [input_fasta] [base_name_for_output] [min_seq_len] > [output.txt]\n\nNOTE: Sequences less than [min_seq_len] will be ignored!\n"

IN=$1
NAME=$2

if test -n "$1" -a "$2" -a "$3"
then
	echo "START:"
	date
	echo -e "\n\nINPUTS:\nInput fasta: $1\nOutput base name: $2\n\n\nDEFINITIONS:\n1. Contig = each fasta sequence of the input file (contigs less than $3 bp are ignored)\n2. N50 = Length of the contig when 50 % of the assembly is achieved after summing the contig lengths sorted in decreasing order\n3. L50 = Serial number of the N50 contig when all the contigs are sorted in decreasing order of their lengths\n\n\nASSEMBLY STATS:"

	# Removing \n and sequences less than [min_seq_len] bases from the input fasta, and retaining fasta header upto first space or tab:
	awk '{RS="\n";ORS=""}{if($0~/>/){printf "\n"$0"\n"}else{printf $0}}' $IN |awk '{print $1}'| awk -v min=$3 'BEGIN {RS = ">" ; ORS = ""} length($2) >= min {print ">"$0}' > $NAME.unwrapped.fasta

	# No. of contigs
	NUM_CONTIGS=`grep -c ">" $NAME.unwrapped.fasta`
	echo "Total number of contigs = $NUM_CONTIGS"

	# Create sequence length distribution of the fasta file (print header until first space)
	sed 's/>/>\n>/' $NAME.unwrapped.fasta |awk '{RS=">";FS="\n"}{print $1"\t"length($2)}' | sort -t$'\t' -k2nr |head -n $NUM_CONTIGS > $NAME.dist
	
	# Rscript
	Rscript --vanilla /home/ssubha/samathmika/scripts/lendist.R $NAME.dist $NAME.dist.pdf 2 > err.log

	# Longest and shortest contigs
	LONGEST_CONTIG=`head -1 $NAME.dist|cut -f2`
	SHORTEST_CONTIG=`tail -1 $NAME.dist|cut -f2`
	echo -e "Longest contig = $LONGEST_CONTIG bp\nShortest contig = $SHORTEST_CONTIG bp"

	# Median contig length
	if [ $(($NUM_CONTIGS % 2)) -eq 0 ]
	then 
		MIDDLE_1=$[$NUM_CONTIGS/2]
		MIDDLE_1_LEN=`head -n $MIDDLE_1 $NAME.dist |tail -n1|cut -f2`
		MIDDLE_2=$[$MIDDLE_1+1]
		MIDDLE_2_LEN=`head -n $MIDDLE_2 $NAME.dist |tail -n1|cut -f2`
		MEDIAN_LEN=$(echo "scale=4;($MIDDLE_1_LEN+$MIDDLE_2_LEN)/2"|bc -l)
	else
		y=$(echo "$NUM_CONTIGS/2"|bc)
		MIDDLE=$[$y+1]
		MEDIAN_LEN=`head -n $MIDDLE $NAME.dist |tail -n1|cut -f2`
	fi
	echo "Median contig length = $MEDIAN_LEN bp"

	# Average contig length
	cut -f2 $NAME.dist|perl -MMath::Round -ne 'chomp(); push(@contigs,$_);$TOTAL_ASSEMBLY+=$_; $AVG_LEN=nearest(0.1111,$TOTAL_ASSEMBLY/@contigs);END{print "Average contig length = $AVG_LEN bp\n\n"}'

	# Total assembly size, N10-N95:
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.1){print "N10 contig length = $L bp\t\tL10 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.2){print "N20 contig length = $L bp\t\tL20 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.3){print "N30 contig length = $L bp\t\tL30 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.4){print "N40 contig length = $L bp\t\tL40 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.5){print "N50 contig length = $L bp\t\tL50 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.6){print "N60 contig length = $L bp\t\tL60 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.7){print "N70 contig length = $L bp\t\tL70 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.8){print "N80 contig length = $L bp\t\tL80 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.9){print "N90 contig length = $L bp\t\tL90 contig count = $count\n";exit;} ;}}'
	cut -f2 $NAME.dist|perl -ne 'chomp(); $count=0; push(@contigs,$_);$TOTAL_ASSEMBLY+=$_;END{foreach(sort{$b<=>$a}@contigs){$sum+=$_;$L=$_;$count++;if($sum>=$TOTAL_ASSEMBLY*0.95){print "N95 contig length = $L bp\t\tL95 contig count = $count\n\nAssembly size = $TOTAL_ASSEMBLY bp\n";exit;} ;}}'

	# Base composition of contigs
	x=$(grep -v ">" $NAME.unwrapped.fasta | awk '{RS="\n";ORS=""}{printf $0}' | wc -m)
	y=`grep -v ">" $NAME.unwrapped.fasta | fold -w1 | egrep -c "(C|G|c|g)"`
	iupac_nuc=`grep -v ">" $NAME.unwrapped.fasta | fold -w1 | egrep -c "(R|r|Y|y|S|s|W|w|K|k|M|m|B|b|D|d|H|h|V|v|N|n)"`
	atgc=`grep -v ">" $NAME.unwrapped.fasta | fold -w1 | egrep -c "(A|a|T|t|G|g|C|c)"`
	n=`grep -v ">" $NAME.unwrapped.fasta | fold -w1 | egrep -c "(N|n)"`
	GC_tot=$(echo "scale=4;$y/$x*100" | bc -l)
	GC=$(echo "scale=4;$y/$atgc*100" | bc -l)
	ATGC=$(echo "scale=4;$atgc/$x*100" | bc -l)
	iupac_nuc_percent=$(echo "scale=4;$iupac_nuc/$x*100" | bc -l)
	n_percent=$(echo "scale=4;$n/$x*100" | bc -l)
	echo "GC % = $GC_tot (of total assembly) and $GC (of ATGC content)"
	echo "Number of ATGCs = $atgc ($ATGC %)"
	echo "Number of IUPAC nucleotides = $iupac_nuc ($iupac_nuc_percent %)"
	echo "Number of Ns = $n ($n_percent %)"
	rm err.log $NAME.unwrapped.fasta

	echo -e "\n\nEND:"
	date
else
	echo -e $USAGE
fi
