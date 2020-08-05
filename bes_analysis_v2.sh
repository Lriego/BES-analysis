#!/bin/bash

############################################################################################
#Version: 6.0                                                                              #
#Update: 2020/08/03                                                                        #
#Authors: JM, CG and LG                                                                    #
#                                                                                          #
#      help:      ./my_pip_v6.0.sh -h	                                                   #
#		                                                                           #
############################################################################################

if [ $# -eq 0 ]; then
	echo "Error, No arguments and options are passed to the script. Please use the -h option to view the help for valid options"
	exit 0
fi

usage() {
	echo -e "\n*****BES grooming and alignment*****\nVersion: 6.0\nUpdate: 2020/08/03\nAuthors: JM, CG and LG"
	echo -e "\nThis program was design to prepare BAC-end Sequences (BESs) from Sanger reads and align them to a reference genome. The process is performed in five eligible stages:"
	echo -e "	1) 'abi2fastq' finds all Sanger reads in ab1 format of a directory and convert them into fastq format."
	echo -e "	2) 'qual' evaluates per base quality of the sequences a gives a report using FastQC."
	echo -e "	3) 'trim' trims and filters low quality BESs using Trimmomatic."
#	echo -e "	4) 'format' formats the fasta genome headers in order to be correctly read by the script."
	echo -e "	4) 'align' aligns them to a reference genome with one of the following methods:"
	echo -e "\t\t-BLAST\n\t\t-MEGA (megablast)\n\t\t-BLAT\n\t\t-NUC (nucmer)\n\t\t-BOW (bowtie)\n\t\t-BOW2 (bowtie2)\n\t\t-BWA (BWA)"
	echo -e "\nAfter the alignmet, the output of the methods that retrieves more than one alignment per sequence goes trough a filtering process in which the longest alignment is selected. In case of a tie, the firts alignment is selected. Then, the script converts any alignment output into gff3 format for easy visualization in the Integrative Genomics Viewer. Also, the script performs by default a screening of the BESs pairs according to their orientation in single-end, paired-end, unpaired-end, opposite, positive and negative. This former option can be disabled (see aditional options)."
	echo -e "\nUsage:\n\thelp:\t\t$0 -h\n\tabi2fastq:\t$0 -s abi2fastq -i INPUT.directory\n\tqual:\t\t$0 -s qual -i INPUT.file\n\ttrim:\t\t$0 -s trim -i INPUT.file\n\talign:\t\t$0 -s align -i INPUT.file -g GENOME.fasta -m METHOD"
	echo -e "\nAdditional options:\n\t-o OUTPUT_DIR\n\t-c HEADCROP (for trim)\n\t-p PHRED_CODING (for trim)\n\t-w WINDOW_SIZE (for trim)\n\t-l LENGTH (for trim)\n\t-t THRESHOLD (for trim)\n\t-e EVALUE (for blastn and megablast)\n\t-q SCREEN_BES (yes/no)"
	echo -e "\nIf input files are already in .fastq format, concatenate your sequences prior implementing the program.\n"
	echo -e "\nNew version modifications:\n\tBowtie and BWA aligners added. Format stage depreciated. Use 'format_genome.R' instead."
}

check() {
	if [ -f $output_file ]; then
		echo "File $output_file already exist, do you wish to overwrite it? [Y/n]"
		read answer
			if [ $answer = "N" ] || [ $answer = "n" ] || [ $answer = "no" ]; then
				echo "Run the program with a different output directory name."
				exit 0
			elif [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
				echo "File $output_file will be overwritten."
			fi
	fi
}

unique() {
	contender=none
	BESname=none
	pivot=0

	while IFS='' read -r line; do
		if [ $BESname != $(echo "$line" | cut -f1) ]; then
			echo "$contender" >> $output_file
			contender=$(echo "$line")
			BESname=$(echo "$line" | cut -f1)
			pivot=$(echo "$line" | cut -f4)
		elif [ $pivot -lt $(echo "$line" | cut -f4) ]; then
			contender=$(echo "$line")
			BESname=$(echo "$line" | cut -f1)
			pivot=$(echo "$line" | cut -f4)
		fi
	done < $OUTPUT/$(basename $output_file .out).all

	echo "$contender" >> $output_file
	sed -i '/none/d' $output_file
}

blast2gff() {
	paste <(cut -f9 $output_file | sed 's/.*[0-9]/&-/') <(cut -f10 $output_file) | sed "s|\t||g" > sum_line.temp

	while IFS='' read -r line; do
		sum=$(echo "$line" | bc)
		if [ $sum -lt 0 ]; then
			sum=+
			echo "$sum" >> orientation.temp
			echo "$line" | sed "s/-/\t/" >> coordenates.temp
		else
			sum=-
			echo "$sum" >> orientation.temp
			paste <(echo "$line" | sed "s/.*-//") <(echo "$line" | sed "s/-.*//") >> coordenates.temp
		fi
	done < sum_line.temp
	paste <(cut -f2 $output_file) <(cut -f1 $output_file | sed 's/.*/./g') <(cut -f1 $output_file) <(cut -f1,2 coordenates.temp) <(cut -f1 $output_file | sed 's/.*/./g') <(cut -f1 orientation.temp) <(cut -f1 $output_file | sed 's/.*/./g') <(cut -f1 $output_file | sed 's/.*/./g') > "$OUTPUT/$(basename "$output_file" .out).gff3"
	rm sum_line.temp orientation.temp coordenates.temp
}

sam2gff() {
	SEQUENCES=$(( $(grep ">" $GENOME | wc -l) +2 ))
	sed "1,"$SEQUENCES"d" $output_file > sam.temp
	cut -f10 sam.temp > bes.temp

	while IFS='' read -r bes; do
		echo ${#bes} >> bes_sizes.temp
	done < bes.temp

	paste <(cut -f4 sam.temp | sed 's/.*[0-9]/&+/') <(cut -f1 bes_sizes.temp) | sed "s|\t||g" > sum_line.temp

	while IFS='' read -r line; do
		sum=$(echo "$line" | bc)
		echo "$sum" >> sum_results.temp
	done < sum_line.temp

	paste <(cut -f3 sam.temp) <(cut -f1 sam.temp | sed 's/.*/./g' sam.temp) <(cut -f1 sam.temp) <(cut -f4 sam.temp) <(cut -f1 sum_results.temp) <(cut -f1 sam.temp | sed 's/.*/./g' sam.temp) <(cut -f2 sam.temp | sed 's/0/+/' | sed 's/16/-/') <(cut -f1 sam.temp | sed 's/.*/./g' sam.temp) <(cut -f1 sam.temp | sed 's/.*/./g' sam.temp) > "$OUTPUT/$(basename "$output_file" .out).gff3"
	rm sam.temp bes.temp bes_sizes.temp sum_line.temp sum_results.temp
}

screening() {
	if [ $SCREEN = "yes" ] || [ $SCREEN = "YES" ] || [ $SCREEN = "y" ] || [ $SCREEN = "Y" ]; then
		NAME=$METHOD
		if [ $NAME = NUC ]; then
			NAME="NUCMER"
		fi
		if [ $NAME = BOW ]; then
			NAME="BOWTIE"
		fi
		if [ $NAME = all ]; then
			NAME="BLAST"
		fi
		if [ $NAME = all.2 ]; then
			NAME="MEGA";
		fi
		if [ $NAME = all.3 ]; then
			NAME="BLAT"
		fi
		if [ $NAME = all.4 ]; then
			NAME="NUCMER"
		fi
		if [ $NAME = all.5 ]; then
			NAME="BOWTIE"
		fi
		echo "$NAME screening"
		./$(dirname $0)/screening.py "$OUTPUT"/"$(basename $output_file .out).gff3"
	fi
}

while getopts ":hs:i:o:p:w:l:m:g:e:t:c:q:" opt; do
	case $opt in
		h)
			usage; exit;;
		i)
			INPUT=$OPTARG;;
		o)
			OUTPUT=$OPTARG;;
		s)
			STAGE=$OPTARG;;
		p)
			PHRED=$OPTARG;;
		w)
			WINDOW=$OPTARG;;
		l)
			LENGTH=$OPTARG;;
		m)
			METHOD=$OPTARG;;
		g)
			GENOME=$OPTARG;;
		e)
			EVALUE=$OPTARG;;
		t)
			SCORE=$OPTARG;;
		c)
			CROP=$OPTARG;;
		q)
			SCREEN=$OPTARG;;
		:)
			echo "Option -$OPTARG requires an argument." >&2; exit 1;;
		\?)
			echo "Invalid option: -$OPTARG" >&2; exit 1;;
		*)
			echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
	esac
done

if [ -z $STAGE ]; then
	echo "Error, option no implemented or option -s (stage) no valid"
	exit 0
fi

if [ -z $INPUT ]; then
	echo "Error, -i (input file) is a mandatory option"
	exit 0
fi

if [ -z $OUTPUT ]; then
	OUTPUT="BES_out"
fi

if [ ! -d $OUTPUT ]; then
	mkdir $OUTPUT
fi

if [ -z $SCREEN ]; then
	SCREEN=yes
fi

#Perform program
while [ TRUE ]; do
case "$STAGE" in
	abi2fastq | all)
		output_file=$OUTPUT/all_bes.fastq
		check

		filename=$(basename "$INPUT")
		extension="${filename##*.}"

		if [[ $STAGE == "all" && $extension == "fastq" ]]; then
			STAGE="all.2"
			continue
		fi

		DIR_seqret=`which seqret`
		if [ -z $DIR_seqret ]; then
			echo "Error, seqret is not installed or is not included in your PATH"
			echo "You can install it with command: 'sudo apt-get install emboss'"
			exit 0
		fi
		mkdir temp
		for file in $(find $INPUT -type f -name '*.ab1')
		do
			seqret -sformat abi -osformat fastq -auto -stdout -sequence $file -outseq temp/"$(basename $file .ab1).fastq"
			sed -i "s|^@1|@$(basename $file .ab1)|" temp/"$(basename $file .ab1).fastq"
		done
		cat temp/*.fastq > $OUTPUT/all_bes.fastq
		rm -r temp
		if [ "$STAGE" == "all" ]; then
			INPUT="$OUTPUT"/"$(basename "$INPUT" .ab1).fastq"
			STAGE="all.2"
		else
			exit 1
		fi
	;;

	qual | all.2)
		output_file="$OUTPUT"/"$(basename "$INPUT" .fastq)_fastqc.html"
		check

		DIR_fastqc=`which fastqc`
		if [ -z $DIR_fastqc ]; then
			echo "Error, FastQC is not installed or is not included in your PATH"
			echo "You can install it with command: 'sudo apt-get install fastqc'"
			exit 0
		fi

		fastqc --noextract -t 4 -o "$OUTPUT" "$INPUT"

		if [ "$STAGE" == "all.2" ]; then
			STAGE="all.3"
		else
			exit 1
		fi
	;;

	trim | all.3)
		output_file="$OUTPUT"/"$(basename "$INPUT" .fastq)_trim.fastq"
		check

		DIR_trimmo=`which TrimmomaticSE`
		if [ -z ${DIR_trimmo} ]; then
			echo "Error, TrimmomaticSE is not installed or is not included in your PATH"
			echo "You can install it with command: 'sudo apt-get install trimmomatic'"
			exit 0
		fi
		if [[ ( -z $PHRED ) || ( $PHRED != "-phred33" ) || ( $PHRED != "-phred64" ) ]]; then
			echo "-p option was not specified or it does not contain a valid value. Default value '-phred33' was assigned"
			PHRED="-phred33"
		fi
		if [ -z $CROP ]; then
			echo "-c was not specified. Default value '91' was assigned"
			CROP=91
		fi
		if [ -z $WINDOW ]; then
			echo "-w option was not specified. Default value '4' was assigned"
			WINDOW=4
		fi
		if [ -z $SCORE ]; then
			echo "-t option was not specified. Default value '20' was assigned"
			SCORE=20
		fi
		if [ -z $LENGTH ]; then
			echo "-l option was not specified. Default value '35' was assigned"
			LENGTH=35
		fi
		TrimmomaticSE $PHRED "$INPUT" "$OUTPUT"/"$(basename "$INPUT" .fastq)_trim.fastq" HEADCROP:"$CROP" LEADING:"$SCORE" TRAILING:"$SCORE" SLIDINGWINDOW:"$WINDOW":"$SCORE" MINLEN:"$LENGTH"
		if [ "$STAGE" == "all.3" ]; then
			INPUT="$OUTPUT"/"$(basename "$INPUT" .fastq)_filt.fastq"
			STAGE="all.4"
		else
			exit 1
		fi
	;;

#	format)
#		if [ -z $GENOME ]; then
#			echo "Error, -g is a mandatory option when formating the headers of the genome. Please, use the -h option to view the help"
#		fi

#		output_file="$OUTPUT"/"$(basename "$GENOME" .fasta).fasta"
#		check

#		sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $GENOME |  sed 's/ /_/g' > temp.fasta
#		for seq in $(grep ">" temp.fasta); do
#			#num=$(grep "$seq" -A1 temp.fasta | sed -n '2~2p' | wc -c)
#			num=$(echo "$(grep "$seq" -A1 temp.fasta | sed -n '2~2p' | wc -c)-1" | bc)
#			echo $seq | sed "s/.*/&|size$num/" >> $OUTPUT/$(basename $GENOME .fasta).fasta
#			grep "$seq" -A1 temp.fasta | sed -n '2~2p' >> $OUTPUT/$(basename $GENOME .fasta).fasta
#		done
#		rm temp.fasta
#		exit 1

#	;;

	align | all.4)
		if [ -z $GENOME ]; then
			echo "Error, -g is a mandatory option when align or all stages are specified. Please, use the -h option to view the help"
		fi

		if [ ! -r $INPUT ]; then
			echo "Error, '$INPUT' is empty or it is not readable"
			exit 0
		fi

		if [ ${INPUT: -5} == fastq ]; then
			DIR_fq2fs=`which fastq_to_fasta`
			if [ -z $DIR_fq2fs ]; then
				echo "Error, fastq_to_fasta is not installed or it is not included in your PATH"
				echo "You can install it with the command: 'sudo apt-get install fastx-toolkit'"
				exit 0
			fi
			fastq_to_fasta -i $INPUT -o "$OUTPUT"/"$(basename "$INPUT" .fastq).fasta"
			INPUT="$OUTPUT"/"$(basename "$INPUT" .fastq).fasta"
		fi

		case "$METHOD" in
			BLAST | all)
				output_file="$OUTPUT"/"$(basename "$INPUT" .fasta)_blast.out"
				if [ -f $output_file ]; then
					echo "$output_file already exists, do you wish to overwrite it? [Y/n]"
					read answer
					if [ $answer = "N" ] || [ $answer = "n" ] || [ $answer = "no" ]; then
						echo "Run the program with a different output directory name."
						exit 0
					elif [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
						echo "File $output_file will be overwritten."
						rm $output_file
					fi
				fi

				DIR_blastn=`which blastn`
				if [ -z $DIR_blastn ]; then
					echo "Error, blastn is not installed or is not included in your PATH"
					echo "You can install it with the command: 'sudo apt-get install ncbi-blast+'"
					exit 0
				fi

				BLAST_DB="$(basename "$GENOME")"

				if [ -r "$GENOME.nin" ]; then
					echo "Database $BLAST_DB.nin for BLAST already exists"
				else
					echo "There is not a data base $BLAST_DB.nin for BLAST. Creating it from $GENOME file"
					#Modificado para nueva version de blast
                        		makeblastdb -dbtype nucl -in "$GENOME"
					#formatdb -i "$GENOME" -p F
				fi

				if [ -z $EVALUE ]; then
					echo "No e-value was specified, e-value of 1e-6 was assigned"
					EVALUE=0.000001
				fi
				echo "BLAST alignment"
				blastn -db "$GENOME" -query "$INPUT" -out "$OUTPUT"/"$(basename "$INPUT" .fasta)_blast.all" -evalue $EVALUE -outfmt 6 -max_target_seqs 1 -num_threads 2
				echo "BLAST filtering"
				unique
				echo "BLAST to gff3"
				blast2gff
				screening

				if [ $METHOD == "all" ]; then
					METHOD=all.2
				else
					echo "Done"
					exit 1
				fi
			;;

			MEGA | all.2)
				output_file="$OUTPUT"/"$(basename "$INPUT" .fasta)_mega.out"
				if [ -f $output_file ]; then
					echo "$output_file already exists, do you wish to overwrite it? [Y/n]"
					read answer
					if [ $answer = "N" ] || [ $answer = "n" ] || [ $answer = "no" ]; then
						echo "Run the program with a different output directory name."
						exit 0
						elif [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
						echo "File $output_file will be overwritten."
						rm $output_file
					fi
				fi

				DIR_mega=`which megablast`
				if [ -z $DIR_mega ]; then
					echo "Error, megablast is not installed or is not included in your PATH"
					echo "You can install it with the command: 'sudo apt-get install ncbi-blast+'"
					exit 0
				fi

				MEGA_DB="$(basename "$GENOME")"
				if [ -r "$GENOME.nin" ]; then
					echo "Database $MEGA_DB.nin for mega-BLAST already exists"
				else echo "There is not a data base $MEGA_DB.nin for mega-BLAST. Creating it from $GENOME file"
					formatdb -i "$GENOME" -p F
				fi

				if [ -z $EVALUE ]; then
					echo "No e-value was specified, e-value of 1e-200 was assigned"
					EVALUE=200
					NEW_EVALUE=$(bc <<< "scale=$EVALUE;1*10^-$EVALUE")
				fi

				if [ $METHOD == all.2 ]; then
					echo "No e-value was specified, e-value of 1e-200 was assigned"
					EVALUE=200
					NEW_EVALUE=$(bc <<< "scale=$EVALUE;1*10^-$EVALUE")
				fi
				echo "MEGA alignment"
				megablast -V F -d "$GENOME" -i "$INPUT" -o mega_file.temp -e "$NEW_EVALUE" -F F -D 3
				sed -i '/#/d' mega_file.temp
				#sed -i '1,4d' mega_file.temp
				sort mega_file.temp > $OUTPUT/$(basename $output_file .out).all
				rm mega_file.temp
				echo "MEGA filtering"
				unique
				echo "MEGA to gff3"
				blast2gff
				screening

				if [ $METHOD == "all.2" ]; then
                                        METHOD=all.3
                                else
					echo "Done"
                                        exit 1
                                fi
			;;

			BLAT | all.3)
				output_file="$OUTPUT"/"$(basename "$INPUT" .fasta)_blat.out"
				check

				DIR_blat=`which blat`
				if [ -z $DIR_blat ]; then
					echo "blat is not istalled or is not included in your PATH"
					echo "You can download it with the command: 'wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/'"
					exit 0
				fi
				echo "BLAT alignment"
				blat "$GENOME" "$INPUT" "$OUTPUT"/"$(basename "$INPUT" .fasta)_blat.all" -noHead -out=blast8
				echo "BLAT filtering"
				unique
				echo "BLAT to gff3"
				blast2gff
				screening

				if [ $METHOD == "all.3" ]; then
                                        METHOD=all.4
                                else
					echo "Done"
                                        exit 1
                                fi
			;;

			NUC | all.4)
				output_file="$OUTPUT"/"$(basename "$INPUT" .fasta)_nuc.out"
				check

				DIR_mummer=`which mummer`
				if [ -z $DIR_mummer ]; then
					echo "Error, mummer is not installed or is not included in your PATH"
					echo "You can install it with the command: 'sudo apt-get install mummer'"
					exit 0
				fi
				echo "NUCMER alignment"
				nucmer -p "$OUTPUT"/"$(basename $INPUT .fasta)_nuc" "$INPUT" "$GENOME"
				show-coords -B "$OUTPUT"/"$(basename $INPUT .fasta)_nuc.delta" > nucmer_file.temp
				paste <(cut -f 6 nucmer_file.temp) <(cut -f 1 nucmer_file.temp) <(cut -f12,13 nucmer_file.temp) <(cut -f14,15 nucmer_file.temp) <(cut -f9,10 nucmer_file.temp) <(cut -f7,8 nucmer_file.temp) <(cut -f14,15 nucmer_file.temp) | sort > "$OUTPUT"/"$(basename "$INPUT" .fasta)_nuc.all"
				rm nucmer_file.temp
				echo "NUCMER filtering"
				unique
				echo "NUCMER to gff3"
				blast2gff
				screening

				if [ $METHOD == "all.4" ]; then
                                        METHOD=all.5
                                else
					echo "Done"
                                        exit 1
                                fi
			;;

			BOW2 | all.5)
				output_file="$OUTPUT"/"$(basename "$INPUT" .fasta)_bow2.out"
				check

				DIR_bowtie2=`which bowtie2`
				if [ -z $DIR_bowtie2 ]; then
					echo "Error, bowtie2 is not installed or is not included in your PATH"
					echo "You can install it with the command: 'sudo apt-get install bowtie2'"
					exit 0
				fi

				#Check if index for bowtie2 exists
				if [ ! -r $(dirname $GENOME)/bt2_indexes/$(basename $GENOME).1.bt2 ]; then
					echo "$OUTPUT/$(basename $GENOME .1.bt2)"
					echo "WARNING!"
					echo "Bowtie2 index of reference genome does not exist"
					echo "Building bowtie2 index..."
					mkdir $(dirname $GENOME)/bt2_indexes
					# Building bowtie2 index of the reference genome
					bowtie2-build -f $GENOME $(dirname $GENOME)/bt2_indexes/$(basename $GENOME) -q
				else
					echo "Bowtie2 index already exists"
				fi

				#Mapping with bowtie2
				echo "BOWTIE2 alignment"
				bowtie2 -f -x $(dirname $GENOME)/bt2_indexes/$(basename $GENOME) -U "$INPUT" -S "$OUTPUT"/"$(basename "$INPUT" .fasta)_bow2.out" --no-unal

				#Convert sam to gff3
				echo "BOWTIE2 to gff3"
				sam2gff
				screening

				if [ $METHOD == "all.5" ]; then
                                        METHOD=all.6
                                else
					echo "Done"
                                        exit 1
                                fi
			;;

			BOW | all.6)
				output_file="$OUTPUT"/"$(basename "$INPUT" .fasta)_bow.out"
				check

				DIR_bowtie=`which bowtie`
				if [ -z $DIR_bowtie ]; then
					echo "Error, bowtie is not installed or is not included in your PATH"
					echo "You can install it with the command: 'sudo apt-get install bowtie'"
					exit 0
				fi

				#Check if index for bowtie exists
				if [ ! -r $(dirname $GENOME)/bt1_indexes/$(basename $GENOME).1.ebwt ]; then
					echo "$OUTPUT/$(basename $GENOME .1.ebwt)"
					echo "WARNING!"
					echo "Bowtie index of reference genome does not exist"
					echo "Building bowtie index..."
					mkdir $(dirname $GENOME)/bt1_indexes
					# Building bowtie index of the reference genome
					bowtie-build -f $GENOME $(dirname $GENOME)/bt1_indexes/$(basename $GENOME) -q
				else
					echo "Bowtie index already exists"
				fi

				#Mapping with bowtie
				echo "BOWTIE alignment"
				# although default for SSPACE is -v 0
				bowtie -v 2 -f $(dirname $GENOME)/bt1_indexes/$(basename $GENOME) "$INPUT" -S "$OUTPUT"/"$(basename "$INPUT" .fasta)_bow.out" --no-unal

				#Convert sam to gff3
				echo "BOWTIE to gff3"
				sam2gff
				screening

				if [ $METHOD == "all.6" ]; then
                                        METHOD=all.7
                                else
					echo "Done"
                                        exit 1
                                fi
			;;

			BWA | all.7)
				output_file="$OUTPUT"/"$(basename "$INPUT" .fasta)_bwa.out"
				check

				DIR_bwa=`which bwa`
				if [ -z $DIR_bwa ]; then
					echo "Error, BWA is not installed or is not included in your PATH"
					echo "You can get it with the command: 'sudo apt-get install bwa'"
					exit 0
				fi

				DIR_sam=`which samtools`
				if [ -z $DIR_sam ]; then
					echo "Error, samtools is not installed or is not included in your PATH"
					echo "You can get it with the command: 'sudo apt-get install samtools'"
					exit 0
				fi

				#Check if index for BWA exists
				if [ ! -r $(dirname $GENOME)/bwa_indexes/$(basename $GENOME).amb ]; then
					echo "$OUTPUT/$(basename $GENOME .amb)"
					echo "BWA index of reference genome does not exist"
					echo "Building BWA index..."
					mkdir $(dirname $GENOME)/bwa_indexes
					bwa index $GENOME
					mv $GENOME.amb $(dirname $GENOME)/bwa_indexes
					mv $GENOME.ann $(dirname $GENOME)/bwa_indexes
					mv $GENOME.bwt $(dirname $GENOME)/bwa_indexes
					mv $GENOME.pac $(dirname $GENOME)/bwa_indexes
					mv $GENOME.sa $(dirname $GENOME)/bwa_indexes

				else
					echo "BWA index already exists"
				fi

				#Mapping with bowtie
				echo "BWA alignment"
				bwa mem $(dirname $GENOME)/bwa_indexes/$GENOME $INPUT > "$OUTPUT"/"$(basename "$INPUT" .fasta)_bwa.all"
				grep "@" "$OUTPUT"/"$(basename "$INPUT" .fasta)_bwa.all" > "$OUTPUT"/"$(basename "$INPUT" .fasta)_bwa.out"
				# Filter unmapped, secondary and supplementary alignments (i.e. 4 + 256 + 2048 = 2308)
				samtools view -F 2308 "$OUTPUT"/"$(basename "$INPUT" .fasta)_bwa.all" >> "$OUTPUT"/"$(basename "$INPUT" .fasta)_bwa.out"

				#Convert sam to gff3
				echo "BWA to gff3"
				sam2gff
				screening

				echo "Done"
				exit 1
			;;

			*)
				echo "Error, -m (method for alignment) is a mandatory option when align or all stages are specified.  Please use the -h option to view the help for valid options"
				exit 0
			;;
		esac
	;;

	*)
		echo "Error, invalid argumento to -s option. Please use the -h option to view the help for valid options"
		exit 0
	;;
esac
done
