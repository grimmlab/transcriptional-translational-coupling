
### Workflow for operon classfication

# Root folder name"
#ROOT_FOLDER_NAME=OperonProject


main(){
	create_folders
	set_variables # -> Never comment this function
	run_fetch_door2_and_ncbi_data
	run_classification_code
	run_concatenate_all_classified_files
	run_modify_concatenate_all_classified_files
	run_occurrence_based_ranking
	run_cooccurrence_based_gene_ranking
	run_cooccurrence_based_functional_ranking
	run_gene_motif_search
	run_extract_sequences_for_phylogenetic_analyses
	run_compress_fasta_for_phylogenetic_analyses
	#run_MSA_analyses
}

#create folder structure

create_folders(){
	echo "Creating folders..."
	for FOLDER in analyses data
		do
			mkdir -p ${ROOT_FOLDER_NAME}/${FOLDER}
		done
	echo "DONE creating folders..."

}

# setting variable path
set_variables(){
	echo "Setting variables for paths..."
	#export ROOT_FOLDER_NAME
	export DATA_FOLDER=$(pwd)/data
	export ANALYSES_FOLDER=$(pwd)/analyses
	export BIN_FOLDER=$(pwd)/bin
	echo "DONE setting variables for paths!"
}


run_fetch_door2_and_ncbi_data(){
	echo "Fetching operon table and genbank using restapi..."
	cd ${BIN_FOLDER}
	unzip ftp_path.zip
	python3 fetch_operon_genbank_restapi.py ${DATA_FOLDER}
	echo "DONE fetching operon table and genbank using restapi!"
}


run_classification_code(){
	echo "Running operon classification"
	for DIR in $(ls ${DATA_FOLDER})
  	do
		FILE1=$(ls ${DATA_FOLDER}/${DIR}/GCA*.txt)
		FILE2=$(ls ${DATA_FOLDER}/${DIR}/*NC_*[0-9].txt)
		#echo $FILE1 $FILE2
	 	cd ${BIN_FOLDER}
	 	python3 oprfile_operon_classification_alldoor2_v2.py ${FILE1} ${FILE2}
		cd ..
   	done
	echo "DONE running operon classification!"
}


run_concatenate_all_classified_files(){
	echo "Running to combine classified files"
   	OP_HEADER=$(cat data/Acaryochloris_marina_MBIC11017/Acaryochloris_marina_MBIC11017_NC_009925_operon_classification_output.txt | head -n 1 | cut -f2-)
   	#printf $OP_HEADER>  ${ANALYSES_FOLDER}/Final_combined_files.txt
   	echo Filename '\t' function '\t' operon '\t' classification '\t' gene_symbol '\t' COG_number >  ${ANALYSES_FOLDER}/Final_combined_files.txt
	for DIR in $(ls ${DATA_FOLDER}/*/*classification_output.txt)
	do
		#echo ${DIR}
      		FILENAME=$(basename $DIR | sed 's/_operon_classification_output.txt//g')
      		#echo $FILENAME
		# grep "both" ${DIR} >> ${ANALYSES_FOLDER}/Final_combined_files.txt
		TAB=`echo 'x' | tr 'x' '\011'`
	      	grep "both" $DIR | while IFS=$TAB read -r col1 col2 col3 col4 col5 col6; do
	      	echo $FILENAME '\t' $col2 '\t' $col3 '\t' $col4 '\t' $col5 '\t' $col6 >> ${ANALYSES_FOLDER}/Final_combined_files.txt;	done
   done
	cat ${ANALYSES_FOLDER}/Final_combined_files.txt | cut -f1 | uniq > ${ANALYSES_FOLDER}/genome_list_containing_both.txt
	echo "DONE running to combine classified files"
}


run_modify_concatenate_all_classified_files(){
	echo "Running to modifying the combined classified files"
	python3 ${BIN_FOLDER}/combined_file_manipulation.py \
	${ANALYSES_FOLDER}/Final_combined_files.txt
	echo "DONE running to modifying the combined classified files"
}


run_occurrence_based_ranking(){
	echo "Running occurrence based ranking"
	mkdir -p ${ANALYSES_FOLDER}/occurrence
	mkdir -p ${ANALYSES_FOLDER}/occurrence/plots
	python3 ${BIN_FOLDER}/occurrence_based_ranking_analysis.py \
			${ANALYSES_FOLDER}/Final_combined_files_corrected.txt \
			${ANALYSES_FOLDER}/occurrence/plots \
			${ANALYSES_FOLDER}/occurrence
	echo "DONE running occurrence based ranking"
}


run_cooccurrence_based_gene_ranking(){
	echo "Running co-occurrence gene based ranking"
	mkdir -p ${ANALYSES_FOLDER}/co-occurrence
	python3 ${BIN_FOLDER}/cooccurrence_gene_based_ranking_analysis.py \
         	${ANALYSES_FOLDER}/Final_combined_files_corrected.txt \
        	${ANALYSES_FOLDER}/co-occurrence
   	echo "DONE running co-occurrence based gene ranking"
}


run_cooccurrence_based_functional_ranking(){
	echo "Running co-occurrence based functional ranking"
	mkdir -p ${ANALYSES_FOLDER}/co-occurrence
   	mkdir  ${ANALYSES_FOLDER}/co-occurrence/top_genome_list
	python3 ${BIN_FOLDER}/cooccurrence_func_based_ranking_analysis.py \
         	${ANALYSES_FOLDER}/Final_combined_files_corrected.txt \
         	${ANALYSES_FOLDER}/co-occurrence
   	echo "DONE running co-occurrence based functional ranking"
}


run_gene_motif_search(){
  	echo "Running to gene motif search "
  	rm ${ANALYSES_FOLDER}/motif_counts_for_genes.txt
	mkdir -p ${ANALYSES_FOLDER}/gene_motif
	echo "">${ANALYSES_FOLDER}/motif_counts.txt
	while read p; do
		echo Running for motif \"${p}\"
		python3 ${BIN_FOLDER}/gene_comb_count.py \
	    		${ANALYSES_FOLDER}/Final_combined_files_corrected.txt \
	    		${ANALYSES_FOLDER}/gene_motif \
	    		"${p}" >> ${ANALYSES_FOLDER}/motif_counts_for_genes.txt
	  	done <${ANALYSES_FOLDER}/gene_motifs_list.txt
			find ${ANALYSES_FOLDER}/gene_motif/ -type f  -ctime -1 |
		while read fname; do
	       		sz=`cat $fname | wc -l`   # Not a UUOC done to get just a line count
	       		[ $sz -lt 2 ] && rm $fname
		done
	echo "DONE running to gene motif search "
}


run_extract_sequences_for_phylogenetic_analyses(){
	echo "Running to extract 16S sequences of genomes"
	rm -rf ${ANALYSES_FOLDER}/gene_motif/16S_sequences
	wget -qO- https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.10.1+-x64-linux.tar.gz | tar -xvz -C ${BIN_FOLDER}
	mkdir -p ${DATA_FOLDER}/16S_DB
	wget -qO- ftp://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz | tar -xvz -C ${DATA_FOLDER}/16S_DB
	mkdir -p ${ANALYSES_FOLDER}/gene_motif/16S_sequences

	for FILE in ${ANALYSES_FOLDER}/gene_motif/*_genome_list.txt
		do
			FILENAME=$(basename ${FILE} | sed -e "s/_genome_list.txt/_genome/")

			#echo "Acidaminococcus fermentans DSM 20731" | grep -F "$(<${ANALYSES_FOLDER}/gene_motif/tmp.txt)"
   			#grep -F $(<analyses/gene_motif/tmp.txt)
			cat ${FILE} | sed 1d | awk -F'\t' '{print $1}'  >  ${ANALYSES_FOLDER}/gene_motif/tmp.txt
			TMPGRP=$(cat ${ANALYSES_FOLDER}/gene_motif/tmp.txt)
			${BIN_FOLDER}/ncbi-blast-2.10.1+/bin/blastdbcmd \
				-db ${DATA_FOLDER}/16S_DB/16S_ribosomal_RNA \
				-entry all -outfmt "%g;;%t" | grep -F "${TMPGRP}" | awk -F";;" '/16S/{print $1}' | ${BIN_FOLDER}/ncbi-blast-2.10.1+/bin/blastdbcmd \
				-db ${DATA_FOLDER}/16S_DB/16S_ribosomal_RNA \
				-entry_batch - \
				-out ${ANALYSES_FOLDER}/gene_motif/16S_sequences/${FILENAME}_sequences.fasta
		done
	rm ${ANALYSES_FOLDER}/gene_motif/tmp.txt
	echo "DONE running to extract 16S sequences of genomes"
}


run_compress_fasta_for_phylogenetic_analyses(){
	echo "Running to compress extract 16S sequences of genomes at genus level"
	mkdir -p  ${ANALYSES_FOLDER}/gene_motif/16S_sequences_genuslevel
	#cd ${ANALYSES_FOLDER}/gene_motif/16S_sequences_genuslevel
	for FILE in ${ANALYSES_FOLDER}/gene_motif/16S_sequences/*.fasta
		do
			echo "running on $(basename ${FILE}) file"
			python3 ${BIN_FOLDER}/compressed_fasta_genuslevel.py ${FILE} ${ANALYSES_FOLDER}/gene_motif/16S_sequences_genuslevel
	done
	echo "Running to compress extract 16S sequences of genomes at genus level"
}


#run_MSA_analyses
run_MSA_analyses(){
	rm mafft-7.471-with-extensions-src.tgz
	wget https://mafft.cbrc.jp/alignment/software/mafft-7.471-with-extensions-src.tgz

	tar xfvz mafft-7.471-with-extensions-src.tgz

   	mkdir -p ${BIN_FOLDER}/mafft

	sed -i "s#PREFIX = /usr/local#PREFIX = ${BIN_FOLDER}/mafft#g" mafft-7.471-with-extensions/core/Makefile

	cd mafft-7.471-with-extensions/core

	make
	make install
}



main
