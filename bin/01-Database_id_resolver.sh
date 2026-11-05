#!/bin/bash

#******************************************************************************
#
#   01-Database_id_resolver.sh
#
#   This script processes miRNA data from three major databases: miRBase,
#   PmiREN, and sRNAanno. It takes raw data from these databases and
#   performs the following tasks:
#
#   1. Filters and prepares the miRNA sequences for both mature and
#      hairpin types from the miRBase database, specifically focusing
#      on plant species.
#
#   2. Retrieves and prepares the species identifiers for each database
#      (miRBase, PmiREN, and sRNAanno), ensuring consistency across all
#      databases. This process avoids identifier mismatches that may
#      hinder accurate miRNA annotation.
#   
#   4. Outputs the final datasets as FASTA files for each database 
#      and for each type of miRNA sequence (mature and hairpin), 
#      with consistent species identifiers across all the databases.
#
#   This is crucial for ensuring the integrity of the annotation process 
#   and providing accurate, reliable data for downstream analyses.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 02/01/2024
#   Version: 1.0
#
#******************************************************************************


### FUNCTIONS

#################################################
#
#   This function contains the usage message
#
#################################################

usage() {
    
    help='''
    usage: miRNAs_annotation.sh [options] ...
    options:
        -h/--help               help
        -i/--inputsps           List with the names of the input species (e.g.
                                Arabidopsis thaliana).
        -m/--mirbase            Path to the FASTA file containing the mature
                                miRNA sequences from miRBase.
        -s/--srnaanno           Path to the GFF3 file containing the mature
                                and hairpin miRNA sequences from sRNAanno.
        -p/--pmiren             Path to the FASTA file containing the mature
                                miRNA sequences from PmiREN.
                                
        -a/--mhairpin           Path to the FASTA file containing the hairpin
                                miRNA sequences from miRBase.
        -n/--phairpin           Path to the FASTA file containing the hairpin
                                miRNA sequences from PmiREN.
        -v/--mirbaseplants      Path to the CSV file containing the plant
                                species and their identifiers used in miRBase.
                                Required to select only the sequences of plant
                                species from miRBase.
    '''
    echo "$help"
    exit 0
}

#####################################################
#
#   This function manages the input of arguments
#   to the program.
#
#####################################################

arguments_management() {

    # Default values
    mirbase=''
    srnaanno=''
    pmiren=''
    mirbase_hairpin=''
    pmiren_hairpin=''

    # Read the options
    TEMP=$(getopt -o h::i:m:s:p:a:n: --long help::,inputsps:,mirbase:,srnaanno:,pmiren:,mhairpin:,phairpin: -- "$@")

    # Check if the arguments are valid
    VALID_ARGUMENTS=$?
    if [ "$VALID_ARGUMENTS" != 0 ];
    then
        usage
    fi

    # Convert options to script arguments
    eval set -- "$TEMP"

    # Extract options and their arguments into variables.
    while true; do
        case "$1" in
            -h|--help)
                usage ;;
            -i|--inputsps)
                input_species="$2"; shift 2 ;;
            -m|--mirbase)
                if [ -n "$2" ] ; then
                    mirbase="$2"
                fi
                shift 2 ;;
            -s|--srnaanno)
                if [ -n "$2" ] ; then
                    srnaanno+=" $2"
                fi
                shift 2 ;;
            -p|--pmiren)
                if [ -n "$2" ] ; then
                    pmiren+=" $2"
                fi
                shift 2 ;;
            -a|--mhairpin)
                if [ -n "$2" ] ; then
                    mirbase_hairpin+=" $2"
                fi
                shift 2 ;;
            -n|--phairpin)
                if [ -n "$2" ] ; then
                    pmiren_hairpin+=" $2"
                fi
                shift 2 ;;
            # -- meands the end of the arguments; drop this, and break out the while loop
            --) shift ; break ;;
            # If invalid options were passed...
            *) echo "Unexpected option: $1 - this should not happen."
                usage ;;
        esac
    done
}


#####################################################
#
#   This function reads the FASTA files of various
#   plant species downloaded from PmiREN and
#   generates a two-column table. In the table,
#   the first column represents the identifiers
#   used in this database to refer to the species,
#   and the second column represents the scientific
#   names of the respective species.
#
#   Arguments:
#       Directory path where the PmiREN database
#           is located
#       Output file path for the identifiers
#           table (CSV)
#       Type of sequences (mature or hairpin)
#
#####################################################

get_PmiREN_species_ids (){

    # Arguments
    local pmiren_files_list="${1}"
    local pmiren_species_id_table="${2}"
    local type="${3}"

    # Iterate through FASTA files
    for file in $pmiren_files_list; do
        # Extract the identifier from the first line of the sequence
        id=$(grep -Eo '^>[^[:space:]]+' "$file" | head -n 1 | cut -c2- | awk -F'-' '{print $1}')

        # Get the file name without the mature.fa or hairpin.fa suffix and perform the split
        species_name=$(basename "$file" | sed "s/_/ /g; s/${type}.fa//" | sed 's/[[:space:]]*$//')

        # Print the result into the table
        echo -e "$id,$species_name" >> "$pmiren_species_id_table"
    done
}


#####################################################
#
#   This function extracts species identifiers and
#   complete organism names from a miRBase FASTA file.
#
#   It extracts the species identifier from the FASTA
#   sequence ID and the organism name from the fields
#   following the miRBase accession number.
#
#   For mature miRNAs (MIMAT), the last field contains
#   the miRNA name and is excluded.
#
#   For precursor miRNAs (MI), the last two fields
#   contain the miRNA name and "stem-loop", and
#   are excluded.
#
#   Duplicate identifier-organism pairs are removed.
#
#   Arguments:
#       FASTA file containing miRBase sequences.
#       Output CSV file for species identifiers.
#       Temporary directory for sorting.
#
#   Output:
#       CSV file with two columns:
#           species_id,species_name
#
#####################################################

get_mirbase_species_ids() {

    local mirbase_fasta_path="${1}"
    local mirbase_species_id_table="${2}"
    local temporary_dir="${3}"

    # Extract species identifiers and complete organism names
    awk '
        /^>/ {
            # Extract species identifier
            split($1, id_parts, "-")
            id = substr(id_parts[1], 2)

            # Identify the last organism-name field
            if ($2 ~ /^MIMAT[0-9]+$/) {
                last = NF - 1
            }
            else if ($2 ~ /^MI[0-9]+$/ && $NF == "stem-loop") {
                last = NF - 2
            }
            else {
                next
            }

            # Reconstruct the complete organism name
            species = ""
            for (i = 3; i <= last; i++) {
                species = species (species == "" ? "" : " ") $i
            }

            # Print valid identifier-organism pairs
            if (id != "" && species != "") {
                print id "," species
            }
        }
    ' "$mirbase_fasta_path" \
        | sort -T "$temporary_dir" -u \
        > "$mirbase_species_id_table"
}


#####################################################
#
#   This function merges and compares two CSV tables 
#   based on specified key columns. It produces two 
#   output files:
#   1. A unified table with species and their final 
#      identifiers.
#   2. A table highlighting species with mismatched 
#      identifiers between the two input tables.
#
#   Arguments:
#       table1: Path to the first input CSV file.
#       table2: Path to the second input CSV file.
#       key1: Key column index in the first table 
#             for the join operation.
#       key2: Key column index in the second table 
#             for the join operation.
#       output_cols: Columns to include in the 
#                    merged output (join format).
#       output_dir_path_species: Path to the output 
#                                file with unified
#                                species.
#       output_dir_path_diff: Path to the output file 
#                             with mismatched species.
#
#####################################################

merge_and_compare_tables(){

    # Arguments
    local table1="${1}"
    local table2="${2}"
    local key1="${3}"
    local key2="${4}"
    local output_cols="${5}"
    local output_dir_path_species="${6}"
    local output_dir_path_diff="${7}"
    local temporary_dir="${8}"

    # Create temporary directory
    mkdir -p $temporary_dir"/tmp_mc"

    # Sort files to be joined    
    LANG=en_EN sort -T $temporary_dir"/tmp_mc" -t',' -f -k "${key1},${key1}" "$table1" > "$temporary_dir/tmp_mc/table1.csv"
    LANG=en_EN sort -T $temporary_dir"/tmp_mc" -t',' -f -k "${key2},${key2}" "$table2" > "$temporary_dir/tmp_mc/table2.csv"

    # Alignment
    LANG=en_EN join -1 $key1 -2 $key2 -t',' -i -o $output_cols -e "NULL" -a 1 -a 2 $temporary_dir"/tmp_mc/table1.csv" $temporary_dir"/tmp_mc/table2.csv" > $temporary_dir"/tmp_mc/pre_results.csv"

    # Merge the columns with the species names into one.
    awk -F',' '{ if ($1 == "NULL") { print $2 "," $3 "," $4 } else { print $1 "," $3 "," $4 } }' $temporary_dir"/tmp_mc/pre_results.csv" > $temporary_dir"/tmp_mc/results.csv"

    # Add a final column with the ultimate identifier that will be assigned to each species (use the mirbase identifier unless it doesn't have one)
    awk 'BEGIN {FS=OFS=","} {print $1, $2, $3, ($2 == "NULL") ? tolower($3) : tolower($2)}' $temporary_dir"/tmp_mc/results.csv" > $output_dir_path_species

    # Verify which species do not have the same identifier in both databases.
    awk -F',' '{ if ($2 != "NULL" && $3 != "NULL" && tolower($2) != tolower($3)) { print } }' $output_dir_path_species > $output_dir_path_diff

    # Delete temporary directory
    rm -r $temporary_dir"/tmp_mc"
}


#####################################################
#
#   This function replaces mismatched identifiers 
#   in a FASTA file with new identifiers based on 
#   a provided mapping table in CSV format. 
#   It updates the identifiers in the FASTA file 
#   and generates a new output file with the changes.
#
#   Arguments:
#       Path to the input FASTA file.
#       Path to the output FASTA file with
#           updated identifiers.
#       Column index in the CSV table containing
#           the current identifiers.
#       Column index in the CSV table containing
#           the new identifiers.
#       Path to the CSV file mapping old
#           identifiers to new identifiers.
#
#####################################################

replace_mismatched_ids(){

    # Arguments
    local fasta_in_path="${1}"
    local fasta_out_path="${2}"
    local database_id_col="${3}"
    local new_id_col="${4}"
    local ids_table_csv="${5}"
    local temporary_dir="${6}"

    # Create temporary directory
    mkdir -p ./tmp_replace

    # Get ids from fasta file
    grep '^>' $fasta_in_path | awk -F'[>-]' '{print $2}' | sort -T "$temporary_dir" | uniq > ./tmp_replace/ids.txt

    # Create a file with new ids
    awk -F, -v database_id_col="$database_id_col" -v new_id_col="$new_id_col" 'NR==FNR { ids[$database_id_col]=$new_id_col; next } { if (tolower($1) in ids) $1=ids[tolower($1)]; print tolower($0)}' "$ids_table_csv" "./tmp_replace/ids.txt" > "./tmp_replace/new_ids.txt"

    # Create old_ids-new_ids table
    paste -d',' ./tmp_replace/ids.txt ./tmp_replace/new_ids.txt > ./tmp_replace/old_new_ids.txt

    # Replace old identifiers with new ones
    while IFS=, read -r old_id new_id; do
        sed -i "s/>${old_id}-/>${new_id}-/" "$fasta_in_path"
    done < ./tmp_replace/old_new_ids.txt

    # Create output fasta file
    cp $fasta_in_path $fasta_out_path

    # Delete temporary directory
    rm -rf ./tmp_replace

}


#####################################################
#
#   This function generates a unique identifier
#   for a given species name by using the first
#   letter of the genus and the first two letters
#   of the specific epithet. If the generated
#   identifier already exists in the provided
#   species ID file, the function iteratively 
#   modifies the identifier until it is unique.
#
#   Arguments:
#       The scientific name of the species 
#           (genus and specific epithet).
#       Path to the CSV file containing
#           existing species identifiers.
#       Column index in the CSV file where
#           the identifiers are stored.
#
#####################################################

change_the_id_using_name(){

    # Arguments
    local species_name="${1}"
    local species_ids_file="${2}"
    local id_col="${3}"

    # Get the two words from the species name
    genus=$(echo "$species_name" | awk '{print $1}')
    specific_epithet=$(echo "$species_name" | awk '{print $2}')

    # Generate the initial code
    new_id=$(echo "${genus:0:1}${specific_epithet:0:2}" | tr '[:upper:]' '[:lower:]')

    # Auto-increment variables
    increment=0
    increment2=0
    sl2_beginning=1

    # Check if the code already exists in the id_table.csv
    while awk -F',' -v col="$id_col" -v code="$new_id" '{ if ($col == code) print }' "$species_ids_file"  | grep -q .; do  
        # Check if the variable 'increment' is equal to the length of the specific epithet
        if [ $increment -eq $(( ${#specific_epithet} - $sl2_beginning )) ]; then
            # Update variables
            let "increment2++"
            increment=0
            sl2_beginning=0
        fi

        # Modify the code if it already exists
        specific_epithet_letter1="${specific_epithet:0+increment2:1}"
        specific_epithet_letter2="${specific_epithet:sl2_beginning+increment:1}"

        # Generate the new code
        new_id=$(echo "${genus:0:1}${specific_epithet_letter1}${specific_epithet_letter2}" | tr '[:upper:]' '[:lower:]')

        # Increment the variable
        let "increment++"
    done

    echo $new_id
}


#####################################################
#
#   This function prepares the sRNAanno
#   database for miRNA annotation by processing
#   GFF3 files from the sRNAanno database. It
#   generates FASTA files containing miRNA hairpin
#   and mature sequences. Additionally, it ensures
#   that the species identifiers in the provided
#   ID file are consistent with previously
#   established identifiers.
#
#   Arguments:
#       Path to the CSV file containing
#           species names and their identifiers.
#       List of paths to GFF3 files
#           from the sRNAanno database with miRNA annotations.
#       Directory where the output FASTA 
#           files will be saved.
#       Path to the output file with 
#           updated species identifiers.
#
#####################################################

prepare_sRNAanno_database(){

    # Arguments
    local species_id_file="${1}"
    local species_gff3_list="${2}"
    local path_dir_out="${3}"
    local path_file_ids_out="${4}"

    # Iterate through the files
    for gff3_file_path in $species_gff3_list; do 

        ## 1. Get the species id
        ################################################################################

        # Get the file name without extension
        file_name_without_extension=$(basename "$gff3_file_path" | cut -d'.' -f1)

        # Get the species name from file name
        species_name="${file_name_without_extension//_/ }"
        
        # Retrieve the line from the IDs file where the species is located and generate a new line
        old_line=$(awk -F',' -v species_name="$species_name" '{ if($1 == species_name) { print $0 }}' "$species_id_file")
        new_line=$(awk -F',' -v species_name="$species_name" '{ if($1 == species_name) { print $0 "," $4 }}' "$species_id_file") # Create a new column using the identifier from the fourth one

        # If the old line exists...
        if [ -n "$old_line" ]; then
            # Generate a new column using the specified identifier for the two previous databases
            sed -i "s|$old_line|$new_line|g" "$species_id_file"

            # Obtain identifier from the original line
            species_id=$(echo "$old_line" | cut -d',' -f4)

        # If it does not exists...
        else
            # Column where the identifiers are located.
            col_id=4
            # Create a new identifier using the species name
            species_id=$(change_the_id_using_name "$species_name" $species_id_file $col_id)

            # Add a new line to the identifier file
            echo "$species_name,NULL,NULL,$species_id,$species_id" >> "$species_id_file"
        fi


        ## 2. Create fasta file using gff3 files obtained from sRNAanno database
        ################################################################################

        # Create output directory
        mkdir -p $path_dir_out

        # Select those rows corresponding to Known miRNAs (mature)
        mature_miRNA_known=$( cat $gff3_file_path | grep -v '^#' | grep Known )

        # Open output files once
        exec 3>>"$path_dir_out/srnaanno_hairpin.fa"
        exec 4>>"$path_dir_out/srnaanno_mature.fa"

        # Read from the variable directly
        while IFS=$'\t' read -r scaffold name type start end point1 strand point2 seq_information; do
            
            # Skip empty records
            [ -z "$seq_information" ] && continue
            
            # Extract ID and sequence using parameter expansion
            id_field="${seq_information%%;seq=*}"
            id="${id_field##*ID=}"
            seq_field="${seq_information##*seq=}"
            
            if [[ $type == 'miRNA_primary_transcript' ]]; then
                # For precursor sequences
                seq_header="${id%%-*}"
                echo ">$species_id-$seq_header" >&3
                echo "$seq_field" >&3
            else
                # For mature miRNA sequences
                seq_header="${id%-*}-${id#*-}"
                echo ">$species_id-$seq_header" >&4
                echo "$seq_field" >&4
            fi
        done <<< "$mature_miRNA_known"

        # Close output file descriptors
        exec 3>&-
        exec 4>&-

    done

    # Add the last column to the ID file in those rows that have not been modified
    awk -F',' '{ if ($5 == "") { print $0",NULL"; } else { print $0; } }' $species_id_file > $path_file_ids_out
    
    # Change column order. Place column 4 in the position of column 5 and vice versa.
    awk 'BEGIN {FS=OFS=","} {temp=$4; $4=$5; $5=temp} 1' $path_file_ids_out > tmp.csv && mv tmp.csv $path_file_ids_out
}


#####################################################
#
#   Resolves duplicated species identifiers.
#   The first species (in table order) keeps its
#   identifier. Any later species with the same one
#   gets a new one using the same strategy as
#   change_the_id_using_name: first letter of the
#   genus + two letters taken from the rest of the
#   name, trying combinations until one is free.
#   Only letters are used (dashes, digits, dots and
#   spaces are ignored). If the letters run out,
#   a number is appended to the original id.
#
#   Arguments:
#       ids_table: CSV table (species,id_m,id_h,final_id)
#       output_table: CSV with the resolved identifiers
#       id_col: column with the identifier (default: 4)
#       changes_table: (optional) species,old_id,new_id
#
#####################################################

resolve_duplicated_species_ids() {

    local ids_table="${1}"
    local output_table="${2}"
    local id_col="${3:-4}"
    local changes_table="${4:-/dev/null}"

    : > "$changes_table"

    awk -F',' -v OFS=',' -v col="$id_col" -v changes="$changes_table" '
        function new_id(name, old_id,    tmp, pool, g, n, i, j, k, cand) {
            # Keep only letters; words separated by spaces
            tmp = name
            gsub(/[^A-Za-z]+/, " ", tmp)
            sub(/^ +/, "", tmp)
            # First letter of the genus
            g = tolower(substr(tmp, 1, 1))
            # Pool of letters: everything after the genus
            sub(/^[^ ]* */, "", tmp)
            gsub(/ /, "", tmp)
            pool = tolower(tmp)
            n = length(pool)

            # Genus letter + 2 letters of the pool, until one is free
            for (i = 1; i <= n; i++) {
                for (j = (i == 1 ? 2 : 1); j <= n; j++) {
                    cand = g substr(pool, i, 1) substr(pool, j, 1)
                    if (!(cand in used)) return cand
                }
            }
            # No letters left: try numbers
            k = 2
            while ((old_id k) in used) k++
            return old_id k
        }
        # Pass 1: every identifier already in use
        NR == FNR { used[$col] = 1; next }
        # Pass 2: first one keeps it, the rest get a new one
        {
            id = $col
            if (!(id in seen)) { seen[id] = 1; print; next }
            nid = new_id($1, id)
            used[nid] = 1
            print $1, id, nid >> changes
            $col = nid
            print
        }
    ' "$ids_table" "$ids_table" > "$output_table"
}


#####################################################
#
#   Creates copies of the PmiREN FASTA files with all
#   species identifiers in the headers replaced by
#   the final one assigned in the identifiers table.
#   It also removes FASTA records with empty or invalid
#   sequences, preserving valid IUPAC nucleotide codes.
#   The original files are not modified.
#
#   Arguments:
#       ids_table: CSV table (species,id_m,id_h,final_id)
#       files_list: list of FASTA files (space-separated)
#       type: mature or hairpin
#       output_dir: directory for the new files
#
#####################################################

update_pmiren_fasta_ids() {

    local ids_table="${1}"
    local files_list="${2}"
    local type="${3}"
    local output_dir="${4}"

    # Column of the table with the original identifier of this type of file
    local old_id_col
    case "$type" in
        mature)  old_id_col=2 ;;
        hairpin) old_id_col=3 ;;
        *) echo "ERROR: type must be 'mature' or 'hairpin'" >&2; return 1 ;;
    esac

    mkdir -p "$output_dir"
    local n_errors=0

    for file in $files_list; do

        # Species name from the file name (same as get_PmiREN_species_ids)
        local species_name
        species_name=$(basename "$file" | sed "s/_/ /g; s/${type}.fa//" | sed 's/[[:space:]]*$//')

        # Original and new identifiers of the species
        local ids old_id new_id
        ids=$(awk -F',' -v sp="$species_name" -v c="$old_id_col" '$1 == sp { print $c "," $4; exit }' "$ids_table")
        if [ -z "$ids" ]; then
            echo "ERROR: '$species_name' ($file) not found in $ids_table" >&2
            n_errors=$((n_errors + 1)); continue
        fi
        old_id="${ids%%,*}"
        new_id="${ids##*,}"

        # Replace identifiers in ALL headers and remove invalid FASTA records
        local output_file="$output_dir/$(basename "$file")"
        awk -v old="$old_id" -v new="$new_id" -v f="$file" '
            function process_record(    p, seq) {
                if (header == "") return

                seq = sequence
                gsub(/[[:space:]]/, "", seq)

                if (seq == "" || seq ~ /[^ACGTURYSWKMBDHVNacgturyswkmbdhvn]/) {
                    print "WARNING: Invalid sequence removed from " f ": " header > "/dev/stderr"
                    return
                }

                p = index(header, "-")
                if (p > 2) {
                    print ">" new substr(header, p)
                } else {
                    print header
                    other++
                }
                printf "%s", sequence
            }

            /^>/ {
                process_record()
                header = $0
                sequence = ""
                next
            }
            { sequence = sequence $0 "\n" }

            END {
                process_record()
                if (other > 0)
                    print "WARNING: " other " headers of " f " do not contain a valid identifier" > "/dev/stderr"
            }
        ' "$file" > "$output_file"

        # Sanity check: output cannot contain more lines than input
        if [ "$(wc -l < "$file")" -lt "$(wc -l < "$output_file")" ]; then
            echo "ERROR: line count differs in $output_file" >&2
            n_errors=$((n_errors + 1))
        fi
    done

    [ "$n_errors" -eq 0 ]
}

### MAIN
main () {
    
    ## 1. PmiREN
    ############################################################################

    # Get arguments
    arguments_management "$@"
    
    # Create a temporary and outpur dir
    workdir=$(pwd)
    mkdir -p "$workdir/tmp"
    mkdir -p tmp
    mkdir -p 01-Mod_databases

    # Species id final file name
    final_ids_file_out=species_ids_db.csv

    # Execute only if pmiren is provided
    if [ -n "$pmiren" ]; then

        # Get PmiREN species ids
        echo "Obtaining identifiers from PmiREN database..."
        get_PmiREN_species_ids "$pmiren" $workdir/tmp/pmiren_species_mature_ids.txt mature
        get_PmiREN_species_ids "$pmiren_hairpin" $workdir/tmp/pmiren_species_hairpin_ids.txt hairpin
        echo "Done!"

        # Merge mature and hairpin
        echo "Comparing identifiers between precursor and mature miRNA files..."
        merge_and_compare_tables $workdir/tmp/pmiren_species_mature_ids.txt $workdir/tmp/pmiren_species_hairpin_ids.txt 2 2 '1.2,2.2,1.1,2.1' $workdir/tmp/01-pmiren_species_ids.csv $workdir/tmp/01-pmiren_species_ids_diff.csv $workdir/tmp
        echo "Done!"

        # Assign new identifiers to species that share the same one
        echo "Resolving duplicated identifiers in PmiREN..."
        resolve_duplicated_species_ids \
            $workdir/tmp/01-pmiren_species_ids.csv \
            $workdir/tmp/01-pmiren_species_ids_unique.csv \
            4 \
            $workdir/tmp/01-pmiren_species_ids_changes.csv
        echo "Done!"

        # Generate PmiREN FASTA files with the new identifiers
        echo "Updating identifiers in PmiREN FASTA files..."
        rm -rf $workdir/pmiren_mod
        update_pmiren_fasta_ids $workdir/tmp/01-pmiren_species_ids_unique.csv "$pmiren" mature $workdir/pmiren_mod || exit 1
        update_pmiren_fasta_ids $workdir/tmp/01-pmiren_species_ids_unique.csv "$pmiren_hairpin" hairpin $workdir/pmiren_mod || exit 1

        # From here on, use the new files
        pmiren=$(find $workdir/pmiren_mod -type f -name "*_mature.fa" | sort | tr '\n' ' ')
        pmiren_hairpin=$(find $workdir/pmiren_mod -type f -name "*_hairpin.fa" | sort | tr '\n' ' ')
        echo "Done!"

        # The new files already use the final identifiers: build a table that
        # matches them (columns 2 and 3 are replaced by column 4)
        awk 'BEGIN {FS=OFS=","} {$2=$4; $3=$4; print}' \
            $workdir/tmp/01-pmiren_species_ids_unique.csv \
            > $workdir/tmp/01-pmiren_species_ids_final.csv

        # Save the path of the file to merge in a variable
        file_to_merge=$workdir/tmp/01-pmiren_species_ids_final.csv
        file_to_merge_diff=$workdir/tmp/01-pmiren_species_ids_diff.csv

        # If only PmiREN has been provided
        if [ -z "$mirbase" ] && [ -z "$srnaanno" ]; then

            # Create the output directory
            mkdir -p 01-Mod_databases/PmiREN

            # Create database file using species files
            cat $pmiren > 01-Mod_databases/PmiREN/pmiren_mature.fa
            cat $pmiren_hairpin > 01-Mod_databases/PmiREN/pmiren_hairpin.fa

            # Modify miR names (MIR -> miR)
            sed -i 's/[Mm][Ii][Rr]/miR/g' 01-Mod_databases/PmiREN/pmiren_hairpin.fa
                
            # Add a NULL column in the positions corresponding to the sRNAanno
            # and miRBase databases
            awk -F, '{OFS=","; print $1,"NULL",$3,"NULL",$3}' $file_to_merge > $final_ids_file_out

            # Sort and add header to IDs file
            sort -T "$workdir/tmp" -t',' -k1,1 $final_ids_file_out > tmp.csv && mv tmp.csv $final_ids_file_out
            sed -i '1i\species,mirbase,pmiren,srnaanno,final_id' $final_ids_file_out
        fi
    fi

    ## 2. miRBase
    ############################################################################

    # Execute only if miRBase is provided
    if [ -n "$mirbase" ]; then

        # Get miRBase species ids
        echo "Obtaining identifiers from miRBase database..."
        get_mirbase_species_ids $mirbase $workdir/tmp/mirbase_mature_species_ids.txt "$workdir/tmp"
        get_mirbase_species_ids $mirbase_hairpin $workdir/tmp/mirbase_hairpin_species_ids.txt "$workdir/tmp"
        echo "Done!"

        # Merge mature and hairpin
        echo "Comparing identifiers between precursor and mature miRNA files..."
        merge_and_compare_tables $workdir/tmp/mirbase_mature_species_ids.txt $workdir/tmp/mirbase_hairpin_species_ids.txt 2 2 '1.2,2.2,1.1,2.1' $workdir/tmp/02-mirbase_species_ids.csv $workdir/tmp/02-mirbase_species_ids_diff.csv $workdir/tmp
        echo "Done!"
        
        # Save the path of the file to merge in a variable
        file_to_merge=$workdir/tmp/02-mirbase_species_ids.csv
        file_to_merge_diff=$workdir/tmp/02-mirbase_species_ids_diff.csv

        # If only miRBase has been provided
        if [ -z "$pmiren" ] && [ -z "$srnaanno" ]; then
            
            # Create the output directory
            mkdir -p 01-Mod_databases/miRBase

            # Create the output mirbase files
            cp $mirbase 01-Mod_databases/miRBase/mirbase_mature.fa
            cp $mirbase_hairpin 01-Mod_databases/miRBase/mirbase_hairpin.fa

            # Modify miR names (MIR -> miR)
            sed -i 's/[Mm][Ii][Rr]/miR/g' 01-Mod_databases/miRBase/mirbase_hairpin.fa
            
            # Add a NULL column in the positions corresponding to the sRNAanno
            # and PmiREN databases
            awk -F, '{OFS=","; print $1,$3,"NULL","NULL",$3}' $file_to_merge > $final_ids_file_out

            # Sort and add header to IDs file
            sort -t',' -k1,1 $final_ids_file_out > tmp.csv && mv tmp.csv $final_ids_file_out
            sed -i '1i\species,mirbase,pmiren,srnaanno,final_id' $final_ids_file_out
        fi

    fi

    ## 3. Merge the tables of mature and precursor miRNAs from each database 
    ############################################################################

    # Execute only if both miRBase and PmiREN have been provided
    if [ -n "$mirbase" ] && [ -n "$pmiren" ]; then

        # Merge miRBase and PmiREN tables
        echo "Comparing identifiers between miRBase and PmiREN databases..."
        merge_and_compare_tables $workdir/tmp/02-mirbase_species_ids.csv $workdir/tmp/01-pmiren_species_ids_final.csv 1 1 '1.1,2.1,1.4,2.4' $workdir/tmp/03-mirbase_pmiren_species_ids.csv $workdir/tmp/03-mirbase_pmiren_species_ids_diff.csv $workdir/tmp
        echo "Done!"
        
        # Save the path of the file to merge in a variable
        file_to_merge=$workdir/tmp/03-mirbase_pmiren_species_ids.csv
        file_to_merge_diff=$workdir/tmp/03-mirbase_pmiren_species_ids_diff.csv

        ## 4. Check if there are species with the same identifier
        ############################################################################

        # Check if there is any duplicate identifier."
        duplicate_ids=$(cut -d ',' -f 4 "$file_to_merge" | sort -T "$workdir/tmp" | uniq -d)
        if [ -n "$duplicate_ids" ]; then

            # Print identifiers
            echo "NOTE: Identifiers have been found that are associated with more than one species"
            echo "IDs:"

            # Iterate through duplicate identifiers
            for id in $duplicate_ids; do

                printf "\t--$id--\n"

                # Get the list of species names (e.g. Arabidopsis_thaliana Cucumis_melo)
                id_line_list=$(awk -F ',' -v OFS=',' -v id="$id" '{gsub(/ /, "_", $1)} $4 == id {print}' "$file_to_merge")

                # Iterate through the lines
                counter=1
                for line in $id_line_list; do  
                    # Get the species name
                    sp_name=$(echo $line | awk -F, '{print $1}' | sed 's/_/ /g')
                    
                    printf "\t\t$sp_name\n"

                    # Run from the second iteration onwards (to keep one species with the original identifier).
                    if [ "$counter" -ne 1 ]; then
                        # Column where the identifiers are located.
                        col_id=4

                        # Create a new id for the species
                        new_id=$(change_the_id_using_name "$sp_name" "$file_to_merge" $col_id)

                        # Build the new line for 'diff' table
                        new_line=$(echo "$line" | awk -F',' -v OFS=',' -v sp_name="$sp_name" -v new_id="$new_id" '{$1=sp_name; $4=new_id; print}')

                        # Add the new identifier to the identifiers table.
                        awk -F ',' -v sp_name="$sp_name" -v new_id="$new_id" 'BEGIN {OFS=","} $1 == sp_name { $4 = new_id } 1' "$file_to_merge" > tmpfile && mv tmpfile "$file_to_merge"
                        
                        # Add the new identifier to the identifiers table (Diff)
                        awk -F ',' -v sp_name="$sp_name" -v new_id="$new_id" 'BEGIN {OFS=","} $1 == sp_name { $4 = new_id } 1' "$file_to_merge_diff" > tmpfile && mv tmpfile "$file_to_merge_diff"
                        
                        # If the new identifier is not in the 'diff' table, add its respective line at the end of the table.
                        awk -F ',' -v sp_name="$sp_name" -v line="$new_line" 'BEGIN {OFS=","} $1 == sp_name {found=1} END {if (!found) print line}' "$file_to_merge_diff" > temp && cat temp >> "$file_to_merge_diff" && rm temp

                    fi

                    # Increment the variable
                    let "counter++"
                done
            done
        fi

        printf '\nNew identifiers assigned to avoid duplications!\n'

        ## 5. Replace identifiers that do not match between the two databases
        ############################################################################

        # Create output directory for filtered database
        output_mirbase=01-Mod_databases/miRBase
        mkdir -p $output_mirbase

        # Create output directory for the PmiREN database
        output_pmiren=01-Mod_databases/PmiREN
        mkdir -p $output_pmiren

        # Create database file using species files
        cat $pmiren > $workdir/tmp/temporal_mature.fa
        cat $pmiren_hairpin > $workdir/tmp/temporal_hairpin.fa

        # Replace identifiers in PmiREN so that both databases match (PmiREN)
        replace_mismatched_ids $workdir/tmp/temporal_mature.fa $output_pmiren/pmiren_mature.fa 3 4 "$file_to_merge_diff" "$workdir/tmp"
        replace_mismatched_ids $workdir/tmp/temporal_hairpin.fa $output_pmiren/pmiren_hairpin.fa 3 4 "$file_to_merge_diff" "$workdir/tmp"

        # Modify miR names (MIR -> miR)
        sed -i 's/[Mm][Ii][Rr]/miR/g' $output_pmiren/pmiren_hairpin.fa

        # Replace identifiers in PmiREN so that both databases match (miRBase)
        replace_mismatched_ids $mirbase $output_mirbase/mirbase_mature.fa 2 4 "$file_to_merge_diff" "$workdir/tmp"
        replace_mismatched_ids $mirbase_hairpin $output_mirbase/mirbase_hairpin.fa 2 4 "$file_to_merge_diff" "$workdir/tmp"

        # Modify miR names (MIR -> miR)
        sed -i 's/[Mm][Ii][Rr]/miR/g' $output_mirbase/mirbase_hairpin.fa

        # Delete temporary directory
        rm -rf $workdir/tmp/temporal_*
            
        # Execute only if srnaanno has not been provided
        if [ -z "$srnaanno" ]; then
            
            # Add a NULL column in the position corresponding to the sRNAanno
            # database.
            awk -F, '{OFS=","; print $1,$2,$3,"NULL",$4}' $file_to_merge > $final_ids_file_out

            # Sort and add header to IDs file
            sort -T "$workdir/tmp" -t',' -k1,1 $final_ids_file_out > tmp.csv && mv tmp.csv $final_ids_file_out
            sed -i '1i\species,mirbase,pmiren,srnaanno,final_id' $final_ids_file_out
        fi
    fi

    # Execute only if srnaanno has been provided
    if [ -n "$srnaanno" ]; then

        # 6. Prepare sRNAanno database 
        ###########################################################################

        # sRNAanno output directory path
        output_srnaanno=01-Mod_databases/sRNAanno

        # Create the output directory
        mkdir -p $output_srnaanno

        # Create an empty file_to_merge if only sRNAanno has been provided
        [ -z "$mirbase" ] && [ -z "$pmiren" ] && file_to_merge="$workdir/tmp/sRNAanno_ids.csv" && touch "$file_to_merge"

        # Check if the file species
        echo "Assigning identifiers to the sRNAanno database..."
        prepare_sRNAanno_database $file_to_merge "$srnaanno" $output_srnaanno $final_ids_file_out
        echo "Done!"

        # miRBase has not been provided, but PmiREN has.
        if [ -z "$mirbase" ] && [ -n "$pmiren" ]; then

            # Create the output directory
            mkdir -p 01-Mod_databases/PmiREN

            # Create database file using species files
            cat $pmiren > 01-Mod_databases/PmiREN/pmiren_mature.fa
            cat $pmiren_hairpin > 01-Mod_databases/PmiREN/pmiren_hairpin.fa

            # Convert the first letter of the identifier to lowercase (in
            # PmiREN, it is uppercase)
            sed -i -E 's/^>([A-Z])/>\L\1/' 01-Mod_databases/PmiREN/pmiren_mature.fa
            sed -i -E 's/^>([A-Z])/>\L\1/' 01-Mod_databases/PmiREN/pmiren_hairpin.fa

            # Modify miR names (MIR -> miR)
            sed -i 's/[Mm][Ii][Rr]/miR/g' 01-Mod_databases/PmiREN/pmiren_hairpin.fa

            # If mirbase has not been provided, remove the third column from the
            # table (which in this case is associated with the mirbase hairpin)
            # and add a column of NAs in the second position.
            awk -F, '{ $3=""; OFS=","; print $1",NULL,"$2","$4","$5 }' $final_ids_file_out > tmp.csv && mv tmp.csv $final_ids_file_out
            
        fi

        # PmiREN has not been provided, but miRBase has.
        if [ -n "$mirbase" ] && [ -z "$pmiren" ]; then

            # Create the output directory
            mkdir -p 01-Mod_databases/miRBase

            # Create the output mirbase files
            cp $mirbase 01-Mod_databases/miRBase/mirbase_mature.fa
            cp $mirbase_hairpin 01-Mod_databases/miRBase/mirbase_hairpin.fa

            # Modify miR names (MIR -> miR)
            sed -i 's/[Mm][Ii][Rr]/miR/g' 01-Mod_databases/miRBase/mirbase_hairpin.fa

            # If pmiren has not been provided, replace the third column of the
            # ID file with NAs (The third column corresponds to this database).
            awk -F, '{ $3=""; OFS=","; print $1","$2",NULL,"$4","$5 }' $final_ids_file_out > tmp.csv && mv tmp.csv $final_ids_file_out

        fi

        # Sort and add header to IDs file
        sort -T "$workdir/tmp" -t',' -k1,1 $final_ids_file_out > tmp.csv && mv tmp.csv $final_ids_file_out
        sed -i '1i\species,mirbase,pmiren,srnaanno,final_id' $final_ids_file_out
    fi

    # 7. Check input species names
    ###########################################################################
    
    # Create an array using the input: ['Glycine max', 'Homo sapiens']
    input_species_pre=$(echo $input_species | tr -d '[]')
    IFS=',' read -ra input_species_array <<< "$input_species_pre"
    for species_name in "${input_species_array[@]}"; do
        # Remove spaces at the beginning and end of the name
        stripped_sp_name=$(echo "$species_name" | sed 's/^ *//; s/ *$//')
        
        # If the species is not in the databases...
        if ! awk -F',' -v name="$stripped_sp_name" '$1 == name' "$final_ids_file_out" | grep -q .; then
            # Get the two words from the species name
            genus=$(echo "$stripped_sp_name" | awk '{print $1}')
            specific_epithet=$(echo "$stripped_sp_name" | awk '{print $2}')

            # Generate an indentifier
            col_ids=5
            new_id=$(echo "${genus:0:1}${specific_epithet:0:2}" | tr '[:upper:]' '[:lower:]')
            species_id=$(change_the_id_using_name "$stripped_sp_name" $final_ids_file_out $col_ids)

            # Add the species
            echo "$stripped_sp_name,NULL,NULL,NULL,$species_id" >> "$final_ids_file_out"
        fi
    done

    #rm -r tmp
}
main "$@"