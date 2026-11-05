process DOWNLOAD_DB {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e4d509f26d12542a7ad16fdb5bc7c700eea8d64be05c54775e7b4c6a693081a/data' :
        'community.wave.seqera.io/library/seqkit_wget:001d8477461a0364' }"
        
    input:
    val species
    val mirbase
    val srnaanno
    val pmiren
    path species_synonyms

    output:
    path "species_ids.csv"                                                                                  , emit: species_ids
    tuple path("01-Mod_databases/miRBase/*_mature.fa")  , path("01-Mod_databases/miRBase/*_hairpin.fa")     , emit: mirbase
    tuple path("01-Mod_databases/sRNAanno/*_mature.fa") , path("01-Mod_databases/sRNAanno/*_hairpin.fa")    , emit: srnaanno
    tuple path("01-Mod_databases/PmiREN/*_mature.fa")   , path("01-Mod_databases/PmiREN/*_hairpin.fa")      , emit: pmiren

    script:
    """
    # 0. FUNCTIONS
    #####################################################

    #####################################################
    #
    #   This function corrects the species names
    #   in miRBase FASTA files using a reference
    #   table that associates species identifiers
    #   with their corresponding scientific names.
    #   It updates the species names in FASTA headers
    #   while preserving sequence identifiers,
    #   accession numbers, annotations, and nucleotide
    #   sequences.
    #
    #   Arguments:
    #       Path to the TSV file containing
    #           species identifiers and scientific names.
    #       Path to the input miRBase FASTA file
    #           containing mature or hairpin sequences.
    #       Path to the output FASTA file where
    #           corrected sequences will be saved.
    #
    #####################################################

    correct_mirbase_species_names() {
        local species_file="\$1"
        local input_file="\$2"
        local output_file="\$3"

        # Correct the species names within the miRBase file
        awk -v species_file="\$species_file" '
        BEGIN {
            # Save the species data into an array
            while (getline < species_file) {
                split(\$0, arr, "\t");
                id2species[arr[1]] = arr[2];
            }
            close(species_file);
        }
        # Only for the headers...
        \$0 ~ /^>/ {
            # Get the species id
            id = substr(\$0, 2, index(\$0, "-") - 2);

            # Check if the species id is in the array
            if (id in id2species) {
                # Get the rest of the elements of the header
                split(\$0, arr, " ");
                first_element = arr[1];
                second_element = arr[2];

                if (arr[length(arr)] == "stem-loop") {
                    last_element = arr[length(arr)-1] " " arr[length(arr)];
                } else {
                    last_element = arr[length(arr)];
                }

                # Get the new species name
                species = id2species[id];
                print first_element" "second_element" "species" "last_element;
            } else {
                print \$0;
            }
            next;
        }
        # Print the sequences
        {
            print \$0;
        }' < "\$input_file" > "\$output_file"
    }


    #####################################################
    #
    #   This function normalizes species names across
    #   miRBase, PmiREN, and sRNAanno databases using
    #   a reference table of canonical scientific names
    #   and their synonyms. It replaces synonymous
    #   species names with their canonical equivalents
    #   while preserving the original sequence data
    #   and species identifiers.
    #
    #   It generates normalized copies of the input
    #   files without modifying the original databases.
    #   For miRBase, it updates species names in FASTA
    #   headers. For PmiREN and sRNAanno, it normalizes
    #   species names in filenames.
    #
    #   Additionally, it normalizes the requested
    #   species names, detects conflicting synonyms
    #   and filename collisions, and generates a
    #   manifest containing the normalized input
    #   paths and species names required by the
    #   database identifier resolver.
    #
    #   Arguments:
    #       Path to the TSV file containing
    #           canonical species names and synonyms.
    #       List of requested species names.
    #       Boolean indicating whether miRBase
    #           database processing is enabled.
    #       Boolean indicating whether PmiREN
    #           database processing is enabled.
    #       Boolean indicating whether sRNAanno
    #           database processing is enabled.
    #       Directory containing the original
    #           miRBase files and species table.
    #       Directory containing the combined
    #           PmiREN mature and star FASTA files.
    #       Directory containing the PmiREN
    #           hairpin FASTA files.
    #       Directory containing the original
    #           sRNAanno GFF3 files.
    #       Directory where the normalized database
    #           files and synonym mapping will be saved.
    #       Path to the output manifest file containing
    #           normalized database arguments and
    #           canonical requested species names.
    #
    #####################################################

    normalize_species_synonyms() {

        local synonyms_file="\$1"
        local requested_species="\$2"
        local use_mirbase="\$3"
        local use_pmiren="\$4"
        local use_srnaanno="\$5"
        local mirbase_dir="\$6"
        local pmiren_mature_dir="\$7"
        local pmiren_hairpin_dir="\$8"
        local srnaanno_dir="\$9"
        local output_dir="\${10}"
        local manifest_file="\${11}"

        local canonical synonyms synonym species_name file output
        local normalized_name species_list
        local mirbase_args='' pmiren_m_args='' pmiren_p_args='' srnaanno_args=''
        local canonical_species_input='[' separator=''

        local -A synonym_map=()
        local -a synonym_list requested_array

        if [ ! -f "\$synonyms_file" ]; then
            echo "ERROR: Species synonyms file not found: \$synonyms_file" >&2
            return 1
        fi

        # Read synonym table
        while IFS=\$'\\t' read -r canonical synonyms || [ -n "\$canonical" ]; do

            canonical=\$(printf '%s' "\$canonical" |
                sed 's/^[[:space:]]*//; s/[[:space:]]*\$//')

            [ "\$canonical" = "canonical_name" ] && continue
            [ -z "\$canonical" ] && continue

            IFS=';' read -ra synonym_list <<< "\$synonyms"

            for synonym in "\${synonym_list[@]}"; do

                synonym=\$(printf '%s' "\$synonym" |
                    sed 's/^[[:space:]]*//; s/[[:space:]]*\$//')

                [ -z "\$synonym" ] && continue

                if [ -n "\${synonym_map[\$synonym]:-}" ] &&
                   [ "\${synonym_map[\$synonym]}" != "\$canonical" ]; then
                    echo "ERROR: Conflicting synonym: \$synonym" >&2
                    return 1
                fi

                synonym_map["\$synonym"]="\$canonical"
            done

        done < "\$synonyms_file"

        mkdir -p "\$output_dir"

        # Export the synonym mapping for AWK
        local mapping_file="\$output_dir/synonyms_map.tsv"
        : > "\$mapping_file"

        for synonym in "\${!synonym_map[@]}"; do
            printf '%s\\t%s\\n' "\$synonym" "\${synonym_map[\$synonym]}" \
                >> "\$mapping_file"
        done

        # Normalize miRBase
        if [ "\$use_mirbase" = true ]; then

            mkdir -p "\$output_dir/miRBase"

            awk -F '\\t' '
                FILENAME == ARGV[1] {
                    mapping[\$1] = \$2
                    next
                }
                {
                    name = (\$2 in mapping) ? mapping[\$2] : \$2
                    print \$1 "\\t" name
                }
            ' "\$mapping_file" "\$mirbase_dir/mirbase_species.txt" \
                > "\$output_dir/miRBase/mirbase_species.txt" || return 1

            # Do not silently combine different source identifiers
            awk -F '\\t' '
                {
                    if (\$2 in ids && ids[\$2] != \$1) {
                        print "ERROR: miRBase synonym collision: " \$2 > "/dev/stderr"
                        error = 1
                    }
                    ids[\$2] = \$1
                }
                END {exit error}
            ' "\$output_dir/miRBase/mirbase_species.txt" || return 1

            correct_mirbase_species_names \
                "\$output_dir/miRBase/mirbase_species.txt" \
                "\$mirbase_dir/mature.fa" \
                "\$output_dir/miRBase/mature_mod.fa" || return 1

            correct_mirbase_species_names \
                "\$output_dir/miRBase/mirbase_species.txt" \
                "\$mirbase_dir/hairpin.fa" \
                "\$output_dir/miRBase/hairpin_mod.fa" || return 1

            mirbase_args="--mirbase \$output_dir/miRBase/mature_mod.fa --mhairpin \$output_dir/miRBase/hairpin_mod.fa"
        fi

        # Normalize PmiREN
        if [ "\$use_pmiren" = true ]; then

            mkdir -p "\$output_dir/PmiREN"

            while IFS= read -r file; do

                species_name=\$(basename "\$file" "_mature.fa" | tr '_' ' ')
                canonical="\${synonym_map[\$species_name]:-\$species_name}"
                normalized_name="\${canonical// /_}_mature.fa"
                output="\$output_dir/PmiREN/\$normalized_name"

                if [ -e "\$output" ]; then
                    echo "ERROR: PmiREN mature collision: \$canonical" >&2
                    return 1
                fi

                cp "\$file" "\$output" || return 1
                pmiren_m_args+="--pmiren \$output "

            done < <(find "\$pmiren_mature_dir" -type f -name '*_mature.fa' | sort)

            while IFS= read -r file; do

                species_name=\$(basename "\$file" "_hairpin.fa" | tr '_' ' ')
                canonical="\${synonym_map[\$species_name]:-\$species_name}"
                normalized_name="\${canonical// /_}_hairpin.fa"
                output="\$output_dir/PmiREN/\$normalized_name"

                if [ -e "\$output" ]; then
                    echo "ERROR: PmiREN hairpin collision: \$canonical" >&2
                    return 1
                fi

                cp "\$file" "\$output" || return 1
                pmiren_p_args+="--phairpin \$output "

            done < <(find "\$pmiren_hairpin_dir" -type f -name '*_hairpin.fa' | sort)
        fi

        # Normalize sRNAanno
        if [ "\$use_srnaanno" = true ]; then

            mkdir -p "\$output_dir/sRNAanno"

            for file in "\$srnaanno_dir"/*.miRNA.gff3; do

                [ -f "\$file" ] || continue

                species_name=\$(basename "\$file" ".miRNA.gff3" | tr '_' ' ')
                canonical="\${synonym_map[\$species_name]:-\$species_name}"
                normalized_name="\${canonical// /_}.miRNA.gff3"
                output="\$output_dir/sRNAanno/\$normalized_name"

                if [ -e "\$output" ]; then
                    echo "ERROR: sRNAanno collision: \$canonical" >&2
                    return 1
                fi

                cp "\$file" "\$output" || return 1
                srnaanno_args+="--srnaanno \$output "

            done
        fi

        # Normalize requested species
        species_list="\${requested_species//[\\[\\]]/}"
        IFS=',' read -ra requested_array <<< "\$species_list"

        for species_name in "\${requested_array[@]}"; do

            species_name=\$(printf '%s' "\$species_name" |
                sed 's/^[[:space:]]*//; s/[[:space:]]*\$//')

            canonical="\${synonym_map[\$species_name]:-\$species_name}"

            canonical_species_input+="\$separator\$canonical"
            separator=','

        done

        canonical_species_input+=']'

        # Write results without modifying caller variables
        printf '%s\\n' \
            "\$mirbase_args" \
            "\$srnaanno_args" \
            "\$pmiren_m_args" \
            "\$pmiren_p_args" \
            "\$canonical_species_input" \
            > "\$manifest_file"
    }

    mirbase_files=''
    srnaanno_files=''
    pmiren_m_files=''
    pmiren_p_files=''

    ## 1. DOWNLOAD DATABASES
    ############################################################################

    # Download miRBase
    if [ ${mirbase} = true ]; then

        # Download miRBase files
        wget --no-check-certificate -P mirbase https://mirbase.org/download/mature.fa
        wget --no-check-certificate -P mirbase https://mirbase.org/download/hairpin.fa

                # Generate species table from miRBase FASTA files
        awk '
            /^>/ {
                # Extract species identifier
                split(\$1, id_parts, "-")
                id = substr(id_parts[1], 2)

                # Identify the last organism-name field
                if (\$2 ~ /^MIMAT[0-9]+\$/) {
                    last = NF - 1
                }
                else if (\$2 ~ /^MI[0-9]+\$/ && \$NF == "stem-loop") {
                    last = NF - 2
                }
                else {
                    next
                }

                # Reconstruct the complete organism name
                species = ""
                for (i = 3; i <= last; i++) {
                    species = species (species == "" ? "" : " ") \$i
                }

                # Store identifier-organism pairs
                if (id != "" && species != "" && !(id in names)) {
                    names[id] = species
                }
            }
            END {
                for (id in names)
                    print id "\t" names[id]
            }
        ' mirbase/hairpin.fa mirbase/mature.fa | sort > mirbase/mirbase_species.txt

        # Check that the species table is not empty
        if [ ! -s mirbase/mirbase_species.txt ]; then
            echo "ERROR: miRBase species table is empty" >&2
            exit 1
        fi

        # Correct species names within the miRBase files
        correct_mirbase_species_names mirbase/mirbase_species.txt mirbase/mature.fa mirbase/mature_mod.fa
        correct_mirbase_species_names mirbase/mirbase_species.txt mirbase/hairpin.fa mirbase/hairpin_mod.fa

        # Create a string for the execution of 01-Database_id_resolver.sh
        mirbase_files='--mirbase mirbase/mature_mod.fa --mhairpin mirbase/hairpin_mod.fa'
    fi

    # Download sRNAanno
    if [ ${srnaanno} = true ]; then
        wget --no-check-certificate http://121.37.229.61:84/sRNAannoDATA/miRNAGFF3/miRNA.gff3.tar.gz
        mkdir -p srnaanno && tar -xzvf miRNA.gff3.tar.gz -C srnaanno
        rm miRNA.gff3.tar.gz

        # Create an string for the execution of 01-Database_id_resolver.sh
        srnaanno_files=''
        for file in srnaanno/*; do
            srnaanno_files+="--srnaanno \$file "
        done
    fi
    
    # Download PmiREN
    if [ ${pmiren} = true ]; then
        pmiren_version="v3.2"

        wget --no-check-certificate -r --no-parent -nH --cut-dirs=1 -P pmiren \\
            --accept "*.fa" --reject "index.html*" --reject-regex '[?]' \\
            https://www.pmiren.com/ftp-download/\$pmiren_version/

        # Combine mature and star sequences into a new mature FASTA
        mkdir -p pmiren_combined

        while IFS= read -r mature_file; do
            star_file="\${mature_file%_mature.fa}_star.fa"
            combined_file="pmiren_combined/\$(basename "\$mature_file")"

            if [ -f "\$star_file" ]; then
                cat "\$mature_file" "\$star_file" > "\$combined_file"
            else
                cp "\$mature_file" "\$combined_file"
            fi

            pmiren_m_files+="--pmiren \$combined_file "
        done < <(find pmiren/\$pmiren_version -type f -name "*_mature.fa" | sort)

        while IFS= read -r hairpin_file; do
            # Rename stemloop files to hairpin
            new_hairpin_file="\${hairpin_file%_stemloop.fa}_hairpin.fa"
            mv "\$hairpin_file" "\$new_hairpin_file"
            pmiren_p_files+="--phairpin \$new_hairpin_file "
        done < <(find pmiren/\$pmiren_version -type f -name "*_stemloop.fa" | sort)
    fi

    ## 1. NORMALIZE SPECIES SYNONYMS
    ############################################################################

    normalize_species_synonyms \
        "${species_synonyms}" \
        "${species}" \
        "${mirbase}" \
        "${pmiren}" \
        "${srnaanno}" \
        "mirbase" \
        "pmiren_combined" \
        "pmiren/\$pmiren_version" \
        "srnaanno" \
        "synonyms_normalized" \
        "synonyms_manifest.txt" || exit 1

    # Read normalization results
    mapfile -t normalized_results < synonyms_manifest.txt

    mirbase_files="\${normalized_results[0]}"
    srnaanno_files="\${normalized_results[1]}"
    pmiren_m_files="\${normalized_results[2]}"
    pmiren_p_files="\${normalized_results[3]}"
    canonical_species_input="\${normalized_results[4]}"

    ## 2. RESOLVE IDENTIFIERS
    ############################################################################

    # Create the output directories
    mkdir -p 01-Mod_databases/miRBase/
    mkdir -p 01-Mod_databases/sRNAanno/
    mkdir -p 01-Mod_databases/PmiREN/

    # Execute the programm
    bash 01-Database_id_resolver.sh \\
        --inputsps "\$canonical_species_input" \\
        \$mirbase_files \\
        \$srnaanno_files \\
        \$pmiren_m_files \\
        \$pmiren_p_files
    
    # Put the entire sequences on the same line
    if [ ${mirbase} = true ]; then
        awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq \$0} END {if (seq) print seq}' 01-Mod_databases/miRBase/mirbase_hairpin.fa > tmp/tmp_s && mv tmp/tmp_s 01-Mod_databases/miRBase/mirbase_hairpin.fa
    else
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_mature.fa
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_hairpin.fa
    fi
    if [ ${srnaanno} = true ]; then
        awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq \$0} END {if (seq) print seq}' 01-Mod_databases/sRNAanno/srnaanno_hairpin.fa > tmp/tmp_s && mv tmp/tmp_s 01-Mod_databases/sRNAanno/srnaanno_hairpin.fa
    else
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_mature.fa
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_hairpin.fa
    fi
    if [ ${pmiren} = true ]; then
        awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq \$0} END {if (seq) print seq}' 01-Mod_databases/PmiREN/pmiren_hairpin.fa > tmp/tmp_s && mv tmp/tmp_s 01-Mod_databases/PmiREN/pmiren_hairpin.fa
    else
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_mature.fa
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_hairpin.fa
    fi

    # Create a new species_ids file
    head -n 1 species_ids_db.csv > species_ids.csv

    # If it is requested that the database be available for each species individually.
    input_species_pre=\$(echo "\$canonical_species_input" | tr -d '[]' | sed 's/, */,/g')
    IFS=',' read -ra input_species_array <<< "\$input_species_pre"
    for species_name in "\${input_species_array[@]}"; do
        # Select the required ids from species_id_db.csv
        awk -F',' -v species="\$species_name" '\$1 == species' species_ids_db.csv >> species_ids.csv
    done
    """

    stub:
    """
    # Create directories for mirbase, srnaanno, and pmiren
    mkdir -p 01-Mod_databases/miRBase/
    mkdir -p 01-Mod_databases/sRNAanno/
    mkdir -p 01-Mod_databases/PmiREN/

    # Simulate miRBase files
    if [ "${mirbase}" == "true" ]; then
        touch 01-Mod_databases/miRBase/mirbase_mature.fa
        touch 01-Mod_databases/miRBase/mirbase_hairpin.fa
    else
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_mature.fa
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_hairpin.fa
    fi

    # Simulate sRNAanno files
    if [ "${srnaanno}" == "true" ]; then
        touch 01-Mod_databases/sRNAanno/srnaanno_mature.fa
        touch 01-Mod_databases/sRNAanno/srnaanno_hairpin.fa
    else
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_mature.fa
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_hairpin.fa
    fi

    # Simulate PmiREN files
    if [ "${pmiren}" == "true" ]; then
        touch 01-Mod_databases/PmiREN/pmiren_mature.fa
        touch 01-Mod_databases/PmiREN/pmiren_hairpin.fa
    else
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_mature.fa
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_hairpin.fa
    fi

    # Create species_id.csv file
    echo "species,mirbase,pmiren,srnaanno,final_id" > species_ids.csv
    input_species_pre=\$(echo "${species}" | tr -d '[]' | sed 's/, */,/g')
    IFS=',' read -ra input_species_array <<< "\$input_species_pre"
    for species_name in "\${input_species_array[@]}"; do
        
        # Create identifier
        identifier=\$(echo \$species_name | awk '{print tolower(substr(\$1,1,1)) substr(\$2,1,2)}')

        # Create row
        if [ "${mirbase}" == "true" ]; then
            id_mirbase=\$identifier
        else
            id_mirbase=NULL
        fi
        if [ "${srnaanno}" == "true" ]; then
            id_srnaanno=\$identifier
        else
            id_srnaanno=NULL
        fi
        if [ "${pmiren}" == "true" ]; then
            id_pmiren=\$identifier
        else
            id_pmiren=NULL
        fi

        # Crete row
        row="\${species_name},\${id_mirbase},\${id_srnaanno},\${id_pmiren},\$identifier"
        echo \$row >> species_ids.csv
    done
    """
}