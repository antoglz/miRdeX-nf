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

    output:
    path "species_ids.csv"                                                                                  , emit: species_ids
    tuple path("01-Mod_databases/miRBase/*_mature.fa")  , path("01-Mod_databases/miRBase/*_hairpin.fa")     , emit: mirbase
    tuple path("01-Mod_databases/sRNAanno/*_mature.fa") , path("01-Mod_databases/sRNAanno/*_hairpin.fa")    , emit: srnaanno
    tuple path("01-Mod_databases/PmiREN/*_mature.fa")   , path("01-Mod_databases/PmiREN/*_hairpin.fa")      , emit: pmiren

    script:
    """
    # FUNCTIONS
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
                split(\$1, arr, "-");
                id = substr(arr[1], 2);
                species = \$3 " " \$4;

                if (!(id in names))
                    names[id] = species;
            }
            END {
                for (id in names)
                    print id "\t" names[id];
            }
        ' mirbase/hairpin.fa mirbase/mature.fa | sort > mirbase/mirbase_species.txt

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

        # Create strings for the execution of 01-Database_id_resolver.sh
        while IFS= read -r mature_file; do
            pmiren_m_files+="--pmiren \$mature_file "
        done < <(find pmiren/\$pmiren_version -type f -name "*_mature.fa" | sort)

        while IFS= read -r hairpin_file; do
            # Rename stemloop files to hairpin
            new_hairpin_file="\${hairpin_file%_stemloop.fa}_hairpin.fa"
            mv "\$hairpin_file" "\$new_hairpin_file"
            pmiren_p_files+="--phairpin \$new_hairpin_file "
        done < <(find pmiren/\$pmiren_version -type f -name "*_stemloop.fa" | sort)
    fi

    ## 1. RESOLVE IDENTIFIERS
    ############################################################################

    # Create the output directories
    mkdir -p 01-Mod_databases/miRBase/
    mkdir -p 01-Mod_databases/sRNAanno/
    mkdir -p 01-Mod_databases/PmiREN/

    # Execute the programm
    bash 01-Database_id_resolver.sh \
        --inputsps "${species}" \
        \$mirbase_files \
        \$srnaanno_files \
        \$pmiren_m_files \
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
    input_species_pre=\$(echo "${species}" | tr -d '[]' | sed 's/, */,/g')
    IFS=',' read -ra input_species_array <<< "\$input_species_pre"
    for species_name in "\${input_species_array[@]}"; do
        # Select the required ids from species_id_db.csv
        grep "\$species_name" species_ids_db.csv >> species_ids.csv
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
