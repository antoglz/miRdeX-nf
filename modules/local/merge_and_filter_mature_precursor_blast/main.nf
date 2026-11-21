process MERGE_AND_FILTER_MATURE_PRECURSOR_BLAST {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), path(mature), path(precursor)

    output:
    tuple val(meta), path("*mpblast.tsv"), emit: mpblast

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash
    
    ### FUNCTIONS
    expand_and_clean_tsv () {
        # Arguments
        local input_file="\${1}"
        local output_file="\${2}"

        # Expand and clean TSV table
        awk -F '\t' '{
            split(\$2, ids, "|");
            for (i in ids) {
                original_id = ids[i];
                id_with_arm = original_id;
                id_clean = original_id;

                # Remove undesired strings
                gsub(/-?(Known|mature|star)-?/, "-", id_with_arm);
                gsub(/^-+|-+\$/, "", id_with_arm);
                gsub(/--+/, "-", id_with_arm);
                
                # Remove duplicated 5p/3p arms
                n = split(id_with_arm, parts, "-")
                delete seen
                id_with_arm = ""
                for (i = 1; i <= n; i++) {
                    if ((parts[i] == "5p" || parts[i] == "3p") && parts[i] in seen)
                        continue
                    seen[parts[i]] = 1
                    id_with_arm = (id_with_arm == "" ? parts[i] : id_with_arm "-" parts[i])
                }

                # Remove undesired strings (5p/3p too)
                gsub(/-?(5p|3p)-?/, "-", id_clean);
                gsub(/-?(Known|mature|star)-?/, "-", id_clean);
                gsub(/[*]/, "", id_clean);
                gsub(/^-+|-+\$/, "", id_clean);
                gsub(/--+/, "-", id_clean);

                print \$1, id_clean, id_with_arm, \$0;
            }
        }' OFS='\t' "\$input_file" > "\$output_file"
    }

    create_key_column (){
            
        # Arguments
        local input_file="\${1}"
        local output_file="\${2}"

        awk -F'\t' '{
            key = tolower(\$1 "|" \$2);
            print key, \$0;
        }' OFS='\t' "\$input_file" > "\$output_file"

    }

    select_most_likely_references(){

        # Arguments
        local input_file="\${1}"
        local output_file="\${2}"

        awk -F'\\t' '
        {
            key = \$14;
            val = \$12;

            if (!(key in min)) {
                min[key] = val;
                lines[key] = \$0;
            } else if (val < min[key]) {
                min[key] = val;
                lines[key] = \$0;
            } else if (val == min[key]) {
                lines[key] = lines[key] "\\n" \$0;
            }
        }
        END {
            for (k in lines) {
                print lines[k];
            }
        }
        '  "\$input_file" > "\$output_file"
    }

    # Select the most probable alignment within the same precursor.
    select_best_alignment_precursor(){
        
        # Arguments
        local input_file="\${1}"
        local output_file="\${2}"

        awk -F'\\t' '
        {
            key = \$1 FS \$16
            val = \$27 + 0
            if (!(key in max) || val > max[key]) {
                max[key] = val
                line[key] = \$0
            }
        }
        END {
            for (k in line) print line[k]
        }
        ' \$input_file > \$output_file
    }

    ### MAIN
    main () {
    
        # Create temporary directory
        workdir=\$(pwd)
        mkdir -p "\$workdir/tmp"

        # Expand and clean TSV tables (Prepare mature miRNA identifiers to find their precursors)
        expand_and_clean_tsv "${mature}" \$workdir/tmp/mature_expanded_cleaned.tsv
        expand_and_clean_tsv "${precursor}" \$workdir/tmp/precursor_expanded_cleaned.tsv

        # Create the key column for the files to be merged
        create_key_column \$workdir/tmp/mature_expanded_cleaned.tsv \$workdir/tmp/mature_with_key.tsv
        create_key_column \$workdir/tmp/precursor_expanded_cleaned.tsv \$workdir/tmp/precursor_with_key.tsv

        # Sort the TSV files to be merged
        sort -T "\$workdir/tmp" -k1,1 tmp/mature_with_key.tsv > \$workdir/tmp/mature_sorted.tsv
        sort -T "\$workdir/tmp" -k1,1 tmp/precursor_with_key.tsv > \$workdir/tmp/precursor_sorted.tsv

        # Inner Join
        join -t  \$'\t' -1 1 -2 1 \$workdir/tmp/mature_sorted.tsv \$workdir/tmp/precursor_sorted.tsv > \$workdir/tmp/joined_table.tsv

        # Remove undesired columns
        cut --complement -f1,3,5,19,21,22,34 \$workdir/tmp/joined_table.tsv > \$workdir/tmp/isomirs_all.tsv
        
        # Select most likely references for each sequence
        select_most_likely_references \$workdir/tmp/isomirs_all.tsv \$workdir/tmp/isomirs_best_alig_within_same_pre.tsv 
        
        # Select the most probable alignment within the same precursor.
        select_best_alignment_precursor \$workdir/tmp/isomirs_best_alig_within_same_pre.tsv "${prefix}.mpblast.tsv"

        # Add empty suffix if the file is empty
        if [ ! -s "${prefix}.mpblast.tsv" ]; then
            mv "${prefix}.mpblast.tsv" "${prefix}.EMPTY_mpblast.tsv"
        fi
    }
    main "\$@"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create output file
    touch "${prefix}.mpblast.tsv"
    """
}