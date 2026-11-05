/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { COUNTS } from '../../../modules/local/counts'
include { RPM    } from '../../../modules/local/rpm'
include { ADD_COUNTS_TO_ISOMIRS_DF as ADD_RAW_COUNTS_TO_ISOMIRS_DF  } from '../../../modules/local/add_counts_to_isomirs_df'
include { ADD_COUNTS_TO_ISOMIRS_DF as ADD_RPM_TO_ISOMIRS_DF         } from '../../../modules/local/add_counts_to_isomirs_df'
include { ISOMIRS_MIRNA_CLASSIFICATION                              } from '../../../modules/local/isomirs_mirna_classification'
include { FILTER_ISOMIRS_BY_ABUNDANCE                               } from '../../../modules/local/filter_isomirs_by_abundance'    

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//

include { summaryToTsv } from '../utils_mirdex_pipeline'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PREPARE_MIRNA_DATABASES } from '../prepare_mirna_databases'
include { ISOMIRS_IDENTIFICATION  } from '../isomirs_identification'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { SEQKIT_FQ2FA } from '../../../modules/nf-core/seqkit/fq2fa/main'
include { SEQKIT_RMDUP } from '../../../modules/nf-core/seqkit/rmdup/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_MATURE    } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_PRECURSOR } from '../../../modules/nf-core/blast/makeblastdb/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MIRNOTE {

    take:
    ch_input                // channel: [[id:val(id), species:val(species), genome:val(genome)], fasta/fastq]
    ch_species              // channel: [ 'species_name', 'species_name2', ... ]
    species_synonyms        // path: TSV file containing species synonyms
    databases               // string: "mirbase,srnaanno,pmiren"
    substitutions           // integer: 1
    five_add                // integer: 0
    three_add               // integer: 3
    ends_modification       // integer: 4
    only_ref_miRNAs         // bool: true o false

    main:
    
    // Create an empty channels
    ch_pipeline_summary  = Channel.empty()
    ch_gff3_files_meta   = Channel.empty()
    ch_versions          = Channel.empty()
    ch_output_gff3       = Channel.empty()

    // Sepecify the filter type based on the parameters
    def filter_type  = params.min_rpm_filt > 0 ? 'RPM' : (params.min_relative_abundance_filt > 0 ? 'Relative_abundance' : null)
    def filter_value = params.min_rpm_filt > 0 ? params.min_rpm_filt : (params.min_relative_abundance_filt > 0 ? params.min_relative_abundance_filt : null)
 
    /*
    ========================================================================================
        1. Prepare miRNA databases
    ========================================================================================
    */

    // Prepare the databases for the annotation
    PREPARE_MIRNA_DATABASES(ch_species, databases, species_synonyms)

    // Use the species name to combine the databases and input channels
    PREPARE_MIRNA_DATABASES.out.dbs
        .map{ meta -> [meta.species, meta]}
        .set{ ch_species_db }

    // Add info to summary channel
    ch_input
        .map{ meta, file -> [meta.species, meta, file]}
        .groupTuple(by:0)
        .join(ch_species_db, by:0, remainder:true)
        .filter { it[1] != null }
        .flatMap { _species, meta_list, _fastq_list, db_info_map -> 

            // Check if the species is in at least one database
            def species_db = db_info_map != null ? 'PASS' : 'REJECT'
            def db = db_info_map != null ? db_info_map.database : 'NA'

            return meta_list.indices.collect { i ->
                def new_meta = meta_list[i].clone()
                return [new_meta.id, new_meta + [species_db:species_db, database:db]]
            }
        }
        .set{ ch_pipeline_summary }

    // Combine databases and input channels
    ch_input
        .map{ meta, file -> [meta.species, meta, file]}
        .combine(ch_species_db, by:0)
        .map{_species, meta1, file, meta2 ->
            def meta_new1 = meta1.clone()
            def meta_new2 = meta2.clone()
            [meta_new1 + meta_new2, file]}
        .set{ ch_input_sp }

    /*
    ========================================================================================
        2. Quantify sequences
    ========================================================================================
    */

    // Check whether a raw counts or RPM threshold is to be used.
    if (filter_type != null) {
        
        // Calculate the raw counts
        COUNTS(ch_input_sp)

        // Save the raw counts inot a channel
        COUNTS.out.raw.set{ ch_raw_counts }

        // Calculate RPM
        RPM(ch_raw_counts)

        // Save the raw counts inot a channel
        RPM.out.rpm.set{ ch_rpm }
    }

    /*
    ========================================================================================
        3. Prepare input files
    ========================================================================================
    */

    // Identify FASTA and FASTQ files in the input channel
    sub_fastq_ch = ch_input_sp.filter{ _meta, file -> file.name ==~ /.*\.(fastq|fq)(\.gz)?$/ }
    sub_fasta_ch = ch_input_sp.filter{ _meta, file -> file.name ==~ /.*\.(fasta|fa)(\.gz)?$/ }

    // Convert FASTQ files to FASTA
    SEQKIT_FQ2FA(sub_fastq_ch)

    // Save the software version
    ch_versions = ch_versions.mix(SEQKIT_FQ2FA.out.versions)

    // Merge the contents of the two previous channels into the same channel
    sub_fasta_ch
        .mix(SEQKIT_FQ2FA.out.fasta)
        .map{ meta, file -> [ meta + [id: "${meta.id}.rmdup", run: meta.id], file]}
        .set{ch_input_sp}

    // Remove duplicates in the input files
    SEQKIT_RMDUP(ch_input_sp)

    // Save the software version
    ch_versions = ch_versions.mix(SEQKIT_RMDUP.out.versions)

    // Save the results in the input channel
    SEQKIT_RMDUP.out.fastx
        .map{ meta, file -> [meta + [id: meta.run], file]}
        .set{ch_input_sp}

    /*
    ========================================================================================
        4. Identify isomiRs
    ========================================================================================
    */

    // Preapare the input channel for ISOMIRS_IDENTIFICATION worflow
    ch_input_sp
        .map{ meta, file -> [meta, file, meta.mature, meta.precursor, meta.genome]}
        .set{ ch_to_identify_isomirs }
    
    // Only reference (canonical) miRNAs
    if (only_ref_miRNAs) {
        // Identify isomiRs
        ISOMIRS_IDENTIFICATION(ch_to_identify_isomirs, true)
    } else {
        // Identify isomiRs
        ISOMIRS_IDENTIFICATION(ch_to_identify_isomirs, false)
    }

    // Save the software version
    ch_versions = ch_versions.mix(ISOMIRS_IDENTIFICATION.out.versions)

    // Get the numnber of sequences identify as potential isomiRs (summary)
    ISOMIRS_IDENTIFICATION.out.iso
        .map { meta, file ->
            def iso = file.readLines()
                        .collect { it.split('\t', -1)[13] }
                        .findAll { it }
                        .toSet()
            [meta.id, iso.size()]
        }
        .join(ch_pipeline_summary, remainder:true)
        .map{ item ->
            // Required items
            def id = item[0]
            def num_iso = item[1]
            def meta = item[2]

            // Fields to add
            def additionalFields = num_iso ? [ num_pot_isomirs: num_iso ] : [ num_pot_isomirs: 'NA']
            
            // Return
            [id, meta + additionalFields]
        }
        .set{ ch_pipeline_summary }

    /*
    ========================================================================================
        5. Add raw counts and RPM to the Blastn isomiRs dataframe
    ========================================================================================
    */

    // Input channel for ISOMIRS_MIRNA_CLASSIFICATION module
    ISOMIRS_IDENTIFICATION.out.iso
        .filter { _meta, file ->
            !file.getName().contains('EMPTY')
        }
        .set{ch_isomirs}

    // Add counts to isomiRs dataframe
    if (filter_type != null) {

        // Prepare the channel for merging
        ch_isomirs
            .map{ meta, file -> [meta.id, meta, file]}
            .set{ ch_isomirs }

        // Prepare the channel for merging
        ch_raw_counts
            .map{ meta, file -> [meta.id, meta, file]}
            .combine(ch_isomirs, by:0)
            .map{ _id, _meta_c, file_c, meta_i, file_i ->
                [meta_i + [id: "${meta_i.id}.rawc", prev_id:meta_i.id], file_i, file_c]
            }
            .set { ch_isomirs_add_raw_counts }
        
        // Add raw counts to isomiRs dataframe
        ADD_RAW_COUNTS_TO_ISOMIRS_DF(ch_isomirs_add_raw_counts)

        // Set the original id and prepare the channel for merging
        ADD_RAW_COUNTS_TO_ISOMIRS_DF.out.isocounts
            .map { meta, file ->
                def updatedMeta = meta + [id: meta.prev_id]
                updatedMeta.remove('prev_id') 
                return [updatedMeta.id, updatedMeta, file]
            }
            .set { ch_isomirs }

        // Prepare the channel for merging
        ch_rpm
            .map{ meta, file -> [meta.id, meta, file]}
            .combine(ch_isomirs, by:0)
            .map{ _id, _meta_c, file_c, meta_i, file_i ->
                [meta_i + [id: "${meta_i.id}.rpmc", prev_id:meta_i.id], file_i, file_c]
            }
            .set { ch_isomirs_add_rpm }

        // Add RPM to isomiRs dataframe
        ADD_RPM_TO_ISOMIRS_DF(ch_isomirs_add_rpm)
        
        // Set the original id and prepare the channel for merging
        ADD_RPM_TO_ISOMIRS_DF.out.isocounts
            .map { meta, file ->
                def updatedMeta = meta + [id: meta.prev_id]
                updatedMeta.remove('prev_id') 
                return [updatedMeta, file]
            }
            .set { ch_isomirs }    
    }

    /*
    ========================================================================================
        6. Classify isomiRs
    ========================================================================================
    */
    
    // Classify isomiRs
    ISOMIRS_MIRNA_CLASSIFICATION(
        ch_isomirs,
        substitutions,
        five_add,
        three_add,
        ends_modification
    )

    // Save the software version
    ch_versions = ch_versions.mix(ISOMIRS_MIRNA_CLASSIFICATION.out.versions)

    // Save the results in the output channel
    ch_output_gff3 = ISOMIRS_MIRNA_CLASSIFICATION.out.gff3

    // Add summary info to the summary channel
    ISOMIRS_MIRNA_CLASSIFICATION.out.sum
        .splitCsv( header: true, sep: '\t' )
        .map{ meta, info -> [meta.id, info]}
        .join(ch_pipeline_summary, remainder:true)
        .map { item ->
            def iso_info = item[1]
            def pip_summary = item[2]

            // Fields to add
            def additionalFields = iso_info ? [
                num_isomirs: iso_info.num_isomirs,
                ref_miRNA: iso_info.ref_miRNA,
                iso_5p: iso_info.iso_5p,
                iso_3p: iso_info.iso_3p,
                iso_add5p: iso_info.iso_add5p,
                iso_add3p: iso_info.iso_add3p,
                iso_snv_seed: iso_info.iso_snv_seed,
                iso_snv_central_offset: iso_info.iso_snv_central_offset,
                iso_snv_central: iso_info.iso_snv_central,
                iso_snv_central_supp: iso_info.iso_snv_central_supp,
                iso_snv: iso_info.iso_snv,
                mixed: iso_info.mixed,
                mixed_shift: iso_info.mixed_shift,
                undefined: iso_info.undefined,
                num_isomirs_filt: 'NA',
                ref_miRNA_filt: 'NA',
                iso_5p_filt: 'NA',
                iso_3p_filt: 'NA',
                iso_add5p_filt: 'NA',
                iso_add3p_filt: 'NA',
                iso_snv_seed_filt: 'NA',
                iso_snv_central_offset_filt: 'NA',
                iso_snv_central_filt: 'NA',
                iso_snv_central_supp_filt: 'NA',
                iso_snv_filt: 'NA',
                mixed_filt: 'NA',
                mixed_shift_filt: 'NA',
                undefined_filt: 'NA'
            ] : [
                num_isomirs: 'NA',
                ref_miRNA: 'NA',
                iso_5p: 'NA',
                iso_3p: 'NA',
                iso_add3p: 'NA',
                iso_add5p: 'NA',
                iso_snv_seed: 'NA',
                iso_snv_central_offset: 'NA',
                iso_snv_central: 'NA',
                iso_snv_central_supp: 'NA',
                iso_snv: 'NA',
                mixed: 'NA',
                mixed_shift: 'NA',
                undefined: 'NA',
                num_isomirs_filt: 'NA',
                ref_miRNA_filt: 'NA',
                iso_5p_filt: 'NA',
                iso_3p_filt: 'NA',
                iso_add5p_filt: 'NA',
                iso_add3p_filt: 'NA',
                iso_snv_seed_filt: 'NA',
                iso_snv_central_offset_filt: 'NA',
                iso_snv_central_filt: 'NA',
                iso_snv_central_supp_filt: 'NA',
                iso_snv_filt: 'NA',
                mixed_filt: 'NA',
                mixed_shift_filt: 'NA',
                undefined_filt: 'NA'
            ]
            // Remove undesired fields for summary channel
            pip_summary.remove('genome')
            pip_summary.remove('single_end')
            
            // Add the additional fields
            [pip_summary.id, pip_summary + additionalFields]
        }
        .set{ ch_pipeline_summary }


    /*
    ========================================================================================
        7. Filter isomiRs by abundance
    ========================================================================================
    */

    if (filter_type != null && !only_ref_miRNAs) {

        // Create a channel with the metadata of the GFF3 files
        ISOMIRS_MIRNA_CLASSIFICATION.out.gff3
            .map{ meta,_file -> [meta.id, meta]}
            .set{ ch_gff3_files_meta }

        // Group gff3 files by species
        ISOMIRS_MIRNA_CLASSIFICATION.out.gff3
            .map { meta, file -> [meta.species, file]}
            .groupTuple(by:0)
            .map{ species, files -> [[id: species], files]}
            .set { ch_gff3_files_by_species }

        // Filter isomiRs by abundance
        FILTER_ISOMIRS_BY_ABUNDANCE(
            ch_gff3_files_by_species,
            filter_type,
            filter_value,
            params.min_samples_filt
        )

        // Prepare the output channel
        FILTER_ISOMIRS_BY_ABUNDANCE.out.filt_gff3
            .flatMap { _meta, files ->
                files.collect { f ->
                    def id = f.getBaseName().replaceFirst(/\.filtered$/, '')
                    [id, f]
                }
            }
            .combine(ch_gff3_files_meta, by:0)
            .map{ _id, file, meta -> [meta, file]}
            .set{ ch_gff3_out }
        
        // Save the results in the output channel
        ch_output_gff3 = ch_gff3_out
        
        // Add summary info to the summary channel
        FILTER_ISOMIRS_BY_ABUNDANCE.out.sum
            .splitCsv( header: true, sep: '\t' )
            .map{ info -> [info.sample, info]}
            .join(ch_pipeline_summary, remainder:true)
            .map { item ->
                def iso_info = item[1]
                def pip_summary = item[2]

                // Fields to add
                def additionalFields = iso_info ? [
                    num_isomirs_filt: iso_info.num_isomirs,
                    ref_miRNA_filt: iso_info.ref_miRNA,
                    iso_5p_filt: iso_info.iso_5p,
                    iso_3p_filt: iso_info.iso_3p,
                    iso_add5p_filt: iso_info.iso_add5p,
                    iso_add3p_filt: iso_info.iso_add3p,
                    iso_snv_seed_filt: iso_info.iso_snv_seed,
                    iso_snv_central_offset_filt: iso_info.iso_snv_central_offset,
                    iso_snv_central_filt: iso_info.iso_snv_central,
                    iso_snv_central_supp_filt: iso_info.iso_snv_central_supp,
                    iso_snv_filt: iso_info.iso_snv,
                    mixed_filt: iso_info.mixed,
                    mixed_shift_filt: iso_info.mixed_shift,
                    undefined_filt: iso_info.undefined
                ] : [
                    num_isomirs_filt: 'NA',
                    ref_miRNA_filt: 'NA',
                    iso_5p_filt: 'NA',
                    iso_3p_filt: 'NA',
                    iso_add5p_filt: 'NA',
                    iso_add3p_filt: 'NA',
                    iso_snv_seed_filt: 'NA',
                    iso_snv_central_offset_filt: 'NA',
                    iso_snv_central_filt: 'NA',
                    iso_snv_central_supp_filt: 'NA',
                    iso_snv_filt: 'NA',
                    mixed_filt: 'NA',
                    mixed_shift_filt: 'NA',
                    undefined_filt: 'NA'
                ]
                
                // Add the additional fields
                [pip_summary.id, pip_summary + additionalFields]
            }
            .set{ ch_pipeline_summary }
    }

    // Remove id from the summary channel
    ch_pipeline_summary
        .map{ _id, info -> info }
        .set{ ch_pipeline_summary }

    emit:
    annotation     = ch_output_gff3                          // channel: [meta, file.gff3]
    summary        = ch_pipeline_summary                     // channel: [id:, num_isomirs:, etc]
    versions       = ch_versions                             // channel: [ path(versions.yml) ]

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/