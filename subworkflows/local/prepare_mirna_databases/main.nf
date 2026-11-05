/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { DOWNLOAD_DB         } from '../../../modules/local/download_db'
include { SPLIT_DB_BY_SPECIES } from '../../../modules/local/split_db_by_species'
include { DNA_RNA_CONVERTER   } from '../../../modules/local/dna_rna_converter'
include { COLLAPSE_DB         } from '../../../modules/local/collapse_db'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { formatDbChannel } from '../utils_mirdex_pipeline/'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW PREPARE_MIRNA_DATABASES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_MIRNA_DATABASES {
    take:
        species             // channel: val([species_name_1, species_name_2, ..., species_name_n])
        databases           // string: e.g. 'mirbase,pmiren,srnaanno'
        species_synonyms    // path: TSV file containing species synonyms

    main:

    // Get the list of databases to be used
    def db_list = databases.split(',').collect { it.trim() }
    
    // Find out which databases will be used.
    def download_mirbase = db_list.contains('mirbase')
    def download_srnaanno = db_list.contains('srnaanno')
    def download_pmiren = db_list.contains('pmiren')

    // Use databases previously downloaded and processed
    if (params.reuse_dbs) {

        // Create a channel for each database using the provided directories
        ch_mirbase = Channel.fromFilePairs("${params.mirbase}/*_{mature,hairpin}.fa", flat: true)
            .map { _id, hairpin, mature -> tuple(mature, hairpin) }
        ch_srnaanno = Channel.fromFilePairs("${params.srnaanno}/*_{mature,hairpin}.fa", flat: true)
            .map { _id, hairpin, mature -> tuple(mature, hairpin) }
        ch_pmiren = Channel.fromFilePairs("${params.pmiren}/*_{mature,hairpin}.fa", flat: true)
            .map { _id, hairpin, mature -> tuple(mature, hairpin) }

        // Create a channel with the species ids file
        ch_species_ids = Channel.fromPath("${params.species_ids}", checkIfExists: true)
    
    } else {

        // Download the required databases
        DOWNLOAD_DB(
            species,
            download_mirbase,
            download_srnaanno,
            download_pmiren,
            species_synonyms
        )
        ch_mirbase = DOWNLOAD_DB.out.mirbase
        ch_srnaanno = DOWNLOAD_DB.out.srnaanno
        ch_pmiren = DOWNLOAD_DB.out.pmiren
        ch_species_ids = DOWNLOAD_DB.out.species_ids
    }

    // Change the format of the db channels
    mirbase_formatted = formatDbChannel(ch_mirbase, 'mirbase')
    srnaanno_formatted = formatDbChannel(ch_srnaanno, 'srnaanno')
    pmiren_formatted = formatDbChannel(ch_pmiren, 'pmiren')

    // Combine the database channels in
    all_databases_ch = mirbase_formatted
            .concat(pmiren_formatted)
            .concat(srnaanno_formatted)

    // Get a list with the species ids
    ch_species_ids
        .splitCsv( header: true, sep: ',' )
        .map{ sps_row -> [sps_row.final_id,sps_row.species]}
        .set{ch_species_ids}
    
    // Convert RNA sequences to DNA
    DNA_RNA_CONVERTER(all_databases_ch, true)

    // Combine databases channel with species_ids channel
    ch_species_ids
        .combine(DNA_RNA_CONVERTER.out.file_conv)
        .map{ species_id, species_name, map, file_db ->
            [[id:"${species_id}_${map.id}", db_id:map.id, species_id: species_id, species: species_name], species_id, file_db]
        }.set{ch_species_ids_db}

    // Create a DB file for each species
    SPLIT_DB_BY_SPECIES(ch_species_ids_db)

    // Collapse the reference databases
    COLLAPSE_DB(SPLIT_DB_BY_SPECIES.out.species_db)

    // Combine the input files channel with the databases channel
    COLLAPSE_DB.out.coldb
        .map{ meta, file -> [meta.species, meta, file] }
        .groupTuple(by:0)
        .map { species_info ->
            // Get the required information from the channel element
            def species_n = species_info[0]
            def db_info = species_info[1]
            def paths = species_info[2]

            //Create a map object with the path of the databases
            def result = [:]
            db_info.eachWithIndex { db, i ->
                def db_name = db.db_id
                def path = paths[i]
                result[db_name] = path
            }
            return [species_id: db_info.species_id[0], species: species_n] + result
        }
        .set{ ch_databases }

    // Select the reference database according to its order of preference
    ch_databases
        .map { meta ->
            def matureFile = null
            def precursorFile = null

            def db = db_list.find { db ->
                def mature = meta["${db}_mature"]
                def precursor = meta["${db}_precursor"]
                if (!mature.getName().contains('EMPTY') && !precursor.getName().contains('EMPTY')) {
                    matureFile = mature
                    precursorFile = precursor
                    return true
                }
                return false
            }

            if (!db) {
                db = 'NULL'
                matureFile = 'NULL'
                precursorFile = 'NULL'
            }

            def result = meta.findAll { key, _value -> 
                key != 'pmiren_mature' &&
                key != 'srnaanno_mature' &&
                key != 'srnaanno_precursor' &&
                key != 'mirbase_precursor' &&
                key != 'pmiren_precursor' &&
                key != 'mirbase_mature'
            }

            result['database'] = db
            result['mature'] = matureFile
            result['precursor'] = precursorFile
            return result
        }
        .filter { meta -> meta.database != 'NULL' }
        .set { ch_selected_database }

    emit:
        dbs = ch_selected_database     // channel: [ id:val(str), db_id:val(str), species_id:val(str), species: val(str), database:val(str), mature: path(File), precursor:path(File) ]
}
