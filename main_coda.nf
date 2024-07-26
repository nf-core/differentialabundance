// include { PROPR_PROPR        } from '../../../modules/nf-core/propr/propr/main'
// include { PROPR_PROPD        } from '../../../modules/nf-core/propr/propd/main'
// include { PROPR_GREA         } from '../../../modules/nf-core/propr/grea/main'
// include { MYGENE             } from '../../../modules/nf-core/mygene/main'
include { EXPERIMENTAL        } from './subworkflows/local/experimental/main.nf'
include { fromSamplesheet } from 'plugin/nf-validation'


// These are local files from my Bachelor Thesis project, I am creating the ch_samples_and_matrix
// manually for testing but it should be be provided by the processing section of nf-core/differentialabundance
Counts_ch = Channel.fromPath(params.matrix)

Sample_ch = Channel.fromPath(params.input)
    .map{ it -> [[id: 'YMC'], it]}

ch_samples_and_matrix = Sample_ch.combine(Counts_ch)

// Convert the samplesheet.csv in a channel with the proper format
ch_tools = Channel.fromSamplesheet('tools')


// TO DO: This should be modified to run one path per default, not all
if (params.pathway == "all") {
    ch_tools
        .set{ ch_tools_single }
} else {
    ch_tools
        .filter{
            it[0]["pathway_name"] == params.pathway // TO DO: change pathway to path also in the tools_samplesheet file
        }
        //.view()
        .set{ ch_tools_single }
}
ch_tools_single.view()

workflow {
    EXPERIMENTAL(ch_samples_and_matrix, ch_tools_single)
    EXPERIMENTAL.out.output.view()
}

