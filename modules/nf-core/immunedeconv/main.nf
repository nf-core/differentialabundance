process IMMUNEDECONV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/r-immunedeconv:2.1.2--e1bb1ea1cf505cb3"

    input:
    tuple val(meta), path(input_file), val(method), val(function)

    output:
    tuple val(meta), path("deconvolution_results.tsv"), emit: deconv_table
    tuple val(meta), path("*png"),                      emit: deconv_plots, optional: true
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run Rscript for deconvolution
    Rscript -e "
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(immunedeconv)
    library(tibble)
    library(readr)

    # Load the TSV file
    gene_expression_matrix <- readr::read_tsv('${input_file}') %>%
        as.data.frame() %>%
        tibble::column_to_rownames('gene_symbol')

    # Generate results
    result <- immunedeconv::${function}(gene_expression_matrix, method = '${method}')

    # Save the result to a CSV file
    readr::write_tsv(result, 'deconvolution_results.tsv')

    # Plot and save results
    # Plot 1: Stacked bar chart
    plot1 <- result %>%
        gather(sample, fraction, -cell_type) %>%
        ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
        geom_bar(stat = 'identity') +
        coord_flip() +
        scale_fill_brewer(palette = 'Paired') +
        scale_x_discrete(limits = rev(levels(result)))

    # Save Plot 1
    ggsave('plot1_stacked_bar_chart.png', plot = plot1, dpi = 300, width = 10, height = 8)

    # Plot 2: Points with facets
    plot2 <- result %>%
        gather(sample, score, -cell_type) %>%
        ggplot(aes(x = sample, y = score, color = cell_type)) +
        geom_point(size = 4) +
        facet_wrap(~cell_type, scales = 'free_x', ncol = 3) +
        scale_color_brewer(palette = 'Paired', guide = FALSE) +
        coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    # Save Plot 2
    ggsave('plot2_points_with_facets.png', plot = plot2, dpi = 300, width = 12, height = 10)
    "

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-immunedeconv: \$(Rscript -e "cat(as.character(packageVersion('immunedeconv')))")
    END_VERSIONS
    """

    stub:
    """
    touch deconvolution_results.tsv
    touch deconvolution_results.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-immunedeconv: \$(Rscript -e "cat(as.character(packageVersion('immunedeconv')))")
    END_VERSIONS
    """
}
