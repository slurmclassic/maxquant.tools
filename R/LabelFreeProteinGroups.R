

prettyToSlotNameMap = list('Peptides' = 'peptides',
                           'Razor...unique.peptides' = 'razor_unique_peptides',
                           'Unique.peptides' = 'unique_peptides',
                           'Identification.type' = 'ident_type',
                           'Intensity' = 'intensity',
                           'LFQ.intensity' = 'lfq_intensity')

#' Load a MaxQuant 'proteinGroups.txt' dataset.
#' @param path Full path to the proteinGroups.txt file.
#' @return An object of class 'MaxQuantProteinGroupsResult'.
#' @export
readMaxQuantProteinMeasurements <- function(path) {
    data <- read.delim(path, stringsAsFactors = F)

    # Find which of the headers are present and calculate the index difference between them.
    possible_first_list_columns = c('Peptides', 'Razor...unique.peptides', 'Unique.peptides', 'Identification.type', 'Intensity', 'LFQ.intensity')
    possible_sample_names = sapply(colnames(data), function(col){
        for (possible_col in possible_first_list_columns ) {
            col = gsub(paste0(possible_col, '.'), '', col)
        }

        return(col)
    })

    names(possible_sample_names) <- NULL
    tabled = table(possible_sample_names)
    sample_names = names(tabled)[tabled > 1]

    present_headers = possible_first_list_columns[sapply(possible_first_list_columns, function(possible_column) {
        return(any(sapply(colnames(data), function(actual_column) {
            return(any(sapply(sample_names, function(s) {
                return(actual_column == paste0(possible_column, '.', s))
                })))

            })))
        })]



    maxquant <- new('MaxQuantProteinGroupsResult')

    for (present_header in present_headers) {

        slot(maxquant, prettyToSlotNameMap[[present_header]]) <- as.matrix(data[,paste0(present_header, '.', sample_names)])

        colnames(slot(maxquant, prettyToSlotNameMap[[present_header]])) <- sample_names
        rownames(slot(maxquant, prettyToSlotNameMap[[present_header]])) <-data$Protein.IDs


    }

    slot(maxquant, 'reverse') <- data$Reverse == '+'
    names(maxquant@reverse) <- data$Protein.IDs
    slot(maxquant, 'potential_contaminant') <- data$Potential.contaminant == '+'
    names(maxquant@potential_contaminant) <- data$Protein.IDs

    slot(maxquant, 'peptide_counts_all') <- sapply(data$Peptide.counts..all., function(c) {
        if (grepl(';', c)) {
            return(max(as.numeric(unlist(strsplit(c, ';')))))
        } else {
            return(as.numeric(c))
        }
    })

    names(maxquant@peptide_counts_all) <- data$Protein.IDs

    slot(maxquant, 'accession') <- data$Protein.IDs
    slot(maxquant, 'samples') <- sample_names
    split = strsplit(path, '/')[[1]]
    slot(maxquant, 'experiment') <- gsub(".txt", "", split[length(split)])

    return(maxquant)
}

setClass("MaxQuantProteinGroupsResult", slots=list(
    accession = "character",
    intensity = "matrix",
    lfq_intensity = "matrix",
    peptides = "matrix",
    peptide_counts_all = "numeric",
    razor_unique_peptides = "matrix",
    unique_peptides = "matrix",
    pvalues = "numeric",
    samples = "character",
    groups = "character",
    ident_type = "matrix",
    reverse = "logical",
    potential_contaminant = "logical",
    experiment = "character"))

setMethod('show', 'MaxQuantProteinGroupsResult', function(object){
    cat('Experiment:', '\t', object@experiment, '\n')
    cat('Samples/Groups:\n')
    print(data.frame(Sample = object@samples, Group = object@groups), row.names=F)
})

mergeMaxQuantProteinResults <-function(mqs, new_experiment_name, session) {
    new_mq = new('MaxQuantProteinGroupsResult')
    slot(new_mq, 'accession') <- unique(unlist(sapply(mqs, function(mq){return(mq@accession)}), recursive = F))
    slot(new_mq, 'samples') <- unique(unlist(sapply(mqs, function(mq){return(mq@samples)}), recursive = F))


    if (all(sapply(mqs, function(mq){return(dim(mq@intensity)[1] > 0)}))) {
        slot(new_mq, 'intensity') <- sapply(new_mq@samples, function(sample) {sapply(new_mq@accession, function(accession){
            val = 0
            for (mq in mqs) {
                if (accession %in% mq@accession && sample %in% mq@samples) {
                    val = mq@intensity[accession, sample]
                    break
                }
            }

            val
        })})
    }

    if (all(sapply(mqs, function(mq){return(dim(mq@lfq_intensity)[1] > 0)}))) {
        slot(new_mq, 'lfq_intensity') <- sapply(new_mq@samples, function(sample) {sapply(new_mq@accession, function(accession){
            val = 0
            for (mq in mqs) {
                if (accession %in% mq@accession && sample %in% mq@samples) {
                    val = mq@lfq_intensity[accession, sample]
                    break
                }
            }

            val
        })})
    }


    slot(new_mq, 'peptide_counts_all') <- sapply(new_mq@accession, function(accession){
        vals = c(NA)
        for (mq in mqs) {
            if (accession %in% mq@accession) {
                vals = c(vals, mq@peptide_counts_all[accession])
            }
        }


        max(as.numeric(vals), na.rm = T)
    })

    slot(new_mq, 'reverse') <- sapply(new_mq@accession, function(accession){
        vals = c(NA)
        for (mq in mqs) {
            if (accession %in% mq@accession) {
                vals = c(vals, mq@reverse[accession])
            }
        }

        any(vals, na.rm = T)
    })

    slot(new_mq, 'potential_contaminant') <- sapply(new_mq@accession, function(accession){
        vals = c()
        for (mq in mqs) {
            if (accession %in% mq@accession) {
                vals = c(vals, mq@potential_contaminant[accession])
            }
        }

        any(vals)
    })

    new_mq
}
