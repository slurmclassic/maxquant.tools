

prettyToSlotNameMap = list('Identification.type' = 'ident_type',
                           'Intensity' = 'intensity',
                           'LFQ.intensity' = 'lfq_intensity')

#' Load a MaxQuant 'peptides.txt' dataset.
#' @param path Full path to the peptides.txt file.
#' @return An object of class 'MaxQuantPeptidesResult'.
#' @export
readMaxQuantPeptideMeasurements <- function(path) {
    data <- read.delim(path, stringsAsFactors = F)

    # Find which of the headers are present and calculate the index difference between them.
    possible_first_list_columns = c('Identification.type', 'Intensity', 'LFQ.intensity', 'Experiment')
    possible_sample_names = sapply(colnames(data), function(col){
        for (possible_col in possible_first_list_columns ) {
            col = gsub(paste0(possible_col, '.'), '', col)
        }

        return(col)
    })

    names(possible_sample_names) <- NULL
    tabled = table(possible_sample_names)
    sample_names = names(tabled)[tabled > 1]

    possible_first_list_columns = c('Identification.type', 'Intensity', 'LFQ.intensity')
    present_headers = possible_first_list_columns[sapply(possible_first_list_columns, function(possible_column) {
        return(any(sapply(colnames(data), function(actual_column) {
            return(any(sapply(sample_names, function(s) {
                return(actual_column == paste0(possible_column, '.', s))
            })))

        })))
    })]

    maxquant <- new('MaxQuantPeptidesResult')
    maxquant@sequence = data$Sequence

    for (present_header in present_headers) {

        slot(maxquant, prettyToSlotNameMap[[present_header]]) <- as.matrix(data[,paste0(present_header, '.', sample_names)])

        colnames(slot(maxquant, prettyToSlotNameMap[[present_header]])) <- sample_names
        rownames(slot(maxquant, prettyToSlotNameMap[[present_header]])) <-data$Sequence


    }

    slot(maxquant, 'reverse') <- data$Reverse == '+'
    names(maxquant@reverse) <- data$Sequence
    slot(maxquant, 'potential_contaminant') <- data$Potential.contaminant == '+'
    names(maxquant@potential_contaminant) <- data$Sequence



    slot(maxquant, 'accession') <- data$Proteins
    slot(maxquant, 'samples') <- sample_names
    split = strsplit(path, '/')[[1]]
    slot(maxquant, 'experiment') <- gsub(".txt", "", split[length(split)])

    return(maxquant)
}

setClass("MaxQuantPeptidesResult", slots=list(
    accession = "character",
    sequence = "character",
    intensity = "matrix",
    lfq_intensity = "matrix",
    pvalues = "numeric",
    samples = "character",
    groups = "character",
    ident_type = "matrix",
    reverse = "logical",
    potential_contaminant = "logical",
    experiment = "character"))

setMethod('show', 'MaxQuantPeptidesResult', function(object){
    cat('Experiment:', '\t', object@experiment, '\n')
    cat('Samples/Groups:\n')
    print(data.frame(Sample = object@samples, Group = object@groups), row.names=F)
})

mergeMaxQuantPeptideResults <-function(mqs, new_experiment_name, session) {
    new_mq = new('MaxQuantPeptidesResult')
    slot(new_mq, 'samples') <- unique(unlist(sapply(mqs, function(mq){return(mq@samples)}), recursive = F))
    slot(new_mq, 'sequence') <- unique(unlist(sapply(mqs, function(mq){return(mq@sequence)}), recursive = F))


    if (all(sapply(mqs, function(mq){return(dim(mq@intensity)[1] > 0)}))) {
        slot(new_mq, 'intensity') <- sapply(new_mq@samples, function(sample) {sapply(new_mq@sequence, function(sequence){
            val = 0
            for (mq in mqs) {
                if (sequence %in% mq@sequence && sample %in% mq@samples) {
                    val = mq@intensity[sequence, sample]
                    break
                }
            }

            val
        })})
    }

    if (all(sapply(mqs, function(mq){return(dim(mq@lfq_intensity)[1] > 0)}))) {
        slot(new_mq, 'lfq_intensity') <- sapply(new_mq@samples, function(sample) {sapply(new_mq@sequence, function(sequence){
            val = 0
            for (mq in mqs) {
                if (sequence %in% mq@sequence && sample %in% mq@samples) {
                    val = mq@lfq_intensity[sequence, sample]
                    break
                }
            }

            val
        })})
    }

    slot(new_mq, 'reverse') <- sapply(new_mq@sequence, function(sequence){
        vals = c(NA)
        for (mq in mqs) {
            if (sequence %in% mq@sequence) {
                vals = c(vals, mq@reverse[sequence])
            }
        }

        any(vals, na.rm = T)
    })

    slot(new_mq, 'potential_contaminant') <- sapply(new_mq@sequence, function(sequence){
        vals = c()
        for (mq in mqs) {
            if (sequence %in% mq@sequence) {
                vals = c(vals, mq@potential_contaminant[sequence])
            }
        }

        any(vals)
    })

    slot(new_mq, 'accession') <- sapply(new_mq@sequence, function(sequence){
        vals = c()
        for (mq in mqs) {
            if (sequence %in% mq@sequence) {
                vals = c(vals, mq@accession[sequence])
            }
        }

        paste(vals, sep=';', collapse = ';')
    })

    new_mq
}
