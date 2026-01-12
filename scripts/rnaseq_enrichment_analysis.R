library(Biobase)
library(SummarizedExperiment)

#########################################
# CREAR ASSAY A PARTIR DE COUNTS MATRIX #
#########################################

# Cargar TSV ignorando la línea de comentario
counts_raw <- read.delim(
  "align/counts/counts_matrix.tsv",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
# Usar la primera columna como rownames (GeneID)
rownames(counts_raw) <- counts_raw[,1]

# Seleccionar solo columnas de conteo (quitando columnas de anotación)
counts_mat <- counts_raw[, -(1:6)]
counts_mat <- as.matrix(counts_mat)

# Limpiar nombres de columnas
colnames(counts_mat) <- gsub("align/bam/", "", colnames(counts_mat))
colnames(counts_mat) <- gsub(".sorted.bam", "", colnames(counts_mat))


###########################################
# CREAR COLDATA -NOMBRES DE LAS MUESTRAS- #
###########################################

samples <- read.delim("samples_baseline_trt.tsv", stringsAsFactors = FALSE)

# Usar un_accession como rownames para coincidir con columnas de counts
rownames(samples) <- samples$run_accession

# Crear columna condition: gut vs probiotic
coldata <- samples[, "condition", drop = FALSE]
colnames(coldata) <- "condition"

coldata$condition <- factor(coldata$condition, levels=c("Baseline", "Treatment"))

# Reordenar para que coincida con columns de counts_mat
coldata <- coldata[colnames(counts_mat), , drop = FALSE]

# Verificación
all(rownames(coldata) == colnames(counts_mat))
table(coldata$condition)


###########################################
# CREAR ROWDATA -FEATURES DE LAS MUESTRAS #
###########################################

library(tidyverse)
library(rtracklayer)

gff_all <- as.data.frame(import("genomes/gut_db_clean.gff"))

tax_map <- gff_all %>%
  filter(type == "region") %>%
  select(seqnames, Dbxref) %>%
  mutate(taxon_id = str_extract(as.character(Dbxref), "(?<=taxon:)[0-9]+")) %>%
  select(seqnames, taxon_id) %>%
  distinct(seqnames, .keep_all = TRUE) 

annot_info <- gff_all %>%
  filter(type == "CDS") %>%
  select(locus_tag, product, protein_id, Name, Ontology_term, Dbxref) %>%
  mutate(
    entrez_id = str_extract(as.character(Dbxref), "(?<=GeneID:)[0-9]+"),
    go_terms = Ontology_term
  ) %>%
  select(-Dbxref, -Ontology_term) %>%
  distinct(locus_tag, .keep_all = TRUE)

# counts_raw es el dataframe original de featureCounts
rowdata_final <- counts_raw[, 2:6] 
rowdata_final$locus_tag <- rownames(counts_raw)

rowdata_final <- rowdata_final %>%
  left_join(tax_map, by = c("Chr" = "seqnames")) %>%
  left_join(annot_info, by = "locus_tag")

anyDuplicated(rowdata_final$locus_tag) == 0

# 6. Sincronización con la matriz de conteos
rownames(rowdata_final) <- rowdata_final$locus_tag
rowdata_final <- rowdata_final[rownames(counts_mat), ]      

# 7. Limpieza de NAs 
rowdata_final$product[is.na(rowdata_final$product)] <- "hypothetical protein"
rowdata_final$taxon_id[is.na(rowdata_final$taxon_id)] <- "Unknown"


##################################
# CREAR EL SUMMARIZED EXPERIMENT #
##################################

se <- SummarizedExperiment(
  assays = list(counts = counts_mat),
  rowData = rowdata_final,
  colData = coldata
)

se

##################################
# NORMALIZACIÓN DE LAS LIBRERÍAS #
##################################

library(edgeR)

# 1. Primero eliminar muestras con 0 counts TOTALES
count_matrix <- assay(se, "counts")
total_counts <- colSums(count_matrix)

samples_to_remove <- which(total_counts == 0)
if (length(samples_to_remove) > 0) {
  cat("Eliminando", length(samples_to_remove), "muestras con 0 counts totales:\n")
  print(colnames(count_matrix)[samples_to_remove])
  se <- se[, -samples_to_remove]
  count_matrix <- assay(se, "counts")
}

# 2. Eliminar NA en condition
elim <- which(is.na(colData(se)$"condition"))
if (length(elim) > 0) {
  se <- se[, -elim]
  count_matrix <- assay(se, "counts")
}


print("Distribución por condición:")
print(table(colData(se)$condition))

# 4. Crear DGEList
dge <- DGEList(counts = count_matrix)

# 5. Filtrar genes de baja expresión (AJUSTA EL 7 si tienes menos muestras)
min_samples <- min(7, floor(ncol(dge) * 0.5))  # 50% de muestras o máximo 7
cat("Filtrando genes con CPM > 0.5 en al menos", min_samples, "muestras\n")

to_keep <- rowSums(cpm(dge) > 0.5) >= min_samples
dge <- dge[to_keep, keep.lib.sizes = FALSE]

print(paste("Genes después de filtrar:", nrow(dge)))

# 6. Normalizar
dge <- calcNormFactors(dge)

#####################################
# ANÁLISIS DE EXPRESIÓN DIFERENCIAL #
#####################################

design<-model.matrix(~colData(se)$condition)
dge<-estimateGLMCommonDisp(dge,design)
fit=glmFit(dge,design=design)
lrt<-glmLRT(fit, coef = 2)
tt<-topTags(lrt, n=Inf, sort.by = "none", adjust.method = "BH")$table

###############################
# REPRESENTACIÓN DE LOS DATOS #
###############################

library(ggplot2)
library(ggrepel) 
library(plotly)
library(dplyr)

res_interactivo <- as.data.frame(tt) %>%
  rownames_to_column("locus_tag") %>%
  left_join(as.data.frame(rowData(se)), by = "locus_tag") %>%
  mutate(
    logP = -log10(PValue),
    diff_state = case_when(
      logFC > 1 & FDR < 0.05 ~ "UP",
      logFC < -1 & FDR < 0.05 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

# 2. Crear el Volcano con Plotly
volcano_plot <- plot_ly(
  data = res_interactivo,
  x = ~logFC,
  y = ~logP,
  type = "scatter",
  mode = "markers",
  color = ~diff_state,
  colors = c("DOWN" = "#1f77b4", "NS" = "#bdbdbd", "UP" = "#d62728"),
  text = ~paste("<b>Gen:</b>", Name, 
                "<br><b>Locus:</b>", locus_tag,
                "<br><b>Producto:</b>", product,
                "<br><b>Taxon:</b>", taxon_id,
                "<br><b>LogFC:</b>", round(logFC, 2),
                "<br><b>FDR:</b>", format.pval(FDR, digits = 3)),
  hoverinfo = "text",
  marker = list(opacity = 0.6, size = 7)
) %>%
  layout(
    title = "Volcano Plot Interactivo: Efecto del Probiótico",
    xaxis = list(title = "Log2 Fold Change"),
    yaxis = list(title = "-log10 P-Value"),
    shapes = list(
      list(type = 'line', y0 = 0, y1 = 1, yref = 'paper', x0 = 1, x1 = 1, line = list(dash = 'dot', width = 1)),
      list(type = 'line', y0 = 0, y1 = 1, yref = 'paper', x0 = -1, x1 = -1, line = list(dash = 'dot', width = 1)),
      list(type = 'line', y0 = -log10(0.05), y1 = -log10(0.05), yref = 'y', x0 = 0, x1 = 1, xref = 'paper', line = list(dash = 'dot', width = 1))
    )
  )
volcano_plot



generate_gene_report_tami <- function(df, fdr_column = "FDR", 
                                      title="", fdr_threshold = 0.05) {
  require(dplyr)
  require(DT)
  require(tami)
  require(htmltools)
  
  report_df <- df %>%
    dplyr::filter(!!sym(fdr_column) < fdr_threshold) %>%
    dplyr::arrange(!!sym(fdr_column)) %>%
    dplyr::mutate(
      Taxonomy_Link = ifelse(!is.na(taxon_id) & taxon_id != "Unknown", 
                             paste0('<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=', taxon_id, '" target="_blank">NCBI Taxonomy</a>'), 
                             NA_character_),
      
      Entrez_Link = if ("entrez_id" %in% colnames(.)) {
        ifelse(!is.na(entrez_id) & entrez_id != "", tami::entrezid2url(entrez_id), NA_character_)
      },
      
      GO_Link = if ("go_terms" %in% colnames(.)) {
        go_id <- gsub(" .*", "", go_terms)
        ifelse(!is.na(go_id) & go_id != "", tami::go2url(go_id), NA_character_)
      }
    ) %>%
    dplyr::select(
      Locus = locus_tag, 
      Symbol = Name, 
      Producto = product, 
      logFC, 
      PValue, 
      fdr_column = !!sym(fdr_column),
      Taxonomy_Link, Entrez_Link, GO_Link
    )
  
    datatable_report <- DT::datatable(
    report_df,
    escape = FALSE,
    extensions = c('Buttons', 'Scroller'),
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf'),
      scrollX = TRUE,
      pageLength = 10,
      columnDefs = list(
        list(
          targets = which(grepl("_Link$", colnames(report_df))) - 1,
          render = JS(
            "
            function(data, type, row) {
              if (data === null || data === undefined || data === 'NA' || data === '') {
                return '';
              }
              return data;
            }
            "
          ),
          className = 'dt-center',
          defaultContent = ''  
        )
      ),
      language = list(
        emptyTable = "No hay datos disponibles" 
      )
    ),
    rownames = FALSE,
    class = 'display nowrap'
  )
  
  htmltools::tagList(
    htmltools::tags$h3(title),
    datatable_report
  )
}

tabla_widget <- generate_gene_report_tami(
  res_interactivo, 
  title = "Detalle de Genes Diferenciales (FDR < 0.05)",
  fdr_threshold = 0.05
)


####################################
#     ENRIQUECIMIENTO FUNCIONAL    #
# ANÁLISIS DE COLECCIONES DE GENES #
####################################

# Universo de genes y set

set0 <- res_interactivo %>%
  filter(FDR < 0.05) %>%
  pull(locus_tag)

universo_genes <- res_interactivo$locus_tag

# Creación de las colecciones GO y taxonómicas

gene_go_map <- as.data.frame(rowData(se)) %>%
  filter(locus_tag %in% universo_genes) %>%
  select(locus_tag, go_terms) %>%
  separate_rows(go_terms, sep = ";|\\,") %>% 
  mutate(go_terms = trimws(go_terms)) %>%
  filter(go_terms != "" & !is.na(go_terms))

# Convertir a lista (GSC)
gsc_go <- split(gene_go_map$locus_tag, gene_go_map$go_terms)

gsc_taxon <- split(as.data.frame(rowData(se))$locus_tag, 
                   as.data.frame(rowData(se))$taxon_id)

ora_go_raw <- tami::ora(set0 = set0, gsc = gsc_go)

# ORA Taxonómico
ora_tax_raw <- tami::ora(set0 = set0, gsc = gsc_taxon)
reporte_tax <- generate_ora_report_meta(ora_tax_raw, 
                                        title = "Análisis Taxonómico", 
                                        analysis_type = "Taxonomy", 
                                        dict = dict_tax)

reporte_go <- generate_ora_report_meta(ora_go_raw, 
                                       title = "Análisis Funcional", 
                                       analysis_type = "GO", 
                                       dict = dict_go)
final_html <- htmltools::tagList(
  htmltools::tags$h1("Análisis de Sobrerepresentación (ORA)"),
  reporte_go,
  htmltools::tags$hr(),
  reporte_tax
)

final_html

# VISUALIZACIÓN

crear_dotplot_interactivo <- function(df_ora, 
                                      dict = NULL, 
                                      titulo = "Análisis de Sobrerepresentación", 
                                      n_top = 20) {
  require(ggplot2)
  require(plotly)
  require(dplyr)
  
  top_terms <- df_ora %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "id") %>%
    filter(rawp < 0.05 & OR > 0) %>%
    mutate(
      logP = -log10(rawp),
      nombre_limpio = if(!is.null(dict)) {
        res <- dict[as.character(id)]
        ifelse(is.na(res), id, res)
      } else {
        sub("GO:\\d+_", "", id) 
      },
      name_short = ifelse(nchar(nombre_limpio) > 40,
                          paste0(substr(nombre_limpio, 1, 37), "..."),
                          nombre_limpio),
      tooltip_text = paste0("Categoría: ", nombre_limpio,
                            "<br>ID: ", id,
                            "<br>p-valor: ", signif(rawp, 3),
                            "<br>Odds Ratio: ", signif(OR, 3))
    ) %>%
    arrange(rawp) %>%
    slice_head(n = n_top)
  
  p <- ggplot(top_terms, aes(x = logP, y = reorder(name_short, logP))) +
    geom_point(aes(size = OR, color = OR, text = tooltip_text)) +
    scale_color_gradient(low = "skyblue", high = "darkred") +
    labs(
      title = titulo,
      x = "-log10(p-valor)",
      y = NULL,
      size = "Odds Ratio",
      color = "Odds Ratio"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold")
    )
  
  ggplotly(p, tooltip = "text")
}

plot_tax <- crear_dotplot_interactivo(ora_tax_raw, 
                                      dict = dict_tax, 
                                      titulo = "Sobrerrepresentación Taxonómica (Especies)",
                                      n_top = 15)

plot_go <- crear_dotplot_interactivo(ora_go_raw, 
                                     dict = dict_go, 
                                     titulo = "Sobrerrepresentación Funcional (GO)",
                                     n_top = 15)

plot_tax
plot_go

