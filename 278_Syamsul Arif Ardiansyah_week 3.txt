# ======================
#SCRIPT ini modifikasi dari script pada modul week 3.
#1. menambahkan script collapsed untuk menghindari gen yang memiliki duplikasi.
#2. menambahkan script untuk output excel.
#3. menambahkan script untuk anaisis GO dan KEGG termasuk visualisasi.
#4. install beberapa aplikasi yang mendukung analisis.
# ======================

# ======================
#MODUL: ANALISIS EKSPRESI GEN KANKER PARU  
# ======================

#Dataset: GSE10072 (Lung Adenocarcinoma vs Normal) 
#Platform: Microarray (Affymetrix GPL96) 
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG)  

# ======================
#PART A. PENGANAR KONSEP  
# ======================
#Analisis ekspresi gen bertujuan untuk membandingkan tingkat ekspresigen  
#antara dua kondisi biologis (misalnya kanker vs normal)  
#Pada modul ini kita menggunakan pendekatan statistik limma (LinearModels 
#for Microarray Data), yang merupakan standar emas untuk data microarray. 

# ======================
#PART B. PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE)  
# ======================
#Apa itu package?  
#Package adalah kumpulan fungsi siap pakai di R 
#Bioinformatika di R sangat bergantung pada package dari CRAN dan Bioconductor  

# ======================
#1. Install BiocManager (manajer paket Bioconductor)  
# ======================

#IF adalah struktur logika : “jika kondisi terpenuhi, lakukan aksi” 
if (!require("BiocManager", quietly = TRUE))  { 
  install.packages("BiocManager")  } 

# ======================
# 2. Install paket Bioconductor (GEOquery & limma)  
# ======================

#GEOquery: mengambil data dari database GEO  
#limma: analisis statistik ekspresi gen  
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)  
#Install annotation package sesuai platform 
#GPL96 = Affymetrix Human Genome U133A  
BiocManager::install("hgu133a.db", ask = FALSE, update = FALSE) 


# ======================
#3. Install paket CRAN untuk visualisasi dan manipulasi data 
# ======================
#phetmap: heatmap ekspresi gen  
#ggplot2: grafik (volcano plot) 
#dplyr: manipulasi tabel data  
install.packages(c("pheatmap", "ggplot2", "dplyr")) 
#umap: grafik (plot UMAP)  
if (!requireNamespace("umap", quietly = TRUE)) { 
  install.packages("umap") } 

# ======================
#4. Library  
# ======================

#library() digunakan agar fungsi di dalam package bisa digunakan  
library(GEOquery) 
library(limma) 
library(pheatmap) 
library(ggplot2) 
library(dplyr) 
library(hgu133a.db) 
library(AnnotationDbi) 
library(umap) 
library(gridExtra)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GEOquery)
library(Biobase)

# ======================
#PART C. PENGAMBILAN DATA DARI GEO  
# ======================
#GEO (Gene Expression Omnibus) adalah database publik milik NCBI 
#getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO 
#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet 
#AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh 
# Load package yang diperlukan
library(GEOquery)
library(Biobase)
gset_list <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)
#ExpressionSet berisi: 
# - exprs() : matriks ekspresi gen 
# - pData() : metadata sampel 
# - fData() : metadata fitur (probe / gen) 

# ======================
#PART D. PRE-PROCESSING DATA EKSPRESI  
# ======================

# exprs(): mengambil matriks ekspresi gen 
# Baris  = probe/gen 
# Kolom  = sampel 
gset_list <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset_list[[1]]
ex <- exprs(gset)

#Mengapa perlu log2 transformasi? 
#Data microarray mentah memiliki rentang nilai sangat besar. 
#Log2 digunakan untuk: 
#1. Menstabilkan varians 
#2. Mendekati asumsi model linear 
#3. Memudahkan interpretasi log fold change 
#quantile(): menghitung nilai kuantil (persentil) 
#as.numeric(): mengubah hasil quantile (yang berupa named vector) 
#menjadi vektor numerik biasa agar mudah dibandingkan 
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)) 
#LogTransform adalah variabel logika (TRUE / FALSE) 
#Operator logika: 
#>  : lebih besar dari 
#| | : OR (atau) 
#&& : AND (dan) 
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) 
#IF statement: 
#Jika LogTransform = TRUE, maka lakukan log2 
if (LogTransform) { 
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA 
  ex[ex <= 0] <- NA 
  ex <- log2(ex) 
} 

#Melakukan Collapsed terhadap probe yang mengkode gen yang sama
#Hasil Collapsed mengeluarkan data gen tana duplikat
# Ambil probe annotation
fdata <- fData(gset)

# Pilih kolom gene symbol yang benar. setelah run script diatas akan muncul data kolom
#pilih nama gene symbol yang sesuai
gene_symbols <- fdata$`Gene symbol`

# Hapus probe tanpa simbol gen
keep <- !is.na(gene_symbols) & gene_symbols != ""
ex_valid <- ex[keep, ]
gene_symbols_valid <- gene_symbols[keep]

# Pastikan matrix numeric
ex_valid <- as.matrix(ex_valid)

# Collapse probe → gene (rata-rata), Gen yang memiliki duplikasi akan diambil rata-rata nya 
ex_gene <- avereps(ex_valid, ID = gene_symbols_valid)

# Cek dimensi dan duplikasi
cat("Dimensi sebelum collapse:", dim(ex), "\n")
cat("Dimensi setelah collapse:", dim(ex_gene), "\n")
any(duplicated(rownames(ex_gene)))

#Hasil harus FALSE menandakan tidak ada gen yang duplikat

# ex_valid sudah numeric matrix, gene_symbols_valid sudah bersih
ex_gene <- avereps(ex_valid, ID = gene_symbols_valid)

# Gunakan SYMBOL sebagai rownames
rownames(ex_gene) <- unique(gene_symbols_valid)


# ======================
#PART E. DEFINISI KELOMPOK SAMPEL 
# ======================

#pData(): metadata sampel 
#source_name_ch1 berisi informasi kondisi biologis sampel 
group_info <- pData(gset)[["source_name_ch1"]] 
#make.names(): mengubah teks menjadi format valid untuk R 
groups <- make.names(group_info) 
#factor(): 
#Mengubah data kategorik menjadi faktor 
#Faktor sangat penting untuk analisis statistik di R 
gset$group <- factor(groups) 
#levels(): melihat kategori unik dalam faktor 
nama_grup <- levels(gset$group) 
print(nama_grup) 

# ======================
#PART F. DESIGN MATRIX (KERANGKA STATISTIK)  
# ======================

#model.matrix(): 
#Membuat matriks desain untuk model linear 
#~0 berarti TANPA intercept (best practice limma) 
design <- model.matrix(~0 + gset$group) 
#colnames(): memberi nama kolom agar mudah dibaca 
colnames(design) <- levels(gset$group) 

#Menentukan perbandingan biologis 
grup_kanker <- nama_grup[1] 
grup_normal <- nama_grup[2] 
contrast_formula <- paste(grup_kanker, "-", grup_normal) 
print(paste("Kontras yang dianalisis:", contrast_formula))

# ======================
#PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA) 
# ======================

# limma dengan ex_gene
fit <- lmFit(ex_gene, design)
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Ambil hasil DEG langsung gene-level
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

# rownames sudah gene-level
topTableResults$SYMBOL <- rownames(topTableResults)

# ======================
#PART H. BOXPLOT DISTRIBUSI NILAI EKSPRESI  
# ======================

# Boxplot digunakan untuk: 
#- Mengecek distribusi nilai ekspresi antar sampel 
#- Melihat apakah ada batch effect 
#- Mengevaluasi apakah normalisasi/log-transform sudah wajar 
#Set warna berdasarkan grup 
group_colors <- as.numeric(gset$group) 
boxplot( 
  ex, 
  col = group_colors, 
  las = 2, 
  outline = FALSE, 
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel", 
  ylab = "Expression Value (log2)" ) 
legend( 
  "topright", 
  legend = levels(gset$group), 
  fill = unique(group_colors), 
  cex = 0.8 ) 

# ======================
#PART I. DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
# ======================

#Density plot menunjukkan sebaran global nilai ekspresi gen 
#Digunakan untuk: 
#- Mengecek efek log-transform 
#- Membandingkan distribusi antar grup 
#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex_gene),
  Group = rep(gset$group, each = nrow(ex_gene)))
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen (Gene Level)",
    x = "Expression Value (log2)",
    y = "Density")

# ======================
#PART J. UMAP (VISUALISASI DIMENSI RENDAH) 
# ======================

#UMAP digunakan untuk: 
#- Mereduksi ribuan gen menjadi 2 dimensi 
#- Melihat pemisahan sampel secara global 
#- Alternatif PCA (lebih sensitif ke struktur lokal) 
#Transpose matriks ekspresi: 
#UMAP bekerja pada OBSERVATION = sampel 
umap_input <- t(ex_gene)
umap_result <- umap(umap_input)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel (Gene Level, No Duplicate)",
    x = "UMAP 1",
    y = "UMAP 2")

# ======================
#PART K. VISUALISASI VOLCANO PLOT  
# ======================

#Volcano plot menggabungkan: 
#- Log fold change (efek biologis) 
#- Signifikansi statistik
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

volcano_data$status <- "NO"

volcano_data$status[
  volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01
] <- "UP"

volcano_data$status[
  volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01
] <- "DOWN"

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Kanker Paru-Paru (No Duplicate Gene)")


# ======================
#PART L. VISUALISASI HEATMAP  
# ======================

#Heatmap digunakan untuk melihat pola ekspresi gen 
#antar sampel berdasarkan gen-gen paling signifikan 
# Urutkan berdasarkan signifikansi
topTableResults <- topTableResults[order(topTableResults$adj.P.Val), ]
# Ambil top 50 gen unik
top50 <- head(topTableResults, 50)
# Ambil expression matrix gene level
mat_heatmap <- ex_gene[top50$SYMBOL, ]
# Gunakan SYMBOL sebagai label
rownames(mat_heatmap) <- top50$SYMBOL
# Hapus NA
mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]
# Hapus varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]
# Annotation
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)
# Plot heatmap
pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes (No Duplicate)"
)

# ======================
#PART M. ANALISIS ENRICHMENT, 
# ======================

#- Gene Ontology 
#- (GO) - KEGG Pathway
#instal terlebi dahulu package ini
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

#data yang diambil dari 50 gen dengan tingkat signifikan tertiggi
top50 <- topTableResults[1:50, ]
gene_symbols <- top50$SYMBOL
gene_df <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)

gene_list <- gene_df$ENTREZID
length(gene_list)

# ======================
#GO
# ======================

ego <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE)

#visualisasi GO
dotplot(ego, showCategory = 10) +
  ggtitle("GO Biological Process Enrichment (Top 50 Genes)")
barplot(ego, showCategory = 10)
cnetplot(ego, showCategory = 5)

# ======================
#KEGG
# ======================

ekegg <- enrichKEGG(
  gene         = gene_list,
  organism     = "hsa",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05
)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

#Visualisasi KEGG
dotplot(ekegg, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment (Top 50 Genes)")


#data top Go dan KEGG
head(as.data.frame(ego))
head(as.data.frame(ekegg))


# ======================
# OUTPUT EXCELL
# ======================
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

# Load library
library(openxlsx)

# ======================
# Top50 DEG
# ======================
# Urutkan berdasarkan adj.P.Val
topTableResults <- topTableResults[order(topTableResults$adj.P.Val), ]
top50_final <- head(topTableResults, 50)

# Pilih kolom yang tersedia di gene-level
top50_export <- top50_final[, c("SYMBOL", "logFC", "AveExpr", "P.Value", "adj.P.Val")]

# ======================
# GO & KEGG enrichment
# ======================
go_table <- as.data.frame(ego)
kegg_table <- as.data.frame(ekegg)

# ======================
# Buat workbook baru
# ======================
wb <- createWorkbook()

# Tambahkan sheet Top50 DEG
addWorksheet(wb, "Top50_DEG")
writeData(wb, sheet = "Top50_DEG", top50_export)

# Tambahkan sheet GO enrichment
addWorksheet(wb, "GO_Enrichment")
writeData(wb, sheet = "GO_Enrichment", go_table)

# Tambahkan sheet KEGG enrichment
addWorksheet(wb, "KEGG_Enrichment")
writeData(wb, sheet = "KEGG_Enrichment", kegg_table)

# ======================
# Save workbook
# ======================
saveWorkbook(wb, file = "GSE10072_GeneLevel_Results.xlsx", overwrite = TRUE)

cat("Export selesai: GSE10072_GeneLevel_Results.xlsx\n")

# ======================
# UPREGULASI & DOWNREGULASI + Export Excel
# ======================

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(enrichplot)

# ----------------------
# Pisahkan gen UP dan DOWN
# ----------------------
up_genes_df <- topTableResults %>%
  filter(logFC > 1 & adj.P.Val < 0.05) %>%
  select(SYMBOL, logFC, adj.P.Val)

down_genes_df <- topTableResults %>%
  filter(logFC < -1 & adj.P.Val < 0.05) %>%
  select(SYMBOL, logFC, adj.P.Val)

cat("Jumlah UP genes:", nrow(up_genes_df), "\n")
cat("Jumlah DOWN genes:", nrow(down_genes_df), "\n")

# Cek duplikasi gen
any(duplicated(topTableResults$SYMBOL))

# ----------------------
# Konversi SYMBOL -> ENTREZID untuk GO & KEGG
# ----------------------
up_df <- bitr(
  up_genes_df$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_df <- bitr(
  down_genes_df$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_list <- up_df$ENTREZID
down_list <- down_df$ENTREZID

cat("UP genes mapped:", length(up_list), "\n")
cat("DOWN genes mapped:", length(down_list), "\n")

# ----------------------
# GO enrichment
# ----------------------
ego_up <- enrichGO(
  gene = up_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = down_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# ----------------------
# KEGG enrichment
# ----------------------
ekegg_up <- enrichKEGG(
  gene = up_list,
  organism = "hsa",
  pvalueCutoff = 0.05
)
ekegg_up <- setReadable(ekegg_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

ekegg_down <- enrichKEGG(
  gene = down_list,
  organism = "hsa",
  pvalueCutoff = 0.05
)
ekegg_down <- setReadable(ekegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# ----------------------
# Export Excel: Up/Down + logFC + adj.P.Val + GO/KEGG results
# ----------------------
wb <- createWorkbook()

# Sheet UP genes
addWorksheet(wb, "Upregulated_Genes")
writeData(wb, "Upregulated_Genes", up_genes_df)

# Sheet DOWN genes
addWorksheet(wb, "Downregulated_Genes")
writeData(wb, "Downregulated_Genes", down_genes_df)

# Sheet GO UP
addWorksheet(wb, "GO_UP")
writeData(wb, "GO_UP", as.data.frame(ego_up))

# Sheet GO DOWN
addWorksheet(wb, "GO_DOWN")
writeData(wb, "GO_DOWN", as.data.frame(ego_down))

# Sheet KEGG UP
addWorksheet(wb, "KEGG_UP")
writeData(wb, "KEGG_UP", as.data.frame(ekegg_up))

# Sheet KEGG DOWN
addWorksheet(wb, "KEGG_DOWN")
writeData(wb, "KEGG_DOWN", as.data.frame(ekegg_down))

# Simpan workbook
saveWorkbook(wb, "GSE10072_UP_DOWN_GO_KEGG.xlsx", overwrite = TRUE)
cat("Excel berhasil dibuat: GSE10072_UP_DOWN_GO_KEGG.xlsx\n")

# ----------------------
# Visualisasi GO & KEGG
# ----------------------
# GO UP
dotplot(ego_up, showCategory = 15) + ggtitle("GO Biological Process - UP Genes")
# GO DOWN
dotplot(ego_down, showCategory = 15) + ggtitle("GO Biological Process - DOWN Genes")
# KEGG UP
dotplot(ekegg_up, showCategory = 10) + ggtitle("KEGG Pathway - UP Genes")
# KEGG DOWN
dotplot(ekegg_down, showCategory = 10) + ggtitle("KEGG Pathway - DOWN Genes")
