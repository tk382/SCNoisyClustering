myheatmap = function(mat, min, max, title){

  palette.gr.marray <- colorRampPalette(c("darkblue", "blue", "white", "red", "brown"))(30)
  heatmap.2(t(as.matrix(mat)),
            trace = "none",
            col = palette.gr.marray,
            Colv = F,
            Rowv = F, labRow = NA, labCol = NA,
            dendrogram = "none",
            key = T,
            rowsep = rowseps,
            breaks = seq(min, max,length=31),
            cexRow = 1,
            symbreaks = T,
            main=title)
}

