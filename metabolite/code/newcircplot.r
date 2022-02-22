#!/usr/bin/Rscript
require(RColorBrewer)
require(circlize)
require(ggplot2)

#circosInput <- read.csv('fsnewedges.csv')
circosInput <- read.csv('metabedges.csv')
#changes <- read.csv('nnnnodes.csv', row.names=1)
elements = (unique(c(circosInput$from, circosInput$to)))

set.seed(1)

gridcol = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
grid.col = gridcol(changes$X2)
names(grid.col) = rownames(changes)

col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
#col = col_fun(circosInput$value)
col = col_fun(circosInput$X0)
names(col) = rownames(circosInput)

#svg('../results/nnsmetab.svg',width=12,height=12)
#svg('../results/fnsmetab.svg',width=12,height=12)
svg('../results/fstmp.svg',width=12,height=12)

circos.par(circle.margin=c(1,1,1,1))
chordDiagram(circosInput,
	col = col,
	#grid.col=grid.col,
	annotationTrack = "grid",
	preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(circosInput))))))

circos.track(track.index = 1, panel.fun = function(x, y) {
	circos.text(CELL_META$xcenter,
		CELL_META$ylim[1],
		CELL_META$sector.index,
		facing = "clockwise",
		niceFacing = TRUE,
		adj = c(0, 0.5))
		},
	bg.border = NA)

circos.clear() 
dev.off()
