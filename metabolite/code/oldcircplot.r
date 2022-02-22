#!/usr/bin/Rscript
require(RColorBrewer)
require(circlize)

#circosInput <- read.csv('../results/gutmetab.csv')
circosInput <- read.csv('../results/oralmetab.csv')
#changes <- read.csv('nnodes.csv', row.names=1)
elements = (unique(c(circosInput$from, circosInput$to)))

set.seed(1)

#gridcol = colorRamp2(range(changes$WEEK04),
	#c("#0000FF", "#FF0000"),
	#transparency = 0.5)

col_fun = colorRamp2(range(circosInput$value),
	c("#0000FF", "#FF0000"),
	transparency = 0.5)

#svg('../results/gmetab.svg',width=12,height=12)
svg('../results/ometab.svg',width=12,height=12)

circos.par(circle.margin=c(1,1,1,1))
#chordDiagram(circosInput, col=col_fun, grid.col = gridcol)
#chordDiagram(circosInput, col=col_fun)
chordDiagram(circosInput,
	col = col_fun,
	annotationTrack = "grid",
	preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(circosInput))))))
#circos.trackHist(changes$X, changes$WEEK04)
#othercol = structure(rep("blue", length(changes$X)), names = changes$X)
#grid.col = c("4" = "red", "3" = "green", "5" = "yellow", othercol)
#col_fun1 = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
#circos.heatmap(changes, col=col_fun1,dend.side = "inside",rownames.side = "outside")
#circos.barplot(changes$X, changes$WEEK04)
#circos.heatmap(changes)
#col_fun1 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
#circos.heatmap(changes, split=changes$X)
#circos.track(ylim = c(-1, 1), panel.fun = function(x, y) {
#    value = changes$WEEK04
#    pos = changes$X
#    circos.barplot(value, pos=pos)
#})
#circos.track(ylim = c(-1, 1), panel.fun = function(x, y) {
#	x = changes$X,
#	y = changes$WEEK04,
#	panel.fun = function(x,y) {
#		circos.barplot(
#			value = y,
#			pos = x+0.3)}

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
