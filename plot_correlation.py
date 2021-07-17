def main(df):
    def get_correlations(df):
        dfcols = pd.DataFrame(columns=df.columns)
        correlations = dfcols.transpose().join(dfcols, how="outer")
        pvalues = dfcols.transpose().join(dfcols, how="outer")
        for ix, r in enumerate(df.columns):
            for jx, c in enumerate(df.columns):
                spear = spearmanr(df[r], df[c])
                correlations[c][r] = spear[0]
                pvalues[c][r] = spear[1]
        return correlations.astype("float"), pvalues.astype("float")

# Calculate and format correlations
    correlations, uncorrected_p_values = get_correlations(numeric_data)
    fcorrelations = correlations.loc[
        msp_data._get_numeric_data().columns,
        metadata._get_numeric_data().columns]
    funcorrected_p_values = uncorrected_p_values.loc[
        msp_data._get_numeric_data().columns,
        metadata._get_numeric_data().columns]
    ffuncorrected_p_values = funcorrected_p_values.loc[
        (fcorrelations.sum(axis=1) != 0),
        (fcorrelations.sum(axis=0) != 0)]
    ffcorrelations = fcorrelations.loc[
        (fcorrelations.sum(axis=1) != 0),
        (fcorrelations.sum(axis=0) != 0)]
    shape = ffuncorrected_p_values.values.shape
    significant_matrix = multipletests(ffuncorrected_p_values.values.flatten())[0].reshape(shape)

# Plot
    g = sns.clustermap(
        ffcorrelations,
        cmap="vlag",
        vmin=-1,
        vmax=1, 
        yticklabels=True,
        xticklabels=True)

    for tick in g.ax_heatmap.get_yticklabels():
        tick.set_rotation(0)

    for i, ix in enumerate(g.dendrogram_row.reordered_ind):
        for j, jx in enumerate(g.dendrogram_col.reordered_ind):
            text = g.ax_heatmap.text(
                j + 0.5,
                i + 0.5,
                "*" if significant_matrix[ix, jx] else "",
                ha="center",
                va="center",
                color="black",
            )
            text.set_fontsize(8)
    return g
