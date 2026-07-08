import os
from pathlib import Path
import numpy as np
import pandas as pd
from plotnine import *
import statsmodels.formula.api as smf

if __name__ == '__main__':
    # Data directory
    wk_dir = Path(f"experiments/metastasis/figures")
    if os.path.exists(wk_dir) == 0:
        os.makedirs(wk_dir, exist_ok=True)

    data_df = pd.read_csv(Path("experiments/metastasis/SuppTable3.csv"), index_col=0)

    n_df = data_df.groupby(['Tumor_ID', 'patient_ID']).size().reset_index(name='n')
    n_df.to_csv(Path("experiments/metastasis/n_df.csv"))

    plot = (
            ggplot(data_df, aes(x='date_len', y='TumorVol_num', color='Damian_annotation', group='Tumor_ID'))
            + geom_point(alpha=0.6, size=2)
            + geom_line(alpha=0.8, size=1)  # Connects points belonging to the exact same replicate
            + facet_wrap('~patient_ID', scales='free_x')  # Separate subplot per patient
            + labs(
        title="Individual PDX Replicate Growth Trajectories by Patient",
        x="Days",
        y="Cubic Root of Tumor Volume",
        color="Metastatic Status"
    )
            + theme_bw()
            + theme(
        figure_size=(12, 8),
        subplots_adjust={'hspace': 0.3, 'wspace': 0.2},
        strip_text=element_text(weight='bold')  # Makes the SA_ID labels stand out
    )
    )


    def compute_r2(data_df):
        model = smf.mixedlm("TumorVol_num ~ date_len * Damian_annotation",
                            data_df,
                            groups=data_df["Tumor_ID"],
                            re_formula="~date_len")
        model_result = model.fit()

        # Fixed effect variance
        fixed_pred = model_result.predict(data_df)
        var_fixed = np.var(fixed_pred)

        # Total Random effect variance (Intercept + Slope)
        random_pred = model_result.fittedvalues - fixed_pred
        var_random = np.var(random_pred)    # Random variance
        var_resid = model_result.scale      # Residual variance

        # R-squared calculations
        marginal_r2 = var_fixed / (var_fixed + var_random + var_resid)
        conditional_r2 = (var_fixed + var_random) / (var_fixed + var_random + var_resid)

        return marginal_r2, conditional_r2


    m_r2, c_r2 = compute_r2(data_df)
    print(f"Marginal R2: {m_r2:.3f}")
    print(f"Conditional R2: {c_r2:.3f}")
    plot.save(wk_dir / "pdx_growth_trajectories.svg", format="svg", verbose=False)
