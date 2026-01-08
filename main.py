import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# --- CONFIGURATION ---
I_CONC = 10.0  
COLORS = {'Control': '#1f77b4', 'Inhibitor': '#d62728'}

def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)

def run_ace_cohesive_project(df):
    results = {}
    
    # Define Layout: 2 rows, 2 columns
    fig = plt.figure(figsize=(12, 10))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 0.6])
    
    ax1 = fig.add_subplot(gs[0, 0]) # Top Left: MM Plot
    ax2 = fig.add_subplot(gs[0, 1]) # Top Right: LB Plot
    ax3 = fig.add_subplot(gs[1, :]) # Bottom: Table (span both columns)

    groups = sorted(df['Condition'].unique(), reverse=True) 

    for group in groups:
        color = COLORS.get(group)
        subset = df[df['Condition'] == group]
        S, V = subset['S_uM'].values, subset['Velocity'].values
        
        # 1. Fit & Statistics
        popt, _ = curve_fit(michaelis_menten, S, V, p0=[max(V), 25])
        vmax, km = popt
        r2 = r2_score(V, michaelis_menten(S, vmax, km))
        results[group] = {'Vmax': vmax, 'Km': km, 'R2': r2}

        # 2. Michaelis-Menten Plot
        s_plot = np.linspace(0, S.max() * 1.1, 100)
        ax1.plot(s_plot, michaelis_menten(s_plot, vmax, km), color=color, lw=3, label=f'{group} Fit')
        ax1.scatter(S, V, color=color, s=60, edgecolors='white', zorder=3, alpha=0.8)
        
        # 3. Lineweaver-Burk Plot
        inv_S, inv_V = 1/S, 1/V
        slope, intercept = np.polyfit(inv_S, inv_V, 1)
        x_ext = np.linspace(-1/km, max(inv_S) * 1.1, 100)
        ax2.plot(x_ext, slope * x_ext + intercept, color=color, lw=2.5, label=f'{group} Linear')
        ax2.scatter(inv_S, inv_V, color=color, s=50, edgecolors='white', zorder=5)

    # 4. Ki Calculation
    ki = I_CONC / ((results['Inhibitor']['Km'] / results['Control']['Km']) - 1)

    # --- STYLING ---
    ax1.set_title('Michaelis-Menten Saturation', fontsize=14, fontweight='bold', pad=10)
    ax1.set_xlabel('[Angiotensin I] (µM)', fontweight='bold')
    ax1.set_ylabel('Velocity (µmol/min)', fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(frameon=True, shadow=True)

    ax2.set_title('Lineweaver-Burk Diagnostic', fontsize=14, fontweight='bold', pad=10)
    ax2.set_xlabel('1/[S] (L/µmol)', fontweight='bold')
    ax2.set_ylabel('1/V (min/µmol)', fontweight='bold')
    ax2.axhline(0, color='black', lw=1.2); ax2.axvline(0, color='black', lw=1.2)
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(frameon=True, shadow=True)

    # --- TABLE (BOTTOM) ---
    ax3.axis('off')
    table_data = [
        ["Parameter", "Control Condition", "Inhibitor Condition"],
        ["Vmax (µmol/min)", f"{results['Control']['Vmax']:.2f}", f"{results['Inhibitor']['Vmax']:.2f}"],
        ["Km (µM)", f"{results['Control']['Km']:.2f}", f"{results['Inhibitor']['Km']:.2f}"],
        ["Goodness of Fit (R²)", f"{results['Control']['R2']:.4f}", f"{results['Inhibitor']['R2']:.4f}"],
        ["Inhibition Mechanism", "COMPETITIVE", f"Ki = {ki:.2f} µM"]
    ]
    
    the_table = ax3.table(cellText=table_data, loc='center', cellLoc='center', colWidths=[0.35, 0.25, 0.25])
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(12)
    the_table.scale(1.2, 4) 

    for (row, col), cell in the_table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold', color='white')
            cell.set_facecolor('#333333')
        if row == 4:
            cell.set_facecolor('#f2f2f2')
            cell.set_text_props(weight='bold')

    plt.suptitle(f"ACE Kinetic Analysis Report\n[Inhibitor Concentration] = {I_CONC} µM", 
                 fontsize=18, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    plt.show()

if __name__ == "__main__":
    data = {
        'S_uM': [5, 10, 25, 50, 100, 200] * 2,
        'Velocity': [18, 30, 52, 68, 80, 89, 7, 13, 28, 45, 62, 75],
        'Condition': ['Control']*6 + ['Inhibitor']*6
    }
    run_ace_cohesive_project(pd.DataFrame(data))