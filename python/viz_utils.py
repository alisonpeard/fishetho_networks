import numpy as np
import matplotlib.pyplot as plt

def invert_dict(dic):
    return {value: key for key, value in dic.items()}

def make_heatmap(df, row1, row2, possible_values, **kwargs):
    """Make a heatmap for value counts of each ranking. 
    
    Clearer than scatter plots for ordinal data.
    Colours indicate frequency a value is recorded."""
    
    n = len(possible_values)
    pv_dict = {index: value for index, value in zip(range(n), possible_values)}
    pv_dict_r = {value: key for key, value in pv_dict.items()}
    
    value_counts = np.zeros((n, n))
    
    for ix, row in df[[row1, row2]].iterrows():
        if (not np.isnan(row[row1])) and (not np.isnan(row[row2])):
            
            idx1 = pv_dict_r[int(row[row1])]
            idx2 = pv_dict_r[int(row[row2])]
            value_counts[idx1, idx2] += 1
        
    return value_counts, pv_dict

def plot_heatmap(measure, col1, col2, values, ax, title, criteria_dict, lims=[0, 20], cbar=True, cmap='viridis'):
    criteria_dict_r = invert_dict(criteria_dict)
    heatmap, pv_dict = make_heatmap(measure, criteria_dict_r[col1], criteria_dict_r[col2], values)
    im = ax.imshow(heatmap, vmin=lims[0], vmax=lims[1], cmap=cmap)
    ax.set_xlabel(col1)
    ax.set_ylabel(col2)
    ax.set_xticks([*pv_dict.keys()])
    ax.set_yticks([*pv_dict.keys()])
    ax.set_xticklabels([*pv_dict.values()])
    ax.set_yticklabels([*pv_dict.values()])
    ax.invert_yaxis()
    ax.set_title(f'Heatmap of value-pairs\nfor {title}')
    if cbar:
        plt.colorbar(im);
    return im


def make_multilevel(likelihood, potential):
    """Combine Likelihood and Potential into a multilevel dataframe."""
    likelihood_t = likelihood.transpose()
    likelihood_t['level'] = ['likelihood'] * 10

    potential_t = potential.transpose()
    potential_t['level'] = ['potential'] * 10
    
    both = pd.concat([likelihood_t, potential_t])
    both = both.groupby(['Criterion', 'level']).agg(lambda x: x)
    both = both.transpose()
    
    return both