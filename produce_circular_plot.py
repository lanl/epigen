# produce_circular_plot.py
#
# This file produces the circular dendrogram and legends from
# the tree written in the newick format (written in get_newick_string.py)
# This is FIGURE 2 in the manuscript.
#
# Input file needed: chr6_newick.txt
#   This file is where the dendrogram for chr 6 is written in the newick format,
#   where labels are epigenetic modifiers.
#
# This file requires ete3, available from https://etetoolkit.org/
#
# To change the chromosome to plot,
# change the variable `newick_filename`
#
# This file should run in less than a minute.
#
# Note: We experienced difficulties producing properly rendered plots
# using macOS; We recommend linux for ete3.
#
# This file will one PDF file:
#     ./results38/png_tree.pdf
# This is the main plot.
#
# This file will output two EPS files:
#     ./results38/legend_solid_lines.eps
#     ./results38/legend.eps
# These are legends for the main plot.
#
#
import random
import re
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, TreeFace, PieChartFace, TextFace, NodeStyle, CircleFace
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import Counter
import itertools


# Rename labels
def replace_tip_label(match):
    label = int(match.group(1))
    return name_map.get(label, label)


# Generate random colors
def get_random_color():
    return '#%02X%02X%02X' % (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))


def layout(node):
    if node.is_leaf():
        label = node.name.split('_')[-1]
        name_face = TextFace(label, fgcolor=label_to_color.get(label, "black"))#, fsize=15)  
        node.add_face(name_face, column=0, position="branch-right")
        return

    ancestor_with_piechart = False
    for ancestor in node.iter_ancestors():
        if ancestor.dist > N_high_importance_color and len(ancestor.get_leaf_names()) >= N_samples:
            ancestor_with_piechart = True
            break

    if not ancestor_with_piechart and node.dist >= N_high_importance_color and len(node.get_leaf_names()) >= N_samples:
        label_counter = Counter([label.split('_')[-1] for label in node.get_leaf_names()])
        
        # Sort label_counter by most common labels
        sorted_items = label_counter.most_common()
        
        # Prepare the values and labels for the pie chart
        if len(sorted_items) > top_N_pies:
            top_5_values = [item[1] for item in sorted_items[:top_N_pies]]
            other_values = sum([item[1] for item in sorted_items[top_N_pies:]])
            values = top_5_values + [other_values]
            labels = [item[0] for item in sorted_items[:top_N_pies]] + ["Other"]
            
            # Ensure the "Other" label has a color assigned
            if "Other" not in label_to_color:
                label_to_color["Other"] = '#000000'
        else:
            values = [item[1] for item in sorted_items]
            labels = [item[0] for item in sorted_items]
        #print(labels)
        labels_visual_list.append(labels)
        # Extract colors based on labels
        colors = [label_to_color[label] for label in labels]
        
        # Normalize to percentages
        total = sum(values)
        percents = [v / total * 100 for v in values]

        piechart = PieChartFace(percents, colors=colors, width=1000, height=1000)
        node.add_face(piechart, column=0, position="float")
        
        # Add label texts below the pie chart
        for label in labels:
            label_face = TextFace(label)#, fsize=40)
            node.add_face(label_face, column=1, position="branch-bottom")


        # Color the branches of the node with the pie chart and all its descendant nodes
        style = NodeStyle()
        style["fgcolor"] = "#0f0f0f"
        style["size"] = 0
        style["vt_line_color"] = "#ff0000"
        style["hz_line_color"] = "#ff0000"
        style["vt_line_width"] = 25
        style["hz_line_width"] = 25
        style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
        style["hz_line_type"] = 0
        
        for desc in node.traverse():
            desc.set_style(style)
    
    elif not ancestor_with_piechart and node.dist >= N_mid_importance_color and node.dist < N_high_importance_color and len(node.get_leaf_names()) >= N_samples:
        green_style = NodeStyle()
        green_style["fgcolor"] = "#009933"
        green_style["vt_line_color"] = "#009933"
        green_style["hz_line_color"] = "#009933"
        green_style["vt_line_width"] = 25
        green_style["hz_line_width"] = 25
        green_style["vt_line_type"] = 0 
        green_style["hz_line_type"] = 0  
        for desc in node.traverse():
            desc.set_style(green_style)
            
    elif not ancestor_with_piechart and node.dist >= N_low_importance_color and node.dist < N_mid_importance_color and len(node.get_leaf_names()) >= N_samples:
        blue_style = NodeStyle()
        blue_style["fgcolor"] = "#0000cc"
        blue_style["vt_line_color"] = "#0000cc"
        blue_style["hz_line_color"] = "#0000cc"
        blue_style["vt_line_width"] = 25
        blue_style["hz_line_width"] = 25
        blue_style["vt_line_type"] = 0 
        blue_style["hz_line_type"] = 0 
        for desc in node.traverse():
            desc.set_style(blue_style)
    else:
        dist_face = TextFace(f"Dist: {node.dist:.2f}")  # Format to 2 decimal places
        node.add_face(dist_face, column=0, position="branch-right")

if __name__ == "__main__":

    random.seed(123)
    newick_filename = './data/chr6_newick.txt'

    # Read metadata
    df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
    df38 = pd.DataFrame(df38)
    df38 = df38.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]
    df38

    # Read chromosome 6 dendrogram saved in Newick format
    # #where labels are integers and branch lengths reflect number of important genes
    # Pull Target names from df38 and rename labels as number_Target
    names = ["{}_{}".format(i, s) for i, s in enumerate(df38['Target'], 0) if i < df38.shape[0]]
    name_map = {i: f"{label}" for i, label in enumerate(names)}
    s = open(newick_filename, 'r').read()

    # Substitute number leaf names with string names using regex
    modified_newick = re.sub(r'(?<=\()(\d+)(?=:)', replace_tip_label, s)
    modified_newick = re.sub(r'(?<=,)(\d+)(?=:)', replace_tip_label, modified_newick)
    print(modified_newick)

    # Read the tree from Newick string and set parameters for visualization (can be changed, see legends below)
    s = modified_newick
    t = Tree(s)
    top_N_pies = 3
    N_samples = 10
    N_low_importance_color = 4
    N_mid_importance_color = 10
    N_high_importance_color = 30

    labels_visual_list = []
    # Extract all unique labels from the tree
    all_labels = set([leaf.name.split('_')[-1] for leaf in t.iter_leaves()])

    # Define available colors
    available_colors = ["red", "green", "blue", "yellow", "purple", "orange", "pink", "cyan", "magenta", "brown"]
    label_to_color = {}  # Will store mapping between labels and colors

    # Populate label_to_color dictionary
    for label in all_labels:
        if label not in label_to_color:
            if available_colors:
                label_to_color[label] = available_colors.pop(0)
            else:
                label_to_color[label] = get_random_color()


    ts = TreeStyle()
    ts.layout_fn = layout
    ts.mode = "c"
    ts.show_leaf_name = False  # Turn off default leaf name display since we're adding custom text faces
    ts.scale = 1000  # Adjust this value to find a good balance
    ts.branch_vertical_margin = 20  # Adjust to space out nodes if necessary
    ts.force_topology = True

    # Add legend to the tree style
    for label, color in label_to_color.items():
        ts.legend.add_face(CircleFace(10, color), column=0)
        ts.legend.add_face(TextFace(" " + label, fgcolor=color), column=1)

    # Save pdf figure (circular plot)
    t.show(tree_style=ts)
    output_png = "./results38/png_tree.pdf"
    t.render(output_png, tree_style=ts)


    # Flatten the list and count frequencies of each label
    flattened_labels = list(itertools.chain(*labels_visual_list))
    label_frequencies = Counter(flattened_labels)

    # Separate single-element and multi-element lists
    single_elements = [item[0] for item in labels_visual_list if len(item) == 1]
    multi_elements = list(set(flattened_labels) - set(single_elements) - {'Other'})

    # Sort single elements by frequency then alphabetically
    sorted_single_elements = sorted(single_elements, key=lambda x: (-label_frequencies[x], x))

    # Sort multi elements alphabetically
    sorted_multi_elements = sorted(multi_elements)

    # Combine lists and add 'Other' at the end
    sorted_labels = sorted_single_elements + sorted_multi_elements
    sorted_labels.append('Other')

    # Get unique labels
    unique_sorted_labels = list(dict.fromkeys(sorted_labels))

    # Create the legend for Target (modifications), corresponds to pie chart colors
    # Assuming 'unique_sorted_labels' and 'label_to_color' are already defined
    # Ensure that the labels are ordered as per 'unique_sorted_labels'
    labels = unique_sorted_labels
    colors = [label_to_color[label] for label in labels]

    fig, ax = plt.subplots(figsize=(6, len(labels) * 0.25))
    ax.axis('off')

    legend_elements = [Line2D([0], [0], marker='o', color='w', label=label, markersize=10, markerfacecolor=color)
                       for label, color in zip(labels, colors)]

    columns = min(3, len(labels))  # Use up to 3 columns, but not more than the number of labels
    ax.legend(handles=legend_elements, loc='center', ncol=columns, title="Epigenetic Modifications", fontsize='small')

    plt.tight_layout()
    plt.savefig("./results38/legend.eps", format='eps')
    plt.show()


    # Create the legend of # important genes (colored lines)
    # Adjusting the order of the legend entries to match the requested sequence: red, green, blue, black
    # Redefining the legend labels and colors in the desired order
    legend_labels = [
        "# of important genes ≥ 30",     # Red
        "10 ≤ # of important genes < 30",  # Green
        "4 ≤ # of important genes < 10",   # Blue
        "# of important genes < 4"        # Black
    ]
    legend_colors = ["#ff0000", "#009933", "#0000cc", "#000000"] # Red, Green, Blue, Black

    fig, ax = plt.subplots(figsize=(6, len(legend_labels) * 0.25))  # Adjust the figsize as needed
    ax.axis('off')

    legend_elements = [Line2D([0], [0], color=color, label=label, linewidth=4)
                       for label, color in zip(legend_labels, legend_colors)]

    columns = min(1, len(legend_labels))

    ax.legend(handles=legend_elements, loc='center', ncol=columns, title="Number of Important Known Genes in Clusters", fontsize='small')
    plt.tight_layout()
    legend_filename_png = "./results38/legend_solid_lines.eps"
    plt.savefig(legend_filename_png, format='eps')
    plt.show()
