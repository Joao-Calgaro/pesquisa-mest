from ete3 import Tree, TreeStyle, NodeStyle
import random

# Create random tree
t = Tree()
t.populate(20)

# ---------- TRANSMISSION RATES ----------
normal_rate = 1.0
green_rate = 10.0  # green subtree has 3x higher transmission rate


# ---------- GREEN SUBTREE ----------
internal_nodes = [
    n for n in t.traverse()
    if not n.is_leaf() and n != t
]

green_root = random.choice(internal_nodes)


# ---------- RANDOMIZE BRANCH LENGTHS ----------
# Branch length = transmission time.
# Higher transmission rate -> shorter average branch length.

for node in t.traverse():
    if node != t:
        if node in green_root.traverse():
            node.dist = random.expovariate(green_rate)
        else:
            node.dist = random.expovariate(normal_rate)


# ---------- RED STYLE ----------
red_style = NodeStyle()
red_style["fgcolor"] = "#B22222"
red_style["size"] = 7
red_style["hz_line_color"] = "#B22222"
red_style["vt_line_color"] = "#B22222"
red_style["hz_line_width"] = 2
red_style["vt_line_width"] = 2

for node in t.traverse():
    node.set_style(red_style)


# ---------- GREEN STYLE ----------
green_style = NodeStyle()
green_style["fgcolor"] = "#228B22"
green_style["size"] = 8
green_style["hz_line_color"] = "#4C794C"
green_style["vt_line_color"] = "#4C794C"
green_style["hz_line_width"] = 3
green_style["vt_line_width"] = 3

for node in green_root.traverse():
    node.set_style(green_style)


# ---------- TREE STYLE ----------
ts = TreeStyle()

ts.show_leaf_name = True
ts.show_scale = True

# Larger leaf labels
for node in t.iter_leaves():
    node.img_style["size"] = 8


# Display
t.show(tree_style=ts)

# Render
t.render(
    "random_tree.pdf",
    tree_style=ts,
    w=180,
    units="mm"
)
