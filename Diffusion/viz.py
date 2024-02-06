import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_cell(cell, ax):
    # Adjusting for the center coordinates and level-dependent size
    half_side_length = cell['width'] / 2
    topLeftX = cell['x'] - half_side_length
    topLeftY = cell['y'] - half_side_length

    rect = patches.Rectangle((topLeftX, topLeftY), cell['width'], cell['width'], linewidth=1, edgecolor='r', facecolor='none')
    if cell['isLeaf']:
        ax.text(cell['x'], cell['y'], str(cell['CID']), ha='center', va='center', fontsize=8, color='blue')

    ax.add_patch(rect)

    if not cell['isLeaf']:
        for child in cell['children']:
            draw_cell(child, ax)

def visualize_quadtree(json_file):
    with open(json_file, 'r') as f:
        quadtree = json.load(f)

    fig, ax = plt.subplots()
    draw_cell(quadtree, ax)
    plt.axis('equal')  # Ensures the plot is square in shape
    plt.show()

visualize_quadtree('quadtree.json')
