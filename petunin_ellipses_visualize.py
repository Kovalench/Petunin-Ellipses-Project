import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def read_data(filename):
    sections = {
        'ScalingFactors:': 'scaling_factors',
        'SquareCenter:': 'square_center',
        'OriginalCenter:': 'original_center',
        'CircleRadii:': 'circle_radii',
        'EllipseRadii:': 'ellipse_radii',
        'OriginalSet1:': 'original_set1',
        'ScaledSet1:': 'scaled_set1',
        'OriginalSet2:': 'original_set2',
        'Set2:': 'set2',
        'OriginalMainLine:': 'original_main_line',
        'MainLine:': 'main_line',
        'OriginalPerpendiculars:': 'original_perpendiculars',
        'Perpendiculars:': 'perpendiculars',
        'ScaledPentagon:': 'scaled_polygon',
        'ScaledPolygon:': 'scaled_polygon', 
        'OriginalPentagon:': 'original_polygon',
        'OriginalPolygon:': 'original_polygon'
    }

    data = {v: [] for v in sections.values()}
    current_section = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line in sections:
                current_section = sections[line]
            elif current_section:
                if current_section in ['scaling_factors', 'square_center', 'original_center']:
                    data[current_section].append(list(map(float, line.split())))
                elif current_section in ['circle_radii', 'ellipse_radii']:
                    data[current_section].append(float(line))
                else:
                    data[current_section].append(list(map(float, line.split())))

    # Determine polygon type based on vertex count
    if data['original_polygon']:
        vertex_count = len(data['original_polygon'])
        data['polygon_type'] = 'quadrilateral' if vertex_count == 4 else 'pentagon'
    
    return data

def show_plot(fig):
    plt.tight_layout()
    plt.show()
    plt.close(fig)

def plot_original(data):
    fig, ax = plt.subplots(figsize=(10, 8))
    title = 'Original Space'
    ax.set_title(title)

    if data['original_polygon']:
        polygon = patches.Polygon(data['original_polygon'], closed=True, 
                                edgecolor='black', fill=False, 
                                label=data['polygon_type'])
        ax.add_patch(polygon)

    if data['original_set1']:
        x1, y1 = zip(*data['original_set1'])
        ax.scatter(x1, y1, c='blue', label=f'Set1 ({len(x1)})')

    if data['original_set2']:
        x2, y2 = zip(*data['original_set2'])
        ax.scatter(x2, y2, c='orange', label=f'Set2 ({len(x2)})')

    if data['original_main_line'] and len(data['original_main_line']) >= 2:
        x_line, y_line = zip(*data['original_main_line'])
        ax.plot(x_line, y_line, 'g-', linewidth=2, label='Main Line')

    for i, p in enumerate(data['original_perpendiculars']):
        if len(p) == 4:
            ax.plot([p[0], p[2]], [p[1], p[3]], 'r-', linewidth=1, 
                   label='Perpendicular' if i == 0 else None)

    ax.set_aspect('equal')
    ax.grid(True)
    ax.legend()
    show_plot(fig)

def plot_scaled(data):
    fig, ax = plt.subplots(figsize=(10, 8))
    title = 'Scaled Space'
    ax.set_title(title)

    cx, cy = data['square_center'][0]
    ax.scatter([cx], [cy], c='black', marker='+', s=100, label='Center')

    for r in data['circle_radii']:
        circle = patches.Circle((cx, cy), r, edgecolor='green', 
                              fill=False, alpha=0.3, label='Circle' if r == data['circle_radii'][0] else None)
        ax.add_patch(circle)

    if data['scaled_set1']:
        x1, y1 = zip(*data['scaled_set1'])
        ax.scatter(x1, y1, c='blue', label=f'Set1 ({len(x1)})')

    if data['set2']:
        x2, y2 = zip(*data['set2'])
        ax.scatter(x2, y2, c='orange', label=f'Set2 ({len(x2)})')

    if data['scaled_polygon']:
        polygon = patches.Polygon(data['scaled_polygon'], closed=True, 
                                edgecolor='black', fill=False,
                                label=data['polygon_type'])
        ax.add_patch(polygon)

    if data['main_line'] and len(data['main_line']) >= 2:
        x_line, y_line = zip(*data['main_line'])
        ax.plot(x_line, y_line, 'g-', linewidth=2, label='Main Line')

    for i, p in enumerate(data['perpendiculars']):
        if len(p) == 4:
            ax.plot([p[0], p[2]], [p[1], p[3]], 'r-', linewidth=1,
                   label='Perpendicular' if i == 0 else None)

    ax.set_aspect('equal')
    ax.grid(True)
    ax.legend()
    show_plot(fig)

def plot_ellipse(data):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_title("Petunin's ellipses")

    ox, oy = data['original_center'][0]
    scale_x, scale_y = data['scaling_factors'][0]

    ax.scatter([ox], [oy], c='black', marker='+', s=100, label='Center')

    angle = 0
    if data['original_main_line'] and len(data['original_main_line']) >= 2:
        x1, y1 = data['original_main_line'][0]
        x2, y2 = data['original_main_line'][1]
        angle = np.degrees(np.arctan2(y2 - y1, x2 - x1))

    for r in data['ellipse_radii']:
        ellipse = patches.Ellipse((ox, oy), width=2*r/scale_x, height=2*r/scale_y, angle=angle,
                                  edgecolor='green', fill=False, alpha=0.3)
        ax.add_patch(ellipse)

    if data['original_set1']:
        x1, y1 = zip(*data['original_set1'])
        ax.scatter(x1, y1, c='blue', label=f'Set1 ({len(x1)})')

    if data['original_set2']:
        x2, y2 = zip(*data['original_set2'])
        ax.scatter(x2, y2, c='orange', label=f'Set2 ({len(x2)})')

    if data['original_polygon']:
        polygon = patches.Polygon(data['original_polygon'], closed=True, 
                                edgecolor='black', fill=False, 
                                label=data['polygon_type'])
        ax.add_patch(polygon)

    if data['original_main_line'] and len(data['original_main_line']) >= 2:
        x_line, y_line = zip(*data['original_main_line'])
        ax.plot(x_line, y_line, 'g-', linewidth=2)

    for p in data['original_perpendiculars']:
        if len(p) == 4:
            ax.plot([p[0], p[2]], [p[1], p[3]], 'r-', linewidth=1)

    ax.set_aspect('equal')
    ax.grid(True)
    ax.legend()
    show_plot(fig)

def main():
    try:
        data = read_data("output_data.txt")
        plot_original(data)
        plot_scaled(data)
        plot_ellipse(data)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
