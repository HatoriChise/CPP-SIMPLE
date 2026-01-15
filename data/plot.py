import matplotlib.pyplot as plt
import json

with open("data/mesh_info.json", "r") as f:
    data = json.load(f)


def plot_grid(grid_data):
    fig, ax = plt.subplots(figsize=(8, 8))

    # 1. 绘制垂直网格线 (x = const)
    for x in grid_data["node_x"]:
        ax.axvline(x, color="gray", linestyle="-", linewidth=0.5, alpha=0.6)

    # 2. 绘制水平网格线 (y = const)
    for y in grid_data["node_y"]:
        ax.axhline(y, color="gray", linestyle="-", linewidth=0.5, alpha=0.6)

    # 3. 绘制胞元中心 (Cell Centers)
    # 使用 meshgrid 生成所有中心点的坐标对
    import numpy as np

    cx, cy = np.meshgrid(grid_data["cell_centers_x"], grid_data["cell_centers_y"])
    ax.scatter(cx, cy, s=2, c="red", marker="o", label="Cell Centers", alpha=0.5)

    # 设置图形属性
    ax.set_aspect("equal")
    ax.set_title(f"Grid Visualization ({data['ncx']}x{data['ncy']})", fontsize=14)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_xlim(min(grid_data["node_x"]), max(grid_data["node_x"]))
    ax.set_ylim(min(grid_data["node_y"]), max(grid_data["node_y"]))
    ax.grid(False)  # 关闭默认网格，因为我们手动画了
    ax.legend(loc="upper right")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_grid(data)
