import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numpy.matrixlib import mat
import matplotlib as mpl
mpl.rcParams['font.size'] = 16
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['text.usetex'] = True  # use latex for all text handling
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['savefig.format'] = 'pdf'  # gives best resolution plots
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['legend.fontsize'] = 16
# print mpl.rcParams.keys()
mpl.rcParams['text.latex.preamble'] = \
    r'\usepackage{siunitx}\usepackage{amsmath}\usepackage{amssymb}'

def model_in_out(x, ident):
    x = 4*x-2
    if ident:
        y1 = x[0]
        y2 = x[1] - x[0]**2
    else:
        y1 = x[0] + x[1]
        y2 = 2 * y1 - 0.5
    return [y1, y2]

def neg_likelihood(x, ident=True):
    x = 4*x-2
    if ident:
        vals = 4 * ((1.-x[0, :])**2 + 1*(x[1, :]-x[0, :]**2)**2)
        vals = np.exp(-vals)
    else:
        vals = 4 * (1 - (x[0, :] + x[1, :])) ** 2 + (2 * (x[0, :] + x[1, :]) - 0.5)** 2
        vals = np.exp(-vals)
    return vals

def contour_plot():
    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    # ax = Axes3D(fig)
    delta = 0.01
    # 生成代表X轴数据的列表
    x = np.arange(0, 1.0, delta)
    # 生成代表Y轴数据的列表
    y = np.arange(0, 1.0, delta)
    # 对x、y数据执行网格化
    X, Y = np.meshgrid(x, y)

    # Calculate Z
    Z1 = neg_likelihood(np.array([X, Y]), ident=True)
    Z2 = neg_likelihood(np.array([X, Y]), ident=False)
    Z = [Z1, Z2]
    for i in range(2):
        pp = axes[i].contourf(X, Y, Z[i],
                    #levels=np.linspace(Z[i].min(),Z[i].max(),20),
                        levels=np.linspace(0, 1, 20), 
                            cmap=matplotlib.cm.coolwarm)
    #     cmap=plt.get_cmap('rainbow'))  # 设置颜色映射
        axes[i].set_xlabel(r'$x_{1}$')
        axes[i].set_ylabel(r'$x_{2}$')
        axes[i].set_title(f'Model-{i+1}')
        axes[i].set_xlim(0, 1)
        axes[i].set_ylim(0, 1)
        plt.colorbar(pp, ax=axes[i])

    plt.savefig("figs/posterior_contour_examples.png", dpi=300)
    plt.show()
# END contour_plot()

def plot_3D():
    ident_bool = [True, False]
    for i in range(len(ident_bool)):
        fig = plt.figure(figsize=(10, 8))
        ax = Axes3D(fig)
        delta = 0.01
        # 生成代表X轴数据的列表
        x = np.arange(0, 1.0, delta)
        # 生成代表Y轴数据的列表
        y = np.arange(0, 1.0, delta)
        # 对x、y数据执行网格化
        X, Y = np.meshgrid(x, y)
        # Calculate Z
        Z = neg_likelihood(np.array([X, Y]), ident=ident_bool[i])
        
        # Plot contours
        ax.plot_surface(X, Y, Z,
            rstride=1,  # rstride（row）指定行的跨度
            cstride=1,  # cstride(column)指定列的跨度 
            cmap=plt.get_cmap('rainbow'))
        ax.set_xlabel(r'$x_{1}$')
        ax.set_ylabel(r'$x_{2}$')
        ax.set_zlabel('Posterior')
        plt.title(f"Model-{i+1}")
        plt.savefig(f"figs/posterior_3D_example_{i+1}.png", dpi=300)
        plt.show()
# END plot_3D()

def plot_xy(ident=True):
    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    delta = 0.01
    # 生成代表X轴数据的列表
    x = np.arange(0, 1.0, delta)
    # 生成代表Y轴数据的列表
    y = np.arange(0, 1.0, delta)
    # 对x、y数据执行网格化
    X, Y = np.meshgrid(x, y)

    # Calculate Z
    Z = model_in_out(np.array([X, Y]), ident=ident)
    # breakpoint()
    for i in range(2):
        pp = axes[i].contourf(X, Y, Z[i],
                    levels=np.linspace(Z[i].min(),Z[i].max(),20),
                        #levels=np.linspace(0, 1, 20), 
                            cmap=matplotlib.cm.coolwarm)
    #     cmap=plt.get_cmap('rainbow'))  # 设置颜色映射
        axes[i].set_xlabel(r'$x_{1}$')
        axes[i].set_ylabel(r'$x_{2}$')
        axes[i].set_title(r'$\hat{d_%i}$'%(i + 1))
        axes[i].set_xlim(0, 1)
        axes[i].set_ylim(0, 1)
        plt.colorbar(pp, ax=axes[i])

    plt.savefig(f"figs/response_surface_Model_{ident}.png", dpi=300)
    plt.show()

if __name__ == '__main__':
    contour_plot()
    plot_3D()
    plot_xy(True)
    plot_xy(False)