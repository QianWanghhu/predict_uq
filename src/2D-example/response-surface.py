import matplotlib
from matplotlib import markers
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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

from example_func import non_identifiable_example, identifiable_example
def model_in_out(x, ident):
    x = 4*x-2
    if ident:
        y1 = x[0]
        y2 = x[1] - x[0]**2
        return [y1, y2]
    else:
        y1 = x[0] + x[1]
        return y1

def neg_likelihood(x, ident=True):
    x = 4*x-2
    if ident:
        vals = 100 * ((1.-x[0, :])**2 + 100*(x[1, :]-x[0, :]**2 - 0.8)**2)
        vals = np.exp(-vals)
    else:
        vals = 4 * (1 - (x[0, :] + x[1, :])) ** 2
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
    # plt.show()
# END contour_plot()

def fittness_landscape(y_true, x_true, x_full, x_red, ident=True):
    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))
    delta = 0.001
    # 生成代表X轴数据的列表
    x = np.arange(0, 1.0, delta)
    # 生成代表Y轴数据的列表
    y = np.arange(0, 1.0, delta)
    # 对x、y数据执行网格化
    X, Y = np.meshgrid(x, y)

    # Calculate Z
    Z = np.array(model_in_out(np.array([X, Y]), ident=ident))
    for ii in range(y_true.shape[0]):
        if ident:
            lsq = np.round(np.sqrt((Z[0] - y_true[ii, 0]) ** 2 + (Z[1] - y_true[ii, 1])**2), 2)
        else:
            lsq = np.round(np.abs(Z - y_true[ii, 0]), 2)
        # breakpoint()
        pp = axes[ii].contourf(X, Y, lsq,
                    levels=np.linspace(lsq.min(),lsq.max(),20),
                        #levels=np.linspace(0, 1, 20), 
                            cmap=matplotlib.cm.coolwarm)
    #     cmap=plt.get_cmap('rainbow'))  # 设置颜色映射
        axes[ii].plot(x_true[ii, 0], x_true[ii, 1], linestyle='', color='r', marker ='o', alpha=0.5, ms=6)
        axes[ii].plot(x_full[ii, 0], x_full[ii, 1], linestyle='', color='k', marker ='o', ms=6)
        axes[ii].plot(x_red[ii, 0], x_red[ii, 1], linestyle='', color='darkgreen', marker ='o', ms=6)
        axes[ii].set_xlabel(r'$x_{1}$')
        axes[ii].hlines(0.5, 0, 1, linestyle='--', linewidth=1.5)
        axes[ii].set_ylabel(r'$x_{2}$')
        axes[ii].set_title(r'$EM-%i$'%(ii + 1))
        axes[ii].set_xlim(0, 1)
        axes[ii].set_ylim(0, 1)
    cbar = plt.colorbar(pp, ax=axes, spacing='uniform')
    cbar.set_ticks([0, 0.5, 1.0, 1.5, 2, 3, 4, 5, 6])
    plt.savefig(f"figs/fittness_landscape_Model_{ident}.png", dpi=300)

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
    for i in range(2):
        pp = axes[i].contourf(X, Y, Z[i],
                    levels=np.linspace(Z[i].min(),Z[i].max(),20),
                        #levels=np.linspace(0, 1, 20), 
                            cmap=matplotlib.cm.coolwarm)
    #     cmap=plt.get_cmap('rainbow'))  # 设置颜色映射
        axes[i].set_xlabel(r'$x_{1}$')
        axes[i].set_ylabel(r'$x_{2}$')
        axes[i].set_title(r'$\hat{O_%i}$'%(i + 1))
        axes[i].set_xlim(0, 1)
        axes[i].set_ylim(0, 1)
        cbar = plt.colorbar(pp, ax=axes[i])
    plt.savefig(f"figs/response_surface_Model_{ident}.png", dpi=300)
    # plt.show()

def solve_analytic(y, ident):
    assert y.shape[1] == 2, "The second dimension of y values should be 2."
    if ident:
        x_analytic = np.zeros_like(y)
        for ii in range(y.shape[0]):
            x_analytic[ii, 0] = (y[ii, 0] + 2) / 4
            x_analytic[ii, 1] = (y[ii, 0]**2 + 2 + y[ii, 1]) / 4
    else:
        x_analytic = np.zeros(shape=(3, 2, 100))
        for ii in range(y.shape[0]):
            x_analytic[ii, 0] = np.linspace(0, 1, num=100)
            x_analytic[ii, 1] = (y[ii, 0] + 4) / 4 - x_analytic[ii, 0]
            # x_analytic[ii, 0] = (y[ii, 0] + 4) / 4
    return x_analytic

if __name__ == '__main__':
    file_names = ['example_1_work_2_pars', 'example_2_work_2_pars']
    file_names_red = ['example_1_work_1_pars', 'example_2_work_1_pars']
    ident_bools = [True, False]
    par_files = ['50', '9', '11', '10'] 
    for ii  in range(2):#len(file_names)
        y = pd.read_csv(f'0119/{file_names[ii]}/example.base.obs.csv', index_col='real_name').values
        x_ies = pd.read_csv(f'0119/{file_names[ii]}/example.{par_files[ii]}.par.csv', index_col='real_name').values/100
        x_ies_red = pd.read_csv(f'0119/{file_names_red[ii]}/example.{par_files[ii+2]}.par.csv', index_col='real_name').values/100
        x_analytic = solve_analytic(y, ident=ident_bools[ii])
        fittness_landscape(y, x_analytic, x_ies, x_ies_red, ident = ident_bools[ii])