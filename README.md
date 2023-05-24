# ReSync: Riemannian Subgradient Synchronization Algorithm

## Disclaimer

The experiment frameowrk is downloaded from https://github.com/ColeWyeth/DESC, which is associated with the paper "Robust Group Synchronization via Quadratic Programming, Yunpeng Shi, Cole Wyeth, and Gilad Lerman, ICML 2022". We upload these dependencies in order to be self-contained. In this experiment framework, we add the following two algorithms:

- SpectrIn.m: generating initialization rotations via spectral initization, which corresponds to our "SpectrIn" (i.e., Algorithm 2).

- ReSync.m: using Riemannian subgradient method (with linearly decay stepsize) to solve our nonconvex nonsmooth optimizaition formulation, which corresponds to our "ReSync" (i.e., Algorithm 1).

## Usage

Download matlab files.

- Open "Demo" folder and then manully add the paths of "Utils", "Models", and "Algorithms" folders to the current work path. 

- run compare_algorithms_demo.m to get a simple and quick comparison results between ReSync and other state-of-the-art methods;

- run convergence_figure.m to get a figure of convergence of our ReSync algorithm under different decay rate of stepsize.

- run compare_algorithms_figure.m to get a figure about comparison results under different problem settings. 

Please cite our work "ReSync: Riemannian Subgradient-based Robust Rotation Synchronization, Huikang Liu, Xiao Li, Anthony Man-Cho So, Preprint 2023." if you use our algorithms "SpectrIn" and/or "ReSync".

Creators:
Huikang Liu & Xiao Li

liuhuikang@sufe.edu.cn & lixiao@cuhk.edu.cn