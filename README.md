# Propagation_MATLAB

If you use this code in published work, please cite paper[1] as listed below.
Summary:

Under the standard power law path loss model, this code demonstrates how a two-tier cellular network can be viewed as a one-tier isotropic network. 

Description: 

Important performance metrics in cellular networks like SINR, SIR etc, can be written as functions of the 'propagation process', which is the ratio of the fading variable and the path-loss/attenuation function. 

Turns out different cellular networks are often equivalent to each other in terms of propagation processes. Hence, two completely different networks will have the same SINR, SIR, etc despite appearing quite different.

The rise of many publications on so-called multi-tier cellular networks motivated the recent paper [1], which has results on 'propagation equivalence' of cellular networks. 

The file PropEquivalenceTwoTier.m examines a HOST231-Hata model propagation model for a two-tier network, and examines the isotropic density of its one-tier equivalent (where spatial averages have been used for its propagation parameters).


[1] B. Błaszczyszyn and H. Keeler, 'Equivalence and comparison of heterogeneous cellular networks', presented at WDN-CN2013, 2013, http://arxiv.org/abs/1306.0772


[2] H.P. Keeler, B. Błaszczyszyn and M. Karray, 'SINR-based coverage probability in cellular networks with arbitrary shadowing ', presented at ISIT, 2013, http://arxiv.org/abs/1301.6491
