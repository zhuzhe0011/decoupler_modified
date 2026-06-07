# decoupler_modified
 
## Overview
On the basis of decoupler, Ucell method is added, and P value calculation is optional.
The usage of this python package completely follows decoupler, please know [decoupler](https://github.com/scverse/decoupler) first.


## Getting started

Fully compatible integration of [Ucell](https://github.com/carmonalab/pyucell) (missing_genes = "skip": missing genes are simply removed) entered decoupler.
```python
import decoupler as dc
adata, net = dc.ds.toy()
dc.mt.ucell(adata, net, tmin=3)
```

Now, when using any method of decoupler, you can choose not to calculate the p value to improve the running speed and avoid possible bugs.
```python
dc.mt.mlm(adata, net, tmin=3, pvalue = False)
```

## Reference

- [decoupler](https://github.com/scverse/decoupler)
- [Ucell](https://github.com/carmonalab/pyucell)

## Citation

> Badia-i-Mompel P., Vélez Santiago J., Braunger J., Geiss C., Dimitrov D.,
Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores R.O.
and Saez-Rodriguez J. 2022. decoupleR: Ensemble of computational methods
to infer biological activities from omics data. Bioinformatics Advances.
<https://doi.org/10.1093/bioadv/vbac016>
