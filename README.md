## Designing Conistent Cortical Surface Features

[Chung et al. 2024](https://github.com/laplcebeltrami/sulcaltree/blob/main/chung.2024.OHBM.pdf) is presented in OHBM 2024.

## Sulcal Pattern Matching with the Wasserstein Distance

The method and codes are published in 

Chen, Z., Das, S., Chung, M.K. 2023, [Sulcal Pattern Matching with the Wasserstein Distance](https://github.com/laplcebeltrami/sulcaltree/blob/main/chen.2023.ISBI.pdf), 
International Symposium in Biomedcial Imaging (ISBI). [Poster version](https://github.com/laplcebeltrami/sulcaltree/blob/main/ISBI2023poster.pdf). The following script will perform the basic sulcal pattern matching.


SCRIPT1_dataPreprocess.m prepares data and performs heat kernel smoothing (Sections 2.1 and 2.2)

SCRIPT2_registration.m performs the Wasserstein disatnce based sulcal pattern matching (Sections 2.3 and 2.4) 

SCRIPT3_validation.m performs the Validation against the Hungarian Algorithm (Section 3.1)





(C) 2022- Zijian Chen, Moo K. Chung
University of Wisconsin-Madison








## 
## SPHARM modeling of sucal pattern

SCRIPT.m generating sulcal tree patterns below published in 

![alt text](https://github.com/laplcebeltrami/sulcaltree/blob/main/resampled.png?raw=true)

Huang, S.-G., Lyu, I., Qiu, A., Chung, M.K. 2020. Fast polynomial approximation of heat kernel convolution on manifolds and its application to brain sulcal and gyral graph pattern analysis, IEEE Transactions on 
Medical Imaging 39:2201-2212  https://pages.stat.wisc.edu/~mchung/papers/huang.2020.TMI.pdf

We can perform the heat kernel smoothing using SPHARM based on  

Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. Weighted Fourier series representation and its application to quantifying the amount of gray matter. Special Issue of  IEEE Transactions on Medical Imaging 26:566-581. 

See https://pages.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/weighted-SPHARM.html for additional codes. The method can resample mesh sizes and reduce the data size further.

![alt text](https://github.com/laplcebeltrami/sulcaltree/blob/main/sulcalpattern.png?raw=true)





(C) 2022- Moo K. Chung
University of Wisconsin-Madison
