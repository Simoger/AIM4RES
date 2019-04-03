# AIM4RES user guide

AIM4RES is a finite-differences forward and inverse anisotropic modeling open source MATLAB library. By considering the resistivity as a tensor, allowing for the anisotropy estimation, AIM4RES provides three model sections from Electrical Resistivity Tomography data: highest resistivity &#961;<sub>1</sub>, lowest resistivity &#961;<sub>3</sub> and the angle of anisotropy &#952;. Details about the theory will be found in the article _AIM4RES, an open-source 2.5D finite differences MATLAB library for anisotropic electrical resistivity modeling_ currently in review.

## Library description

The library is open source and freely available for any purpose. It contains 4 folders:
* _functions_: this folder contains all the needed functions for AIM4RES proper functioning. Each function contains a header describing its purpose, the needed input and its outputs. A description of the 3 major structures __param__, __XYZ__ and __Inv__ are thoroughly explained below.

* _scripts_: this folder contains the executable example scripts.


* _data_: this folder contains the data used the scripts previously described, and the output data resulting from the scripts. These latter allow for direct parameters exploration or drawing without previous calculation (for the inversion).

* _figures_: this folder contains the .fig files obtain from the previous scripts, and correspond to the figures presented in the article.

All scripts are executable form the folder __./aim4res__

### Prerequisites

Having a working MATLAB distribution, installed on any operating system. No installation is needed for this library and it is directly usable in a MATLAB environment.

## What AIM4RES has inside

### Key functions

AIM4RES is based on some key functions, which functioning is described in the accompanying article:

* *calcul_u_S_anis* drives forward modeling as well as sensitivity calculation,  

* *matrix_coeff_anis* is the function building the capacitance matrix thoroughly described in the article and is called by *calcul_u_S_anis*,

* *gauss_newton_inversion_anis* drives the inverse modeling. This function needs the two structures __param__ and __XYZ__ and return the __Inv__ structure.

### Key parameters

Three structures are needed in the modeling process:

* __param__: parameters structure. Some of its items are computed by the different functions. The others are user filled out:
  * **flag.geo_factor**: geometric factor associated to each measurements. Needed for the apparent resistivity calculation used in the inverse modeling,
  * **cell_size**: determine the block size. Assuming a raw block is 1 m<sup>2</sup>, **cell_size** slices it dividing its horizontal dimension by **cell_size(1)** and its vertical dimension by **cell_size(2)**,
  * **nb_pad_bloc**, **nb_raff**,  **nb_surr** and **fact** are described in *grille2d_elect.m* header,
  * **rho.xx**, **rho.zz**, **rho.xz** and **rho.yy**: synthetic model resistivity components (&#961;<sub>xx</sub>, &#961;<sub>zz</sub>, &#961;<sub>xz</sub> and &#961;<sub>yy</sub>),
  * **flag.inv.p**: determines wether inverse modeling computes only the two &#961;<sub>xx</sub> and &#961;<sub>zz</sub> components (**flag.inv.p** = 2) ; or the three &#961;<sub>xx</sub>, &#961;<sub>zz</sub> and &#961;<sub>xz</sub> components (**flag.inv.p** = 3),
  * **data_nb**: &#961;<sub>xx</sub>, &#961;<sub>zz</sub> or &#961;<sub>xz</sub> components size,
  * **K**: geometric factor,
  * **MEAS.Res**: resistance data (U/I), synthetically computed or measured on the field,
  * **anis_init**: homogeneous initial anisotropy considered in the inversion,
  * **const_ind**, **const_vect** and **const_TrueFalse**: respectively constrained model cells indices, constrained model cells resistivity value and application (or not) of the constraints,
  * **const.***: other constraints parameters (see *script_preparation_inversion.m*),
  * **inv.***: inverse modeling features
    * **invparam** = 'log resistivity', 'resistivity' or 'resistance' according to the desired inversion. 'log resistivity' is highly recommended,
    * **fct_reg** = 'flatness' or 'smoothness', depending on wether first derivative or second derivative is prefered as regularization, respectively,
    * **appl_fct_reg** = 'model' or 'model perturbation'
    * **dataweight**: see *Calc_data_weight.m* header,
    * **alx**, **alz**, **als**: see *CtC_anis.m* header
    * **BETA**: regularization coefficient. Usually higher for the anisotropic problem than for the isotropic one,
    * **rms_model**
    * **tol**: tolerance threshold iteration criterion,
    * **maxit**: maximum amount of iterations,
    * **weight**: weights applied on the regularization matrix (true of false),
    * **weightMinit** and **weightMaxit**: respectively first and last iteration for the weigthing application,
    * **weightFun**: 'Distance weighting' or 'Sensitivity weighting' (see Li and Oldenburg, 1996)
    * **alp**: Armijo coefficient.


* __XYZ__: geometric structures
  * **surface_electrode**, **borehole1_electrode** and **borehole2_electrode**: respectively surface, first borehole and second borehole admissible coordinates. Can contain unused coordinates but __*MUST*__ contain all the used coordinates (x: first column, z: second column),
  * **MEAS.***: quadrupoles coordinates
    * **C1**: current electrode 1 coordinates. n<sup>th</sup> row = n<sup>th</sup> measure (x: first column, z: second column)
    * **C2**:current electrode 2 coordinates,
    * **P1**:potential electrode 1 coordinates,
    * **P2**:potential electrode 2 coordinates,
  * **areas**: area formed b the quadrupoles. n<sup>th</sup> row = n<sup>th</sup> measure.


* __Inv__: inverse modeling products
  * **rho.***: &#961;<sub>1</sub>, &#961;<sub>2</sub>, &#961;<sub>xx</sub>, &#961;<sub>zz</sub>, &#961;<sub>xz</sub> and &#952; inverted sections at each iteration,
  * **D**: data weighting matrix (see *Calc_data_weight.m* header),
  * **d_cal**: apparent resistivity computed on the inverted model at each iteration
  * **rho_app_pos_index**: indices of the considered measures at each iteration. Inverse modeling considers logarithmic values of apparent resistivity, **rho_app_pos_index** only keep positive apparent resistivities,
  * **rms**: root mean square value at each iteration,
  * **beta_weighting**: weighting coefficient stored at each iteration,
  * **sens_weighting**: weighting sections stored at each iteration,
  * **Ki2**: &#967;<sup>2</sup> value stored at each iteration,
  * **CTC**: regularization matrix stored at each iteration

## Running the tests

Ready for use scripts have been added to AIM4RES library. They reproduce the figures presented in the article. These scripts can be edited for any use to adapt any ERT data.

Scripts are to be executed from AIM4RES root folder:

  * *script_forward_validation.m*: forward modeling considering the two synthetic models described in the article. The script waits for user input to compute the two different forward modeling proposed.

  * *script_inversion.m*: inverse modeling considering the synthetic model described in the article. The script waits for user input to determine if synthetic data have to be computed beforehand from *script_preparation_inversion.m*, or if already existing MATLAB file is to be considered ('synt_data_for_inversion.mat').

  * *script_drawing.m*: script producing the figures presented in the article.


## Authors

* **Simon GERNEZ**

* **Abderrezak BOUCHEDDA**

Members of the [Laboratoire d'Interprétation et Acquisition des Mesures en Géosciences](https://github.com/groupeLIAMG) of INRS-ETE University, Quebec, Canada.

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE. See LICENSE.txt for details.
