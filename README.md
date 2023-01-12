# ILFHC
The Inertial Lift Force Helper Class provides a minimal interface to stripped research data relating to particle migration in curved microfluidic ducts.

Please ensure you cite our Journal of Fluid Mechanics paper (https://doi.org/10.1017/jfm.2019.323) if you use this code/data. 
Note that interface makes use of the non-dimensional scaling described within the paper.
This code is provided under an MIT license (see the included LICENSE.txt or refer to https://opensource.org/licenses/MIT). 
I would also appreciate it if you contact me and to let me know if you use this code/data.
Please also don't hesitate to contact me if you have any questions/queries.

Brendan Harding, 2019.



2023 update:

The folder MDN-python contains data and a helper class associated with the results published in the Journal of Fluid Mechanics paper "Inertial focusing of spherical particles in curved microfluidic ducts at moderate Dean numbers" (accepted on 4 Jan 2023).
This data has been pre-processed a little differently than was used in the paper. 
This has been done primarily to facilitate ease of access and efficiency of the helper class (the original data is several GB in size but has been reduced here to just ~100MB).
A consequence of this pre-processing is that results derived from this class/data may have very minor quantitative differences with the results in the paper.
However, we don't expect any significant qualitative differences (if you notice such a difference it is likely a bug in the code and you should get in touch).
The interface currently consists of a minumum building of interpolants to query migration velocities according to a number of input parameters.
I intend to develop this further in collaboration with users as the need arises, so please get in touch if you intend to use this class/data in your research and I will be happy to help update the interface to suit your needs.
