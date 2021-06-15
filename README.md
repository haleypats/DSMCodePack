# DSMCodePack
This repository provides a set of python functions that facilitate the calculation of dose-surface maps (DSMs) from radiotherapy treatment planning data. These functions can be used to create both of the two most common types of DSMs: (1) planar DSMs, as described in the methods of Buettner et al, and (2) non-planar DSMs, as described in the methods of Hoogeman et al.

## Dose-Surface Maps
Radiotherapy research typically uses dose-volume histograms (DVHs), which condense 3D dose distributions within a region of interest (ROI) into a 2D representation, as their main source of dosimetric data for dose-outcome studies. However, DVHs lack spatial information, making it difficult to identify specific sub regions or dose-spread patterns in organs that can be indicative of higher risk of adverse events. Dose-surface maps (DSMs) preserve this spatial information and make it possible to identify spatio-dosimetric features on the surface of any ROI. 

DSMs are most commonly used for hollow, tubular organs (like the esophagus and rectum) to characterize the dose to the actual tissue wall and ignore dose to the filling contents. They can be thought of as tubes “painted” with dose information that have been cut open to allow them to lay flat as a 2D dose map, similar to a piece of radiochromic film irradiated rolled up then unfurled. 

While relatively simple in concept, different groups have developed different methods to calculate them, with the largest methodological difference being how groups define the planes along which they sample dose to the ROI. The simpler method, used by Buettner *et al*, is to merely use the axial slices as they are typically drawn in a treatment planning system and ignore any curvature of the organ. The more complex method accounts for the organ’s curvature and uses sampling planes perpendicular to the central path of the organ, and has been employed by researchers like Hoogeman *et al*. This repository facilitates the calculation of both types of DSMs and will now on refer to them as **planar** (Buettner method) and **non-planar** (Hoogemen method) DSMs respectively. 

## DSM calculation methods
<img width="696" alt="Screen Shot 2021-06-15 at 9 10 51 AM" src="https://user-images.githubusercontent.com/14338664/122058470-a4a73f00-cdb9-11eb-84c3-2ee42b6389d9.png">

#### Sampling style:
Two different in-plane sampling styles exist in the literature, equidistant and equiangular. Currently our code only supports the dominant equiangular approach.

<img width="440" alt="Screen Shot 2021-06-15 at 9 11 33 AM" src="https://user-images.githubusercontent.com/14338664/122058544-b852a580-cdb9-11eb-89ea-a1d26748fca2.png">

#### Cutting style:
While cut location often varies based on the ROI being investigated (eg. posterior/anterior), some groups opt to define the cut location as the **direction**-most point on a slice of the ROI, others instead used the point directly **direction** to the centroid of the slice. Our code uses the latter, centroid-based approach as it has been reported to be more robust against irregularly shaped slices. As it stands right now the code slices at the point posterior to the centroid, but this can be adjusted by the user in the source code if desired.

<img width="456" alt="Screen Shot 2021-06-15 at 9 12 00 AM" src="https://user-images.githubusercontent.com/14338664/122058597-c86a8500-cdb9-11eb-8150-9426aed88696.png">

#### Slice Collision Correction (non-planar only):
When using the non-planar DSM approach it is possible to define sampling planes that intersect with one another, which can lead to the same voxel on the ROI’s surface being attributed to two different voxels in a DSM. In the event of plane collisions we use the method described by Witztum *et al* to correct the positioning of the colliding slices. In summary, the collision correction function draws connecting rays between the two slices sandwiching the collision at set equiangular points. It then defines new sampling planes for the colliding slices using points a set fractional distance along the connecting rays (eg. ⅓ the distance along each ray) to create the new plane equation.

<img width="572" alt="Screen Shot 2021-06-15 at 9 12 43 AM" src="https://user-images.githubusercontent.com/14338664/122058701-e2a46300-cdb9-11eb-8fb8-179b508b6f11.png">

### Dependencies
In addition to the common modules **numpy, scipy,** and **matplotlib**, this repository also requires **pyvista, pydicom,** and **shapely** to run. Please note that pyvista is only compatible with python versions 3.5 and greater.

### Open-source use conditions
We ask that any publications, conference presentations, or similar disseminations of research performed using our code base credit this repository.

### References
* Buettner *et al* (2009) DOI: 10.1088/0031-9155/54/17/005
* Hoogeman *et al* (2004) DOI: 10.1016/j.radonc.2003.11.015
* Witztum *et al* (2006) DOI: 10.1118/1.4964790
