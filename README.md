# DIC-Simulator
This simulator is a Matlab-based GUI that assesses how “well” a particular user-provided image correlates given 4 different simulated displacement fields. The purpose of the simulator is to allow the user to optimize the speckle and intensity patterns in an image, without having to run experimentally-based calibration tests, but rather have Matlab simulate them. The DIC simulator performs the following steps:

1. User loads DIC image of interest

2. Matlab maps the intensities in the uploaded image (undeformed configuration) into a set of 4 different analytically predetermined deformation fields

3. A standard DIC algorithm correlates each image pair (deformed-undeformed config.) per given deformation field and returns the displacement matrices u1 and u2

4. The Matlab-GUI plots the determined displacement field as color contours, and as histograms of the error difference between the analytically-prescribed and DIC correlated displacement fields
