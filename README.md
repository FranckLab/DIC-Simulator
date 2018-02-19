## Important pages
* [Download latest version v1.1!](https://github.com/FranckLab/DIC-Simulator/releases)
* [FAQ](https://github.com/FranckLab/DIC-Simulator/blob/master/README.md#faq)
* [Questions/Issues](https://github.com/FranckLab/DIC-Simulator/issues)
* [Cite](https://github.com/FranckLab/DIC-Simulator/blob/master/README.md#cite)
* [Franck Lab](http://franck.engin.brown.edu)

## Purpose
The Digital Image Correlation (DIC) simulator is a Matlab-based GUI that assesses how “well” a particular user-provided image correlates given 4 different simulated displacement fields. The purpose of the simulator is to allow the user to optimize the speckle and intensity patterns in an image, without having to run experimentally-based calibration tests, but rather have Matlab simulate them. The DIC simulator performs the following steps:

1. User loads DIC image of interest

2. Matlab maps the intensities in the uploaded image (undeformed configuration) into a set of 4 different analytically predetermined deformation fields

3. A standard DIC algorithm correlates each image pair (deformed-undeformed config.) per given deformation field and returns the displacement matrices u1 and u2

4. The Matlab-GUI plots the determined displacement field as color contours, and as histograms of the error difference between the analytically-prescribed and DIC correlated displacement fields

**Health Warning!**

The algorithm expects the speckle pattern to cover the entire field-of-view of the image and will either error out or provide eroneous result if this is not the case.

## FAQ
**The GUI looks squashed and button don't appear properly, what should I do?**

The GUI in newer versions of Matlab and high resolution screens, can sometimes get squashed with various buttons and icons not appearing properly. The easiest way to solve this issue is to resize the GUI layout manually on your computer. Use the command `guide DIC_simulator.fig` and resize the buttons.


**What is the recommended minimum size of the input image stack?**

We recommend that the input image stack should have at least have 96 pixels in each dimension. 


## Cite
If used please cite:
[Estrada, J., Franck, C.,”Intuitive Interface for the Quantitative Evaluation of Speckle Patterns for Use in Digital Image and Volume Correlation Techniques”, J. Applied Mechanics, 2015.](https://appliedmechanics.asmedigitalcollection.asme.org/article.aspx?articleid=2336782)

```bibtex
@article{estrada2015intuitive,
  title={Intuitive Interface for the Quantitative Evaluation of Speckle Patterns for Use in Digital Image and Volume Correlation Techniques},
  author={Estrada, Jonathan B and Franck, Christian},
  journal={Journal of Applied Mechanics},
  volume={82},
  number={9},
  pages={095001},
  year={2015},
  publisher={American Society of Mechanical Engineers}
}
```

## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/DIC-Simulator/blob/master/README.md#faq) and [Questions/Issues](https://github.com/FranckLab/DIC-Simulator/issues). Add a new question if similar issue hasn't been reported. We shall help you at the earliest. The author's contact information can be found at [Franck Lab](http://franck.engin.brown.edu).
