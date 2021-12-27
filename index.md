<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
<style>
* {
  box-sizing: border-box;
}

.column {
  float: left;
  width: 50%;
  padding: 5px;
}

/* Clearfix (clear floats) */
.row::after {
  content: "";
  clear: both;
  display: table;
}
</style>

## Overview
Label placement in maps is a very challenging task that is critical for the overall map quality.
Most previous work focused on designing and implementing fully automatic solutions, but the
resulting visual and aesthetic quality has not reached the same level of sophistication that skilled
human cartographers achieve. 
We investigate a different strategy that combines the strengths of
humans and algorithms. 
In our proposed method, first an initial labeling is computed that has
many well-placed labels but is not claiming to be perfect. 
Instead it serves as a starting point for
an expert user who can then interactively and locally modify the labeling where necessary. In an
iterative human-in-the-loop process alternating between user modifications and local algorithmic
updates and refinements the labeling can be tuned to the user’s needs.
We demonstrate our approach by performing different possible modification steps in a sample
workflow with a prototypical interactive labeling editor. Further, we report computational
performance results from a simulation experiment in QGIS, which investigates the differences
between exact and heuristic algorithms for semi-automatic map labeling. 
To that end, we compare several alternatives for recomputing the labeling after local modifications and updates,
as a major ingredient for an interactive labeling editor.


This is the QGIS labeling integration model of the open source project Hi-Labeling---an interactive labeling project developed at the  Institute of Logic and Computation, TU Wien, Vienna, Austria.
So far our framework contains various labeling algoritghms, intrgated into the QGIS labeling engine.



<iframe width="420" height="315" src="https://youtu.be/embed/1yXoEF6kSD8" frameborder="0" allowfullscreen></iframe>


## Licence

This project is under MIT licence. 
## Reference
If you want to know more about our implemented algorithms, please refer to our paper:<br>
**Exploring Semi-Automatic Map Labeling**<br>
Fabian Klute, Guangping Li, Raphael Löffler, Martin Nöllenburg, Manuela Schmidt<br>
Advances in Geographic Information Systems (SIGSPATIAL'19) (eds. Farnoush Banaei Kashani, Goce Trajcevski etc.), pages 13–22. ACM, 2019.<br>
[[bibtex]](https://www.ac.tuwien.ac.at/publications/fs-esl-19?file=../../publications/noellenburg-ac-web.bib) [[pdf]](https://arxiv.org/abs/1910.07799) [[doi]](https://dx.doi.org/10.1145/3347146.3359359)


Please acknowledge our work if you publish your result using our algorithms or code.


## Support
Please write us an [Email](mailto:guangping@ac.tuwien.ac.at) if you have questions.

We are glad to get any comments and error reports.

A random instance generator is available upon request.
## Acknowledgments
DynaMIS is part of the project ["Human-Centered Algorithm Engineering: Graph and Map Visualization"](https://www.ac.tuwien.ac.at/research/humalgo/) supported by the Austrian Science Fund (FWF) under Grant P31119.
