# KNN-Classifier
[Document link](https://docs.google.com/document/d/1DlprCeO_yCAFl605km7NrIExPAnbm1X1axNf4igZ-3I/edit?usp=sharing)
[Stability Paper By Adam Dickey](https://www-proquest-com.proxy.uchicago.edu/docview/914436445)

Tasks:

[x] 1. Compare population activity irrespective of tongue kinematics.
- To support the claim that overall, firing is more sparce for the majority of units. This is measured per trial and averaged across all trials per unit, and the comparison can be depicted with four scatter subplots.
- Needs refinement, spesific to phases of gape cycle as control

[ ] 1.1 Compare population activity with respect to gape cycle
- Compare duration of sequence, rythmic chewing, and contacts

[x] 2. Increase k value (k = 111?) for Type A and Type B
- Since votes are based on activity across units, the number of units in a cortical regions does not matter.

[x] 3. Run Type A classification on 3 or 5 (should be odd) most frequent contact patterns.
- The heavy skew in contact patterns might mean the high performances we have so far do not accurately refect the ability of the cortex to predict contact patterns.

[x] 4. Sort grouped horizontal bar graphs by highest (combined) frequency and abridge graph.
- Most contact patterns occurs very rarely. Either combine the rare ones as 'other' group or exclude that data and leave a note that the rest do not total over e.g. 1% occurance.

[ ] 5. Refactor original contact events code
- Features missing: Gape cycle type, trajectory, mean position, angle, general stats

[ ] 6. Train and predict classifiers using data from different days.
- First, train on control and predict on nerve block. The same level of performance suggests the cortext responds the same way to LP contact.
- Second, train on nerve block and predict on control, which is expected to be worse. This is a control case mostly. If 1 is the same and 2 is worse, the firing pattern is nonlinear.
- If it's not worse, then try training using combined and predicting on combined/control/nerve block
- If it is, no step 3

[ ] 7. Correlation analysis between contact with palate and teeth in sets of three.
- Reviewers might point out that classifier perf that excludes teeth might include neural data explained by contact with teeth. If contact with teeth are highly correlated, we can eliminate this concern.
- PC analysis for x-y-z of teeth (sets of three) close to palate region.
- Overlay average 3D position when contact with both palate and teeth to show they are not variable.

[ ] 8. Type B Classifier
- The script control flow correctly produces classifiers for each marker. However we should examine whether the absence of 'contact events' in markers 2 3 5 and 6 is expected.

[ ] 9. Save and compare 'spiketable' for control and nerve block for Type A and Type B.
- Each classifier type needs a supplemental figure that characterizes the changes in firing rate relative to contact. Fig. 3 of the poster attempts something similar.

[ ] 10. Visualize area around tongue markers
- We need to know whether we should interpolate marker positions to better represent the tongue and its contact with oral structures, is 8mm enough?
- Including information pertaining to the size of each oral structure, like distance between extreme verticies in the mesh, would also be helpful

[ ] 11. Determine how oral structure size biases touch frequencies
- Mid and Pos palate regions have nearly twice the locators than Ant. what justifications do we have for choosing them where they are?
- We can eliminate the confound by normalizing FR to the frequency of contact, as a weight.

[ ] 12. Create any contact vs no contact comparison figure for each condition for each cortical region
-Should take the form of 8 double bar graphs. Start with palate to explain classifier results, then move to teeth.

Type A Classifier, 7 most frequent contact patterns using all tongue markers, k = 1

Ry20190228 (Control)
Ry20190227 (All nerve blocks)
    
Ye20190508 (No Lingual Nerve)
Ye20190509 (Control)
Ye20190510 (All nerve blocks)

Type B Classifier, contact based on tongue marker

Electrode Array locations:
M1U - Rostral Orofacial M1
M1F - Caudal Orofacial M1(more proprioceptive-dependant activity)
S1U - area 1,2
S1F - area 3a/3b(3a is first level processing of for sensory information, Toda, 2002)

Utah arrays @ 1.5mm depth to layers IV and V of neocortex
FM arrays @ 1.0mm depth(?)

KNN Classifier for neuronal spiking
Outline of all the data and options for classifier construction:

Data/features that we suspect are associate with feeding/mouth somatosensation:
- neural activity in four cortical regions over the entire duration of feeding trials, including transport, chewing, swallowing: M1, S1, M1FMA, S1FMA. Channels are fixed but can have 0 to 3 units each. Most have 1.
-Also EMG data? Need to ask about this.

Trained with and used to predict different kinematic behavior:

- Contact event locations on the palate (e.g. 13, 14), and their trajectories
    -Two days for Rocky only, one nerve block, one without. How to include more?
    -Three types of interventions: nerve block with lingual nerve, nerve block without     lingual nerve, and control
- Binary contact between palate and any part of the tongue
- Need to refine teeth contact + not all tongue contact can be inferred, only proximity of markers to oral structures.

We can draw different conclusions based on performance quality (low or high, confidence in prediction) for different combinations of neural activity (which regions/channels/units we exclude and our reasons for doing so) predicting kinematic behavior (contact events, gape, tongue position and deformation, etc.) under the presence of selective knockout of sensation (three types of nerve blocks)

Need to be able to generalize out of single day trials to more. Trained classifiers have little value if they correspond to only one monkey on one day.


M1U Ry20190227 (All nerve blocks) 190 units

M1F Ry20190227 (All nerve blocks) 36 units

S1U Ry20190227 (All nerve blocks) 33 units

S1F Ry20190227 (All nerve blocks) 27 units


Nerve block is:

Inferior alveolar nerve @ mandibular foramen

[Optional] Lingual nerve @ 

Buccal nerve(cheeky!) @ 

Greater palatine nerve @ greater palatine foramen

11/20/2020 meeting

- Run the scripts 10-20 more times to produce confidence intervals
- Run scripts with Yosemite’s duplet
- Could the number of neurons in M1F explain the difference in performance? No, says our analysis.
- JD: there’s lots of cross talk between M1 and S1 could explain the result. Also, a change in kinematics between nerve block and no nerve block.
- Look at receptive fields of neurons, i.e. single unit analysis for receptive fields, 
1. find if a single neuron is highly sensitive as a function of touch events. 
2. correct this to gape cycle 
3. When during contact events/gape cycle does the classifier perform better or worse.

Project Milestones:

SfN, Jan 11th-13th, 2021

Cosyne, Feb 23rd-26th, 2021
