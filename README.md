# KNN-Classifier
[Document link](https://docs.google.com/document/d/1DlprCeO_yCAFl605km7NrIExPAnbm1X1axNf4igZ-3I/edit?usp=sharing)

Timeline:

By Mon, Dec 28th:
Part I:
10 repetitions for the KNN classifier on each of the following areas:

Ry20190228 (Control)

    M1U 11/10 (Done)

    M1F 0/10 (Running)

    S1U

    S1F

Ry20190227 (All nerve blocks)

    M1U

    M1F

    S1U

    S1F
    
Triplet

Ye20190508(No lingual nerve)

    M1U

    M1F

    S1U

    S1F

Ye20190509 (Control)

    M1U

    M1F

    S1U

    S1F
    
Ye20190510 (All nerve blocks)

    M1U

    M1F

    S1U

    S1F
    
Part II:
Some sensitivity critera
    
Mon, Jan 4th: Internal Soft Deadline for Materials

Wed, Jan 6th: Official Hard Deadline for Materials

[ ] One page pdf poster (Mandatory)
[ ] submit your bio
[ ] review presenter and photography/recording guidelines
[ ] acknowledge your presentation time
[ ] review your abstract's hashtags
[ ] 20-min audio accompaniment (Optional)
[ ] QR Code for additional resources (Optional)


Mon, Jan 11th @ 10am CT: Presentation

Barry: only item to consider IF IF you have room, is to insert in parentheses a p value and stats test used to reflect 'better' and /or 'improved"'

Electrode Array locations:
M1U - Rostral Orofacial M1
M1F - Caudal Orofacial M1(more proprioceptive-dependant activity)
S1U - area 1,2
S1F - area 3a/3b(3a is first level processing of )

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

On 10/28/2020, we have the following kinematic data sets complete:
    
Duplet

Ry20190227 (All nerve blocks)

Ry20190228 (Control)
    
Triplet

Ye20190508(No lingual nerve, going through the pipeline, won’t get it in time)

Ye20190509 (Control)
    
Ye20190510 (All nerve blocks)

Sorted NEV -> Eliminate low SNR units (<2) -> Each region

Completed:

Binary palate region contact with the following regions:

M1U Ry20190227 (All nerve blocks) 190 units

M1F Ry20190227 (All nerve blocks) 36 units

S1U Ry20190227 (All nerve blocks) 33 units

S1F Ry20190227 (All nerve blocks) 27 units

To do, by order of importance:

S1U Ry20190227 (All nerve blocks)

S1U Ry20190228 (No nerve blocks)

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
