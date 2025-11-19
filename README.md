### intralaminarCoherenceAnalysis
Brock Carlson - brock.m.carlson@vanderbilt.edu
MaierLab, Vanderbilt University.

This repository produces the figures for the PNAS (submission).
Raw data is available from the authors upon request (2.5 TB).
Epoched data is available for download with this repository at Zenodo ________ (200 MB). 
Data was recorded from V1 in two awake macaques.


## EXPERIMENTAL MODEL AND STUDY PARTICIPANT DETAILS
Animals were trained for a passive fixation task using positive reinforcement. Once both animals performed
this task satisfactorily (90% trial completion rate), we began electrophysiology recordings. We recorded
32 V1 penetrations from two adult male monkeys (Macaca radiata and macac mulatta) using laminar microelectrodes. Probes had 31 electrode contacts spaced 100um apart (Plexon UProbe).
10 recording sessions were performed in subject ‘‘B’’ (male, age 13) and 21 recording session were performed in subject ‘‘J’’ (male,
age 15). During the experimental period, both animals received their daily fluid ration in the form of juice as
a reward for the behavioral paradigm. All procedures were approved by the Institutional Animal Care and
Use Committee at Vanderbilt University and followed regulations by the National Institutions of Health as
well as the Association for the Assessment and Accreditation of Laboratory Animal Care.
## METHOD DETAILS
# Recording sites and surgical procedures
Recording chambers and headposts were implanted on each monkey in a series of surgeries. Both chambers and headposts were made of custom-designed MRI-compatible plastic (Schmiedt et al., 2014). Sterile
surgical conditions were used for all surgeries. Vital signs were monitored continuously and isoflurane
anesthesia (1.5-2.0%) was maintained. A craniotomy removed the skull over primary visual cortex’s perifoveal visual-field representation. Self-curing denture acrylic (Lang Dental Manufacturing, Wheeling, IL) and
transcranial ceramic screws (Thomas Recording, Gießen, Germany) attached the headpost and recording
chamber to the skull during surgery. The recording chamber was placed over the position of the craniotomy. Each monkey was given antibiotics and analgesics for postsurgical care.
# Visual stimulation
Animals sat, head-fixed, in their personalized, custom-designed experimental chairs. Their chairs were
placed in front of a ViewPIXX monitor running at either 144Hz (resolution 1,289 x 1024 pixels).
Before the start of the experiments, the luminance of the CRT monitors was measured at
17 brightness increments. The CRT monitors were then linearized (‘‘gamma-corrected’’) between the minimum and maximum for each gun (red, green, and blue) using a spectroradiometer (Photoresearch, Syracuse, NY). 
Between the animal and the monitor, a custom-designed dual cold-mirrored stereoscope split
the monitor into a left and right eye view (Carmel, Arcaro, Kastner, & Hasson, 2010; Leopold, Maier, & Logothetis, 2003). 
Images on the right-hand side of the monitor were viewed with the right eye. A black, nonreflective septum divided the two visual fields. 20.5 to 34.5 pixels per degree of visual angle (dva) resulted from
monitor positioning at 46cm-57cm viewing distance. The monocular visual fields were considered to be
aligned when gaze position united across all calibration points on both sides of the monitor (Carlson et al., 2023; Cox et al.,
2019; Dougherty, Cox, Westerberg, & Maier, 2019; Dougherty, Schmid, & Maier, 2019; Westerberg, Cox,
Dougherty, & Maier, 2019). Fusion was facilitated by placing an oval aperture at the edge of each halfscreen.
# Behavioral task
Subjects needed to fixate within 1 dva of a fixation cross for the duration of the trial. If the subject fixated for
the entire trial (1.6s), a juice reward was given. No behavioral responses were required. Infrared light-sensitive
cameras measured gaze position though infrared-transparent (cold) mirrors.28,190 EyeLink II and EyeLink 1000 eye-tracking systems were used to follow the subject’s eye position throughout the experiment.
The temporal resolution of the eye-tracking systems was sampled at 500Hz, giving a 2ms resolution

# 1) laminarLabeling.m
This script plots the LFP, CSD, PSD, Coherence and MUA of our 32 successful V1 penetrations.
The results of this analysis are stored in officialLaminarAssignments.xlsx
- V1: 29 penetrations had clear layer IV/V boundary dynamics. 27 penetrations covered all laminar of V1.
- V2: 10 penetrations had clear V2 laminar presentations.
- 10 penetrations in J occured simultaneously - potentially allowing for furure analysis of horizontal interlaminar communication.



# 2) createTrialTriggeredLFPandMUA.m
This script demonstrates how the the rawData is epoched and sorted by different visual stimulation condition.
A dependency for this code is the laminar assignments
The following 20 unique conditions were sorted for each session
