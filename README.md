# SFND_Radar_Target_Generation_and_Detection
Udacity Sensor Fusion nanodegree radar target generation final project

1. Implementation steps for the 2D CFAR process.
The 2D CFAR process consists of creating a sliding window through the dataset.  For each point that is far from the edges of the dataset, we take a patch of size (2 * training_r + 2 * guard_r + 1) x (2 * training_d + 2 * guard_d + 1) centered on that point; we mask out the central (2 * training_r + 1) x (2 * training_r + 1) points, and compute the mean noise level on the rest.  Note that we convert the RDM (in dB) to power in order to compute the average noise level, and convert that back to dB for storage.

2. Selection of Training, Guard cells and offset.
The number of training cells was set at 5 in both directions.  The number of guard cells was set as 3 in the range direction, and 5 in the velocity direction, because there seemed to be more ambiguity about velocity and using too few guard cells in that direction resulted in two separate velocities detected instead of a single batch.  The signal/noise threshold was set at 10dB somewhat arbitrarily, since there was no explicit noise in the dataset, but the main signal seemed to be about 15dB higher than the other frequency content.

3. Steps taken to suppress the non-thresholded cells at the edges.
The 0/1 cfar matrix was initialized as all-zero, and only filled in for the region where we had the full training and guard data.