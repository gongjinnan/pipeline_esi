This toolbox is a pipeline for calculating the inverse and direct problems of MEEG data. It's origins are custom modifications of various matlab based EEG processing toolboxes, mainly Brainstorm, Meth, OpenMEEG and SMP. Currently under development. The main function is pipeline_brainstorm_fs.m and it's description is: 

INPUT:

% -MriFile : Mri Freesufer File Path

% -electrodPos : Electrode Position Matrix

% -nVertices : Number of vertices for the cortex

% -resamplingMethod : Resampling Method options: 'reducepatch' || 'iso2mesh'

% -erodeFactor : Parameter for convexifing the scalp surface, ranges from 0 (min) to 3 (max)(default: 1)

% -fillFactor : Parameter for filling holes in the scalp surface ranges from 0 (min) to 3 (max)(default: 2)

% -headVertices : Number of vertices on the estimated scalp surface(default: 1922)

% -isEEG : 1 for EEG or 0 for MEG (default: 1)

% -conductivity : parameter modelling the conductivity of the skull surface, relative to the conductivities of the scalp and the cortex surfaces (default: 0.0125)

% -NoiseReg : Noise Covariance Regularization (default:0.1)

% -SnrFixed :Signal to noise ratio 1/lambda (default:3)

% -showFigure :Produce Graphics


OUTPUT

%- Gain  :  The calculated LeadField matrix, obtained using the Boundary Element Method implemented by OpenMEEG

%-InverseSolution   : The Imaging Kernel produced by the sLoreta method. Multiply this matrix with the data EEG to obtain the inverse solution.
