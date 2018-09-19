# PIPELINE_ESI

README
INPUT:

% -MriFile               : Mri Freesufer File Path

% -electrodPos           : Electrode Position Matrix

% -nVertices             : Number of vertices for the cortex

% -resamplingMethod      : Resampling Method options: 'reducepatch' ||  'iso2mesh'

% -erodeFactor           : Parameter for convexifing the scalp surface, ranges from 0 (min) to 3 (max)(default: 1)

% -fillFactor            : Parameter for filling holes in the scalp surface ranges from 0 (min) to 3 (max)(default: 2)

% -headVertices          : Number of vertices on the estimated scalp surface(default: 1922)

% -isEEG                 : 1 for EEG or 0 for MEG (default: 1)

% -conductivity          : parameter modelling the conductivity of the skull surface, relative to the conductivities of the scalp and the cortex surfaces (default: 0.0125)

% -NoiseReg              : Noise Covariance Regularization (default:0.1)

% -SnrFixed              :Signal to noise ratio 1/lambda (default:3)

% -showFigure            :Produce Graphics

%End Pipeline ESI Readme
