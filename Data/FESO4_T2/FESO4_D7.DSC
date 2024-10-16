#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	7501
XMIN	1699.800000
XWID	3000.000000
*
*	Documentational Text:
*
TITL	'FESO4_D7'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    10/16/24
TIME    15:31:10
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.32
A1SW    0.30004
MWFQ    9.647725e+09
MWPW    0.0003162
AVGS    1
SPTP    0.008
RCAG    40
RCHM    1
B0MA    0.0004
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    7501
RCTC    0.00128
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     ezcw, 1.0

MetalIonsDef       True
OrgRadicalsDef     False

.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        3200.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     1699.800 G
SweepDirection     Up
SweepWidth         3000.4 G
WidthTM            200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.647725 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        101.8
MarkerPos          0.0000 %
Power              0.3162 mW
PowerAtten         25.00 dB
QValue             3317

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         1
NbScansDone        1
NbScansToDo        1
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       10

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           1.21e-08
ConvTime           8.00 ms
DModDetectSCT      First
DualDetect         OFF
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               40 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             4.000 G
ModFreq            100.00 kHz
ModPhase           0.0
ModResol           10.000
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0.993287,0,0,0.00240032,0,0,-0.00984044,0,0,8.8985e-05,0,0,-4.75288e-06,0,0,-1.05279e-06,0,0,2.9393e-07,0,0,2.70257e-09,0,0,-7.36288e-10,0
QuadMode           False
ResCenter          125 mm
ResLength          25 mm
Resolution         7501
Resonator          1
SctNorm            True
SctRevision        Allegro
SpuExtension       False
SweepTime          60.01 s
TimeConst          1.28 ms
TimeExp            False
TuneCaps           1

*
************************************************************
