Sat Sep 18 09:48:30 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/A2/dat34.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m34.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1300.19355738295        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4606E+02 -1.0676E+02 -9.7191E+00 -1.2094E+02  9.3839E+01 -4.0473E+01 -1.1212E+01 -5.8029E+00  1.0226E+01 -5.3077E+01
            -5.5804E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1457.22729961883        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.1758E-01  1.1053E+00  1.0150E+00  1.0712E+00  9.8119E-01  1.0280E+00  9.6461E-01  9.8129E-01  8.0100E-01  1.0639E+00
             1.9156E+00
 PARAMETER:  1.3987E-02  2.0014E-01  1.1487E-01  1.6878E-01  8.1016E-02  1.2760E-01  6.3973E-02  8.1110E-02 -1.2190E-01  1.6198E-01
             7.5004E-01
 GRADIENT:  -8.8205E+01  3.2772E+01 -6.9482E+00  5.7611E+01  4.5271E+00 -2.0410E+01 -6.8698E+00  1.1733E+00 -3.6279E+00 -3.3247E+00
             6.2735E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1458.47101880981        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.2828E-01  9.4562E-01  1.1151E+00  1.1828E+00  9.4543E-01  1.0484E+00  1.2538E+00  8.3353E-01  6.9317E-01  1.0312E+00
             1.9323E+00
 PARAMETER:  2.5582E-02  4.4086E-02  2.0896E-01  2.6790E-01  4.3885E-02  1.4724E-01  3.2617E-01 -8.2086E-02 -2.6648E-01  1.3077E-01
             7.5869E-01
 GRADIENT:  -5.9217E+01  4.3556E+01  5.4295E+00  7.6788E+01 -7.8745E+00 -9.4638E+00 -8.1728E-01 -1.1209E+00 -5.1817E+00 -6.5736E+00
            -4.2276E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1461.49552067981        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.5826E-01  9.4530E-01  8.6324E-01  1.1246E+00  8.3745E-01  1.0690E+00  1.1960E+00  4.1751E-01  7.3358E-01  9.9948E-01
             1.9092E+00
 PARAMETER:  5.7368E-02  4.3752E-02 -4.7059E-02  2.1739E-01 -7.7393E-02  1.6673E-01  2.7899E-01 -7.7344E-01 -2.0981E-01  9.9481E-02
             7.4666E-01
 GRADIENT:   2.5727E+00 -3.9858E-01 -1.0382E+00 -9.1636E-01  1.2202E+00  3.5286E-01 -3.1569E-01  2.7028E-01 -1.9565E-01  3.7469E-01
             6.5671E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1462.12871824848        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.5899E-01  7.4353E-01  8.6381E-01  1.2456E+00  7.5696E-01  1.0694E+00  1.4423E+00  9.2441E-02  6.7940E-01  9.7553E-01
             1.9027E+00
 PARAMETER:  5.8121E-02 -1.9635E-01 -4.6401E-02  3.1961E-01 -1.7844E-01  1.6713E-01  4.6627E-01 -2.2812E+00 -2.8654E-01  7.5230E-02
             7.4326E-01
 GRADIENT:   5.9321E+00  4.7008E+00  1.8636E+00  8.5302E+00 -3.6984E+00  9.0405E-01 -2.8066E-01  2.2076E-02 -2.7314E+00 -1.9336E-01
            -1.6693E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1462.42562491573        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  9.5655E-01  6.3465E-01  8.9940E-01  1.3166E+00  7.4005E-01  1.0677E+00  1.6211E+00  2.4462E-02  6.6508E-01  9.8775E-01
             1.8999E+00
 PARAMETER:  5.5578E-02 -3.5468E-01 -6.0317E-03  3.7502E-01 -2.0104E-01  1.6547E-01  5.8310E-01 -3.6106E+00 -3.0785E-01  8.7677E-02
             7.4181E-01
 GRADIENT:   3.5378E+00  7.2848E+00  5.2987E-01  2.1779E+01 -2.8779E+00  6.1965E-01  1.5147E-01  1.9869E-03 -2.9053E+00 -5.2020E-01
            -2.0879E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1462.57885948623        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      454
 NPARAMETR:  9.5456E-01  5.9896E-01  9.1329E-01  1.3310E+00  7.3935E-01  1.0655E+00  1.6758E+00  1.2335E-02  6.6921E-01  9.8866E-01
             1.9000E+00
 PARAMETER:  5.3498E-02 -4.1257E-01  9.3008E-03  3.8592E-01 -2.0198E-01  1.6342E-01  6.1627E-01 -4.2953E+00 -3.0165E-01  8.8596E-02
             7.4185E-01
 GRADIENT:   1.0791E+00  1.7440E+00 -1.2863E+00  7.3340E+00  1.9973E+00  1.5411E-01  1.8110E-01  5.4429E-04 -1.1455E+00 -1.2684E+00
            -1.4280E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1463.23196792660        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      633
 NPARAMETR:  9.5484E-01  4.1697E-01  1.0848E+00  1.4604E+00  7.7298E-01  1.0673E+00  2.1308E+00  1.0000E-02  6.4973E-01  1.0748E+00
             1.9127E+00
 PARAMETER:  5.3787E-02 -7.7473E-01  1.8136E-01  4.7870E-01 -1.5750E-01  1.6510E-01  8.5651E-01 -8.2243E+00 -3.3120E-01  1.7215E-01
             7.4851E-01
 GRADIENT:  -8.5292E-01  4.0640E+00 -7.0064E-01  1.4754E+01  5.4905E-01  6.3636E-01  4.7571E-01  0.0000E+00 -1.3386E+00 -1.2736E+00
            -4.2991E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1463.89899530372        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      810
 NPARAMETR:  9.5314E-01  2.0056E-01  1.1375E+00  1.5870E+00  7.3878E-01  1.0636E+00  3.1880E+00  1.0000E-02  6.3755E-01  1.1131E+00
             1.9056E+00
 PARAMETER:  5.2004E-02 -1.5067E+00  2.2882E-01  5.6186E-01 -2.0276E-01  1.6165E-01  1.2594E+00 -1.8007E+01 -3.5012E-01  2.0714E-01
             7.4478E-01
 GRADIENT:   3.2620E+00  2.2918E+00  2.0277E+00  1.2076E+01 -4.7491E+00  4.9711E-01  1.1412E+00  0.0000E+00 -8.5304E-01 -4.1318E-01
            -2.8685E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1464.42647799488        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      987
 NPARAMETR:  9.4861E-01  9.0234E-02  1.2697E+00  1.6661E+00  7.6848E-01  1.0572E+00  4.4302E+00  1.0000E-02  6.3091E-01  1.1946E+00
             1.9032E+00
 PARAMETER:  4.7245E-02 -2.3053E+00  3.3878E-01  6.1050E-01 -1.6334E-01  1.5560E-01  1.5884E+00 -2.9626E+01 -3.6059E-01  2.7779E-01
             7.4356E-01
 GRADIENT:  -1.7885E+00 -1.7002E-01  2.9644E+00  3.1226E+01 -8.3464E+00 -1.1944E+00 -5.6477E+00  0.0000E+00  2.4589E+00  3.3947E+00
            -6.9312E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1467.15247660766        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1166
 NPARAMETR:  9.4909E-01  1.0000E-02  1.3420E+00  1.7048E+00  7.8887E-01  1.0607E+00  1.1704E+01  1.0000E-02  6.0433E-01  1.1900E+00
             1.9158E+00
 PARAMETER:  4.7744E-02 -4.5665E+00  3.9414E-01  6.3347E-01 -1.3715E-01  1.5892E-01  2.5600E+00 -6.4389E+01 -4.0364E-01  2.7398E-01
             7.5012E-01
 GRADIENT:   2.3417E+00  0.0000E+00 -2.3479E+00  1.1346E+01  4.2943E+00  6.2296E-01 -4.4339E+00  0.0000E+00  5.9107E-01  1.1046E+00
             3.8140E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1467.42487816655        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1344
 NPARAMETR:  9.4869E-01  1.0000E-02  1.2654E+00  1.6950E+00  7.5712E-01  1.0597E+00  1.2236E+01  1.0000E-02  6.0784E-01  1.1615E+00
             1.9095E+00
 PARAMETER:  4.7322E-02 -4.6760E+00  3.3536E-01  6.2770E-01 -1.7823E-01  1.5801E-01  2.6044E+00 -6.6165E+01 -3.9784E-01  2.4973E-01
             7.4682E-01
 GRADIENT:   1.7930E+00  0.0000E+00 -1.3938E+00  3.7358E+00  2.6594E+00  3.0040E-01 -2.3508E+00  0.0000E+00  4.4807E-01  2.9573E-01
             2.0198E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1467.44322678631        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1506
 NPARAMETR:  9.4779E-01  1.0000E-02  1.2243E+00  1.6902E+00  7.3818E-01  1.0589E+00  1.2240E+01  1.0000E-02  6.1140E-01  1.1485E+00
             1.9054E+00
 PARAMETER:  4.6379E-02 -4.6762E+00  3.0235E-01  6.2487E-01 -2.0356E-01  1.5723E-01  2.6047E+00 -6.6196E+01 -3.9200E-01  2.3848E-01
             7.4471E-01
 GRADIENT:   2.7112E-02  0.0000E+00 -1.3999E-01 -1.0872E+00  2.6313E-01 -1.5929E-03  7.1901E-01  0.0000E+00 -4.5013E-01 -1.6778E-01
            -7.9399E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1506
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9849E-04  1.0570E-02 -3.4641E-05 -1.6958E-02 -2.7880E-02
 SE:             2.9419E-02  6.3613E-03  1.2490E-04  2.6409E-02  2.2887E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8919E-01  9.6594E-02  7.8152E-01  5.2079E-01  2.2318E-01

 ETASHRINKSD(%)  1.4434E+00  7.8689E+01  9.9582E+01  1.1525E+01  2.3324E+01
 ETASHRINKVR(%)  2.8659E+00  9.5458E+01  9.9998E+01  2.1722E+01  4.1209E+01
 EBVSHRINKSD(%)  1.2846E+00  8.4139E+01  9.9551E+01  1.0535E+01  2.1414E+01
 EBVSHRINKVR(%)  2.5527E+00  9.7484E+01  9.9998E+01  1.9961E+01  3.8242E+01
 RELATIVEINF(%)  9.7309E+01  1.9894E+00  1.3343E-04  3.6366E+01  4.1438E+00
 EPSSHRINKSD(%)  3.5356E+01
 EPSSHRINKVR(%)  5.8211E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1467.4432267863099     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -732.29240022257170     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.54
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1467.443       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.48E-01  1.00E-02  1.22E+00  1.69E+00  7.38E-01  1.06E+00  1.22E+01  1.00E-02  6.11E-01  1.15E+00  1.91E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.11E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.94E+00  0.00E+00  1.61E+06
 
 TH 4
+       -4.86E+01  0.00E+00 -2.79E+02  1.98E+05
 
 TH 5
+        2.33E+00  0.00E+00 -3.97E+06  7.19E+02  9.78E+06
 
 TH 6
+        6.49E+00  0.00E+00  6.31E+00  1.03E+01  3.19E+01  1.72E+02
 
 TH 7
+        6.41E-01  0.00E+00  6.72E+00  3.68E+01 -2.42E+01 -6.86E-01  2.17E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.90E+01  0.00E+00 -1.14E+03 -9.72E+02  4.11E+03  1.02E+02  2.64E+01  0.00E+00  3.84E+06
 
 TH10
+       -5.01E-01  0.00E+00  2.18E+06 -1.95E+02 -5.37E+06  9.15E+00  4.98E+00  0.00E+00 -1.10E+03  2.94E+06
 
 TH11
+       -9.61E+00  0.00E+00  4.20E+05 -7.07E+01 -1.04E+06 -2.33E+00  1.55E+00  0.00E+00 -2.82E+02  5.68E+05  1.10E+05
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       24.457
Stop Time:
Sat Sep 18 09:48:56 CDT 2021
