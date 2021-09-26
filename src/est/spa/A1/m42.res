Sat Sep 25 08:03:04 CDT 2021
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
$DATA ../../../../data/spa/A1/dat42.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1346.23188919429        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.1081E+01  4.2497E+01 -4.4031E+00  5.3454E+01  8.8491E+01  5.0906E+00 -2.1351E+01  2.6382E+00 -3.9991E+01 -3.6140E+01
            -6.1966E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1523.34871821904        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0520E+00  9.3698E-01  1.0346E+00  1.0448E+00  9.2133E-01  9.6143E-01  1.0108E+00  9.2829E-01  1.1095E+00  9.3434E-01
             2.0348E+00
 PARAMETER:  1.5069E-01  3.4910E-02  1.3400E-01  1.4382E-01  1.8059E-02  6.0671E-02  1.1073E-01  2.5587E-02  2.0389E-01  3.2083E-02
             8.1039E-01
 GRADIENT:   9.7199E+00  2.9369E+01  4.1636E+00  3.1380E+01 -1.1963E+01 -3.7106E+00  3.7614E+00  7.0234E+00  5.5923E+00  7.8156E+00
             1.4119E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1528.62050692582        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0428E+00  7.7596E-01  7.2321E-01  1.1312E+00  7.1680E-01  9.6931E-01  8.3781E-01  4.3688E-01  1.0519E+00  6.4878E-01
             2.0127E+00
 PARAMETER:  1.4187E-01 -1.5366E-01 -2.2406E-01  2.2328E-01 -2.3295E-01  6.8832E-02 -7.6961E-02 -7.2810E-01  1.5061E-01 -3.3266E-01
             7.9950E-01
 GRADIENT:  -1.8564E+01  2.5198E+01 -1.9680E+01  6.9985E+01  3.0048E+01 -3.1446E+00 -6.3944E+00  1.7152E+00 -2.5259E+00 -4.5155E+00
             7.7234E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1537.99584404936        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0509E+00  5.2246E-01  4.6624E-01  1.1854E+00  4.6481E-01  9.8991E-01  1.6056E+00  6.4131E-02  8.9131E-01  5.5159E-01
             1.8469E+00
 PARAMETER:  1.4964E-01 -5.4921E-01 -6.6306E-01  2.7011E-01 -6.6614E-01  8.9860E-02  5.7351E-01 -2.6468E+00 -1.5064E-02 -4.9496E-01
             7.1350E-01
 GRADIENT:  -1.4028E+00  3.3171E+01  2.4539E+01  4.9152E+01 -4.2709E+01  2.4033E+00  2.6558E+00  8.5027E-02 -2.7609E+00  2.1766E+00
            -4.0554E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1538.90826838865        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      321
 NPARAMETR:  1.0534E+00  4.7664E-01  3.7319E-01  1.1465E+00  3.9938E-01  9.9055E-01  1.6224E+00  3.7596E-02  8.8520E-01  4.7675E-01
             1.8276E+00
 PARAMETER:  1.5203E-01 -6.4100E-01 -8.8568E-01  2.3669E-01 -8.1783E-01  9.0506E-02  5.8393E-01 -3.1808E+00 -2.1937E-02 -6.4077E-01
             7.0298E-01
 GRADIENT:  -1.7954E+01  1.4421E+01  8.0979E+00  1.9646E+01 -2.9424E+01 -1.9585E-01  1.1892E+00  3.0391E-02 -2.7488E+00  9.8704E-01
            -1.3000E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1540.50908841089        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      497
 NPARAMETR:  1.0539E+00  3.6729E-01  5.4348E-01  1.2743E+00  4.7875E-01  9.8120E-01  1.9105E+00  2.7417E-02  8.9144E-01  6.4652E-01
             1.8600E+00
 PARAMETER:  1.5250E-01 -9.0161E-01 -5.0975E-01  3.4237E-01 -6.3657E-01  8.1023E-02  7.4736E-01 -3.4966E+00 -1.4922E-02 -3.3615E-01
             7.2056E-01
 GRADIENT:  -4.1380E+00  1.0943E+01  1.4434E+01  1.8454E+01 -2.5957E+01  4.3024E-01 -3.4283E-01  1.5824E-02  1.1651E+00  3.7089E+00
            -4.4455E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1541.99115651338        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      672
 NPARAMETR:  1.0509E+00  2.0679E-01  5.5011E-01  1.3411E+00  4.5708E-01  9.7721E-01  2.6890E+00  1.0000E-02  8.7171E-01  6.4639E-01
             1.8557E+00
 PARAMETER:  1.4968E-01 -1.4761E+00 -4.9764E-01  3.9349E-01 -6.8290E-01  7.6944E-02  1.0892E+00 -5.1281E+00 -3.7304E-02 -3.3635E-01
             7.1824E-01
 GRADIENT:   1.9878E+00 -3.5866E-01  6.1815E+00  7.4722E+00 -8.5892E+00  5.4762E-01 -2.9498E+00  0.0000E+00  2.8133E+00 -1.3375E-01
             9.4648E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1542.33853219818        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.0468E+00  1.2237E-01  5.3540E-01  1.3729E+00  4.3624E-01  9.6937E-01  3.8787E+00  1.0000E-02  8.5180E-01  6.4131E-01
             1.8384E+00
 PARAMETER:  1.4570E-01 -2.0007E+00 -5.2475E-01  4.1691E-01 -7.2955E-01  6.8889E-02  1.4555E+00 -6.8047E+00 -6.0402E-02 -3.4423E-01
             7.0889E-01
 GRADIENT:   1.0240E+00  3.0390E+00  1.2710E+01  1.3584E+01 -1.9219E+01 -1.8026E+00  2.1293E+00  0.0000E+00 -2.7131E+00  1.0553E+00
            -1.6410E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1542.63610566114        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  1.0477E+00  9.9893E-02  5.0111E-01  1.3589E+00  4.1671E-01  9.7396E-01  4.3192E+00  1.0000E-02  8.5811E-01  6.1426E-01
             1.8396E+00
 PARAMETER:  1.4655E-01 -2.2037E+00 -5.9093E-01  4.0668E-01 -7.7537E-01  7.3616E-02  1.5631E+00 -7.5675E+00 -5.3022E-02 -3.8734E-01
             7.0953E-01
 GRADIENT:   7.8929E+00  1.6567E+01 -1.3931E+01 -1.9606E+01  2.0416E+01  2.6899E-01  2.7483E+01  0.0000E+00 -1.3179E+01 -7.4203E+00
            -5.2538E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1542.66136946150        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  1.0455E+00  9.6361E-02  5.0795E-01  1.3631E+00  4.2061E-01  9.7510E-01  4.4064E+00  1.0000E-02  8.5234E-01  6.1645E-01
             1.8389E+00
 PARAMETER:  1.4445E-01 -2.2397E+00 -5.7738E-01  4.0976E-01 -7.6605E-01  7.4785E-02  1.5831E+00 -7.6539E+00 -5.9771E-02 -3.8379E-01
             7.0916E-01
 GRADIENT:   7.1903E-01  2.2310E+00 -1.9821E+00 -1.5550E+00  2.9704E+00  2.3770E-01  3.5460E+00  0.0000E+00 -1.9675E+00 -1.0423E+00
            -8.6081E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1542.66554334147        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1394
 NPARAMETR:  1.0451E+00  9.2220E-02  5.0781E-01  1.3639E+00  4.1999E-01  9.7473E-01  4.5193E+00  1.0000E-02  8.5100E-01  6.1693E-01
             1.8390E+00
 PARAMETER:  1.4419E-01 -2.2833E+00 -5.7770E-01  4.1033E-01 -7.6742E-01  7.4435E-02  1.6086E+00 -7.7997E+00 -6.1445E-02 -3.8317E-01
             7.0912E-01
 GRADIENT:   2.0012E-01  2.8239E+02 -1.1135E+00 -2.3863E+00  8.4530E+02  2.3185E-02  4.0003E+02  0.0000E+00 -7.9884E-01 -3.2377E-01
            -3.3753E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1394
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6328E-04  3.5772E-02 -1.5782E-04 -1.7927E-02  3.3651E-03
 SE:             2.9465E-02  1.3936E-02  2.3599E-04  2.7106E-02  2.0889E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7392E-01  1.0262E-02  5.0366E-01  5.0837E-01  8.7202E-01

 ETASHRINKSD(%)  1.2890E+00  5.3313E+01  9.9209E+01  9.1927E+00  3.0020E+01
 ETASHRINKVR(%)  2.5614E+00  7.8203E+01  9.9994E+01  1.7540E+01  5.1028E+01
 EBVSHRINKSD(%)  1.3872E+00  6.5897E+01  9.9117E+01  7.2425E+00  2.6230E+01
 EBVSHRINKVR(%)  2.7552E+00  8.8370E+01  9.9992E+01  1.3961E+01  4.5580E+01
 RELATIVEINF(%)  9.6308E+01  5.7203E+00  3.5833E-04  3.6541E+01  2.5348E+00
 EPSSHRINKSD(%)  3.7812E+01
 EPSSHRINKVR(%)  6.1327E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1542.6655433414680     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -807.51471677772986     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.22
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1542.666       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  9.22E-02  5.08E-01  1.36E+00  4.20E-01  9.75E-01  4.52E+00  1.00E-02  8.51E-01  6.17E-01  1.84E+00
 


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
+        1.14E+03
 
 TH 2
+       -9.21E+02  3.63E+06
 
 TH 3
+       -1.13E+03  7.86E+03  1.89E+06
 
 TH 4
+       -1.03E+03  4.43E+03  9.89E+05  5.19E+05
 
 TH 5
+       -6.69E+02 -1.64E+04 -1.72E+06 -8.97E+05  1.57E+06
 
 TH 6
+        2.04E+01 -2.39E+02 -2.17E+02 -2.12E+02 -1.79E+02  2.15E+02
 
 TH 7
+       -2.24E+01  1.05E+05  1.91E+02  1.07E+02 -4.03E+02 -6.87E+00  3.05E+03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.33E+02  4.53E+03  6.22E+03  6.36E+03  3.85E+03 -1.06E+02  1.21E+02  0.00E+00  3.10E+03
 
 TH10
+       -3.65E+02  2.29E+03  4.68E+03  1.23E+06  2.06E+03 -7.28E+01  5.73E+01  0.00E+00  1.90E+03  1.45E+03
 
 TH11
+       -1.09E+02  7.16E+02  4.25E+05  2.23E+05  6.25E+02 -1.88E+01  1.89E+01  0.00E+00  5.28E+02  3.89E+02  1.74E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.713
Stop Time:
Sat Sep 25 08:03:32 CDT 2021
