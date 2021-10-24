Sun Oct 24 03:52:47 CDT 2021
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
$DATA ../../../../data/SD4/TD1/dat74.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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
$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1682.56623932299        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1742E+02 -3.1180E+01 -7.3488E+00 -2.4526E+01 -7.6344E+00  5.1504E+01 -8.4436E+00  7.7112E+00  3.2292E+00  8.8892E+00
             2.0624E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1687.83726271043        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9797E-01  1.0684E+00  1.0646E+00  1.0234E+00  1.0547E+00  9.5473E-01  1.0549E+00  9.6842E-01  9.9423E-01  9.7760E-01
             9.3880E-01
 PARAMETER:  9.7971E-02  1.6613E-01  1.6261E-01  1.2316E-01  1.5324E-01  5.3669E-02  1.5349E-01  6.7910E-02  9.4217E-02  7.7346E-02
             3.6852E-02
 GRADIENT:   1.0286E+00  1.3389E+00  4.6272E+00 -7.9240E-01 -1.0856E+01 -4.8972E+00 -7.4836E+00 -1.1965E-02 -3.2375E-01 -5.1179E+00
            -1.0428E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1688.54182396789        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  9.9975E-01  1.0909E+00  1.1003E+00  1.0123E+00  1.0988E+00  9.6315E-01  1.1498E+00  9.5667E-01  9.6958E-01  1.0431E+00
             9.5262E-01
 PARAMETER:  9.9746E-02  1.8701E-01  1.9561E-01  1.1220E-01  1.9424E-01  6.2449E-02  2.3962E-01  5.5708E-02  6.9105E-02  1.4219E-01
             5.1457E-02
 GRADIENT:   5.5813E+00  2.7078E+00  2.3923E+00 -3.3342E-01 -6.5719E-01 -1.2630E+00  2.0631E-02 -1.2639E+00 -1.1005E+00 -9.1902E-01
            -4.0316E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1688.75396458303        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.9503E-01  8.8578E-01  1.2530E+00  1.1481E+00  1.0720E+00  9.6550E-01  1.3087E+00  1.0503E+00  9.0279E-01  1.0437E+00
             9.6284E-01
 PARAMETER:  9.5019E-02 -2.1286E-02  3.2551E-01  2.3815E-01  1.6950E-01  6.4895E-02  3.6902E-01  1.4904E-01 -2.2657E-03  1.4276E-01
             6.2133E-02
 GRADIENT:  -1.0239E+00  3.0634E+00  8.8738E-01  2.2178E+00 -1.3662E+00  5.6969E-01  1.8392E-01 -2.3517E-01 -6.2954E-01 -6.0041E-01
             1.0520E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1688.78741876626        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.9447E-01  7.8355E-01  1.3319E+00  1.2196E+00  1.0614E+00  9.6448E-01  1.4160E+00  1.0850E+00  8.7291E-01  1.0487E+00
             9.6226E-01
 PARAMETER:  9.4458E-02 -1.4391E-01  3.8659E-01  2.9852E-01  1.5957E-01  6.3829E-02  4.4780E-01  1.8159E-01 -3.5920E-02  1.4755E-01
             6.1525E-02
 GRADIENT:   4.8695E-01  6.4448E+00  2.0137E+00  9.3634E+00 -2.8144E+00  6.1931E-01  3.6425E-01 -6.3506E-01 -8.2514E-01 -7.3532E-01
             6.2402E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1688.86815151775        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  9.9312E-01  6.4115E-01  1.4618E+00  1.3138E+00  1.0601E+00  9.6123E-01  1.5792E+00  1.1715E+00  8.4241E-01  1.0710E+00
             9.5937E-01
 PARAMETER:  9.3101E-02 -3.4449E-01  4.7967E-01  3.7296E-01  1.5834E-01  6.0460E-02  5.5692E-01  2.5826E-01 -7.1484E-02  1.6859E-01
             5.8526E-02
 GRADIENT:   1.8195E+00  5.6671E+00  1.7807E+00  1.0543E+01 -2.6408E+00  1.7439E-02  2.8646E-01 -5.1186E-01 -3.3090E-01 -1.4750E-01
            -5.9446E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1688.93375003771        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  9.9111E-01  4.7892E-01  1.6115E+00  1.4208E+00  1.0596E+00  9.5787E-01  1.8401E+00  1.2828E+00  8.1092E-01  1.0939E+00
             9.5738E-01
 PARAMETER:  9.1067E-02 -6.3623E-01  5.7717E-01  4.5122E-01  1.5790E-01  5.6958E-02  7.0983E-01  3.4905E-01 -1.0959E-01  1.8977E-01
             5.6443E-02
 GRADIENT:   2.4255E+00  4.6715E+00  1.2953E+00  1.1081E+01 -2.2658E+00 -4.9399E-01  1.2018E-01 -3.4545E-01  7.3644E-02  4.3690E-01
            -1.4767E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1688.96058916716        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.8926E-01  3.6133E-01  1.7248E+00  1.4992E+00  1.0594E+00  9.5595E-01  2.1324E+00  1.3776E+00  7.8812E-01  1.1049E+00
             9.5725E-01
 PARAMETER:  8.9197E-02 -9.1797E-01  6.4510E-01  5.0491E-01  1.5769E-01  5.4950E-02  8.5726E-01  4.2032E-01 -1.3811E-01  1.9976E-01
             5.6309E-02
 GRADIENT:   1.9718E+00  4.1677E+00  1.0637E+00  1.2717E+01 -2.3971E+00 -6.1674E-01 -6.2900E-02 -2.6934E-01 -1.7242E-01  5.8660E-01
            -1.5647E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1688.96972520767        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  9.8800E-01  2.9060E-01  1.7971E+00  1.5464E+00  1.0599E+00  9.5483E-01  2.3924E+00  1.4413E+00  7.7408E-01  1.1091E+00
             9.5747E-01
 PARAMETER:  8.7926E-02 -1.1358E+00  6.8617E-01  5.3595E-01  1.5814E-01  5.3780E-02  9.7231E-01  4.6555E-01 -1.5608E-01  2.0353E-01
             5.6542E-02
 GRADIENT:   1.3572E+00  3.6433E+00  9.5105E-01  1.3920E+01 -2.5406E+00 -6.4469E-01 -3.7123E-01 -2.4632E-01 -5.6465E-01  6.0857E-01
            -1.5092E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1689.11132475293        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1585
 NPARAMETR:  9.8515E-01  1.5587E-01  1.8737E+00  1.6363E+00  1.0402E+00  9.5223E-01  3.4848E+00  1.5106E+00  7.4460E-01  1.0872E+00
             9.5817E-01
 PARAMETER:  8.5034E-02 -1.7587E+00  7.2791E-01  5.9245E-01  1.3945E-01  5.1052E-02  1.3484E+00  5.1254E-01 -1.9491E-01  1.8356E-01
             5.7268E-02
 GRADIENT:  -7.9640E-01  3.6312E+00  1.2387E+00  2.4264E+01 -4.1835E+00 -8.2817E-01  1.3050E-01 -3.5964E-01 -3.6616E+00  2.7145E-01
            -1.6301E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1689.41656902894        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1763
 NPARAMETR:  9.8286E-01  7.3091E-02  1.9634E+00  1.6929E+00  1.0406E+00  9.5107E-01  4.9813E+00  1.6052E+00  7.4487E-01  1.0840E+00
             9.6148E-01
 PARAMETER:  8.2707E-02 -2.5161E+00  7.7466E-01  6.2641E-01  1.3979E-01  4.9832E-02  1.7057E+00  5.7327E-01 -1.9455E-01  1.8067E-01
             6.0719E-02
 GRADIENT:  -3.6663E+00  8.6789E-01  1.0694E+00  3.2515E+01 -6.3121E+00 -7.9444E-01 -2.5759E+00  8.8154E-01  3.8958E+00  1.0020E+00
             9.8738E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1690.11186500903        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1942            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8347E-01  3.7408E-02  1.9368E+00  1.6974E+00  1.0291E+00  9.5185E-01  6.5037E+00  1.5734E+00  7.3400E-01  1.0702E+00
             9.6144E-01
 PARAMETER:  8.3333E-02 -3.1859E+00  7.6103E-01  6.2908E-01  1.2872E-01  5.0649E-02  1.9724E+00  5.5324E-01 -2.0925E-01  1.6787E-01
             6.0676E-02
 GRADIENT:   4.4217E+02  9.4688E+00  9.6237E+00  1.2669E+03  7.0266E+00  3.9728E+01  3.4333E+01  2.2807E+00  2.9221E+01  1.6079E+00
             8.2606E-01

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1690.29549877921        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     2074
 NPARAMETR:  9.8416E-01  3.7835E-02  1.9350E+00  1.6975E+00  1.0283E+00  9.5222E-01  6.6943E+00  1.5711E+00  7.3187E-01  1.0693E+00
             9.6161E-01
 PARAMETER:  8.3032E-02 -3.1658E+00  7.5988E-01  6.2799E-01  1.2876E-01  5.0774E-02  2.0062E+00  5.5165E-01 -2.1277E-01  1.6590E-01
             6.0705E-02
 GRADIENT:  -2.1190E+00  1.6270E+01 -1.2586E-01 -6.3436E+01  9.0665E-01 -1.8156E-01  2.3109E+01 -2.0180E-02 -1.2989E+02 -1.7408E-01
            -1.0912E-01
 NUMSIGDIG:         1.8         2.3         3.3         2.5         1.9         2.3         2.4         3.5         2.3         1.9
                    2.6

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2074
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2150E-03  7.2488E-03 -3.5950E-02 -1.6102E-02 -4.8078E-02
 SE:             2.9849E-02  6.8939E-03  1.8350E-02  2.8802E-02  1.9891E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6753E-01  2.9304E-01  5.0102E-02  5.7612E-01  1.5647E-02

 ETASHRINKSD(%)  2.4291E-03  7.6905E+01  3.8525E+01  3.5113E+00  3.3362E+01
 ETASHRINKVR(%)  4.8581E-03  9.4666E+01  6.2208E+01  6.8993E+00  5.5594E+01
 EBVSHRINKSD(%)  4.2719E-01  8.3199E+01  4.1908E+01  3.5293E+00  2.9390E+01
 EBVSHRINKVR(%)  8.5256E-01  9.7177E+01  6.6253E+01  6.9340E+00  5.0142E+01
 RELATIVEINF(%)  9.8883E+01  1.0040E+00  9.8010E+00  3.4862E+01  1.3993E+01
 EPSSHRINKSD(%)  4.5048E+01
 EPSSHRINKVR(%)  6.9803E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1690.2954987792104     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -955.14467221547227     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1690.295       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  3.82E-02  1.93E+00  1.70E+00  1.03E+00  9.52E-01  6.73E+00  1.57E+00  7.31E-01  1.07E+00  9.61E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       63.781
Stop Time:
Sun Oct 24 03:53:00 CDT 2021
