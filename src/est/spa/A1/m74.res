Wed Sep 29 12:19:52 CDT 2021
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
$DATA ../../../../data/spa/A1/dat74.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1543.92500999939        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2575E+02 -1.6566E+01  2.3330E+01 -3.8745E+01 -1.5115E+00  5.2200E+01 -1.3809E+01 -4.0388E+00 -1.5184E+01  1.6318E+00
            -2.0163E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1577.25536616988        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0158E+00  1.0759E+00  9.8392E-01  1.0101E+00  1.0058E+00  9.0593E-01  1.1116E+00  9.6972E-01  1.0254E+00  8.8331E-01
             1.5576E+00
 PARAMETER:  1.1569E-01  1.7317E-01  8.3791E-02  1.1000E-01  1.0582E-01  1.2068E-03  2.0581E-01  6.9255E-02  1.2505E-01 -2.4084E-02
             5.4315E-01
 GRADIENT:   2.4805E+02  2.4269E+01  1.5221E+01  5.1612E+00 -1.6972E+01 -2.7124E+00  2.5311E+00  9.2074E-01  4.6534E+00  3.4224E+00
             4.0871E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1579.78472979466        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      180
 NPARAMETR:  9.9989E-01  8.9240E-01  6.2671E-01  1.1068E+00  7.1034E-01  9.1562E-01  1.4082E+00  7.6159E-01  8.8599E-01  4.6506E-01
             1.5000E+00
 PARAMETER:  9.9891E-02 -1.3839E-02 -3.6727E-01  2.0151E-01 -2.4202E-01  1.1843E-02  4.4235E-01 -1.7235E-01 -2.1048E-02 -6.6560E-01
             5.0545E-01
 GRADIENT:   6.1948E+00  1.4142E+01 -9.5392E+00  2.1336E+01  3.7302E+00 -1.2224E+01  8.4489E+00  3.4837E+00  9.3492E-01 -1.2434E+00
             2.6504E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1581.99358071419        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      357
 NPARAMETR:  9.8723E-01  8.6610E-01  6.8038E-01  1.1242E+00  7.5088E-01  9.3886E-01  1.3779E+00  5.1358E-01  9.0026E-01  6.5472E-01
             1.4050E+00
 PARAMETER:  8.7147E-02 -4.3756E-02 -2.8510E-01  2.1708E-01 -1.8651E-01  3.6915E-02  4.2058E-01 -5.6636E-01 -5.0761E-03 -3.2355E-01
             4.4001E-01
 GRADIENT:  -2.1602E+01 -9.2436E-01 -1.7649E+01  1.7565E+01  2.5253E+01 -2.8425E+00  1.3933E+00 -5.9750E-01 -1.0643E-01 -1.3830E+00
             1.2128E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1583.08382941609        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  9.9414E-01  6.6230E-01  8.2205E-01  1.2622E+00  7.4608E-01  9.4171E-01  1.6334E+00  6.9294E-01  8.5516E-01  7.0570E-01
             1.3972E+00
 PARAMETER:  9.4127E-02 -3.1203E-01 -9.5948E-02  3.3285E-01 -1.9292E-01  3.9944E-02  5.9068E-01 -2.6681E-01 -5.6462E-02 -2.4857E-01
             4.3446E-01
 GRADIENT:   3.7629E+00  8.6710E+00  4.2244E+00  1.3009E+01 -9.6368E+00 -1.1745E-01 -1.2416E+00  2.9922E-01 -1.8883E+00  9.3781E-01
            -3.5144E-04

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1583.78851600660        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  9.9019E-01  4.5299E-01  9.1918E-01  1.3827E+00  7.3564E-01  9.4073E-01  2.0928E+00  7.7277E-01  8.2512E-01  7.2029E-01
             1.3939E+00
 PARAMETER:  9.0145E-02 -6.9188E-01  1.5730E-02  4.2403E-01 -2.0701E-01  3.8897E-02  8.3852E-01 -1.5777E-01 -9.2231E-02 -2.2810E-01
             4.3212E-01
 GRADIENT:   4.8276E+00  3.5245E+00  8.4714E+00 -1.8586E+00 -1.0180E+01  1.2718E+00 -1.1992E-01 -4.1507E-02 -6.9031E-01 -5.9437E-01
            -1.3486E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1584.37271706866        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.8404E-01  2.4086E-01  9.4580E-01  1.5252E+00  6.9506E-01  9.3950E-01  2.9251E+00  7.6256E-01  8.0634E-01  7.6121E-01
             1.3916E+00
 PARAMETER:  8.3916E-02 -1.3235E+00  4.4272E-02  5.2215E-01 -2.6376E-01  3.7588E-02  1.1733E+00 -1.7108E-01 -1.1525E-01 -1.7285E-01
             4.3049E-01
 GRADIENT:   1.0012E+00  6.2855E+00  9.1587E+00  3.8247E+01 -1.7666E+01  2.6228E+00 -6.8119E-01 -4.9613E-01  6.2185E-01  8.1872E-01
            -1.6581E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1585.73216974489        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.7942E-01  9.4744E-02  9.8314E-01  1.5880E+00  6.8515E-01  9.3071E-01  4.4343E+00  8.5478E-01  7.8803E-01  7.4123E-01
             1.3957E+00
 PARAMETER:  7.9205E-02 -2.2566E+00  8.2992E-02  5.6245E-01 -2.7811E-01  2.8190E-02  1.5894E+00 -5.6914E-02 -1.3821E-01 -1.9945E-01
             4.3340E-01
 GRADIENT:  -8.7407E-01 -7.5877E-01  4.5589E+00  3.4046E+00 -7.3235E+00  4.6023E-01 -3.6320E+00  1.4703E-01  3.9961E+00  1.2360E-01
             1.1157E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1586.04348781879        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.7921E-01  5.6079E-02  1.0252E+00  1.6176E+00  6.9961E-01  9.2654E-01  5.6264E+00  9.0252E-01  7.7818E-01  7.3311E-01
             1.3929E+00
 PARAMETER:  7.8992E-02 -2.7810E+00  1.2487E-01  5.8097E-01 -2.5723E-01  2.3703E-02  1.8275E+00 -2.5589E-03 -1.5079E-01 -2.1045E-01
             4.3136E-01
 GRADIENT:   1.2049E+00 -2.5484E+00  3.4458E+00  1.8110E+01 -3.2603E+00 -7.0853E-01 -6.5988E+00  6.7388E-02  5.5717E+00 -3.1272E-01
            -6.7502E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1586.15682639545        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1416
 NPARAMETR:  9.7854E-01  5.7119E-02  9.8671E-01  1.6051E+00  6.8242E-01  9.2915E-01  5.5358E+00  8.6211E-01  7.7776E-01  7.3608E-01
             1.3925E+00
 PARAMETER:  7.8308E-02 -2.7626E+00  8.6617E-02  5.7316E-01 -2.8212E-01  2.6517E-02  1.8112E+00 -4.8376E-02 -1.5134E-01 -2.0642E-01
             4.3112E-01
 GRADIENT:  -5.2015E-01 -1.8167E+00  1.2335E+00  2.2810E+00 -1.3416E+00  2.3777E-01 -4.6411E+00  4.4337E-01  2.9987E+00  6.2045E-01
             3.3885E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1586.45718596024        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1594
 NPARAMETR:  9.7863E-01  4.2219E-02  9.5998E-01  1.6124E+00  6.6590E-01  9.2743E-01  6.3454E+00  8.1274E-01  7.7957E-01  7.1860E-01
             1.3927E+00
 PARAMETER:  7.8393E-02 -3.0649E+00  5.9161E-02  5.7771E-01 -3.0661E-01  2.4663E-02  1.9477E+00 -1.0734E-01 -1.4901E-01 -2.3045E-01
             4.3122E-01
 GRADIENT:   3.6008E-01 -2.4975E+00  5.4178E+00  1.5704E+01 -5.7650E+00 -3.5096E-01 -5.4730E+00 -7.6817E-01  5.6419E+00 -6.2949E-01
            -5.4610E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1586.45769590488        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1769
 NPARAMETR:  9.7857E-01  4.3280E-02  9.6063E-01  1.6106E+00  6.6664E-01  9.2790E-01  6.2704E+00  8.2082E-01  7.7856E-01  7.2009E-01
             1.3923E+00
 PARAMETER:  7.8334E-02 -3.0401E+00  5.9834E-02  5.7662E-01 -3.0550E-01  2.5169E-02  1.9358E+00 -9.7448E-02 -1.5031E-01 -2.2838E-01
             4.3096E-01
 GRADIENT:   2.0105E-01 -2.3150E+00  4.2561E+00  1.2207E+01 -4.6764E+00 -1.7405E-01 -5.1207E+00 -4.0693E-01  4.8322E+00 -2.7396E-01
            -2.8480E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1586.45994428527        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1944
 NPARAMETR:  9.7849E-01  4.4286E-02  9.6036E-01  1.6084E+00  6.6705E-01  9.2844E-01  6.2029E+00  8.2835E-01  7.7723E-01  7.2203E-01
             1.3919E+00
 PARAMETER:  7.8254E-02 -3.0171E+00  5.9555E-02  5.7527E-01 -3.0489E-01  2.5748E-02  1.9250E+00 -8.8317E-02 -1.5201E-01 -2.2569E-01
             4.3067E-01
 GRADIENT:   6.4036E-02 -1.2163E+00  2.0904E+00  5.9028E+00 -2.3420E+00 -4.8785E-02 -2.7410E+00 -1.5144E-01  2.4912E+00 -6.8902E-02
            -8.9932E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1586.63388442182        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2129
 NPARAMETR:  9.7791E-01  3.6821E-02  9.5350E-01  1.6037E+00  6.6407E-01  9.2855E-01  6.6453E+00  8.2692E-01  7.6543E-01  7.1541E-01
             1.3889E+00
 PARAMETER:  7.7665E-02 -3.2017E+00  5.2380E-02  5.7231E-01 -3.0937E-01  2.5874E-02  1.9939E+00 -9.0044E-02 -1.6731E-01 -2.3490E-01
             4.2850E-01
 GRADIENT:  -8.0815E-01 -2.7368E+00  4.2787E-01 -5.2591E+00  1.9337E+00  2.4853E-01 -4.8342E+00  1.6248E-01 -1.5981E-01 -1.4409E-02
            -9.0402E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1586.69450088756        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2305
 NPARAMETR:  9.7740E-01  3.7178E-02  9.4782E-01  1.6040E+00  6.5936E-01  9.2845E-01  6.6777E+00  8.2034E-01  7.7238E-01  7.2046E-01
             1.3900E+00
 PARAMETER:  7.7143E-02 -3.1920E+00  4.6408E-02  5.7247E-01 -3.1649E-01  2.5757E-02  1.9988E+00 -9.8041E-02 -1.5828E-01 -2.2787E-01
             4.2933E-01
 GRADIENT:  -2.1473E+00 -1.2619E+00  3.0563E+00 -4.9055E+00 -3.3621E+00  5.9055E-02 -1.7835E+00  2.5159E-01  8.2032E-01  6.0748E-01
            -2.9395E-02

0ITERATION NO.:   73    OBJECTIVE VALUE:  -1586.70373930566        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     2401
 NPARAMETR:  9.7834E-01  3.6996E-02  9.4299E-01  1.6064E+00  6.5743E-01  9.2856E-01  6.6545E+00  8.1898E-01  7.7674E-01  7.1631E-01
             1.3899E+00
 PARAMETER:  7.7100E-02 -3.1890E+00  4.1014E-02  5.7253E-01 -3.1850E-01  2.5240E-02  1.9994E+00 -9.9493E-02 -1.5310E-01 -2.3392E-01
             4.2898E-01
 GRADIENT:  -2.1829E+00  1.8961E+01 -7.3187E+02 -1.1684E+02  2.2836E+02 -2.2382E-01  2.5594E+01  4.9478E-02 -4.7411E+02 -1.3886E-01
            -2.3590E-01
 NUMSIGDIG:         1.8         2.4         2.3         2.4         2.3         1.9         2.4         2.5         2.3         2.7
                    3.0

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2401
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6426E-03  1.9259E-02 -1.5332E-02 -1.2370E-02 -1.9382E-02
 SE:             2.9670E-02  1.0045E-02  1.6218E-02  2.7874E-02  1.9044E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5585E-01  5.5202E-02  3.4446E-01  6.5720E-01  3.0881E-01

 ETASHRINKSD(%)  6.0118E-01  6.6348E+01  4.5669E+01  6.6177E+00  3.6199E+01
 ETASHRINKVR(%)  1.1987E+00  8.8676E+01  7.0481E+01  1.2797E+01  5.9294E+01
 EBVSHRINKSD(%)  8.9875E-01  7.8260E+01  4.5906E+01  5.3307E+00  3.4742E+01
 EBVSHRINKVR(%)  1.7894E+00  9.5274E+01  7.0738E+01  1.0377E+01  5.7414E+01
 RELATIVEINF(%)  9.7965E+01  2.8498E+00  2.7696E+00  5.2785E+01  3.9474E+00
 EPSSHRINKSD(%)  4.1174E+01
 EPSSHRINKVR(%)  6.5395E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1586.7037393056635     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -851.55291274192530     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1586.704       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.77E-01  3.73E-02  9.43E-01  1.60E+00  6.58E-01  9.28E-01  6.68E+00  8.19E-01  7.76E-01  7.16E-01  1.39E+00
 


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
+        1.33E+03
 
 TH 2
+       -7.43E+01  1.25E+05
 
 TH 3
+        7.84E+01  3.01E+02  2.06E+05
 
 TH 4
+       -1.22E+01  7.47E+02 -9.03E+01  2.73E+03
 
 TH 5
+       -1.38E+01 -8.32E+02 -8.13E+02 -1.53E+01  4.34E+04
 
 TH 6
+       -6.95E+00  3.33E+01 -1.05E+02 -1.08E+01  3.52E+01  2.38E+02
 
 TH 7
+       -2.23E-01 -2.28E+01  1.94E+00  4.43E+00 -4.79E+00  3.72E-01  1.03E+01
 
 TH 8
+       -2.45E+01  4.72E+01  2.37E+05 -1.06E+01  5.77E+01  3.34E+01  5.67E-01  1.10E+02
 
 TH 9
+        5.57E+01  4.29E+02 -3.56E+02 -9.29E+01  5.36E+02 -6.92E+01  4.14E+00 -1.06E+02  1.29E+05
 
 TH10
+       -6.43E+01  4.03E+01  1.16E+05 -1.86E+01 -1.93E+01  9.01E+01  6.30E-01  1.33E+05 -1.21E+02  6.53E+04
 
 TH11
+       -2.47E+01  2.44E+00  3.25E+04 -1.26E+01  5.16E+00  2.14E+01  1.17E-01  6.57E+01 -2.23E+01  1.83E+04  5.26E+03
 
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
 #CPUT: Total CPU Time in Seconds,       37.671
Stop Time:
Wed Sep 29 12:20:31 CDT 2021
