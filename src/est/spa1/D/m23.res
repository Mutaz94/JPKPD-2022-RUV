Thu Sep 30 02:48:51 CDT 2021
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
$DATA ../../../../data/spa1/D/dat23.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   3285.69646996215        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7024E+02  8.0311E+01 -1.2010E+02 -4.4628E+00  2.9200E+02 -8.3428E+02 -3.1688E+02 -2.0154E+01 -7.4512E+02 -3.9296E+02
            -8.3839E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1091.48919631296        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3714E+00  8.5371E-01  1.9630E+00  2.0729E+00  1.2013E+00  2.2315E+00  1.7127E+00  9.2336E-01  4.6277E+00  1.4987E+00
             5.8378E+00
 PARAMETER:  4.1580E-01 -5.8162E-02  7.7445E-01  8.2893E-01  2.8339E-01  9.0266E-01  6.3809E-01  2.0268E-02  1.6321E+00  5.0460E-01
             1.8644E+00
 GRADIENT:   1.1089E+02 -1.3972E+01 -2.5174E+01  7.6129E+01 -4.3626E+01  5.6467E+01  1.2992E+01  3.3992E+00  1.3858E+02  2.1595E+01
             2.2271E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1119.09437520609        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      220
 NPARAMETR:  1.9106E+00  1.2867E+00  6.6855E+00  1.4291E+00  1.3858E+00  3.5523E+00  1.9244E+00  9.0301E-01  4.5345E+00  1.4103E+00
             4.8139E+00
 PARAMETER:  7.4743E-01  3.5207E-01  1.9999E+00  4.5702E-01  4.2626E-01  1.3676E+00  7.5461E-01 -2.0187E-03  1.6117E+00  4.4382E-01
             1.6715E+00
 GRADIENT:   8.2165E+01 -2.4868E+00 -1.9306E+00  3.2303E+01 -5.4680E+01  7.5279E+01  1.9464E+01  3.6147E-01  5.3089E+01  2.7487E+01
             9.7490E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1194.43570252742        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      402
 NPARAMETR:  1.1789E+00  1.4325E+00  5.8108E+00  8.8087E-01  4.1714E+00  2.2530E+00  8.2322E-01  1.6853E+00  4.1518E+00  2.9390E+00
             4.6218E+00
 PARAMETER:  2.6457E-01  4.5944E-01  1.8597E+00 -2.6846E-02  1.5283E+00  9.1228E-01 -9.4537E-02  6.2196E-01  1.5235E+00  1.1781E+00
             1.6308E+00
 GRADIENT:  -5.3761E+00  2.0725E+01  5.9640E+00  6.2789E+00 -5.6023E-01 -4.7355E+00  3.5821E+00 -1.1601E+00 -1.1000E+00  5.7611E-01
             5.8876E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1208.31248608151        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  1.1678E+00  9.3864E-01  1.5636E+00  1.0895E+00  6.0654E+00  1.9871E+00  1.9444E-01  1.5412E+00  3.0924E+00  4.1494E+00
             4.3434E+00
 PARAMETER:  2.5512E-01  3.6678E-02  5.4698E-01  1.8572E-01  1.9026E+00  7.8670E-01 -1.5376E+00  5.3259E-01  1.2290E+00  1.5230E+00
             1.5687E+00
 GRADIENT:   1.5100E+01  3.0501E+00 -8.0955E+00 -4.7461E+00 -1.1407E+01 -3.7129E+01  3.1729E-01  4.6068E+00  1.6499E+01  1.2333E+00
            -1.4546E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1221.06701636399        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  1.0592E+00  5.7371E-01  8.7048E-01  1.3124E+00  1.1605E+01  2.0887E+00  5.0121E-02  1.1773E-01  2.1050E+00  6.9998E+00
             4.3530E+00
 PARAMETER:  1.5755E-01 -4.5563E-01 -3.8714E-02  3.7186E-01  2.5514E+00  8.3653E-01 -2.8933E+00 -2.0393E+00  8.4434E-01  2.0459E+00
             1.5709E+00
 GRADIENT:  -8.0948E+00  2.1127E+01  5.0732E+00  1.5690E+01 -1.2964E+01  3.9720E+00  1.0415E-02 -4.6666E-02 -1.3165E+01  1.4863E+01
            -3.2556E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1225.00637885143        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  1.0591E+00  4.7328E-01  7.0114E-01  1.3466E+00  1.3241E+01  2.1094E+00  2.7258E-02  6.8355E-02  1.9415E+00  7.4519E+00
             4.4464E+00
 PARAMETER:  1.5739E-01 -6.4806E-01 -2.5505E-01  3.9761E-01  2.6833E+00  8.4640E-01 -3.5024E+00 -2.5830E+00  7.6344E-01  2.1085E+00
             1.5921E+00
 GRADIENT:   2.0775E+00  2.8174E+01 -3.0195E+00  2.6028E+01 -9.3412E+00  1.4557E+01  2.5051E-03 -1.6760E-02 -1.4550E+01  1.4926E+01
            -1.0304E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1240.59087577783        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  1.0573E+00  2.9177E-01  6.8093E-01  1.3513E+00  1.3288E+01  2.0938E+00  1.0000E-02  1.2467E+00  1.9188E+00  7.4185E+00
             4.4321E+00
 PARAMETER:  1.5574E-01 -1.1318E+00 -2.8429E-01  4.0105E-01  2.6868E+00  8.3899E-01 -4.6293E+00  3.2051E-01  7.5168E-01  2.1040E+00
             1.5889E+00
 GRADIENT:   5.3560E+01  8.4305E+00  2.2361E+01  3.9070E+01 -1.0688E+00  1.0545E+02  0.0000E+00 -7.4496E+00  5.0115E+01  1.6345E-01
             3.3203E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1248.08868848253        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  1.0568E+00  1.5761E-01  6.8087E-01  1.3497E+00  1.3304E+01  2.0820E+00  1.0000E-02  2.7049E+00  1.9125E+00  7.4077E+00
             4.4234E+00
 PARAMETER:  1.5527E-01 -1.7477E+00 -2.8438E-01  3.9986E-01  2.6881E+00  8.3334E-01 -4.6474E+00  1.0951E+00  7.4841E-01  2.1025E+00
             1.5869E+00
 GRADIENT:   2.8424E+01 -7.4324E-01  1.3501E+01 -8.5632E+00 -4.8084E-01  3.9748E+01  0.0000E+00 -3.0855E+00  7.1098E+01 -6.4763E-03
             1.9271E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1255.28479924794        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1343
 NPARAMETR:  1.0405E+00  1.8392E-02  6.7767E-01  1.3322E+00  1.3574E+01  1.7452E+00  1.0000E-02  3.5393E+00  1.6696E+00  7.2838E+00
             4.0174E+00
 PARAMETER:  1.3966E-01 -3.8959E+00 -2.8910E-01  3.8682E-01  2.7082E+00  6.5687E-01 -5.4854E+00  1.3639E+00  6.1256E-01  2.0857E+00
             1.4906E+00
 GRADIENT:   5.9922E+01 -6.7917E-01  1.1867E+01 -1.5795E+01  2.3973E-01 -4.0105E-01  0.0000E+00  4.7817E+00  7.3394E+01 -2.9037E-02
            -9.7600E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1266.97212444075        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  1.0243E+00  1.1807E-01  6.7486E-01  1.3422E+00  1.3665E+01  1.6044E+00  1.0000E-02  2.9981E+00  1.4568E+00  7.2605E+00
             4.4236E+00
 PARAMETER:  1.2397E-01 -2.0364E+00 -2.9325E-01  3.9433E-01  2.7149E+00  5.7277E-01 -5.9194E+00  1.1980E+00  4.7624E-01  2.0824E+00
             1.5870E+00
 GRADIENT:   4.0305E+01  4.7145E+00  2.2544E+01  1.2349E+01 -1.9344E-01 -3.2366E+01  0.0000E+00 -1.3993E+01  3.9800E+01 -1.2234E-01
             1.1679E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1286.91385046508        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1657
 NPARAMETR:  8.8050E-01  6.5856E-02  4.1675E-01  1.2014E+00  1.6427E+01  1.9168E+00  1.0000E-02  3.5439E+00  7.6173E-01  3.0165E+00
             4.3638E+00
 PARAMETER: -2.7262E-02 -2.6203E+00 -7.7526E-01  2.8350E-01  2.8989E+00  7.5066E-01 -6.3130E+00  1.3652E+00 -1.7216E-01  1.2041E+00
             1.5733E+00
 GRADIENT:   1.4761E+01  7.0004E+00  1.4059E+01  1.7181E+02  6.4537E-01  1.2016E+02  0.0000E+00  4.4878E+01 -1.6774E+01 -1.0027E-02
             6.2062E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1306.13344712018        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     1788
 NPARAMETR:  8.7090E-01  3.9344E-02  3.4730E-01  1.0357E+00  2.0326E+01  1.4439E+00  1.0000E-02  3.1652E+00  7.8914E-01  5.7764E+00
             4.3704E+00
 PARAMETER: -3.8224E-02 -3.1354E+00 -9.5757E-01  1.3506E-01  3.1119E+00  4.6738E-01 -6.3130E+00  1.2522E+00 -1.3682E-01  1.8538E+00
             1.5749E+00
 GRADIENT:   5.2951E+01  6.2244E-03  5.0873E+01  1.0064E+00  1.5476E-01  1.1522E+01  0.0000E+00  3.9882E+01  1.9220E+01 -2.7628E-02
             2.4556E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1314.44741280652        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1884
 NPARAMETR:  7.5576E-01  1.3707E-02  1.6101E-01  8.2547E-01  1.8054E+01  1.1549E+00  1.0000E-02  2.0673E+00  6.5190E-01  6.0399E+00
             4.3330E+00
 PARAMETER: -1.8004E-01 -4.1898E+00 -1.7263E+00 -9.1800E-02  2.9934E+00  2.4399E-01 -6.3130E+00  8.2626E-01 -3.2786E-01  1.8984E+00
             1.5663E+00
 GRADIENT:   6.0231E+00  4.8965E-02  1.3894E+01  2.5985E+01  1.6348E-01 -1.1168E+02  0.0000E+00 -1.0905E+01  1.3831E+01 -7.6806E-02
            -5.7111E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1339.33118150142        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2064
 NPARAMETR:  6.4503E-01  1.0000E-02  8.0446E-02  5.9529E-01  1.3285E+01  1.4694E+00  1.0000E-02  1.7041E+00  2.5949E-01  6.8856E+00
             4.3266E+00
 PARAMETER: -3.3847E-01 -5.3471E+00 -2.4202E+00 -4.1870E-01  2.6866E+00  4.8484E-01 -6.3130E+00  6.3303E-01 -1.2491E+00  2.0294E+00
             1.5648E+00
 GRADIENT:   5.1360E+00  0.0000E+00 -2.4582E+01  4.3826E+01  1.3993E+00  2.1181E+00  0.0000E+00 -4.0056E-01  1.0928E+00 -2.0399E+00
            -4.3521E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1340.22502697929        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2248
 NPARAMETR:  6.3998E-01  1.0000E-02  8.1656E-02  5.8338E-01  1.3381E+01  1.4634E+00  1.0000E-02  1.7159E+00  1.5551E-01  6.9071E+00
             4.3389E+00
 PARAMETER: -3.4632E-01 -5.3513E+00 -2.4052E+00 -4.3891E-01  2.6938E+00  4.8079E-01 -6.3130E+00  6.3993E-01 -1.7610E+00  2.0326E+00
             1.5676E+00
 GRADIENT:   5.3082E+00  0.0000E+00  4.1116E+00 -6.9901E+00  9.3980E-01  1.0414E+00  0.0000E+00 -1.7078E+00  3.7568E-01 -1.5937E+00
             1.1739E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1340.35076561908        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2423
 NPARAMETR:  6.3681E-01  1.0000E-02  8.1985E-02  5.8500E-01  1.3427E+01  1.4588E+00  1.0000E-02  1.7623E+00  8.0443E-02  6.8971E+00
             4.3341E+00
 PARAMETER: -3.5128E-01 -5.3513E+00 -2.4012E+00 -4.3615E-01  2.6972E+00  4.7759E-01 -6.3130E+00  6.6661E-01 -2.4202E+00  2.0311E+00
             1.5665E+00
 GRADIENT:  -3.1762E-01  0.0000E+00  1.2073E+00  2.8193E-01 -1.0664E+00 -2.0101E-01  0.0000E+00  7.3670E-01  7.1383E-02  7.7511E-01
            -3.8311E-01

0ITERATION NO.:   84    OBJECTIVE VALUE:  -1340.38810155711        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     2560
 NPARAMETR:  6.3812E-01  1.0000E-02  8.1932E-02  5.8463E-01  1.3556E+01  1.4616E+00  1.0000E-02  1.7554E+00  3.2554E-02  6.8671E+00
             4.3377E+00
 PARAMETER: -3.4921E-01 -5.3513E+00 -2.4022E+00 -4.3670E-01  2.7063E+00  4.7954E-01 -6.3130E+00  6.6256E-01 -3.2919E+00  2.0271E+00
             1.5671E+00
 GRADIENT:   7.3818E-01  0.0000E+00 -3.6048E+01  2.2050E+02 -3.4591E+01 -6.0025E-03  0.0000E+00 -1.4671E+02  1.4321E-02  4.6736E+01
            -6.1437E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2560
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2171E-03 -3.4755E-06 -5.0245E-03 -1.9335E-03 -7.6564E-03
 SE:             2.9558E-02  1.3123E-06  2.7352E-02  1.0563E-03  4.6454E-03
 N:                     100         100         100         100         100

 P VAL.:         9.1333E-01  8.0884E-03  8.5425E-01  6.7190E-02  9.9321E-02

 ETASHRINKSD(%)  9.7658E-01  9.9996E+01  8.3686E+00  9.6461E+01  8.4437E+01
 ETASHRINKVR(%)  1.9436E+00  1.0000E+02  1.6037E+01  9.9875E+01  9.7578E+01
 EBVSHRINKSD(%)  1.3388E+00  9.9996E+01  6.8663E+00  9.6578E+01  8.9886E+01
 EBVSHRINKVR(%)  2.6597E+00  1.0000E+02  1.3261E+01  9.9883E+01  9.8977E+01
 RELATIVEINF(%)  3.0681E+01  4.4721E-09  1.1655E+00  7.5732E-04  9.5267E-01
 EPSSHRINKSD(%)  1.7568E+01
 EPSSHRINKVR(%)  3.2050E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1340.3881015571119     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -421.44956835243920     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.13
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1340.388       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.38E-01  1.00E-02  8.19E-02  5.85E-01  1.35E+01  1.46E+00  1.00E-02  1.76E+00  3.36E-02  6.87E+00  4.34E+00
 


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
+        1.25E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.82E+02  0.00E+00  1.09E+05
 
 TH 4
+       -6.30E+02  0.00E+00 -5.95E+04  4.05E+04
 
 TH 5
+        2.24E-01  0.00E+00  2.95E-01  1.70E+00  1.79E+00
 
 TH 6
+        5.23E+00  0.00E+00  4.38E+01 -1.82E+01  1.69E-02  8.85E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        7.63E+00  0.00E+00 -9.54E+01 -8.25E+03  1.31E-01  1.05E+00  0.00E+00  1.85E+03
 
 TH 9
+       -4.58E+00  0.00E+00  1.00E+01 -6.24E+00 -3.00E-02 -4.28E-01  0.00E+00  1.56E+00  1.48E+01
 
 TH10
+       -1.57E-01  0.00E+00 -2.97E+00 -3.62E+00 -1.45E-01 -1.57E-02  0.00E+00 -3.05E-01  6.70E-02  1.30E+01
 
 TH11
+       -1.15E+01  0.00E+00  1.91E+03 -1.42E+03  3.03E-02  9.46E-01  0.00E+00  3.96E+00  1.73E-01 -1.07E-01  8.35E+01
 
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
 #CPUT: Total CPU Time in Seconds,       52.102
Stop Time:
Thu Sep 30 02:49:44 CDT 2021
