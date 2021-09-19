Sat Sep 18 08:30:46 CDT 2021
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
$DATA ../../../../data/spa/B/dat44.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1579.18185789349        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9466E+02  3.1647E+01  1.3784E+01  4.7642E+01  2.1742E+01 -2.0577E+01  1.5401E+00 -1.5285E+01  1.0008E+01 -1.1104E+01
            -3.5263E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1584.07887076903        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9570E-01  9.7344E-01  8.9410E-01  1.0079E+00  9.0918E-01  1.0450E+00  9.8207E-01  1.1365E+00  9.0479E-01  9.9117E-01
             1.0722E+00
 PARAMETER:  9.5690E-02  7.3077E-02 -1.1942E-02  1.0789E-01  4.7886E-03  1.4399E-01  8.1905E-02  2.2799E-01 -4.9280E-05  9.1134E-02
             1.6970E-01
 GRADIENT:   1.6419E+02  3.2020E+01  2.9736E+00  4.2924E+01 -1.8497E+01  1.0356E+00 -3.3079E+00  9.0237E-01 -4.8558E+00  9.1093E+00
            -8.7054E-02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1586.06928081141        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.8074E-01  9.8710E-01  8.1698E-01  9.8850E-01  8.7471E-01  1.0397E+00  1.1957E+00  1.1211E+00  8.4488E-01  8.4713E-01
             1.0758E+00
 PARAMETER:  8.0549E-02  8.7015E-02 -1.0215E-01  8.8437E-02 -3.3859E-02  1.3898E-01  2.7874E-01  2.1432E-01 -6.8563E-02 -6.5900E-02
             1.7308E-01
 GRADIENT:   1.3305E+02  3.0140E+01  3.0477E+00  2.9992E+01 -8.0458E+00  2.0759E+00  7.5395E+00  2.0712E+00 -1.7821E+00 -1.3567E+00
             2.8694E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1589.14099120980        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.2386E-01  9.0521E-01  8.2267E-01  1.0194E+00  8.5112E-01  1.0123E+00  1.1577E+00  1.0188E+00  8.5800E-01  8.5470E-01
             1.0636E+00
 PARAMETER:  2.0810E-02  4.0703E-04 -9.5204E-02  1.1919E-01 -6.1198E-02  1.1221E-01  2.4639E-01  1.1860E-01 -5.3152E-02 -5.7002E-02
             1.6170E-01
 GRADIENT:  -2.3327E+01  1.8074E+00  4.0392E-01 -1.7063E+00 -1.6597E+00 -9.9081E+00 -2.8579E-01 -2.2053E-01 -5.9155E-01 -5.2608E-01
            -2.6935E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1589.82371306849        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  9.3211E-01  7.0240E-01  1.0311E+00  1.1615E+00  8.7456E-01  1.0393E+00  1.3405E+00  1.1584E+00  8.0334E-01  9.1281E-01
             1.0714E+00
 PARAMETER:  2.9694E-02 -2.5325E-01  1.3064E-01  2.4974E-01 -3.4037E-02  1.3853E-01  3.9306E-01  2.4704E-01 -1.1898E-01  8.7683E-03
             1.6893E-01
 GRADIENT:   2.3521E+00  4.0104E+00 -5.2632E-01  8.2894E+00  2.4334E-01  2.3120E+00 -4.7346E-01 -2.0250E-01 -7.7206E-01 -5.0832E-01
             7.7405E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1590.31750474200        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  9.2654E-01  4.1546E-01  1.3225E+00  1.3581E+00  8.9270E-01  1.0235E+00  1.6702E+00  1.3558E+00  7.5782E-01  9.7672E-01
             1.0719E+00
 PARAMETER:  2.3701E-02 -7.7836E-01  3.7953E-01  4.0612E-01 -1.3510E-02  1.2325E-01  6.1293E-01  4.0439E-01 -1.7731E-01  7.6442E-02
             1.6946E-01
 GRADIENT:   3.9663E-01  7.3750E+00  5.7250E+00  2.5063E+01 -9.3706E+00 -1.7063E+00  8.3904E-01 -1.9902E+00 -3.8820E-01 -4.9329E-01
            -3.7748E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1590.81167101478        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  9.2411E-01  2.9860E-01  1.4499E+00  1.4251E+00  9.1059E-01  1.0256E+00  1.5175E+00  1.5006E+00  7.5425E-01  1.0115E+00
             1.0712E+00
 PARAMETER:  2.1071E-02 -1.1086E+00  4.7149E-01  4.5422E-01  6.3377E-03  1.2523E-01  5.1707E-01  5.0589E-01 -1.8203E-01  1.1146E-01
             1.6882E-01
 GRADIENT:  -6.5846E-01  7.7721E-01 -8.6738E-02  3.5557E+00 -7.6641E-01 -2.6696E-01  5.4312E-01  4.7120E-01  1.1309E+00  5.8251E-01
             6.3691E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1590.86480941639        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  9.2372E-01  2.4443E-01  1.4925E+00  1.4584E+00  9.0898E-01  1.0262E+00  1.3488E+00  1.5424E+00  7.4536E-01  1.0107E+00
             1.0716E+00
 PARAMETER:  2.0649E-02 -1.3088E+00  5.0045E-01  4.7734E-01  4.5666E-03  1.2582E-01  3.9919E-01  5.3335E-01 -1.9389E-01  1.1061E-01
             1.6918E-01
 GRADIENT:   3.5722E-01 -1.2490E-02  2.5285E-02  1.7914E+00 -4.0245E-01  1.9281E-01  3.1233E-01  1.5804E-01  4.0820E-01  1.3863E-01
             1.1460E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1591.41590851833        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1141
 NPARAMETR:  9.2581E-01  3.7469E-01  1.4652E+00  1.3790E+00  9.3529E-01  1.0273E+00  8.9071E-02  1.5045E+00  8.0269E-01  1.0332E+00
             1.0728E+00
 PARAMETER:  2.2917E-02 -8.8165E-01  4.8201E-01  4.2134E-01  3.3098E-02  1.2694E-01 -2.3183E+00  5.0849E-01 -1.1979E-01  1.3270E-01
             1.7024E-01
 GRADIENT:   4.7498E-01  1.9927E+00  3.1505E+00  1.1115E+01 -5.1110E+00 -9.9221E-02  1.0411E-02 -7.7686E-01 -6.1576E-01 -2.9206E-01
             2.9677E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1591.64909488203        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1317
 NPARAMETR:  9.2759E-01  5.2645E-01  1.3851E+00  1.2790E+00  9.6214E-01  1.0297E+00  1.0000E-02  1.4687E+00  8.6504E-01  1.0523E+00
             1.0712E+00
 PARAMETER:  2.4840E-02 -5.4161E-01  4.2580E-01  3.4607E-01  6.1407E-02  1.2931E-01 -6.1860E+00  4.8435E-01 -4.4974E-02  1.5101E-01
             1.6873E-01
 GRADIENT:  -8.9771E-02  1.2410E+00 -7.4357E-01  5.4602E+00  4.0439E-01  7.5933E-02  0.0000E+00 -2.2591E-02 -1.1712E+00 -4.3087E-03
            -2.6867E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1591.69247934917        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1492
 NPARAMETR:  9.2861E-01  6.1384E-01  1.3757E+00  1.2211E+00  9.8725E-01  1.0306E+00  1.0000E-02  1.4974E+00  9.0869E-01  1.0666E+00
             1.0719E+00
 PARAMETER:  2.5929E-02 -3.8802E-01  4.1893E-01  2.9972E-01  8.7170E-02  1.3013E-01 -9.0850E+00  5.0374E-01  4.2484E-03  1.6446E-01
             1.6947E-01
 GRADIENT:   6.8351E-02  2.4210E-01  1.5071E-01  5.6888E-01 -1.7195E-01  1.6578E-02  0.0000E+00 -5.6467E-02 -4.0821E-02 -2.8831E-02
            -2.4115E-02

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1591.69262088521        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1619
 NPARAMETR:  9.2859E-01  6.1545E-01  1.3740E+00  1.2197E+00  9.8742E-01  1.0306E+00  1.0000E-02  1.4974E+00  9.0962E-01  1.0667E+00
             1.0719E+00
 PARAMETER:  2.5917E-02 -3.8539E-01  4.1775E-01  2.9861E-01  8.7339E-02  1.3011E-01 -9.1552E+00  5.0374E-01  5.2731E-03  1.6459E-01
             1.6947E-01
 GRADIENT:   4.8731E-03  2.2049E-03  3.3434E-03  1.7502E-03 -6.2141E-03  2.3348E-03  0.0000E+00 -3.3298E-04  7.5526E-04 -2.2701E-03
             1.9380E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1619
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.2142E-04 -3.5332E-04 -4.0177E-02 -4.3749E-03 -4.1443E-02
 SE:             2.9849E-02  1.2498E-04  1.8095E-02  2.9083E-02  2.0942E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7537E-01  4.6994E-03  2.6392E-02  8.8043E-01  4.7822E-02

 ETASHRINKSD(%)  9.1226E-04  9.9581E+01  3.9381E+01  2.5667E+00  2.9841E+01
 ETASHRINKVR(%)  1.8245E-03  9.9998E+01  6.3254E+01  5.0676E+00  5.0778E+01
 EBVSHRINKSD(%)  4.6083E-01  9.9605E+01  4.3557E+01  3.0120E+00  2.6418E+01
 EBVSHRINKVR(%)  9.1954E-01  9.9998E+01  6.8142E+01  5.9332E+00  4.5857E+01
 RELATIVEINF(%)  9.8450E+01  1.0722E-04  8.7245E+00  7.5823E+00  1.3087E+01
 EPSSHRINKSD(%)  4.5177E+01
 EPSSHRINKVR(%)  6.9945E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1591.6926208852135     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -856.54179432147532     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1591.693       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.29E-01  6.15E-01  1.37E+00  1.22E+00  9.87E-01  1.03E+00  1.00E-02  1.50E+00  9.10E-01  1.07E+00  1.07E+00
 


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
+        1.20E+03
 
 TH 2
+       -2.08E+01  4.90E+02
 
 TH 3
+        2.72E+00  6.32E+01  9.03E+01
 
 TH 4
+       -1.36E+01  5.58E+02 -1.24E+01  8.50E+02
 
 TH 5
+        4.06E+00 -2.07E+02 -1.65E+02 -4.74E+01  5.63E+02
 
 TH 6
+       -9.24E-01 -4.06E+00  4.52E-01 -3.00E+00 -2.05E+00  1.85E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -4.07E-01 -1.52E+01 -2.35E+01 -2.97E+00 -2.12E+00 -4.07E-01  0.00E+00  2.23E+01
 
 TH 9
+        1.06E+00 -1.04E+02  4.22E+00  1.67E+00 -1.01E+00 -2.70E+00  0.00E+00 -4.52E-01  2.18E+02
 
 TH10
+       -5.07E-02  4.07E+00 -4.14E+00  3.06E-01 -7.05E+01  1.10E+00  0.00E+00  1.21E+01 -1.33E+00  6.06E+01
 
 TH11
+       -1.03E+01 -1.61E+01 -7.20E+00 -1.01E+01 -7.57E+00  3.39E+00  0.00E+00  4.49E+00  1.33E+01  1.08E+01  1.78E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.712
Stop Time:
Sat Sep 18 08:31:12 CDT 2021
