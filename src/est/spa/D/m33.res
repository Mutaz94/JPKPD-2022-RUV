Sat Sep 25 14:14:55 CDT 2021
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
$DATA ../../../../data/spa/D/dat33.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12267.2199480101        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4818E+02  1.4025E+02 -2.2554E-01 -5.6549E+01  1.7830E+02 -1.6539E+03 -5.8103E+02 -6.1417E+01 -1.0649E+03 -4.4396E+02
            -2.3969E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -619.145550481289        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3294E+00  1.2108E+00  9.6600E-01  1.7684E+00  1.1194E+00  1.9131E+00  1.1712E+00  9.7364E-01  1.2425E+00  1.1557E+00
             1.4378E+01
 PARAMETER:  3.8469E-01  2.9125E-01  6.5404E-02  6.7010E-01  2.1278E-01  7.4871E-01  2.5801E-01  7.3284E-02  3.1714E-01  2.4471E-01
             2.7657E+00
 GRADIENT:  -1.7337E+01  3.0391E+01  1.0178E-01  4.9836E+01 -1.5194E+01  2.7542E+01 -1.9587E+00  4.3069E+00  7.0713E+00  4.1718E+00
             1.3633E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -635.885176126624        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.3226E+00  1.1718E+00  1.3437E+00  1.7766E+00  2.3587E+00  1.6592E+00  2.0907E+00  2.7804E-01  1.1987E+00  5.2170E+00
             1.3178E+01
 PARAMETER:  3.7958E-01  2.5852E-01  3.9542E-01  6.7473E-01  9.5812E-01  6.0635E-01  8.3751E-01 -1.1800E+00  2.8121E-01  1.7519E+00
             2.6786E+00
 GRADIENT:   7.5884E+00  2.7732E+01 -3.0984E+00  6.0003E+01 -1.2948E+01 -1.3315E+01  6.5133E-01  1.2177E-01  7.7668E+00  1.1921E+01
             9.7254E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -662.495736463285        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1163E+00  4.8833E-01  2.1065E+00  1.7301E+00  3.5270E+00  1.6463E+00  2.2200E+00  1.2540E-01  6.4820E-01  7.7018E+00
             1.1185E+01
 PARAMETER:  2.1002E-01 -6.1677E-01  8.4505E-01  6.4820E-01  1.3605E+00  5.9855E-01  8.9752E-01 -1.9763E+00 -3.3355E-01  2.1414E+00
             2.5146E+00
 GRADIENT:  -1.8620E+01  9.9741E+00  1.1644E+01  2.9961E+01 -6.6397E+00  1.5344E+01 -5.2010E-01 -3.9384E-04 -7.0191E+00  2.2833E+00
             2.6556E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -693.271613593307        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0175E+00  2.8462E-01  5.9511E-01  1.4587E+00  4.6747E+00  1.4480E+00  1.5354E+00  1.0000E-02  3.9157E-01  6.0609E+00
             1.0126E+01
 PARAMETER:  1.1734E-01 -1.1566E+00 -4.1901E-01  4.7755E-01  1.6422E+00  4.7019E-01  5.2878E-01 -6.3437E+00 -8.3759E-01  1.9019E+00
             2.4151E+00
 GRADIENT:  -6.5622E-01  1.2659E+01  9.1520E+00 -2.8498E+01 -2.0260E+01 -1.8077E+01  6.4842E-01  0.0000E+00  4.0783E+00  4.6948E-01
            -2.4776E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -735.771960505940        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  8.9877E-01  3.9160E-02  2.0102E-01  1.1490E+00  5.5905E+00  1.6654E+00  5.7738E-01  1.0000E-02  4.0118E-02  2.0093E+00
             9.9330E+00
 PARAMETER: -6.7278E-03 -3.1401E+00 -1.5044E+00  2.3893E-01  1.8211E+00  6.1007E-01 -4.4926E-01 -1.5321E+01 -3.1159E+00  7.9778E-01
             2.3959E+00
 GRADIENT:   4.1308E+01 -1.0236E-01 -1.1296E+01  7.6947E+01  1.0670E+01  4.8519E+00  3.9407E-03  0.0000E+00 -3.1091E-03 -1.3580E+00
            -4.2307E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -760.603126506144        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  4.8385E-01  1.0000E-02  3.2304E-02  3.7588E-01  6.0460E+00  1.4873E+00  3.9980E-02  1.0000E-02  1.0000E-02  1.6821E-01
             9.7763E+00
 PARAMETER: -6.2598E-01 -6.0491E+00 -3.3326E+00 -8.7848E-01  1.8994E+00  4.9694E-01 -3.1194E+00 -2.9968E+01 -7.4495E+00 -1.6825E+00
             2.3800E+00
 GRADIENT:  -9.9654E+00  0.0000E+00 -8.5982E+00  2.5085E+01  4.3541E+00 -1.4450E+00  3.1273E-05  0.0000E+00  0.0000E+00  8.7604E-02
            -1.5215E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -761.049270880101        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      619
 NPARAMETR:  5.2145E-01  1.0000E-02  3.7838E-02  4.1897E-01  5.7091E+00  1.5097E+00  5.0620E-02  1.0000E-02  1.0000E-02  2.0067E-01
             9.8350E+00
 PARAMETER: -5.5114E-01 -5.8024E+00 -3.1744E+00 -7.6995E-01  1.8421E+00  5.1194E-01 -2.8834E+00 -2.8656E+01 -7.0699E+00 -1.5061E+00
             2.3859E+00
 GRADIENT:   2.0819E-01  0.0000E+00 -2.3602E-01 -1.3719E-01 -1.3904E-01  5.1780E-02  4.3128E-05  0.0000E+00  0.0000E+00  1.7912E-01
             7.2531E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -761.131382408127        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  5.1791E-01  1.0000E-02  3.7538E-02  4.1677E-01  5.6603E+00  1.5060E+00  5.0220E-02  1.0000E-02  1.0000E-02  4.1594E-02
             9.8148E+00
 PARAMETER: -5.5796E-01 -5.8203E+00 -3.1824E+00 -7.7521E-01  1.8335E+00  5.0946E-01 -2.8913E+00 -2.8711E+01 -7.0874E+00 -3.0798E+00
             2.3839E+00
 GRADIENT:   2.9620E-01  0.0000E+00  5.8213E+00  3.0774E-01 -6.2645E-02 -1.4433E-01  5.1459E-05  0.0000E+00  0.0000E+00  9.2076E-03
             1.2119E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -761.138726295465        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      853
 NPARAMETR:  5.1869E-01  1.0000E-02  3.7456E-02  4.1657E-01  5.6493E+00  1.5081E+00  5.0217E-02  1.0000E-02  1.0000E-02  3.8069E-02
             9.8178E+00
 PARAMETER: -5.5644E-01 -5.8203E+00 -3.1846E+00 -7.7570E-01  1.8315E+00  5.1086E-01 -2.8914E+00 -2.8711E+01 -7.0874E+00 -3.1684E+00
             2.3842E+00
 GRADIENT:  -1.6336E+00  0.0000E+00  3.2812E-01 -7.6526E-02 -7.6882E-02 -3.6251E-01  5.2209E-05  0.0000E+00  0.0000E+00  7.8555E-03
            -9.1395E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -761.145550921044        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1029
 NPARAMETR:  5.2022E-01  1.0000E-02  3.7567E-02  4.1736E-01  5.6568E+00  1.5099E+00  5.0097E-02  1.0000E-02  1.0000E-02  1.0000E-02
             9.8309E+00
 PARAMETER: -5.5350E-01 -5.8203E+00 -3.1816E+00 -7.7380E-01  1.8329E+00  5.1204E-01 -2.8938E+00 -2.8711E+01 -7.0874E+00 -4.5488E+00
             2.3855E+00
 GRADIENT:  -8.7817E-02  0.0000E+00  1.0764E-01 -4.8106E-01  1.6803E-01  4.5562E-02  5.0668E-05  0.0000E+00  0.0000E+00  0.0000E+00
            -1.3980E-01

0ITERATION NO.:   54    OBJECTIVE VALUE:  -761.148751204949        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1157
 NPARAMETR:  5.2437E-01  1.0000E-02  3.8374E-02  4.2395E-01  5.6243E+00  1.5119E+00  4.9830E-02  1.0000E-02  1.0000E-02  1.0000E-02
             9.8350E+00
 PARAMETER: -5.4557E-01 -5.8203E+00 -3.1604E+00 -7.5813E-01  1.8271E+00  5.1338E-01 -2.8991E+00 -2.8711E+01 -7.0874E+00 -6.2439E+00
             2.3860E+00
 GRADIENT:   1.5388E-02  0.0000E+00 -4.1436E-03 -3.3036E-03 -4.5320E-03 -9.4164E-04  4.8357E-05  0.0000E+00  0.0000E+00  0.0000E+00
             1.7236E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1157
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4993E-03  1.1613E-05  1.1793E-04 -2.5514E-04  8.0770E-06
 SE:             2.8822E-02  2.2151E-05  2.7553E-04  3.7492E-04  6.7855E-05
 N:                     100         100         100         100         100

 P VAL.:         9.3090E-01  6.0009E-01  6.6864E-01  4.9618E-01  9.0525E-01

 ETASHRINKSD(%)  3.4422E+00  9.9926E+01  9.9077E+01  9.8744E+01  9.9773E+01
 ETASHRINKVR(%)  6.7659E+00  1.0000E+02  9.9991E+01  9.9984E+01  9.9999E+01
 EBVSHRINKSD(%)  3.5880E+00  9.9912E+01  9.9077E+01  9.8746E+01  9.9723E+01
 EBVSHRINKVR(%)  7.0473E+00  1.0000E+02  9.9991E+01  9.9984E+01  9.9999E+01
 RELATIVEINF(%)  8.0534E+00  5.7278E-06  6.6148E-05  1.1989E-04  6.5654E-05
 EPSSHRINKSD(%)  7.4548E+00
 EPSSHRINKVR(%)  1.4354E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -761.14875120494946     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -25.997924641211284     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -761.149       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.24E-01  1.00E-02  3.84E-02  4.24E-01  5.62E+00  1.51E+00  4.98E-02  1.00E-02  1.00E-02  1.00E-02  9.84E+00
 


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
+        1.62E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.62E+03  0.00E+00  3.68E+05
 
 TH 4
+       -2.58E+02  0.00E+00 -4.47E+04  6.11E+03
 
 TH 5
+        8.32E+00  0.00E+00 -5.23E+02  6.71E+01  2.09E+00
 
 TH 6
+        1.07E+00  0.00E+00  2.97E+02 -5.99E+01  5.16E-01  7.27E+01
 
 TH 7
+        1.98E-01  0.00E+00  5.84E-01 -1.67E-01 -6.81E-03  1.38E-02 -3.63E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.96E+01  0.00E+00  2.53E+02 -2.30E+01 -4.38E-01  9.98E-01 -5.74E-03  0.00E+00  0.00E+00  0.00E+00  4.66E+00
 
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
 #CPUT: Total CPU Time in Seconds,       19.965
Stop Time:
Sat Sep 25 14:15:19 CDT 2021
