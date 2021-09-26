Sat Sep 25 09:32:08 CDT 2021
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
$DATA ../../../../data/spa/A3/dat83.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   728.503554108023        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0188E+02  1.1398E+02  1.7474E+02 -8.1821E+01  3.9459E+01 -6.2190E+00 -5.9038E+01 -6.9882E+01 -1.4977E+02 -1.8832E+02
            -4.2627E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1178.48479033986        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0319E+00  9.4831E-01  8.8193E-01  1.1934E+00  1.2704E+00  8.1015E-01  8.2486E-01  9.7520E-01  7.3458E-01  9.0739E-01
             5.4629E+00
 PARAMETER:  1.3137E-01  4.6923E-02 -2.5639E-02  2.7678E-01  3.3932E-01 -1.1053E-01 -9.2536E-02  7.4891E-02 -2.0845E-01  2.8193E-03
             1.7980E+00
 GRADIENT:  -1.3324E+00 -3.5626E+01 -5.0894E+01 -1.4632E+01  2.7861E+01 -3.7757E+01  8.6350E+00  6.9981E+00  1.6084E+01  7.9924E+00
             1.5087E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1199.75607843548        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0224E+00  1.1320E+00  1.8384E+00  1.1743E+00  3.4730E+00  9.5205E-01  1.3099E-01  1.0492E+00  6.5374E-01  3.8959E+00
             5.0264E+00
 PARAMETER:  1.2214E-01  2.2396E-01  7.0891E-01  2.6071E-01  1.3450E+00  5.0858E-02 -1.9327E+00  1.4805E-01 -3.2504E-01  1.4599E+00
             1.7147E+00
 GRADIENT:   4.5994E+00  4.1107E+01 -8.1424E+00  8.6051E+01  2.3043E+00  8.2112E+00  1.6895E-01  1.1959E-01  1.1044E+01  9.2929E+00
             6.4130E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1212.29214174113        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0042E+00  1.0772E+00  2.0175E+00  1.1253E+00  2.6976E+00  9.0959E-01  1.8117E-01  3.2384E+00  4.8796E-01  1.9070E+00
             4.6672E+00
 PARAMETER:  1.0422E-01  1.7437E-01  8.0188E-01  2.1804E-01  1.0924E+00  5.2406E-03 -1.6083E+00  1.2751E+00 -6.1752E-01  7.4554E-01
             1.6406E+00
 GRADIENT:  -5.3643E+00 -7.5760E+00 -1.0134E+01 -2.8805E+00 -3.9726E+00 -6.8578E+00  9.8537E-02  6.2331E+00  2.7683E+00  1.7684E+00
             7.3407E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1213.04227466856        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      356
 NPARAMETR:  1.0096E+00  1.0892E+00  2.0796E+00  1.1293E+00  2.7831E+00  9.4480E-01  2.0678E-01  3.2721E+00  4.0192E-01  1.6010E+00
             4.6459E+00
 PARAMETER:  1.0957E-01  1.8542E-01  8.3217E-01  2.2161E-01  1.1236E+00  4.3213E-02 -1.4761E+00  1.2854E+00 -8.1150E-01  5.7060E-01
             1.6360E+00
 GRADIENT:   5.2091E+00  1.5745E+00 -4.8736E+00  5.7811E+00  3.9077E-01  4.1281E+00 -9.1350E-02 -1.4081E+00 -9.8618E-02  3.2810E-01
            -1.9216E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1213.22639754984        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      426
 NPARAMETR:  1.0070E+00  1.1053E+00  2.0796E+00  1.1214E+00  2.7831E+00  9.3001E-01  4.5295E-01  3.2724E+00  3.2179E-01  1.0889E+00
             4.6451E+00
 PARAMETER:  1.0698E-01  2.0010E-01  8.3216E-01  2.1458E-01  1.1236E+00  2.7441E-02 -6.9196E-01  1.2855E+00 -1.0339E+00  1.8513E-01
             1.6358E+00
 GRADIENT:  -2.6297E+00  4.1361E+00 -3.2011E+00  9.1117E+00  1.8537E+00 -8.7684E-01 -9.5579E-02 -3.2512E+00 -6.4050E-02  6.7100E-02
             1.5840E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1213.25891632049        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0078E+00  1.1222E+00  2.0795E+00  1.1081E+00  2.7828E+00  9.3191E-01  5.3996E-01  3.2737E+00  2.4932E-01  7.0507E-01
             4.6428E+00
 PARAMETER:  1.0779E-01  2.1527E-01  8.3211E-01  2.0262E-01  1.1234E+00  2.9476E-02 -5.1627E-01  1.2859E+00 -1.2890E+00 -2.4945E-01
             1.6353E+00
 GRADIENT:  -1.1799E-01  1.3933E+00 -2.0786E+00  4.2195E+00  2.7271E+00  5.2235E-02 -4.8065E-04 -4.3161E+00 -1.3971E-01  1.8515E-02
             1.4501E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1213.29602743681        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      569
 NPARAMETR:  1.0086E+00  1.2103E+00  2.0790E+00  1.0482E+00  2.7813E+00  9.3117E-01  5.5388E-01  3.2794E+00  1.7436E-01  2.7268E-01
             4.6329E+00
 PARAMETER:  1.0852E-01  2.9088E-01  8.3188E-01  1.4710E-01  1.1229E+00  2.8691E-02 -4.9080E-01  1.2877E+00 -1.6466E+00 -1.1994E+00
             1.6332E+00
 GRADIENT:   4.4594E-01  3.9206E-01 -4.8175E-01  1.1348E+00  3.8945E+00 -4.8601E-02 -1.1175E-01 -5.5473E+00 -1.3566E-01  5.9158E-03
            -1.2758E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1213.37473306059        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      641
 NPARAMETR:  1.0089E+00  1.2449E+00  2.0788E+00  1.0273E+00  2.7803E+00  9.3149E-01  4.3690E-01  3.2827E+00  3.3522E-01  5.0354E-01
             4.6297E+00
 PARAMETER:  1.0889E-01  3.1905E-01  8.3177E-01  1.2692E-01  1.1226E+00  2.9034E-02 -7.2806E-01  1.2887E+00 -9.9296E-01 -5.8609E-01
             1.6325E+00
 GRADIENT:  -3.9484E-01  1.2285E+00  1.5085E-01  1.0641E+00  4.3327E+00 -7.8338E-02 -2.8263E-03 -5.8588E+00  2.3306E-02  2.6346E-02
             4.8928E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1213.38203493550        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.0089E+00  1.2563E+00  2.0787E+00  1.0181E+00  2.7798E+00  9.3162E-01  4.5461E-01  3.2850E+00  3.2056E-01  3.8111E-01
             4.6273E+00
 PARAMETER:  1.0886E-01  3.2813E-01  8.3172E-01  1.1796E-01  1.1224E+00  2.9170E-02 -6.8832E-01  1.2894E+00 -1.0377E+00 -8.6467E-01
             1.6320E+00
 GRADIENT:   1.6244E-02 -7.0213E-01  3.0900E-01 -1.5425E+00  4.5161E+00  1.0858E-01  6.3915E-02 -5.9297E+00  2.5574E-02  1.6084E-02
             3.6791E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1213.65305182793        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  1.0085E+00  1.3185E+00  2.0791E+00  9.7377E-01  2.7280E+00  9.3293E-01  3.9604E-01  3.5272E+00  3.8230E-01  1.0000E-02
             4.5363E+00
 PARAMETER:  1.0849E-01  3.7652E-01  8.3193E-01  7.3425E-02  1.1036E+00  3.0578E-02 -8.2624E-01  1.3605E+00 -8.6156E-01 -3.4542E+01
             1.6121E+00
 GRADIENT:   1.0531E+00 -7.8743E-01 -2.6546E+00  3.6198E-02  4.1164E+00  3.1768E-02 -3.7689E-01 -3.5008E+00 -3.2235E-01  0.0000E+00
            -1.9982E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1213.69021964620        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      863
 NPARAMETR:  1.0087E+00  1.3338E+00  2.0791E+00  9.6541E-01  2.7247E+00  9.3256E-01  4.5479E-01  3.5431E+00  3.7278E-01  1.0000E-02
             4.5313E+00
 PARAMETER:  1.0863E-01  3.8803E-01  8.3192E-01  6.4797E-02  1.1024E+00  3.0179E-02 -6.8791E-01  1.3650E+00 -8.8676E-01 -3.6633E+01
             1.6110E+00
 GRADIENT:  -3.6476E-02  9.3691E-02 -2.1859E+00  5.2322E-01  4.3313E+00  3.7830E-03  7.1490E-02 -3.8614E+00  3.7613E-02  0.0000E+00
            -1.7772E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1214.73251546589        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  1.0123E+00  1.5714E+00  2.0784E+00  8.0472E-01  2.5960E+00  9.2898E-01  3.7331E-01  4.2048E+00  2.2927E-01  1.0000E-02
             4.6688E+00
 PARAMETER:  1.1223E-01  5.5197E-01  8.3161E-01 -1.1726E-01  1.0540E+00  2.6337E-02 -8.8536E-01  1.5362E+00 -1.3728E+00 -1.0800E+02
             1.6409E+00
 GRADIENT:  -1.2932E-01  9.4654E-01 -1.7923E+00  5.0912E-01  4.3487E+00 -5.6325E-01 -6.8500E-01 -3.2999E+00 -2.1657E-01  0.0000E+00
             3.3900E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1215.01813130093        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1015
 NPARAMETR:  1.0130E+00  1.6686E+00  2.0782E+00  7.4813E-01  2.5473E+00  9.2766E-01  4.3086E-01  4.5025E+00  1.7428E-01  1.0000E-02
             4.6747E+00
 PARAMETER:  1.1293E-01  6.1198E-01  8.3148E-01 -1.9017E-01  1.0350E+00  2.4913E-02 -7.4198E-01  1.6046E+00 -1.6471E+00 -1.3777E+02
             1.6422E+00
 GRADIENT:  -5.4120E+00  1.8536E+01 -6.6028E+00  1.6605E+01  1.3443E+00 -1.1315E+00  1.6922E-01  1.5047E+00 -4.8906E-02  0.0000E+00
             4.7823E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1215.27643179561        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  1.0124E+00  1.7033E+00  2.0779E+00  7.1595E-01  2.5171E+00  9.2886E-01  4.4158E-01  4.7040E+00  1.6961E-01  1.0000E-02
             4.6452E+00
 PARAMETER:  1.1233E-01  6.3256E-01  8.3137E-01 -2.3415E-01  1.0231E+00  2.6206E-02 -7.1738E-01  1.6484E+00 -1.6743E+00 -1.5723E+02
             1.6358E+00
 GRADIENT:  -2.3596E+00  2.0563E+01 -2.5700E+01  3.8144E+01 -8.3001E+00 -3.3834E-01  2.5485E-01  2.0213E+01 -5.0567E-02  0.0000E+00
            -1.5389E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1216.77575553314        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  1.0135E+00  1.8902E+00  2.0770E+00  5.9511E-01  2.3941E+00  9.2912E-01  3.9627E-01  5.6398E+00  1.4865E-01  1.0000E-02
             4.6142E+00
 PARAMETER:  1.1338E-01  7.3669E-01  8.3091E-01 -4.1900E-01  9.7300E-01  2.6478E-02 -8.2566E-01  1.8298E+00 -1.8062E+00 -2.3589E+02
             1.6291E+00
 GRADIENT:  -3.4509E+00  1.6821E+02 -1.9390E+02  2.9653E+02 -1.2755E+02 -2.3794E+00 -7.2153E+00  1.6215E+02 -7.2912E-01  0.0000E+00
            -1.5669E+02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1217.34551616694        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  1.0156E+00  1.9683E+00  2.0764E+00  5.3696E-01  2.3245E+00  9.3451E-01  4.6045E-01  6.2846E+00  1.2192E-01  1.0000E-02
             4.5455E+00
 PARAMETER:  1.1552E-01  7.7719E-01  8.3066E-01 -5.2184E-01  9.4352E-01  3.2271E-02 -6.7556E-01  1.9381E+00 -2.0044E+00 -2.8415E+02
             1.6141E+00
 GRADIENT:  -7.2175E-01  1.9426E+01 -8.1163E+00  1.2832E+01 -1.7901E+00  1.3306E+00  2.2567E+00  7.4379E-01  3.5937E-03  0.0000E+00
            -5.3619E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1217.34785004511        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:     1368
 NPARAMETR:  1.0155E+00  1.9677E+00  2.0769E+00  5.3688E-01  2.3250E+00  9.2884E-01  4.6052E-01  6.2832E+00  1.1537E-01  1.0000E-02
             4.5471E+00
 PARAMETER:  1.1536E-01  7.7684E-01  8.3087E-01 -5.2197E-01  9.4373E-01  2.6176E-02 -6.7539E-01  1.9379E+00 -2.0596E+00 -2.8415E+02
             1.6145E+00
 GRADIENT:  -1.3540E+00  1.2251E+01  2.0275E+00 -2.0854E+00  5.6704E+00 -4.6957E-01  2.6959E+00 -8.1337E+00  2.2143E-02  0.0000E+00
             3.1329E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1217.36863902315        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1444
 NPARAMETR:  1.0142E+00  1.9654E+00  2.0770E+00  5.3686E-01  2.3246E+00  9.3291E-01  4.6051E-01  6.2927E+00  8.1112E-02  1.0000E-02
             4.5464E+00
 PARAMETER:  1.1412E-01  7.7570E-01  8.3093E-01 -5.2202E-01  9.4353E-01  3.0557E-02 -6.7542E-01  1.9394E+00 -2.4119E+00 -2.8415E+02
             1.6143E+00
 GRADIENT:  -3.3202E+00  1.3879E+01 -4.2718E+00  6.0835E+00  1.2111E+00  9.2907E-01  2.4297E+00 -2.5734E+00  5.4008E-03  0.0000E+00
            -1.7226E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1217.41112855895        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1537
 NPARAMETR:  1.0143E+00  1.9531E+00  2.0784E+00  5.3664E-01  2.3227E+00  9.2846E-01  4.6028E-01  6.3268E+00  2.0341E-02  1.0000E-02
             4.5496E+00
 PARAMETER:  1.1416E-01  7.6944E-01  8.3158E-01 -5.2243E-01  9.4274E-01  2.5777E-02 -6.7592E-01  1.9448E+00 -3.7951E+00 -2.8415E+02
             1.6150E+00
 GRADIENT:  -2.2465E+00 -8.7451E+00  4.5479E-01 -7.5946E+00  5.7053E+00 -1.0700E-01  2.7613E+00 -6.8828E+00  6.4834E-04  0.0000E+00
             2.7512E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1217.41579730899        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1717
 NPARAMETR:  1.0154E+00  1.9549E+00  2.0784E+00  5.3666E-01  2.3229E+00  9.2909E-01  4.6030E-01  6.3241E+00  2.4662E-02  1.0000E-02
             4.5497E+00
 PARAMETER:  1.1523E-01  7.7035E-01  8.3158E-01 -5.2238E-01  9.4280E-01  2.6446E-02 -6.7588E-01  1.9444E+00 -3.6025E+00 -2.8415E+02
             1.6151E+00
 GRADIENT:   2.9275E-01  6.6136E+00 -1.9724E+01  2.2588E+01 -9.1742E+00 -1.6214E-01  1.9187E+00  1.0836E+01 -6.1624E-04  0.0000E+00
            -1.4245E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1217.45798646992        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1898
 NPARAMETR:  1.0161E+00  1.9587E+00  2.0898E+00  5.3623E-01  2.3181E+00  9.2840E-01  4.5943E-01  6.4002E+00  1.0391E-02  1.0000E-02
             4.5906E+00
 PARAMETER:  1.1598E-01  7.7226E-01  8.3708E-01 -5.2320E-01  9.4073E-01  2.5709E-02 -6.7777E-01  1.9563E+00 -4.4668E+00 -2.8415E+02
             1.6240E+00
 GRADIENT:  -7.4093E-01 -5.7332E+00  1.4888E+00 -7.7257E+00  6.0782E+00  1.5968E-02  3.1907E+00 -7.6178E+00  2.1906E-04  0.0000E+00
             1.1972E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1217.61444241377        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     1993
 NPARAMETR:  1.0134E+00  1.9489E+00  2.0900E+00  5.3620E-01  2.3182E+00  9.2822E-01  4.1077E-01  6.3986E+00  1.0505E-02  1.0000E-02
             4.5912E+00
 PARAMETER:  1.1334E-01  7.6726E-01  8.3718E-01 -5.2325E-01  9.4080E-01  2.5510E-02 -7.8972E-01  1.9561E+00 -4.4559E+00 -2.8415E+02
             1.6241E+00
 GRADIENT:  -1.0123E+00 -1.0517E+01  1.3122E+00 -9.5250E+00  6.3171E+00  9.1415E-02  8.4864E-01 -7.0245E+00  6.8481E-05  0.0000E+00
             8.1548E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1217.62599570424        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     2063
 NPARAMETR:  1.0140E+00  1.9491E+00  2.0901E+00  5.3620E-01  2.3181E+00  9.2835E-01  4.0476E-01  6.3981E+00  1.0906E-02  1.0000E-02
             4.5912E+00
 PARAMETER:  1.1394E-01  7.6737E-01  8.3722E-01 -5.2324E-01  9.4077E-01  2.5655E-02 -8.0446E-01  1.9560E+00 -4.4184E+00 -2.8415E+02
             1.6241E+00
 GRADIENT:   5.5683E-01 -7.6600E+00 -2.0354E+00 -4.3183E+00  3.7790E+00  6.9606E-02  4.5700E-01 -4.0858E+00 -9.6651E-06  0.0000E+00
             4.4058E+00

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1217.62722762376        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     2177
 NPARAMETR:  1.0141E+00  1.9491E+00  2.0902E+00  5.3620E-01  2.3181E+00  9.2876E-01  4.0308E-01  6.3977E+00  1.1212E-02  1.0000E-02
             4.5912E+00
 PARAMETER:  1.1403E-01  7.6736E-01  8.3724E-01 -5.2324E-01  9.4074E-01  2.6093E-02 -8.0861E-01  1.9559E+00 -4.3907E+00 -2.8415E+02
             1.6242E+00
 GRADIENT:  -9.9147E-01  1.2686E+00 -1.8487E+01  1.9720E+01 -8.4666E+00 -5.7032E-02 -2.2392E-01  1.0011E+01 -3.2663E-04  0.0000E+00
            -1.0930E+01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1217.63318997959        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2354
 NPARAMETR:  1.0157E+00  1.9525E+00  2.0905E+00  5.3623E-01  2.3176E+00  9.2910E-01  4.1019E-01  6.3940E+00  1.5691E-02  1.0000E-02
             4.5926E+00
 PARAMETER:  1.1560E-01  7.6911E-01  8.3742E-01 -5.2319E-01  9.4052E-01  2.6462E-02 -7.9114E-01  1.9554E+00 -4.0547E+00 -2.8415E+02
             1.6245E+00
 GRADIENT:   1.3099E+00 -2.4028E+00 -9.1356E+00  7.0347E+00 -1.7149E+00  6.6127E-02  3.4014E-01  1.7512E+00 -2.8755E-04  0.0000E+00
            -2.3446E+00

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1217.64666228168        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2531
 NPARAMETR:  1.0155E+00  1.9567E+00  2.0924E+00  5.3631E-01  2.3155E+00  9.2959E-01  3.9558E-01  6.3803E+00  1.3039E-01  1.0000E-02
             4.5975E+00
 PARAMETER:  1.1540E-01  7.7126E-01  8.3833E-01 -5.2303E-01  9.3962E-01  2.6989E-02 -8.2741E-01  1.9532E+00 -1.9373E+00 -2.8415E+02
             1.6255E+00
 GRADIENT:  -5.4581E-01 -3.4084E+00  1.9038E+00 -6.4167E+00  5.9572E+00  3.9263E-02  2.0945E-01 -7.9723E+00 -2.5789E-03  0.0000E+00
             5.9412E+00

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1217.64905518026        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2707
 NPARAMETR:  1.0155E+00  1.9552E+00  2.0926E+00  5.3629E-01  2.3156E+00  9.2948E-01  3.9385E-01  6.3839E+00  1.6797E-01  1.0000E-02
             4.5968E+00
 PARAMETER:  1.1540E-01  7.7048E-01  8.3840E-01 -5.2309E-01  9.3967E-01  2.6868E-02 -8.3178E-01  1.9538E+00 -1.6840E+00 -2.8415E+02
             1.6254E+00
 GRADIENT:  -6.4633E-02 -5.9174E+00  2.1980E+00 -7.7760E+00  6.3404E+00  7.0070E-02  2.6076E-01 -8.2046E+00 -2.2305E-04  0.0000E+00
             6.3482E+00

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1217.69077792478        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2891
 NPARAMETR:  1.0151E+00  1.9526E+00  2.1047E+00  5.3579E-01  2.3101E+00  9.2785E-01  4.0163E-01  6.4357E+00  5.1226E-02  1.0000E-02
             4.5967E+00
 PARAMETER:  1.1497E-01  7.6918E-01  8.4415E-01 -5.2402E-01  9.3729E-01  2.5111E-02 -8.1222E-01  1.9619E+00 -2.8715E+00 -2.8415E+02
             1.6253E+00
 GRADIENT:  -8.3337E-02  3.5050E-01 -1.1787E+01  1.1367E+01 -3.8375E+00 -4.1221E-01 -4.8021E-02  4.2432E+00 -4.7020E-03  0.0000E+00
            -4.6213E+00

0ITERATION NO.:  145    OBJECTIVE VALUE:  -1218.42831408787        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3068
 NPARAMETR:  1.0164E+00  1.9580E+00  2.1926E+00  5.3332E-01  2.2637E+00  9.3070E-01  3.3113E-01  6.6628E+00  4.7181E-01  1.0000E-02
             4.6276E+00
 PARAMETER:  1.1630E-01  7.7191E-01  8.8507E-01 -5.2864E-01  9.1699E-01  2.8185E-02 -1.0052E+00  1.9965E+00 -6.5118E-01 -2.8415E+02
             1.6320E+00
 GRADIENT:   1.0794E+00 -3.1548E+00  2.0377E+00 -7.4208E+00  5.5909E+00  5.3091E-01 -2.4002E-01 -8.1168E+00  8.3872E-02  0.0000E+00
             1.1284E+01

0ITERATION NO.:  150    OBJECTIVE VALUE:  -1218.58487293923        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3245
 NPARAMETR:  1.0135E+00  1.9492E+00  2.2273E+00  5.3189E-01  2.2507E+00  9.2252E-01  3.6421E-01  6.8166E+00  3.3055E-01  1.0000E-02
             4.6181E+00
 PARAMETER:  1.1337E-01  7.6740E-01  9.0077E-01 -5.3132E-01  9.1122E-01  1.9350E-02 -9.1002E-01  2.0194E+00 -1.0070E+00 -2.8415E+02
             1.6300E+00
 GRADIENT:  -2.5871E+00 -1.3536E+01 -1.2778E+00 -7.4311E+00  4.0210E+00 -1.8572E+00  7.9541E-02 -4.8551E+00  1.3414E-02  0.0000E+00
             8.9677E+00

0ITERATION NO.:  155    OBJECTIVE VALUE:  -1218.64615562076        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3427
 NPARAMETR:  1.0136E+00  1.9496E+00  2.2493E+00  5.3124E-01  2.2404E+00  9.2270E-01  3.5699E-01  6.8796E+00  3.2921E-01  1.0000E-02
             4.6219E+00
 PARAMETER:  1.1350E-01  7.6762E-01  9.1063E-01 -5.3254E-01  9.0664E-01  1.9548E-02 -9.3005E-01  2.0286E+00 -1.0110E+00 -2.8415E+02
             1.6308E+00
 GRADIENT:  -4.0195E-01  8.5824E+00 -3.6364E+01  4.4341E+01 -2.2526E+01 -2.0800E+00 -6.8456E-01  2.7958E+01 -3.0645E-01  0.0000E+00
            -2.2066E+01

0ITERATION NO.:  156    OBJECTIVE VALUE:  -1218.64615562076        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:     3458
 NPARAMETR:  1.0135E+00  1.9482E+00  2.2514E+00  5.3097E-01  2.2424E+00  9.2279E-01  3.5732E-01  6.8664E+00  3.2955E-01  1.0000E-02
             4.6289E+00
 PARAMETER:  1.1350E-01  7.6762E-01  9.1063E-01 -5.3254E-01  9.0664E-01  1.9548E-02 -9.3005E-01  2.0286E+00 -1.0110E+00 -2.8415E+02
             1.6308E+00
 GRADIENT:   2.3651E+03  3.3473E+02 -1.4708E+02  4.9376E+02 -2.9065E+02 -1.3456E+03 -1.4477E+02  1.2531E+02 -1.3303E+02  0.0000E+00
            -1.5337E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3458
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.4587E-03 -7.3783E-03 -9.0271E-03 -2.3005E-03 -3.7941E-05
 SE:             2.8552E-02  1.2463E-02  7.4634E-03  4.6205E-03  2.9283E-05
 N:                     100         100         100         100         100

 P VAL.:         8.4838E-01  5.5383E-01  2.2646E-01  6.1857E-01  1.9509E-01

 ETASHRINKSD(%)  4.3457E+00  5.8248E+01  7.4997E+01  8.4521E+01  9.9902E+01
 ETASHRINKVR(%)  8.5025E+00  8.2568E+01  9.3748E+01  9.7604E+01  1.0000E+02
 EBVSHRINKSD(%)  4.7842E+00  5.8741E+01  8.6407E+01  8.4899E+01  9.9846E+01
 EBVSHRINKVR(%)  9.3395E+00  8.2977E+01  9.8152E+01  9.7720E+01  1.0000E+02
 RELATIVEINF(%)  8.4988E+01  3.8731E-01  1.3236E+00  5.3238E-02  1.3175E-04
 EPSSHRINKSD(%)  1.3812E+01
 EPSSHRINKVR(%)  2.5716E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1218.6461556207576     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -483.49532905701938     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.46
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.19
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1218.646       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.95E+00  2.25E+00  5.31E-01  2.24E+00  9.23E-01  3.57E-01  6.88E+00  3.29E-01  1.00E-02  4.62E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.14E+06
 
 TH 2
+        7.48E+02  2.61E-01
 
 TH 3
+       -2.53E+02 -8.85E-02  2.99E-02
 
 TH 4
+        1.43E+05  5.00E+01 -1.69E+01  9.58E+03
 
 TH 5
+       -2.57E+02 -9.00E-02  3.04E-02 -1.72E+01  3.10E-02
 
 TH 6
+       -4.14E+06 -1.45E+03  4.90E+02 -2.77E+05  4.98E+02  8.01E+06
 
 TH 7
+       -8.85E+04 -3.09E+01  1.05E+01 -5.92E+03  1.06E+01  1.71E+05  3.66E+03
 
 TH 8
+       -1.16E+01 -4.07E-03  1.38E-03 -7.79E-01  1.40E-03  2.25E+01  4.82E-01  6.33E-05
 
 TH 9
+       -8.78E+04 -3.07E+01  1.04E+01 -5.88E+03  1.06E+01  1.70E+05  3.63E+03  4.78E-01  3.61E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        4.45E+00  1.56E-03 -5.27E-04  2.98E-01 -5.36E-04 -8.61E+00 -1.84E-01 -2.42E-05 -1.83E-01  0.00E+00  9.27E-06
 
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
+        5.07E+06
 
 TH 2
+       -1.87E+02  3.00E+04
 
 TH 3
+        2.23E+01  1.92E+02  1.59E+04
 
 TH 4
+       -3.81E+02 -9.00E+02  6.05E+02  8.34E+05
 
 TH 5
+        3.19E+01  1.80E+02 -8.27E+01  6.28E+02  1.62E+04
 
 TH 6
+        3.47E+02 -2.13E+01 -3.07E+00 -1.43E+01 -3.10E+00  7.89E+06
 
 TH 7
+        4.81E+02 -1.97E+01  4.39E-01 -3.47E+01  1.70E+00 -7.03E+01  6.09E+05
 
 TH 8
+       -3.00E+00 -2.60E+01  8.96E+00 -5.90E+01  8.82E+00  3.59E-01 -4.65E-01  3.43E+02
 
 TH 9
+        6.66E+02  1.60E+00 -1.97E+00  1.14E+01 -2.27E+00 -9.60E+01 -1.89E+02  6.55E-02  6.06E+05
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -9.79E+00  3.54E+01 -1.99E+01  1.08E+02 -1.82E+01  3.57E+00  1.51E+01  4.80E+00  1.93E+00  0.00E+00  1.19E+03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.37E+06
 
 TH 2
+        1.05E+05  8.52E+03
 
 TH 3
+       -7.68E+04 -5.92E+03  4.32E+03
 
 TH 4
+        5.57E+05  4.35E+04 -3.13E+04  2.28E+05
 
 TH 5
+       -7.78E+04 -6.01E+03  4.38E+03 -3.17E+04  4.44E+03
 
 TH 6
+       -1.70E+06 -1.31E+05  9.56E+04 -6.93E+05  9.69E+04  2.12E+06
 
 TH 7
+       -4.73E+05 -3.65E+04  2.66E+04 -1.93E+05  2.69E+04  5.88E+05  1.64E+05
 
 TH 8
+        1.13E+04  8.70E+02 -6.36E+02  4.61E+03 -6.44E+02 -1.41E+04 -3.91E+03  9.37E+01
 
 TH 9
+       -4.72E+05 -3.64E+04  2.65E+04 -1.92E+05  2.69E+04  5.87E+05  1.63E+05 -3.90E+03  1.63E+05
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.07E+04 -1.64E+03  1.16E+03 -8.51E+03  1.18E+03  2.58E+04  7.17E+03 -1.71E+02  7.15E+03  0.00E+00  3.73E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       50.730
Stop Time:
Sat Sep 25 09:33:00 CDT 2021
