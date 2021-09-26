Sat Sep 25 01:40:10 CDT 2021
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
$DATA ../../../../data/int/SL2/dat86.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      996
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      896
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1880.85426726791        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.5112E+01  8.5691E+01  1.3504E+02 -1.1460E+02  1.0851E+02  5.2048E+00 -8.1247E+01 -5.9089E+02 -1.3887E+02 -4.6254E+01
            -3.1848E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3077.60916501369        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0257E+00  9.9548E-01  1.2724E+00  1.0409E+00  1.0625E+00  9.5123E-01  1.0381E+00  1.0286E+00  9.2511E-01  1.0198E+00
             1.9656E+00
 PARAMETER:  1.2537E-01  9.5474E-02  3.4094E-01  1.4006E-01  1.6063E-01  5.0001E-02  1.3744E-01  1.2822E-01  2.2153E-02  1.1963E-01
             7.7579E-01
 GRADIENT:  -3.1984E+01 -3.0670E+01  5.6348E+00 -5.1668E+00  2.2893E+01 -1.3187E+01 -4.4566E+00 -4.7617E+00 -1.0103E+01 -7.8253E+00
            -1.2351E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3078.91501248478        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0251E+00  1.1614E+00  1.3197E+00  9.6760E-01  1.1846E+00  9.3384E-01  8.5671E-01  9.0600E-01  9.6286E-01  1.1775E+00
             2.0012E+00
 PARAMETER:  1.2474E-01  2.4963E-01  3.7737E-01  6.7068E-02  2.6940E-01  3.1550E-02 -5.4659E-02  1.2827E-03  6.2151E-02  2.6338E-01
             7.9377E-01
 GRADIENT:  -3.7949E+01  1.0862E+01  8.7835E+00  1.4624E+01  1.9141E+01 -2.1109E+01 -5.8665E+00 -7.3117E+00 -1.1304E+01 -2.3068E+00
            -9.7811E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3083.45258004263        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.0412E+00  1.1079E+00  1.2599E+00  9.9044E-01  1.1132E+00  9.8211E-01  9.3784E-01  1.1698E+00  9.8572E-01  1.1203E+00
             2.0729E+00
 PARAMETER:  1.4036E-01  2.0247E-01  3.3107E-01  9.0392E-02  2.0721E-01  8.1953E-02  3.5824E-02  2.5684E-01  8.5616E-02  2.1360E-01
             8.2895E-01
 GRADIENT:   2.0315E+00  7.5087E-01 -1.3635E+00  4.8945E-01 -9.0003E-01  4.1358E-01  2.5491E-01  1.2348E+00  1.3402E-01  2.9507E-01
             1.9541E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3083.64517845454        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      391
 NPARAMETR:  1.0479E+00  1.1628E+00  1.3403E+00  9.6238E-01  1.1767E+00  9.8481E-01  8.9557E-01  1.2739E+00  1.0025E+00  1.1578E+00
             2.0736E+00
 PARAMETER:  1.4676E-01  2.5087E-01  3.9293E-01  6.1653E-02  2.6268E-01  8.4689E-02 -1.0297E-02  3.4205E-01  1.0251E-01  2.4649E-01
             8.2929E-01
 GRADIENT:   6.9793E-01 -4.1661E-01 -4.9194E-01 -5.2812E-01  7.1764E-01  1.6679E-01  1.5486E-02 -1.3140E-02  9.1874E-02 -8.6331E-02
            -1.6988E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3087.32243831873        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      570
 NPARAMETR:  1.0489E+00  1.2095E+00  2.1006E+00  9.5621E-01  1.3796E+00  9.8398E-01  8.7891E-01  2.7684E+00  9.6849E-01  1.1834E+00
             2.0536E+00
 PARAMETER:  1.4770E-01  2.9021E-01  8.4223E-01  5.5226E-02  4.2183E-01  8.3848E-02 -2.9072E-02  1.1183E+00  6.7984E-02  2.6836E-01
             8.1959E-01
 GRADIENT:   3.6967E+00  6.5805E+00 -1.0115E+01 -8.0769E+00  1.2077E+01 -7.6269E-01  3.9540E+00  6.0974E+00  3.8900E+00 -2.3727E+00
             1.7197E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3090.37449160253        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  1.0436E+00  1.0403E+00  2.4316E+00  1.0690E+00  1.2865E+00  9.8957E-01  9.6747E-01  2.7672E+00  8.7018E-01  1.1096E+00
             2.0264E+00
 PARAMETER:  1.4270E-01  1.3950E-01  9.8856E-01  1.6674E-01  3.5191E-01  8.9520E-02  6.6927E-02  1.1179E+00 -3.9051E-02  2.0397E-01
             8.0626E-01
 GRADIENT:  -6.1039E+00  8.8082E+00 -8.1471E-01 -6.3570E+00 -1.2286E+01  1.5595E+00  7.2294E-01  1.7928E+00  2.2962E-01  3.0520E+00
            -1.0660E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3090.74468375921        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  1.0376E+00  1.0198E+00  2.4315E+00  1.0801E+00  1.2905E+00  9.8030E-01  9.5393E-01  2.7671E+00  8.7003E-01  1.0756E+00
             2.0265E+00
 PARAMETER:  1.3689E-01  1.1963E-01  9.8853E-01  1.7704E-01  3.5503E-01  8.0099E-02  5.2831E-02  1.1178E+00 -3.9225E-02  1.7285E-01
             8.0630E-01
 GRADIENT:  -3.3469E+00  3.7086E+00 -3.1806E+00 -2.9782E+00  5.7346E+00 -8.4558E-01 -1.6274E-01  3.2261E-01  2.3374E-01 -1.4256E+00
             2.6900E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3090.98866692651        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      942
 NPARAMETR:  1.0396E+00  9.5013E-01  2.4315E+00  1.1228E+00  1.2466E+00  9.8469E-01  1.0137E+00  2.7670E+00  8.4055E-01  1.0341E+00
             2.0265E+00
 PARAMETER:  1.3888E-01  4.8841E-02  9.8851E-01  2.1586E-01  3.2044E-01  8.4576E-02  1.1361E-01  1.1178E+00 -7.3696E-02  1.3350E-01
             8.0631E-01
 GRADIENT:   1.9400E+00 -2.5415E+00 -5.6322E+00 -9.9088E-01  5.0776E-01  1.1096E+00 -9.2827E-02  4.8730E-01  2.0812E-01  1.1320E-01
             7.0334E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3090.98898935086        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  1.0394E+00  9.5462E-01  2.4315E+00  1.1204E+00  1.2489E+00  9.8368E-01  1.0101E+00  2.7670E+00  8.4196E-01  1.0369E+00
             2.0265E+00
 PARAMETER:  1.3861E-01  5.3555E-02  9.8851E-01  2.1367E-01  3.2225E-01  8.3546E-02  1.1007E-01  1.1178E+00 -7.2020E-02  1.3625E-01
             8.0630E-01
 GRADIENT:   1.2180E+00 -1.6995E+00 -5.4297E+00 -6.6812E-01  3.6686E-01  6.9768E-01 -7.6371E-02  4.8893E-01  1.3481E-01  8.1117E-02
             6.7877E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3091.05909839422        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1190
 NPARAMETR:  1.0466E+00  9.5984E-01  2.4316E+00  1.1212E+00  1.2548E+00  9.8481E-01  1.0104E+00  2.7671E+00  8.4204E-01  1.0422E+00
             2.0264E+00
 PARAMETER:  1.4557E-01  5.9009E-02  9.8855E-01  2.1438E-01  3.2702E-01  8.4695E-02  1.1030E-01  1.1178E+00 -7.1928E-02  1.4137E-01
             8.0627E-01
 GRADIENT:   4.2751E-01 -9.6875E-02 -7.5373E+00 -1.0518E-02  1.0939E-01 -1.2370E-01 -4.8649E-02 -5.2745E-01  9.7457E-04  5.4489E-02
             4.4749E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3091.12262978788        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1327
 NPARAMETR:  1.0461E+00  9.5945E-01  2.4535E+00  1.1213E+00  1.2546E+00  9.8500E-01  1.0115E+00  2.7697E+00  8.4178E-01  1.0414E+00
             2.0254E+00
 PARAMETER:  1.4506E-01  5.8604E-02  9.9751E-01  2.1445E-01  3.2682E-01  8.4887E-02  1.1145E-01  1.1187E+00 -7.2239E-02  1.4055E-01
             8.0577E-01
 GRADIENT:   1.6462E+01  1.9429E+00 -4.5937E+00  5.4633E+00  1.4708E+00  1.2804E+00  2.6340E-02  2.0454E-01  2.7205E-01  1.1944E-01
             4.9645E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3091.12868109641        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1492
 NPARAMETR:  1.0462E+00  9.5904E-01  2.4551E+00  1.1217E+00  1.2559E+00  9.8507E-01  1.0128E+00  2.7700E+00  8.4139E-01  1.0417E+00
             2.0255E+00
 PARAMETER:  1.4520E-01  5.8182E-02  9.9818E-01  2.1486E-01  3.2789E-01  8.4961E-02  1.1269E-01  1.1188E+00 -7.2696E-02  1.4081E-01
             8.0582E-01
 GRADIENT:  -3.7115E-01  2.0210E-01 -6.2944E+00 -6.9343E-01 -5.5700E-01 -3.2184E-02 -7.5626E-02 -7.9314E-01  1.1525E-01 -4.4908E-02
             3.6081E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3091.12992935209        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1670
 NPARAMETR:  1.0464E+00  9.5622E-01  2.4552E+00  1.1237E+00  1.2554E+00  9.8511E-01  1.0177E+00  2.7700E+00  8.3923E-01  1.0395E+00
             2.0255E+00
 PARAMETER:  1.4534E-01  5.5234E-02  9.9819E-01  2.1666E-01  3.2743E-01  8.5000E-02  1.1756E-01  1.1188E+00 -7.5270E-02  1.3874E-01
             8.0582E-01
 GRADIENT:  -3.4244E-02 -1.6740E-02 -6.6484E+00 -1.0405E-01  3.0319E-01 -1.5738E-02  2.9373E-03 -8.0165E-01  2.5395E-02 -9.6454E-02
             3.6943E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3091.13313206694        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1827
 NPARAMETR:  1.0463E+00  9.5497E-01  2.4546E+00  1.1245E+00  1.2543E+00  9.8508E-01  1.0187E+00  2.7789E+00  8.3859E-01  1.0392E+00
             2.0236E+00
 PARAMETER:  1.4522E-01  5.3928E-02  9.9795E-01  2.1732E-01  3.2661E-01  8.4965E-02  1.1854E-01  1.1221E+00 -7.6030E-02  1.3845E-01
             8.0490E-01
 GRADIENT:  -2.5517E-01 -3.5092E-02 -6.8084E+00 -1.3892E-01  3.6532E-03 -3.8115E-02  1.1966E-03 -4.6539E-01 -2.9459E-02 -1.7100E-02
             2.0698E+00

0ITERATION NO.:   71    OBJECTIVE VALUE:  -3091.13313206694        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     1854
 NPARAMETR:  1.0463E+00  9.5499E-01  2.4541E+00  1.1245E+00  1.2543E+00  9.8513E-01  1.0187E+00  2.7783E+00  8.3868E-01  1.0393E+00
             2.0240E+00
 PARAMETER:  1.4522E-01  5.3928E-02  9.9795E-01  2.1732E-01  3.2661E-01  8.4965E-02  1.1854E-01  1.1221E+00 -7.6030E-02  1.3845E-01
             8.0490E-01
 GRADIENT:  -1.8325E-01 -3.6080E-02  1.8584E+03 -1.3775E-01  5.7032E+03 -2.8251E-02  7.7870E-04  1.6625E+03 -2.4932E-02 -1.3778E-02
            -2.3203E+03
 NUMSIGDIG:         2.9         3.4         3.3         3.5         3.3         2.8         4.0         3.3         2.6         2.9
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1854
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4250E-03 -2.1221E-02 -2.1870E-02  8.7626E-03 -3.7145E-02
 SE:             2.9639E-02  1.8417E-02  2.1999E-02  2.4854E-02  2.1520E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6165E-01  2.4922E-01  3.2014E-01  7.2442E-01  8.4342E-02

 ETASHRINKSD(%)  7.0669E-01  3.8302E+01  2.6302E+01  1.6736E+01  2.7904E+01
 ETASHRINKVR(%)  1.4084E+00  6.1933E+01  4.5686E+01  3.0672E+01  4.8021E+01
 EBVSHRINKSD(%)  9.8497E-01  3.8707E+01  2.7522E+01  1.7943E+01  2.6587E+01
 EBVSHRINKVR(%)  1.9602E+00  6.2432E+01  4.7469E+01  3.2666E+01  4.6105E+01
 RELATIVEINF(%)  9.8010E+01  9.6510E+00  3.1272E+01  1.9611E+01  2.8022E+01
 EPSSHRINKSD(%)  1.8843E+01
 EPSSHRINKVR(%)  3.4135E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          896
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1646.7378515027735     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3091.1331320669356     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1444.3952805641620     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    48.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3091.133       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  9.55E-01  2.45E+00  1.12E+00  1.25E+00  9.85E-01  1.02E+00  2.78E+00  8.39E-01  1.04E+00  2.02E+00
 


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
+        1.03E+03
 
 TH 2
+       -4.62E+00  5.11E+07
 
 TH 3
+        5.71E+01  2.48E+02  7.74E+04
 
 TH 4
+       -1.23E+01  4.55E+02 -1.22E+02  8.26E+02
 
 TH 5
+        3.36E+02 -1.19E+07 -1.01E+03 -4.65E+06  2.77E+06
 
 TH 6
+        6.49E+00 -9.95E+00  2.37E+01 -3.95E+00  1.40E+02  2.00E+02
 
 TH 7
+        5.01E+00  1.35E+01  9.33E+01 -5.66E+00  5.61E+02 -1.66E+00  3.05E+01
 
 TH 8
+        4.47E+01  1.82E+02  7.93E+01 -8.54E+01 -7.68E+02  1.87E+01  7.58E+01  4.80E+04
 
 TH 9
+        7.18E-01 -2.11E+01 -4.42E+01  2.78E+01 -2.63E+02  2.55E-01  3.26E+01 -3.53E+01  1.35E+02
 
 TH10
+       -7.06E-01 -3.39E+07 -2.00E+02  1.05E+01  7.90E+06 -1.19E+00  4.46E+00 -1.56E+02  5.52E+00  2.25E+07
 
 TH11
+       -9.77E+01 -3.67E+02 -5.80E+02  1.46E+02  1.46E+03 -3.27E+01 -1.37E+02 -1.16E+02  7.70E+01  3.12E+02  1.77E+05
 
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
 #CPUT: Total CPU Time in Seconds,       63.270
Stop Time:
Sat Sep 25 01:41:15 CDT 2021
