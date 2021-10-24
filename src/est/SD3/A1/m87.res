Sat Oct 23 21:58:04 CDT 2021
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
$DATA ../../../../data/SD3/A1/dat87.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1455.36155623331        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1276E+02 -1.5723E+01  4.8949E+01  8.7593E+00  1.7046E+02  4.4645E+01 -5.5422E+01 -2.9652E+02 -1.0184E+02 -3.8160E+01
            -8.3291E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1830.92315944610        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.6934E-01  9.9327E-01  1.0182E+00  1.0471E+00  8.7636E-01  9.7372E-01  1.0884E+00  9.4381E-01  1.0890E+00  1.0051E+00
             1.7038E+00
 PARAMETER:  6.8864E-02  9.3251E-02  1.1805E-01  1.4599E-01 -3.1975E-02  7.3371E-02  1.8470E-01  4.2170E-02  1.8523E-01  1.0508E-01
             6.3286E-01
 GRADIENT:   3.0514E+01  4.3191E+01  1.4088E+01  4.3623E+01 -2.4715E+01  1.9723E+00 -9.2464E+00  1.3656E+01 -1.0654E+00  4.9014E+00
            -7.9173E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1838.37344032859        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      203
 NPARAMETR:  9.8126E-01  6.5078E-01  7.1610E-01  1.2850E+00  6.4442E-01  9.9656E-01  1.5026E+00  8.5818E-01  1.0610E+00  6.9503E-01
             1.6717E+00
 PARAMETER:  8.1081E-02 -3.2959E-01 -2.3394E-01  3.5079E-01 -3.3941E-01  9.6555E-02  5.0721E-01 -5.2941E-02  1.5921E-01 -2.6379E-01
             6.1385E-01
 GRADIENT:  -9.4444E+01  1.9882E+01 -5.1514E+01  1.0428E+02  8.1370E+01 -7.9450E+00 -1.7743E+00  2.3904E+01  7.5309E+00 -5.5631E+00
            -8.7419E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1856.05311732750        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      381
 NPARAMETR:  1.0203E+00  4.6918E-01  7.5939E-01  1.3535E+00  5.8928E-01  9.9805E-01  1.7666E+00  5.9412E-01  9.3536E-01  6.5381E-01
             1.9827E+00
 PARAMETER:  1.2005E-01 -6.5678E-01 -1.7524E-01  4.0269E-01 -4.2886E-01  9.8046E-02  6.6908E-01 -4.2068E-01  3.3174E-02 -3.2493E-01
             7.8447E-01
 GRADIENT:  -9.5104E+00  2.2188E+01  1.2318E+01  5.4969E+01 -1.5084E+01 -6.6764E-01 -3.7705E+00  1.0974E+01 -9.7433E+00 -6.8593E-01
             6.1092E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1870.41048507549        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  1.0142E+00  1.9250E-01  6.2942E-01  1.4244E+00  4.7840E-01  9.8869E-01  3.2462E+00  2.4161E-01  9.2018E-01  6.4221E-01
             1.8141E+00
 PARAMETER:  1.1412E-01 -1.5476E+00 -3.6296E-01  4.5373E-01 -6.3732E-01  8.8626E-02  1.2775E+00 -1.3204E+00  1.6817E-02 -3.4284E-01
             6.9558E-01
 GRADIENT:  -1.9609E+00  2.7178E+00 -5.4713E+00  1.4584E+01  6.5916E+00 -2.0250E+00  3.3317E+00  1.6673E+00 -4.6916E+00  6.2218E-01
             1.1872E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1871.01028842652        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  1.0125E+00  1.3164E-01  6.1051E-01  1.4341E+00  4.5876E-01  9.9342E-01  3.7495E+00  1.4984E-01  9.2909E-01  6.3638E-01
             1.8154E+00
 PARAMETER:  1.1242E-01 -1.9277E+00 -3.9346E-01  4.6056E-01 -6.7924E-01  9.3401E-02  1.4216E+00 -1.7982E+00  2.6450E-02 -3.5196E-01
             6.9629E-01
 GRADIENT:  -6.5529E-01 -4.2374E-01  1.6977E+00 -2.9554E+00 -2.8190E+00  7.6474E-02 -5.4654E-02  5.6366E-01  5.3880E-01  6.6600E-02
             8.0572E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1871.05168106307        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      910            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0120E+00  1.1440E-01  5.9639E-01  1.4379E+00  4.4947E-01  9.9483E-01  3.9698E+00  1.1716E-01  9.3025E-01  6.2883E-01
             1.8148E+00
 PARAMETER:  1.1190E-01 -2.0680E+00 -4.1685E-01  4.6319E-01 -6.9969E-01  9.4813E-02  1.4787E+00 -2.0442E+00  2.7699E-02 -3.6389E-01
             6.9600E-01
 GRADIENT:   1.3953E+02  7.8331E+00  1.5854E+01  2.7536E+02  5.0174E+01  1.3940E+01  1.1756E+01  4.3976E-01  9.6733E+00  2.6836E+00
             8.4611E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1871.20268214367        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1087
 NPARAMETR:  1.0122E+00  1.1427E-01  5.9703E-01  1.4393E+00  4.5032E-01  9.9381E-01  3.9580E+00  3.2474E-02  9.2741E-01  6.3010E-01
             1.8165E+00
 PARAMETER:  1.1212E-01 -2.0692E+00 -4.1578E-01  4.6414E-01 -6.9779E-01  9.3791E-02  1.4757E+00 -3.3273E+00  2.4638E-02 -3.6187E-01
             6.9693E-01
 GRADIENT:  -1.1655E-02 -1.6178E+00  1.2401E+00  3.9569E+00 -3.4316E+00  2.7518E-01 -3.1843E+00  2.7182E-02  2.4645E-01  2.7333E-01
             1.4296E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1871.23459201745        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  1.0125E+00  1.1799E-01  6.0735E-01  1.4420E+00  4.5575E-01  9.9339E-01  3.9535E+00  1.0000E-02  9.3009E-01  6.4432E-01
             1.8144E+00
 PARAMETER:  1.1244E-01 -2.0372E+00 -3.9865E-01  4.6606E-01 -6.8581E-01  9.3363E-02  1.4746E+00 -1.0933E+01  2.7525E-02 -3.3956E-01
             6.9578E-01
 GRADIENT:   6.4801E-01  1.1684E-01 -2.4749E-01  1.5743E+00 -1.7570E+00  6.4396E-02  5.7839E-01  0.0000E+00  9.7829E-01  9.2663E-01
             3.1854E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1871.23961815052        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1440
 NPARAMETR:  1.0125E+00  1.1979E-01  6.1258E-01  1.4428E+00  4.5896E-01  9.9330E-01  3.9252E+00  1.0000E-02  9.2316E-01  6.4550E-01
             1.8125E+00
 PARAMETER:  1.1239E-01 -2.0220E+00 -3.9008E-01  4.6656E-01 -6.7879E-01  9.3276E-02  1.4674E+00 -1.3357E+01  2.0046E-02 -3.3773E-01
             6.9470E-01
 GRADIENT:   3.9697E-01 -2.0165E-01 -5.3708E-01  1.3855E+00 -3.0893E-01  5.3890E-02  8.2766E-02  0.0000E+00 -1.2303E+00  7.0178E-01
            -7.3722E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1871.24866081473        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1615
 NPARAMETR:  1.0123E+00  1.2334E-01  6.1959E-01  1.4423E+00  4.6319E-01  9.9302E-01  3.8821E+00  1.0000E-02  9.2589E-01  6.4517E-01
             1.8145E+00
 PARAMETER:  1.1218E-01 -1.9928E+00 -3.7870E-01  4.6627E-01 -6.6963E-01  9.2997E-02  1.4564E+00 -1.6255E+01  2.3002E-02 -3.3824E-01
             6.9583E-01
 GRADIENT:  -3.4900E-01 -3.4894E-01  1.3434E-01 -1.0496E+00  9.9509E-01 -3.5725E-02  1.9180E-01  0.0000E+00  6.6179E-02 -7.2001E-02
             3.3023E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1871.26368134344        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1794
 NPARAMETR:  1.0145E+00  1.3764E-01  6.2261E-01  1.4421E+00  4.6630E-01  9.9419E-01  3.6911E+00  1.0000E-02  9.2696E-01  6.4330E-01
             1.8211E+00
 PARAMETER:  1.1442E-01 -1.8831E+00 -3.7384E-01  4.6608E-01 -6.6294E-01  9.4169E-02  1.4059E+00 -1.9308E+01  2.4155E-02 -3.4114E-01
             6.9946E-01
 GRADIENT:   2.9684E+00 -8.5638E-01  3.1553E+00  8.7672E+00 -3.5823E+00  4.1439E-01 -1.7901E+00  0.0000E+00  4.2952E-01  2.6770E-02
             3.4354E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1871.30869326015        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1970
 NPARAMETR:  1.0148E+00  1.4431E-01  6.2272E-01  1.4387E+00  4.6755E-01  9.9412E-01  3.6314E+00  1.0000E-02  9.2905E-01  6.4259E-01
             1.8176E+00
 PARAMETER:  1.1474E-01 -1.8358E+00 -3.7366E-01  4.6371E-01 -6.6025E-01  9.4104E-02  1.3896E+00 -1.8122E+01  2.6409E-02 -3.4225E-01
             6.9754E-01
 GRADIENT:   3.3329E+00 -1.4584E-01  5.0403E-01  6.8902E+00  6.7140E-02  3.4994E-01  1.0811E-01  0.0000E+00  3.7968E-01 -7.0364E-01
             1.5235E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1871.36120623156        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2147
 NPARAMETR:  1.0135E+00  1.4819E-01  6.1987E-01  1.4312E+00  4.6669E-01  9.9323E-01  3.5783E+00  1.0000E-02  9.2962E-01  6.4567E-01
             1.8135E+00
 PARAMETER:  1.1343E-01 -1.8092E+00 -3.7824E-01  4.5848E-01 -6.6209E-01  9.3207E-02  1.3749E+00 -1.5477E+01  2.7016E-02 -3.3746E-01
             6.9528E-01
 GRADIENT:   2.9520E-01 -2.0219E-01  1.3763E-01 -9.3508E-01  1.0965E+00  2.7042E-02  6.0048E-01  0.0000E+00 -4.6685E-02 -9.2939E-02
             2.3948E-01

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1871.36120623156        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2169
 NPARAMETR:  1.0135E+00  1.4819E-01  6.1987E-01  1.4312E+00  4.6669E-01  9.9323E-01  3.5783E+00  1.0000E-02  9.2962E-01  6.4567E-01
             1.8135E+00
 PARAMETER:  1.1343E-01 -1.8092E+00 -3.7824E-01  4.5848E-01 -6.6209E-01  9.3207E-02  1.3749E+00 -1.5477E+01  2.7016E-02 -3.3746E-01
             6.9528E-01
 GRADIENT:   2.9520E-01 -2.0219E-01  1.3763E-01 -9.3508E-01  1.0965E+00  2.7042E-02  6.0048E-01  0.0000E+00 -4.6685E-02 -9.2939E-02
             2.3948E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2169
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5092E-03  3.2141E-02 -1.9468E-04 -1.8967E-02  4.0901E-03
 SE:             2.9532E-02  1.4922E-02  2.1206E-04  2.7656E-02  2.0661E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5924E-01  3.1244E-02  3.5860E-01  4.9283E-01  8.4307E-01

 ETASHRINKSD(%)  1.0638E+00  5.0010E+01  9.9290E+01  7.3487E+00  3.0784E+01
 ETASHRINKVR(%)  2.1164E+00  7.5010E+01  9.9995E+01  1.4157E+01  5.2091E+01
 EBVSHRINKSD(%)  1.0883E+00  6.3152E+01  9.9162E+01  5.7380E+00  2.6657E+01
 EBVSHRINKVR(%)  2.1647E+00  8.6422E+01  9.9993E+01  1.1147E+01  4.6208E+01
 RELATIVEINF(%)  9.7092E+01  5.2250E+00  5.2580E-04  4.1033E+01  4.0408E+00
 EPSSHRINKSD(%)  2.9448E+01
 EPSSHRINKVR(%)  5.0224E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1871.3612062315565     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -952.42267302688379     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1871.361       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.48E-01  6.20E-01  1.43E+00  4.67E-01  9.93E-01  3.58E+00  1.00E-02  9.30E-01  6.46E-01  1.81E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,      182.980
Stop Time:
Sat Oct 23 21:58:30 CDT 2021
