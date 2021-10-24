Sun Oct 24 03:38:43 CDT 2021
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
$DATA ../../../../data/SD4/SL3/dat90.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m90.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1581.84752146051        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6492E+02 -9.1307E+00 -3.9922E+01  5.4669E+01  9.2927E+01  2.8408E+01 -2.0933E+00  4.6950E+00 -3.9489E+00 -8.5785E+00
            -1.0223E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1592.88353622642        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.7955E-01  9.8921E-01  1.0259E+00  9.8581E-01  9.6245E-01  9.4704E-01  9.9410E-01  9.5320E-01  1.0272E+00  9.3401E-01
             1.2266E+00
 PARAMETER:  7.9342E-02  8.9154E-02  1.2555E-01  8.5711E-02  6.1725E-02  4.5584E-02  9.4087E-02  5.2068E-02  1.2686E-01  3.1735E-02
             3.0426E-01
 GRADIENT:   2.7692E+02 -2.1889E+01 -1.0794E+01 -7.7532E+00  3.5998E+01 -2.5899E+00 -8.4429E-01  2.5540E+00  1.3823E+00 -5.5166E+00
             2.6744E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1594.04225273212        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      246
 NPARAMETR:  9.8292E-01  9.7453E-01  1.0474E+00  1.0483E+00  9.2667E-01  1.0744E+00  1.0170E+00  8.5594E-01  1.0308E+00  9.7365E-01
             1.2913E+00
 PARAMETER:  8.2775E-02  7.4202E-02  1.4630E-01  1.4712E-01  2.3845E-02  1.7177E-01  1.1683E-01 -5.5558E-02  1.3031E-01  7.3296E-02
             3.5563E-01
 GRADIENT:   2.7952E+00  1.5638E+01  8.9028E+00  1.3107E+01 -2.3763E+01  2.1378E+01 -1.9061E-02  1.8367E-01  1.6861E+00  2.9606E+00
             2.2485E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1595.58691402104        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      422
 NPARAMETR:  9.8047E-01  9.9796E-01  1.0036E+00  1.0196E+00  9.3553E-01  1.0124E+00  1.0240E+00  8.1029E-01  1.0347E+00  9.6356E-01
             1.2275E+00
 PARAMETER:  8.0276E-02  9.7961E-02  1.0354E-01  1.1939E-01  3.3360E-02  1.1230E-01  1.2370E-01 -1.1036E-01  1.3414E-01  6.2878E-02
             3.0500E-01
 GRADIENT:  -1.7519E+00 -2.6151E-01  1.7205E-01  1.6553E-01  1.2088E-01 -1.0793E+00 -3.9556E-02 -1.6374E-01 -1.1128E-03 -3.1592E-01
             1.6369E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1595.65404657953        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      598
 NPARAMETR:  9.8332E-01  1.1927E+00  8.7160E-01  8.9549E-01  9.6016E-01  1.0178E+00  9.3336E-01  7.3180E-01  1.1286E+00  9.5939E-01
             1.2181E+00
 PARAMETER:  8.3181E-02  2.7620E-01 -3.7428E-02 -1.0381E-02  5.9343E-02  1.1769E-01  3.1030E-02 -2.1225E-01  2.2096E-01  5.8544E-02
             2.9727E-01
 GRADIENT:   1.6087E+00  6.4671E+00  1.7278E+00  5.0227E+00 -4.3527E+00  3.8023E-01 -9.3212E-02  5.6360E-02 -3.3172E-01 -1.6142E-02
            -1.7878E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1595.70421337185        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      774
 NPARAMETR:  9.8226E-01  1.3054E+00  7.9138E-01  8.1909E-01  9.7886E-01  1.0155E+00  8.9019E-01  6.5595E-01  1.1951E+00  9.5695E-01
             1.2264E+00
 PARAMETER:  8.2098E-02  3.6650E-01 -1.3398E-01 -9.9563E-02  7.8633E-02  1.1542E-01 -1.6322E-02 -3.2166E-01  2.7819E-01  5.5992E-02
             3.0412E-01
 GRADIENT:  -2.0412E+00  5.3811E+00  1.5618E+00  3.2040E+00 -3.6687E+00 -6.8995E-01 -7.4774E-02  4.1891E-02 -2.6298E-01 -1.0511E-01
             9.7779E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1595.73630202813        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      951
 NPARAMETR:  9.8388E-01  1.4297E+00  6.9147E-01  7.3655E-01  9.9621E-01  1.0194E+00  8.5324E-01  5.2012E-01  1.2734E+00  9.5598E-01
             1.2244E+00
 PARAMETER:  8.3753E-02  4.5747E-01 -2.6893E-01 -2.0578E-01  9.6205E-02  1.1922E-01 -5.8716E-02 -5.5369E-01  3.4170E-01  5.4978E-02
             3.0245E-01
 GRADIENT:   3.3188E-01  6.7984E+00  7.4150E-01  4.3285E+00 -2.8549E+00  5.5927E-01  4.5179E-02  8.9650E-02 -4.0509E-01  2.9361E-01
             3.4769E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1595.75312413995        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1131
 NPARAMETR:  9.8441E-01  1.4893E+00  6.4139E-01  6.8589E-01  1.0095E+00  1.0171E+00  8.3379E-01  3.6492E-01  1.3260E+00  9.5274E-01
             1.2232E+00
 PARAMETER:  8.4286E-02  4.9828E-01 -3.4412E-01 -2.7703E-01  1.0949E-01  1.1693E-01 -8.1773E-02 -9.0807E-01  3.8220E-01  5.1588E-02
             3.0146E-01
 GRADIENT:   1.3938E+00 -7.1807E+00  1.4782E+00 -5.3046E+00  3.3092E-01 -3.5018E-01 -2.4726E-01 -5.7621E-02 -6.1136E-01 -6.7551E-01
            -3.2752E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1595.77449240153        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1306
 NPARAMETR:  9.8377E-01  1.4980E+00  6.3504E-01  6.8614E-01  1.0090E+00  1.0178E+00  8.3409E-01  3.6840E-01  1.3320E+00  9.5393E-01
             1.2233E+00
 PARAMETER:  8.3633E-02  5.0413E-01 -3.5407E-01 -2.7667E-01  1.0898E-01  1.1761E-01 -8.1408E-02 -8.9859E-01  3.8669E-01  5.2837E-02
             3.0154E-01
 GRADIENT:  -2.3603E-01  3.5347E-01  4.6015E-01  9.8129E-01 -1.0402E-01 -1.4636E-01  1.7623E-01  6.5138E-03  2.4045E-01 -1.5943E-01
            -9.8309E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1595.77692808030        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1482
 NPARAMETR:  9.8381E-01  1.5323E+00  6.0019E-01  6.6355E-01  1.0081E+00  1.0168E+00  8.2831E-01  2.4906E-01  1.3594E+00  9.4214E-01
             1.2224E+00
 PARAMETER:  8.3675E-02  5.2680E-01 -4.1050E-01 -3.1014E-01  1.0803E-01  1.1665E-01 -8.8373E-02 -1.2901E+00  4.0706E-01  4.0398E-02
             3.0078E-01
 GRADIENT:  -4.1154E-01  3.1231E+00  9.6055E-01  1.9789E+00 -1.7885E+00 -5.8873E-01  5.3782E-01  1.8680E-02  7.6671E-01 -7.6906E-01
            -4.9732E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1595.77702002337        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1660
 NPARAMETR:  9.8384E-01  1.5440E+00  5.8911E-01  6.5566E-01  1.0084E+00  1.0167E+00  8.2555E-01  2.1080E-01  1.3683E+00  9.3948E-01
             1.2222E+00
 PARAMETER:  8.3712E-02  5.3440E-01 -4.2914E-01 -3.2211E-01  1.0840E-01  1.1659E-01 -9.1702E-02 -1.4568E+00  4.1358E-01  3.7576E-02
             3.0063E-01
 GRADIENT:  -4.0435E-01  3.5014E+00  8.6888E-01  2.0924E+00 -2.1330E+00 -6.2764E-01  5.8450E-01  2.4378E-02  8.2273E-01 -8.2846E-01
            -5.3154E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1595.79828930556        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1841
 NPARAMETR:  9.8397E-01  1.5438E+00  5.8619E-01  6.5360E-01  1.0103E+00  1.0182E+00  8.2301E-01  1.3041E-01  1.3642E+00  9.4716E-01
             1.2231E+00
 PARAMETER:  8.3841E-02  5.3427E-01 -4.3412E-01 -3.2527E-01  1.1024E-01  1.1799E-01 -9.4788E-02 -1.9371E+00  4.1058E-01  4.5712E-02
             3.0137E-01
 GRADIENT:  -1.4782E-01 -1.3374E+00 -1.0117E-01  7.5829E-01  7.5205E-01 -6.1906E-02  6.4240E-02  6.2546E-03  7.2005E-02 -9.2498E-03
             6.2220E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1595.80146729722        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2020             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8474E-01  1.5504E+00  5.7860E-01  6.4719E-01  1.0114E+00  1.0186E+00  8.1969E-01  3.7710E-02  1.3709E+00  9.4617E-01
             1.2230E+00
 PARAMETER:  8.4620E-02  5.3854E-01 -4.4715E-01 -3.3511E-01  1.1136E-01  1.1840E-01 -9.8833E-02 -3.1778E+00  4.1543E-01  4.4671E-02
             3.0131E-01
 GRADIENT:   2.8527E+02  3.0357E+02  2.5364E+00  6.4556E+01  7.7465E+00  3.1532E+01  3.9378E+00  3.4468E-03  1.4895E+01  4.4838E-01
             2.5116E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1595.80363215402        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2202             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8455E-01  1.5514E+00  5.7857E-01  6.4739E-01  1.0103E+00  1.0186E+00  8.2050E-01  1.0000E-02  1.3708E+00  9.4590E-01
             1.2228E+00
 PARAMETER:  8.4432E-02  5.3916E-01 -4.4720E-01 -3.3481E-01  1.1026E-01  1.1838E-01 -9.7838E-02 -5.7317E+00  4.1543E-01  4.4378E-02
             3.0111E-01
 GRADIENT:   2.8497E+02  3.0639E+02  3.2378E+00  6.4811E+01  5.7790E+00  3.1551E+01  4.1096E+00  0.0000E+00  1.4904E+01  5.0457E-01
             2.3893E+00

0ITERATION NO.:   67    OBJECTIVE VALUE:  -1595.80363215402        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:     2260
 NPARAMETR:  9.8455E-01  1.5514E+00  5.7857E-01  6.4739E-01  1.0103E+00  1.0186E+00  8.2050E-01  1.0000E-02  1.3708E+00  9.4590E-01
             1.2228E+00
 PARAMETER:  8.4432E-02  5.3916E-01 -4.4720E-01 -3.3481E-01  1.1026E-01  1.1838E-01 -9.7838E-02 -5.7317E+00  4.1543E-01  4.4378E-02
             3.0111E-01
 GRADIENT:   1.0992E+00 -2.2470E+00  1.7909E-01 -2.8861E-01 -3.6970E-02  8.7808E-02 -3.4270E-03  0.0000E+00  4.5737E-02  2.3454E-02
            -1.9654E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2260
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6695E-04 -2.7394E-02 -2.7195E-04  1.8863E-02 -3.2228E-02
 SE:             2.9781E-02  2.2204E-02  1.1514E-04  2.3570E-02  2.2100E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9553E-01  2.1729E-01  1.8187E-02  4.2354E-01  1.4478E-01

 ETASHRINKSD(%)  2.3066E-01  2.5615E+01  9.9614E+01  2.1038E+01  2.5961E+01
 ETASHRINKVR(%)  4.6078E-01  4.4668E+01  9.9999E+01  3.7649E+01  4.5182E+01
 EBVSHRINKSD(%)  6.1127E-01  2.4957E+01  9.9670E+01  2.2139E+01  2.5044E+01
 EBVSHRINKVR(%)  1.2188E+00  4.3685E+01  9.9999E+01  3.9376E+01  4.3816E+01
 RELATIVEINF(%)  9.8717E+01  3.4925E+00  1.3911E-04  4.2947E+00  9.1581E+00
 EPSSHRINKSD(%)  4.2838E+01
 EPSSHRINKVR(%)  6.7325E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1595.8036321540176     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -860.65280559027940     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1595.804       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.55E+00  5.79E-01  6.47E-01  1.01E+00  1.02E+00  8.21E-01  1.00E-02  1.37E+00  9.46E-01  1.22E+00
 


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
 #CPUT: Total CPU Time in Seconds,       64.745
Stop Time:
Sun Oct 24 03:38:56 CDT 2021
