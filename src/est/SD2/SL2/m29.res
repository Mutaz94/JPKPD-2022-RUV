Sat Oct 23 18:53:19 CDT 2021
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
$DATA ../../../../data/SD2/SL2/dat29.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2756.26969209809        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7985E+02 -6.6760E+01 -5.7442E+00  9.1773E+01  1.2947E+02  1.7570E+01 -5.2633E+01 -1.4245E+01 -2.2810E+01 -1.2014E+01
            -2.9660E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2792.23520008699        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8870E-01  1.1078E+00  8.2490E-01  9.1691E-01  9.4424E-01  1.1583E+00  1.3145E+00  1.0249E+00  1.1477E+00  9.8463E-01
             1.2695E+00
 PARAMETER:  8.8635E-02  2.0239E-01 -9.2495E-02  1.3255E-02  4.2623E-02  2.4694E-01  3.7349E-01  1.2455E-01  2.3773E-01  8.4507E-02
             3.3861E-01
 GRADIENT:   2.7914E+02  4.2142E+01 -2.9940E+01  2.4019E+01  3.1870E+01  1.0769E+02  3.1878E+01  8.9148E+00  2.2419E+01  1.2584E+01
             1.0406E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2798.98287180680        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      264
 NPARAMETR:  9.9584E-01  1.3703E+00  9.2295E-01  8.4043E-01  1.0782E+00  1.1115E+00  1.2098E+00  9.9163E-01  1.2144E+00  1.1055E+00
             1.2418E+00
 PARAMETER:  9.5834E-02  4.1500E-01  1.9825E-02 -7.3843E-02  1.7533E-01  2.0567E-01  2.9043E-01  9.1595E-02  2.9427E-01  2.0029E-01
             3.1655E-01
 GRADIENT:   1.2640E+01  3.7712E+01  1.3902E+00  2.3587E+01 -2.1642E+01  1.0833E+01  1.0617E+01  2.4814E+00  8.0560E+00  1.5700E+00
             6.2398E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2801.83839922248        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  9.8603E-01  1.3914E+00  9.1442E-01  7.9391E-01  1.1271E+00  1.0828E+00  1.1248E+00  9.4302E-01  1.2007E+00  1.1535E+00
             1.1931E+00
 PARAMETER:  8.5936E-02  4.3033E-01  1.0531E-02 -1.3079E-01  2.1963E-01  1.7953E-01  2.1759E-01  4.1336E-02  2.8293E-01  2.4281E-01
             2.7652E-01
 GRADIENT:  -4.5444E+00  1.2755E+00 -3.9512E-01 -9.3318E-01 -1.1627E+00  7.5152E-01  8.7167E-01  5.0991E-01  8.2853E-01  2.2015E-01
             2.9321E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2801.88580797644        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      615
 NPARAMETR:  9.8857E-01  1.3704E+00  9.1672E-01  8.0647E-01  1.1155E+00  1.0806E+00  1.1298E+00  9.0302E-01  1.1900E+00  1.1416E+00
             1.1913E+00
 PARAMETER:  8.8501E-02  4.1513E-01  1.3052E-02 -1.1508E-01  2.0928E-01  1.7749E-01  2.2207E-01 -2.0086E-03  2.7393E-01  2.3242E-01
             2.7502E-01
 GRADIENT:   2.7651E-01 -1.0163E-02  1.2390E-02  4.1581E-02  1.3339E-01 -2.3583E-02  1.6969E-02 -3.4301E-02 -4.5833E-02 -3.9336E-02
            -1.3935E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2801.88675566218        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  9.8909E-01  1.3695E+00  9.1710E-01  8.0667E-01  1.1151E+00  1.0817E+00  1.1298E+00  9.0373E-01  1.1905E+00  1.1414E+00
             1.1913E+00
 PARAMETER:  8.9029E-02  4.1442E-01  1.3461E-02 -1.1484E-01  2.0893E-01  1.7849E-01  2.2204E-01 -1.2248E-03  2.7436E-01  2.3222E-01
             2.7501E-01
 GRADIENT:   1.2818E+00 -4.4015E-01  7.3970E-02 -3.2633E-01  4.8274E-02  3.7366E-01 -3.9023E-03 -2.8722E-02  3.5996E-02  9.8565E-03
            -1.0227E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2801.88701695779        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  9.8899E-01  1.3690E+00  9.1724E-01  8.0705E-01  1.1147E+00  1.0815E+00  1.1300E+00  9.0434E-01  1.1903E+00  1.1410E+00
             1.1913E+00
 PARAMETER:  8.8930E-02  4.1409E-01  1.3612E-02 -1.1437E-01  2.0858E-01  1.7835E-01  2.2224E-01 -5.4930E-04  2.7417E-01  2.3189E-01
             2.7501E-01
 GRADIENT:   1.0943E+00 -3.2823E-01  8.2190E-02 -2.3968E-01 -3.1282E-02  3.1836E-01 -2.3456E-03 -2.4044E-02  3.2454E-02  1.1163E-02
            -8.2109E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2801.88732443156        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1162
 NPARAMETR:  9.8923E-01  1.3680E+00  9.1759E-01  8.0742E-01  1.1142E+00  1.0819E+00  1.1304E+00  9.0397E-01  1.1901E+00  1.1405E+00
             1.1912E+00
 PARAMETER:  8.9172E-02  4.1331E-01  1.3995E-02 -1.1391E-01  2.0813E-01  1.7870E-01  2.2259E-01 -9.5724E-04  2.7404E-01  2.3148E-01
             2.7497E-01
 GRADIENT:   1.5561E+00 -6.1143E-01  1.5550E-01 -5.0623E-01 -9.7861E-02  4.5542E-01 -3.1608E-03 -3.3861E-02  4.0432E-02  2.4647E-02
            -1.1922E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2801.88754833971        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1343
 NPARAMETR:  9.8909E-01  1.3676E+00  9.1779E-01  8.0781E-01  1.1138E+00  1.0817E+00  1.1306E+00  9.0401E-01  1.1898E+00  1.1401E+00
             1.1912E+00
 PARAMETER:  8.9028E-02  4.1305E-01  1.4217E-02 -1.1343E-01  2.0779E-01  1.7849E-01  2.2277E-01 -9.1326E-04  2.7378E-01  2.3114E-01
             2.7498E-01
 GRADIENT:   1.2828E+00 -3.9524E-01  2.1138E-01 -3.8141E-01 -2.2822E-01  3.7436E-01 -9.7059E-03 -3.9434E-02  2.0765E-02  1.2551E-02
            -1.0395E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2801.88779268837        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1526
 NPARAMETR:  9.8913E-01  1.3670E+00  9.1788E-01  8.0820E-01  1.1135E+00  1.0817E+00  1.1309E+00  9.0406E-01  1.1896E+00  1.1397E+00
             1.1912E+00
 PARAMETER:  8.9069E-02  4.1260E-01  1.4309E-02 -1.1295E-01  2.0748E-01  1.7855E-01  2.2302E-01 -8.5406E-04  2.7364E-01  2.3075E-01
             2.7495E-01
 GRADIENT:   1.3608E+00 -4.2913E-01  1.8480E-01 -3.2629E-01 -1.7084E-01  3.9682E-01 -5.8429E-03 -3.9320E-02  2.7387E-02  2.8977E-03
            -1.3685E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2801.88797142872        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1706
 NPARAMETR:  9.8909E-01  1.3666E+00  9.1789E-01  8.0845E-01  1.1132E+00  1.0817E+00  1.1311E+00  9.0480E-01  1.1895E+00  1.1394E+00
             1.1912E+00
 PARAMETER:  8.9028E-02  4.1232E-01  1.4320E-02 -1.1263E-01  2.0725E-01  1.7849E-01  2.2320E-01 -4.3473E-05  2.7352E-01  2.3054E-01
             2.7497E-01
 GRADIENT:   1.2811E+00 -4.3054E-01  1.4378E-01 -2.8105E-01 -1.6926E-01  3.7395E-01  6.6169E-03 -2.7970E-02  3.3646E-02  1.5033E-02
            -7.9711E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -2801.88806035502        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1802
 NPARAMETR:  9.8909E-01  1.3664E+00  9.1769E-01  8.0854E-01  1.1133E+00  1.0816E+00  1.1312E+00  9.0661E-01  1.1896E+00  1.1394E+00
             1.1912E+00
 PARAMETER:  8.9027E-02  4.1220E-01  1.4099E-02 -1.1252E-01  2.0729E-01  1.7849E-01  2.2332E-01  1.9565E-03  2.7361E-01  2.3054E-01
             2.7499E-01
 GRADIENT:   1.2795E+00 -5.7480E-01 -3.3622E-02 -1.8614E-01  5.2651E-02  3.7354E-01  4.2005E-02  4.8119E-03  7.7545E-02  4.0307E-02
             2.1670E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1802
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.5675E-04 -2.1511E-02 -2.0211E-02  1.7549E-02 -2.6081E-02
 SE:             2.9836E-02  2.3449E-02  1.1193E-02  2.4264E-02  2.3857E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8779E-01  3.5896E-01  7.0950E-02  4.6954E-01  2.7430E-01

 ETASHRINKSD(%)  4.5302E-02  2.1443E+01  6.2504E+01  1.8711E+01  2.0076E+01
 ETASHRINKVR(%)  9.0584E-02  3.8288E+01  8.5940E+01  3.3922E+01  3.6122E+01
 EBVSHRINKSD(%)  3.6721E-01  2.1334E+01  6.5981E+01  2.1104E+01  1.7784E+01
 EBVSHRINKVR(%)  7.3307E-01  3.8116E+01  8.8427E+01  3.7754E+01  3.2405E+01
 RELATIVEINF(%)  9.9263E+01  1.4958E+01  6.9121E+00  1.5680E+01  2.1483E+01
 EPSSHRINKSD(%)  2.4627E+01
 EPSSHRINKVR(%)  4.3189E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2801.8880603550233     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1515.3741138684816     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2801.888       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.89E-01  1.37E+00  9.18E-01  8.09E-01  1.11E+00  1.08E+00  1.13E+00  9.07E-01  1.19E+00  1.14E+00  1.19E+00
 


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
 #CPUT: Total CPU Time in Seconds,      171.516
Stop Time:
Sat Oct 23 18:53:43 CDT 2021
